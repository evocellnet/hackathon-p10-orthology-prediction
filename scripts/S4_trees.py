#!/usr/bin/env python
"""S4 benchmark: gene-tree inference methods for orthology on OrthoBench DEV RefOGs
(docs/S4_trees.md). Benchmark-only code (reads RefOG trees; never imported by the method).

Sets (units): each DEV RefOG alone ('single', 34) and each DEV multi-RefOG homology component
('comp', 4: RefOG pairs that are deep paralogs; reference = the two curated trees joined at
their S5 roots, so cross-RefOG pairs are paralogs by construction).
Scopes: 'c3di' = set genes with a 3Di string (all arms, same proteins); 'full' = all set
genes (arms without 3Di only).
Arms (single-threaded):
  a_nj      NJ on -ln S, S = pipeline norm_max score (run_v2 ScoreTable), missing -> max+1
            (ancestra.localtree.nj, read-only reuse; what localtree.py does, without outgroup)
  b_ft      FAMSA + FastTree -lg
  b_iq      FAMSA + IQ-TREE LG+G4 -fast (AA-only control for d_part)
  c_ft3di   AA-MSA columns mapped to 3Di + FastTree -trans (Foldseek-derived 3Di matrix)
  c_iq3di   same alignment + IQ-TREE 3DI_foldseek+G4 -fast
  d_part    AA+3Di partitioned IQ-TREE (LG+G4 | 3DI_foldseek+G4, proportional edges) -fast
  d_ftcat   FastTree -lg on concatenated AA+3Di (cheap proxy)
  e_fident  Foldseek 3Di+AA all-pairs alignment, d = -ln(fident), NJ  (FoldTree-style)
  e_fsbits  same alignments, d = -ln(bits / max(self_a, self_b)), NJ
  ref       the reference tree itself (positive control)
Rooting (all arms + reference): LCA reconciliation with the 12-species guide tree,
min duplications, then losses, then balance. Labels: species overlap.
Truth: 'matched' = reference tree pruned to the same leaves, same rooting/labelling;
'S5' = work/S5/refog_pair_labels.parquet (curated tree rooted with its outgroups).

    python scripts/S4_trees.py --workers 6
"""
import argparse, json, os, resource, shutil, subprocess, sys, time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import load_proteomes, read_gene_table  # noqa: E402
from ancestra.localtree import nj  # noqa: E402
from ancestra.merge import leaf_scores  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from ancestra import structure as STR  # noqa: E402
from ancestra.treebench import (SpeciesTree, UTree, nrf, pair_labels, quartet_distance,  # noqa: E402
                                root_overlap, root_reconcile)

RUN = f"{ROOT}/work/core_v0_ob/run_v2"
OB = f"{ROOT}/data/orthobench/Open_Orthobench"
W = f"{ROOT}/work/S4/sets"
OUT = f"{ROOT}/results/S4_trees"
MODEL_PAML = f"{ROOT}/work/S4/3DI_foldseek.paml"
MODEL_FT = f"{ROOT}/work/S4/3DI_foldseek.fasttree.tsv"
ARMS_ALL = ["a_nj", "b_ft", "b_iq", "c_ft3di", "c_iq3di", "d_part", "d_ftcat", "e_fident", "e_fsbits"]
ARMS_FULL = ["a_nj", "b_ft", "b_iq"]
VERT = {"Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Rattus_norvegicus", "Canis_familiaris",
        "Monodelphis_domestica", "Gallus_gallus", "Danio_rerio", "Tetraodon_nigroviridis"}
MAMMAL = {"Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Rattus_norvegicus", "Canis_familiaris",
          "Monodelphis_domestica"}

G = {}  # worker globals


def short(sp):
    return sp.split(".")[0]


def pair_class(a, b):
    if a in VERT and b in VERT:
        return "mammal" if (a in MAMMAL and b in MAMMAL) else "vert"
    if a in VERT or b in VERT:
        return "vert_invert"
    return "invert"


def cpu_children():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def cpu_self():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def run(cmd, cwd=None):
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))}: {r.stderr.decode()[-500:]}")


def write_fa(path, recs):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n{s}\n")


def read_fa(path):
    out, n, buf = {}, None, []
    for l in open(path):
        l = l.rstrip()
        if l.startswith(">"):
            if n is not None:
                out[n] = "".join(buf)
            n, buf = l[1:].split()[0], []
        else:
            buf.append(l)
    if n is not None:
        out[n] = "".join(buf)
    return out


def nj_tree(D, names):
    edges, n = nj(D)
    return UTree.from_edges(edges, n, names)


def dist_from_sim(S):
    """-ln S clipped like localtree (S in [1e-3, 1]); NaN -> max + 1."""
    d = -np.log(np.clip(S, 1e-3, 1.0))
    off = ~np.eye(len(S), dtype=bool)
    fin = np.isfinite(d) & off
    dmax = d[fin].max() if fin.any() else 1.0
    d = np.where(np.isfinite(d), d, dmax + 1.0)
    np.fill_diagonal(d, 0.0)
    return d


# ----------------------------------------------------------------------------
# worker: all arms for one (set, scope)
# ----------------------------------------------------------------------------

def do_task(task):
    key, names, seqs, ss, S_seq, arms = (task[k] for k in ("key", "names", "seqs", "ss", "S", "arms"))
    wd = os.path.join(W, key)
    os.makedirs(wd, exist_ok=True)
    res = {}
    ids = [f"g{i}" for i in range(len(names))]
    back = dict(zip(ids, names))

    def timed(fn):
        c0, s0, t0 = cpu_children(), cpu_self(), time.time()
        tree = fn()
        return tree, (cpu_children() - c0) + (cpu_self() - s0), time.time() - t0

    def a_nj():
        return nj_tree(dist_from_sim(S_seq), names)
    res["a_nj"] = timed(a_nj)

    aln_path = os.path.join(wd, "aa.aln")

    def msa():
        write_fa(os.path.join(wd, "aa.fa"), zip(ids, seqs))
        run(["famsa", "-t", "1", os.path.join(wd, "aa.fa"), aln_path])
        return read_fa(aln_path)
    aln, msa_cpu, msa_wall = timed(msa)
    res["_msa"] = (None, msa_cpu, msa_wall)

    def fasttree(path, extra):
        out = path + ".ft.nwk"
        with open(out, "w") as fh:
            r = subprocess.run(["FastTree", "-quiet", "-nopr", "-nosupport", *extra, path], stdout=fh,
                               stderr=subprocess.PIPE)
        if r.returncode != 0:
            raise RuntimeError(r.stderr.decode()[-500:])
        return UTree.from_newick(open(out).read(), rename=back.get)

    def iqtree(path, model=None, part=None, tag="iq"):
        pre = os.path.join(wd, tag)
        cmd = ["iqtree3", "-s", path, "-fast", "-T", "1", "--seed", "11", "--quiet", "--redo",
               "--prefix", pre]
        cmd += ["-p", part] if part else ["-m", model]
        run(cmd)
        return UTree.from_newick(open(pre + ".treefile").read(), rename=back.get)

    if "b_ft" in arms:
        res["b_ft"] = timed(lambda: fasttree(aln_path, ["-lg"]))
    if "b_iq" in arms:
        res["b_iq"] = timed(lambda: iqtree(aln_path, model="LG+G4", tag="b_iq"))
    if ss is not None:
        # 3Di alignment: AA-MSA columns mapped to 3Di characters
        tdi = {}
        for i, n in zip(ids, names):
            a, s, k, row = aln[i], ss[names.index(n)], 0, []
            for ch in a:
                if ch == "-":
                    row.append("-")
                else:
                    row.append(s[k])
                    k += 1
            assert k == len(s), (n, k, len(s))
            tdi[i] = "".join(row)
        p3 = os.path.join(wd, "3di.aln")
        write_fa(p3, tdi.items())
        pc = os.path.join(wd, "cat.aln")
        write_fa(pc, ((i, aln[i] + tdi[i]) for i in ids))
        L = len(aln[ids[0]])
        part = os.path.join(wd, "part.nex")
        with open(part, "w") as fh:
            fh.write(f"#nexus\nbegin sets;\n  charset aa = 1-{L};\n  charset tdi = {L + 1}-{2 * L};\n"
                     f"  charpartition mine = LG+G4:aa, {MODEL_PAML}+G4:tdi;\nend;\n")
        res["c_ft3di"] = timed(lambda: fasttree(p3, ["-trans", MODEL_FT]))
        res["c_iq3di"] = timed(lambda: iqtree(p3, model=f"{MODEL_PAML}+G4", tag="c_iq3di"))
        res["d_part"] = timed(lambda: iqtree(pc, part=part, tag="d_part"))
        res["d_ftcat"] = timed(lambda: fasttree(pc, ["-lg"]))

        def fs_all():
            fwd = os.path.join(wd, "fs")
            shutil.rmtree(fwd, ignore_errors=True)
            os.makedirs(fwd)
            write_fa(os.path.join(fwd, "aa.fasta"), ((str(k), s) for k, s in enumerate(seqs)))
            write_fa(os.path.join(fwd, "ss.fasta"), ((str(k), s) for k, s in enumerate(ss)))
            db = os.path.join(fwd, "db")
            run(["foldseek", "base:createdb", os.path.join(fwd, "aa.fasta"), db, "--shuffle", "0", "-v", "1"])
            run(["foldseek", "base:createdb", os.path.join(fwd, "ss.fasta"), db + "_sstmp", "--shuffle", "0", "-v", "1"])
            for suf in ("", ".index", ".dbtype"):
                os.replace(db + "_sstmp" + suf, db + "_ss" + suf)
            m = len(seqs)
            I, J = np.triu_indices(m, 0)
            h = STR.align_pairs(db, {k: k for k in range(m)}, I, J, os.path.join(fwd, "aln"), threads=1)
            return h
        h, fs_cpu, fs_wall = timed(fs_all)
        m = len(seqs)
        F = np.full((m, m), np.nan)
        B = np.full((m, m), np.nan)
        for q, t, b, fi in zip(h["q"], h["t"], h["bits"], h["pident"]):
            # keep the best-scoring alignment per pair
            if not np.isfinite(B[q, t]) or b > B[q, t]:
                B[q, t] = B[t, q] = b
                F[q, t] = F[t, q] = fi
        selfb = np.diag(B).copy()
        selfb = np.where(np.isfinite(selfb), selfb, np.nanmax(B, axis=1))
        Sb = B / np.maximum(selfb[:, None], selfb[None, :])
        t0 = cpu_self()
        res["e_fident"] = (nj_tree(dist_from_sim(F), names), fs_cpu + cpu_self() - t0, fs_wall)
        t0 = cpu_self()
        res["e_fsbits"] = (nj_tree(dist_from_sim(Sb), names), fs_cpu + cpu_self() - t0, fs_wall)
        res["_fs_pairs"] = (None, float(len(h["q"])), float(m * (m + 1) / 2))
    # MSA cost is part of every MSA-based arm
    out = {}
    for arm, (tree, cpu, wall) in res.items():
        if arm.startswith("_"):
            out[arm] = dict(cpu=cpu, wall=wall)
            continue
        if arm[0] in "bcd":
            cpu, wall = cpu + msa_cpu, wall + msa_wall
        with open(os.path.join(wd, f"{arm}.nwk"), "w") as fh:
            fh.write(to_newick(tree) + "\n")
        out[arm] = dict(cpu=cpu, wall=wall, newick=to_newick(tree))
    return key, out


def to_newick(t: UTree):
    leaves = t.leaves
    if not leaves:
        return ";"
    start = next(iter(t.adj[leaves[0]]))
    def rec(u, p):
        kids = [v for v in t.adj[u] if v != p]
        if not kids:
            return f"{t.name[u]}:{t.adj[u][p]:.6f}"
        s = "(" + ",".join(rec(v, u) for v in kids) + ")"
        return s + (f":{t.adj[u][p]:.6f}" if p is not None else "")
    return rec(start, None) + ";"


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", nargs="*", default=None, help="restrict to these set keys (debug)")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(W, exist_ok=True)
    t_start = time.time()
    split = json.load(open(f"{ROOT}/splits/orthobench_split.json"))
    dev = split["dev"]
    gt = read_gene_table(f"{RUN}/stages/input/genes.tsv")
    gp = load_proteomes(f"{OB}/BENCHMARKS/Input")
    assert gp.names == gt.names, "proteome order differs from the run_v2 gene table"
    idx = gt.index()
    sp_of = {n: short(gt.species[s]) for n, s in zip(gt.names, gt.species_of)}
    ss_all = STR.load_3di(f"{ROOT}/data/structure/orthobench_3di_best.fa", gt)
    stree = SpeciesTree(open(f"{RUN}/stages/input/guide_tree.nwk").read(), rename=short)
    st = ScoreTable.load(f"{RUN}/stages/scores/scores.npz")
    S = leaf_scores(gt, st, "norm_max")
    n = gt.n
    key_all = st.q.astype(np.int64) * n + st.t.astype(np.int64)
    o = np.argsort(key_all)
    key_s, S_s = key_all[o], S[o]
    del st, key_all, S

    def leaf_name(nm):
        parts = nm.split("_")
        for L in range(len(parts), 0, -1):
            for i in range(0, len(parts) - L + 1):
                g = "_".join(parts[i:i + L])
                if g in idx:
                    return g
        return nm

    def species_of_name(nm):
        return sp_of.get(nm) or "_".join(nm.split("_")[:2])

    # reference trees (curated), RefOG members
    refs, members = {}, {}
    for r in dev:
        genes = {l.strip() for l in open(f"{OB}/BENCHMARKS/RefOGs/{r}.txt") if l.strip()}
        t = UTree.from_newick(open(f"{OB}/Supporting_Data/Data_for_RefOGs/trees_tight/{r}.tre").read(),
                              rename=leaf_name)
        members[r] = sorted(genes & t.leaf_names())
        refs[r] = t
    comps = [c for c in split["multi_refog_components"] if all(x in dev for x in c)]
    sets = {r: [r] for r in dev}
    for c in comps:
        sets["+".join(c)] = c

    def set_reference(refogs, keep):
        if len(refogs) == 1:
            return refs[refogs[0]].prune(keep)
        assert len(refogs) == 2
        parts = []
        for r in refogs:
            t = refs[r]
            ch, root, info = root_overlap(t, species_of_name)
            a0, b0 = info["root_edge"]
            side = set(t._side(a0, b0)) & set(keep)
            tp = t.prune(set(keep) & t.leaf_names())
            edge = None
            for x, y in tp.edges():
                sx = set(tp._side(x, y))
                if sx == side or sx == tp.leaf_names() - side:
                    edge = (x, y)
                    break
            if edge is None:  # root fell on a pendant to pruned leaves: reconcile root
                _, _, inf2 = root_reconcile(tp, species_of_name, stree)
                edge = inf2["root_edge"]
            parts.append((tp, edge))
        return UTree.join(parts[0][0], parts[0][1], parts[1][0], parts[1][1])

    # S5 truth
    s5 = pd.read_parquet(f"{ROOT}/work/S5/refog_pair_labels.parquet")
    s5_lab = {tuple(sorted((x, y))): int(l) for x, y, l in zip(s5.gene_a, s5.gene_b, s5.label)}

    tasks, meta = [], {}
    for sk, refogs in sets.items():
        genes = sorted({g for r in refogs for g in members[r]})
        for scope in ("full", "c3di"):
            gs = genes if scope == "full" else [g for g in genes if g in ss_all]
            if len(gs) < 4:
                continue
            key = f"{sk}__{scope}"
            if a.only and key not in a.only and sk not in a.only:
                continue
            ga = np.array([idx[g] for g in gs], dtype=np.int64)
            m = len(gs)
            I, J = np.triu_indices(m, 1)
            kk = np.minimum(ga[I], ga[J]) * n + np.maximum(ga[I], ga[J])
            pos = np.minimum(np.searchsorted(key_s, kk), len(key_s) - 1)
            found = key_s[pos] == kk
            Sm = np.full((m, m), np.nan)
            Sm[I, J] = np.where(found, S_s[pos], np.nan)
            Sm[J, I] = Sm[I, J]
            np.fill_diagonal(Sm, 1.0)
            tasks.append(dict(key=key, names=gs, seqs=[gp.seqs[i] for i in ga],
                              ss=[ss_all[g] for g in gs] if scope == "c3di" else None, S=Sm,
                              arms=ARMS_ALL if scope == "c3di" else ARMS_FULL))
            meta[key] = dict(set=sk, kind="single" if len(refogs) == 1 else "comp", scope=scope,
                             refogs=refogs, n=m, n_set=len(genes), seq_pairs_missing=int((~found).sum()),
                             n_pairs=int(len(kk)))
    print(f"{len(tasks)} tasks; prep {time.time() - t_start:.0f}s", flush=True)
    tasks.sort(key=lambda t: -len(t["names"]))
    results = {}
    with Pool(a.workers) as pool:
        for key, out in pool.imap_unordered(do_task, tasks):
            results[key] = out
            print(key, {k: round(v["cpu"], 2) for k, v in out.items()}, flush=True)

    # ------------------------------------------------------------------ evaluation
    rows = []
    for key, out in results.items():
        mt = meta[key]
        gs = next(t["names"] for t in tasks if t["key"] == key)
        ref = set_reference(mt["refogs"], gs)
        rch, rr, rinfo = root_reconcile(ref, species_of_name, stree)
        truth_m = pair_labels(ref, rch, rr, species_of_name)
        arm_trees = {arm: UTree.from_newick(v["newick"]) for arm, v in out.items() if "newick" in v}
        arm_trees["ref"] = ref
        for arm, t in arm_trees.items():
            r_nrf, rf, mx, nl = nrf(t, ref)
            qd, nq = quartet_distance(t, ref)
            ch, rt, info = root_reconcile(t, species_of_name, stree)
            lab = pair_labels(t, ch, rt, species_of_name)
            row = dict(key=key, **{k: mt[k] for k in ("set", "kind", "scope", "n", "n_set", "seq_pairs_missing")},
                       arm=arm, nrf=r_nrf, rf=rf, rf_max=mx, qd=qd, n_quartets=nq, dups=info["dups"],
                       losses=info["losses"], cpu=out.get(arm, {}).get("cpu", np.nan),
                       wall=out.get(arm, {}).get("wall", np.nan))
            for tname, truth in (("m", truth_m), ("s5", s5_lab)):
                cnt = defaultdict(int)
                for p, y in lab.items():
                    if p not in truth:
                        continue
                    cls = pair_class(species_of_name(p[0]), species_of_name(p[1]))
                    z = truth[p]
                    for c in ("all", cls):
                        cnt[f"{tname}_{c}_n"] += 1
                        cnt[f"{tname}_{c}_ok"] += int(y == z)
                    cnt[f"{tname}_tp"] += int(y == 1 and z == 1)
                    cnt[f"{tname}_fp"] += int(y == 1 and z == 0)
                    cnt[f"{tname}_fn"] += int(y == 0 and z == 1)
                    cnt[f"{tname}_tn"] += int(y == 0 and z == 0)
                row.update(cnt)
            rows.append(row)
    df = pd.DataFrame(rows).fillna({c: 0 for c in rows[0] if c.endswith(("_n", "_ok", "_tp", "_fp", "_fn", "_tn"))})
    for c in df.columns:
        if c.endswith(("_n", "_ok", "_tp", "_fp", "_fn", "_tn")):
            df[c] = df[c].fillna(0).astype(int)
    df.to_csv(f"{OUT}/per_set.tsv", sep="\t", index=False)
    extra = {k: {a2: v for a2, v in out.items() if a2.startswith("_")} for k, out in results.items()}
    json.dump(dict(meta=meta, aux=extra, wall_total_s=time.time() - t_start, workers=a.workers),
              open(f"{OUT}/run_meta.json", "w"), indent=1)
    print(f"done {time.time() - t_start:.0f}s; rows {len(df)}")


if __name__ == "__main__":
    main()
