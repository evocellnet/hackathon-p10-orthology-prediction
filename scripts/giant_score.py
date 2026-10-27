"""Incremental S2 scoring for re-partitioned giant S0 components (G1/G3) and cross-piece rep edges (G2).

final families of an arm = the S2 cross-family rescue (`scoring.join_families`, min cov 0.5, e <= 1e-5,
joins capped at rescue_max_join = 7,500 as in the pipeline) applied to the new S0 partition, using the cross-family
hits saved by the original S2 run (pairs of different new S0 families among them). Only pairs that are new (same new family, different old S2-final
family) are scored: DIAMOND block self-search (same flags as S2: --sensitive -k0 --max-hsps 1 -b4 -c1, e<=10)
of the changed families packed into <= 100k-gene blocks; families <= 1000 genes also get the MMseqs2 SW
rescue of missed pairs (S2 hybrid rule). Scored once over the union of all arms' changed families
(G1 at 15k refines G3 at 20k: same Leiden ladder and seed), then each arm keeps its own pairs and the
S2 cap rule (families > 5000: top 2 per other species + 20 same-species per query) is re-applied.

    python scripts/giant_score.py {qfo,ob} THREADS
"""
import json
import os
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import load_proteomes, read_gene_table, write_fasta  # noqa
from ancestra.scoring import (ScoreTable, ScoringConfig, _iter_m8, _orient_dedup, _run, apply_cap,  # noqa
                              join_families, pack_blocks)

CFG = {"qfo": dict(run="work/rescue/qfo/base", prot="data/inputs/qfo2022/proteomes", out="work/giant/qfo"),
       "ob": dict(run="work/core_v0_ob/run_det", prot="data/orthobench/Open_Orthobench/BENCHMARKS/Input",
                  out="work/giant/ob")}


def union_labels(*labs):
    n = len(labs[0])
    offs, rows, cols = 0, [], []
    off = []
    for l in labs:
        off.append(offs)
        offs += int(l.max()) + 1
    # graph over label nodes; each gene links its labels
    for i in range(1, len(labs)):
        rows.append(labs[0] + off[0]); cols.append(labs[i] + off[i])
    r, c = np.concatenate(rows), np.concatenate(cols)
    g = coo_matrix((np.ones(len(r)), (r, c)), shape=(offs, offs))
    _, cc = connected_components(g, directed=False)
    out = cc[labs[0] + off[0]]
    _, out = np.unique(out, return_inverse=True)
    return out.astype(np.int64)


def main():
    bench, threads = sys.argv[1], int(sys.argv[2])
    C = CFG[bench]
    os.makedirs(C["out"], exist_ok=True)
    log = open(os.path.join(C["out"], "score_commands.log"), "a")
    t0 = time.time()
    gt = read_gene_table(os.path.join(C["run"], "stages/input/genes.tsv"))
    full = load_proteomes(C["prot"])
    assert full.names == gt.names
    gt.seqs = full.seqs
    fold = np.load(os.path.join(C["run"], "stages/scores/families_final.npy")).astype(np.int64)
    arms = {}
    z = np.load(os.path.join(C["run"], "stages/scores/cross_family_hits.npz"))
    cross = {k: z[k] for k in z.files}
    sc = ScoringConfig()
    for g in (15000, 20000):
        s0 = np.load(os.path.join(C["out"], f"s0_g{g}.npy"))
        _, s0 = np.unique(s0, return_inverse=True)
        arms[g], nj = join_families(s0.astype(np.int64), cross, sc.rescue_cov, sc.rescue_evalue, sc.rescue_max_join)
        arms[g] = arms[g].astype(np.int64)
        print(g, "S2 rescue joins", nj, flush=True)
    # scoring groups: union of both arms' families (G1 nests in G3 up to reassignment)
    grp = union_labels(arms[15000], arms[20000])
    # changed groups: hold genes of >= 2 old final families
    first = {}
    changed = np.zeros(int(grp.max()) + 1, bool)
    o = np.lexsort((fold, grp))
    gs, fs = grp[o], fold[o]
    diff = np.r_[False, (gs[1:] == gs[:-1]) & (fs[1:] != fs[:-1])]
    changed[np.unique(gs[diff])] = True
    sub = np.where(changed[grp], grp, -1)
    genes_ch = np.flatnonzero(sub >= 0)
    stats = dict(n_changed_groups=int(changed.sum()), n_genes_changed=int(len(genes_ch)),
                 max_group=int(np.bincount(grp).max()))
    print(stats, flush=True)
    cfg = ScoringConfig(threads=threads)
    lab = np.where(sub >= 0, sub, int(grp.max()) + 1 + np.arange(gt.n))   # unchanged genes -> singletons
    blocks = pack_blocks(lab, cfg.block_genes)
    stats["n_blocks"] = len(blocks)
    wd = os.path.join(C["out"], "tmp")
    os.makedirs(wd, exist_ok=True)
    keep_cols = []
    td = time.time()
    n_raw = 0
    for bi, genes in enumerate(blocks):
        fa = os.path.join(wd, f"b{bi}.fa")
        write_fasta(gt, fa, genes)
        _run(["diamond", "makedb", "--in", fa, "-d", fa + ".dmnd", "--threads", threads, "--quiet"], log)
        out = os.path.join(wd, f"b{bi}.m8")
        _run(["diamond", "blastp", "-q", fa, "-d", fa + ".dmnd", "-o", out, "-k", "0", "--max-hsps", "1",
              "-b", cfg.diamond_b, "-c", "1", "-e", cfg.evalue, "--threads", threads, "--quiet", "--outfmt", "6",
              "qseqid", "sseqid", "bitscore", "qcovhsp", "scovhsp", "pident", "evalue", cfg.diamond_sens], log)
        for c in _iter_m8(out, 7):
            q, t = c[0].astype(np.int64), c[1].astype(np.int64)
            n_raw += len(q)
            k = (lab[q] == lab[t]) & (fold[q] != fold[t])
            keep_cols.append([x[k] for x in c])
        for f in (fa, fa + ".dmnd", out):
            os.remove(f)
        print(f"block {bi + 1}/{len(blocks)} {len(genes)} genes, {time.time() - td:.0f}s", flush=True)
    stats["diamond_s"] = round(time.time() - td, 1)
    stats["n_raw_hits"] = n_raw
    c = [np.concatenate([x[i] for x in keep_cols]) for i in range(7)]
    c = _orient_dedup([c[0].astype(np.int64), c[1].astype(np.int64)] + c[2:])
    c[3], c[4], c[5] = c[3] / 100.0, c[4] / 100.0, c[5] / 100.0
    stats["n_new_pairs_raw"] = int(len(c[0]))
    np.savez(os.path.join(C["out"], "new_pairs.npz"), q=c[0].astype(np.int32), t=c[1].astype(np.int32),
             bits=c[2].astype(np.float32), qcov=c[3].astype(np.float32), tcov=c[4].astype(np.float32),
             pident=c[5].astype(np.float32), evalue=c[6].astype(np.float32))
    st = ScoreTable.load(os.path.join(C["run"], "stages/scores/scores.npz"))
    for g, name in ((15000, "G1"), (20000, "G3")):
        fam = arms[g]
        m = fam[c[0]] == fam[c[1]]
        cols = {k: np.concatenate([getattr(st, k), v[m].astype(getattr(st, k).dtype)])
                for k, v in zip(("q", "t", "bits", "qcov", "tcov", "pident", "evalue"), c)}
        u = ScoreTable(**cols, self_bits=st.self_bits)
        s2 = {}
        u = apply_cap(u, fam, gt.species_of, cfg, s2)
        u.save(os.path.join(C["out"], f"scores_{name}.npz"))
        np.save(os.path.join(C["out"], f"families_{name}.npy"), fam)
        sizes = np.bincount(fam)
        stats[name] = dict(n_pairs=int(len(u)), n_pairs_old=int(len(st)), n_new_kept=int(m.sum()),
                           dropped_by_cap=s2.get("n_pairs_dropped_by_cap"), n_families=int((sizes > 0).sum()),
                           max_family=int(sizes.max()), n_fam_gt5000=int((sizes > 5000).sum()),
                           n_fam_gt20000=int((sizes > 20000).sum()))
        print(name, stats[name], flush=True)
    stats["wall_s"] = round(time.time() - t0, 1)
    json.dump(stats, open(os.path.join(C["out"], "score_stats.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
