#!/usr/bin/env python
"""L2-network step (a) on TreeFam-A DEV (359 families; splits/treefamA_split.json): interolog and
abundance-free standalone signal on QfO 2022 species with dense own-species networks (yeast, fission
yeast, Arabidopsis, human, mouse, fly, worm). Extension of the pre-registered step (a) (docs/L2_network.md
§8 deviation 1): identical score, channels, controls and gate; only the benchmark differs.

Truth: TreeFam-A reference relations of dev families (S = ortholog, D = paralog), QfO 2022
ReconciledTrees (scripts/make_treefamA_split.py). Candidates = all reference pairs between two species.
Sequence score: DIAMOND --sensitive all-vs-all over the dev-family proteins (max bits both directions).
Anchors: 1:1 ortholog pairs of the frozen QfO accurate run (work/qfo2022_det/cand2_accurate/pairs.tsv);
anchors in the S0 family of a or b ignored. Network: STRING v12 own channels mapped via UniProt_AC aliases.
"""
import gzip, json, os, subprocess, sys, time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import L2_network_standalone as S  # noqa: E402

MAIN = S.MAIN
sys.path.insert(0, f"{MAIN}/scripts")
from make_treefamA_split import entry_to_acc, load_relations  # noqa: E402

RUN = f"{MAIN}/work/qfo2022_det/cand2_accurate"
SD = f"{MAIN}/data/extra/l2_network/string"
W = f"{MAIN}/work/l2_network/treefam"
STRING_TAX = {"9606": 9606, "10090": 10090, "3702": 3702, "7227": 7227, "6239": 6239, "9031": 9031,
              "7719": 7719, "7955": 7955, "8364": 8364, "559292": 4932, "284812": 284812, "9598": 9598,
              "10116": 10116}
CHANNELS = {k: S.CHANNELS[k] for k in ("exp150", "exp400", "coexp400", "own", "db400_CIRCULAR")}


def dev_labels():
    sp = json.load(open(f"{MAIN}/splits/treefamA_split.json"))
    g2f = {g: f for f, d in sp["families"].items() for g in d["members"]}
    sub = {f: d["subset"] for f, d in sp["families"].items()}
    inv = entry_to_acc()
    rows = []
    for x, y, e in load_relations():
        a, b = next(iter(inv[x])), next(iter(inv[y]))
        f = g2f.get(a)
        if f is None or sub[f] != "dev" or g2f.get(b) != f:
            continue
        rows.append((a, b, f, int(e == "S")))
    return pd.DataFrame(rows, columns=["gene_a", "gene_b", "refog_a", "label"])


def string_map(taxid, accs):
    al = pd.read_csv(f"{SD}/{taxid}.protein.aliases.v12.0.txt.gz", sep="\t", header=0, names=["sid", "alias", "source"],
                     dtype=str)
    al = al[al.source.str.contains("UniProt_AC") & al.alias.isin(accs)]
    m = {}
    for s_, a in zip(al.sid, al.alias):
        m.setdefault(s_, set()).add(a)
    return {k: sorted(v) for k, v in m.items()}


def edges(taxid, m):
    d = pd.read_csv(f"{SD}/{taxid}.protein.links.full.v12.0.txt.gz", sep=" ",
                    usecols=["protein1", "protein2", "experiments", "coexpression", "database"])
    d = d[(d[["experiments", "coexpression", "database"]] > 0).any(axis=1) & (d.protein1 < d.protein2)]
    a = d.protein1.map(m); b = d.protein2.map(m); ok = a.notna() & b.notna()
    x = d[ok].assign(a=a[ok], b=b[ok]).explode("a").explode("b")
    x = x[x.a != x.b]
    lo = np.where(x.a < x.b, x.a, x.b); hi = np.where(x.a < x.b, x.b, x.a)
    return pd.DataFrame({"a": lo, "b": hi, "exp": x.experiments.values, "coexp": x.coexpression.values,
                         "db": x.database.values}).groupby(["a", "b"], as_index=False).max()


def seq_scores(genes, sp_of):
    fa = f"{W}/dev.fa"; out = f"{W}/dev_allvsall.tsv"
    if not os.path.exists(out):
        want = set(genes)
        with open(fa, "w") as fo:
            for sp in sorted(set(sp_of[g] for g in genes)):
                keep = False
                for line in open(f"{MAIN}/data/inputs/qfo2022/proteomes/{sp}.fasta"):
                    if line.startswith(">"):
                        acc = line.split("|")[1] if "|" in line else line[1:].split()[0]
                        keep = acc in want
                        if keep:
                            fo.write(f">{acc}\n")
                    elif keep:
                        fo.write(line)
        subprocess.run(["diamond", "makedb", "--in", fa, "-d", f"{W}/dev", "--quiet"], check=True)
        subprocess.run(["diamond", "blastp", "-q", fa, "-d", f"{W}/dev", "-o", out, "--sensitive", "-k", "1000",
                        "-e", "1e-3", "-p", "4", "--quiet", "-f", "6", "qseqid", "sseqid", "bitscore"], check=True)
    h = pd.read_csv(out, sep="\t", header=None, names=["q", "t", "bits"])
    h = h[h.q != h.t]
    lo = np.where(h.q < h.t, h.q, h.t); hi = np.where(h.q < h.t, h.t, h.q)
    return pd.Series(h.bits.values).groupby([lo, hi]).max()


def main():
    t0 = time.time()
    os.makedirs(W, exist_ok=True)
    g = pd.read_csv(f"{RUN}/stages/input/genes.tsv", sep="\t")
    sp_of = dict(zip(g.protein, g.species))
    fam = dict(zip(g.protein, np.load(f"{RUN}/stages/scores/families_final.npy")))
    lab = dev_labels()
    lab["sp_a"] = lab.gene_a.map(sp_of); lab["sp_b"] = lab.gene_b.map(sp_of)
    lab = lab[lab.sp_a.notna() & lab.sp_b.notna() & (lab.sp_a != lab.sp_b)].reset_index(drop=True)
    lab["refog_b"] = lab.refog_a; lab["kind"] = "within"
    canon = seq_scores(sorted(set(lab.gene_a) | set(lab.gene_b)), sp_of)
    lo = np.where(lab.gene_a < lab.gene_b, lab.gene_a, lab.gene_b); hi = np.where(lab.gene_a < lab.gene_b, lab.gene_b, lab.gene_a)
    lab["seq"] = canon.reindex(pd.MultiIndex.from_arrays([lo, hi])).fillna(0.0).values
    sps = sorted(set(lab.sp_a) | set(lab.sp_b))
    # networks over whole proteomes (neighbours and anchors live outside the dev families)
    nets, cov = {}, {}
    for sp in sps:
        tx = STRING_TAX.get(sp.split("_")[1])
        if tx is None or not os.path.exists(f"{SD}/{tx}.protein.links.full.v12.0.txt.gz"):
            cov[sp] = None; continue
        accs = set(g.protein[g.species == sp])
        m = string_map(tx, accs)
        e = edges(tx, m)
        mapped = {a for v in m.values() for a in v}
        devg = set(lab.gene_a[lab.sp_a == sp]) | set(lab.gene_b[lab.sp_b == sp])
        cov[sp] = dict(string_taxid=tx, n_proteome=len(accs), mapped_frac=round(len(mapped) / len(accs), 4),
                       dev_genes=len(devg), dev_mapped_frac=round(len(devg & mapped) / max(1, len(devg)), 4))
        for ch, fn in CHANNELS.items():
            x = e[fn(e)]
            nodes = set(x.a) | set(x.b)
            cov[sp][f"nodes_{ch}"] = len(nodes)
            if len(nodes) >= S.MIN_NODES:
                nets[(sp, ch)] = (x.a.values, x.b.values, nodes)
        print(sp, cov[sp], flush=True)
    # anchors: 1:1 pairs of the frozen QfO accurate run (only among species we need)
    keep = set(sps)
    rows = []
    with open(f"{RUN}/pairs.tsv") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            a, b = line.split("\t", 2)[:2]
            sa, sb = sp_of.get(a), sp_of.get(b)
            if sa in keep and sb in keep:
                rows.append((a, b, sa, sb) if sa < sb else (b, a, sb, sa))
    p = pd.DataFrame(rows, columns=["a", "b", "sa", "sb"]).drop_duplicates(["a", "b"])
    anc = {}
    for (sa, sb), d in p.groupby(["sa", "sb"]):
        ca = d.a.value_counts(); cb = d.b.value_counts()
        d = d[d.a.map(ca).eq(1) & d.b.map(cb).eq(1)]
        anc[(sa, sb)] = dict(zip(d.a, d.b))
    adj, adj_shuf = {}, {}
    for (sp, ch), (a, b, nodes) in nets.items():
        adj[(sp, ch)] = S.adj_from(a, b)
        for s in S.SEEDS:
            adj_shuf[(sp, ch, s)] = S.adj_from(a, b, S.permutation(nodes, s * 1000 + sps.index(sp)))
    tabs, per_pair = [], {}
    for i, A in enumerate(sps):
        for B in sps[i + 1:]:
            c = lab[((lab.sp_a == A) & (lab.sp_b == B)) | ((lab.sp_a == B) & (lab.sp_b == A))].copy()
            if len(c) == 0:
                continue
            sw = (c.sp_a != A).values
            for x, y in (("gene_a", "gene_b"), ("sp_a", "sp_b")):
                c.loc[sw, [x, y]] = c.loc[sw, [y, x]].values
            c = c.drop_duplicates(["gene_a", "gene_b"]).reset_index(drop=True)
            M = anc.get((A, B), {})
            pp = dict(n=int(len(c)), n_orth=int(c.label.sum()), n_anchors=len(M))
            for ch in CHANNELS:
                ok = (A, ch) in adj and (B, ch) in adj
                if ok:
                    cnt, cos, na_, nb_ = S.interolog(c, adj[(A, ch)], adj[(B, ch)], M, fam)
                else:
                    cnt = cos = np.zeros(len(c)); na_ = nb_ = np.zeros(len(c), int)
                c[f"{ch}_count"] = cnt; c[f"{ch}_cos"] = cos; c[f"{ch}_eval"] = (na_ > 0) & (nb_ > 0)
                for s in S.SEEDS:
                    if ok:
                        k1, k2, _, _ = S.interolog(c, adj_shuf[(A, ch, s)], adj_shuf[(B, ch, s)], M, fam)
                    else:
                        k1 = k2 = np.zeros(len(c))
                    c[f"shuf{s}_{ch}_count"] = k1; c[f"shuf{s}_{ch}_cos"] = k2
                pp[ch] = dict(evaluable=round(float(c[f"{ch}_eval"].mean()), 4),
                              orth_count_gt0=round(float((cnt[c.label.values == 1] > 0).mean()), 4) if c.label.any() else None,
                              par_count_gt0=round(float((cnt[c.label.values == 0] > 0).mean()), 4) if (c.label == 0).any() else None)
            c["pair"] = f"{A}:{B}"
            per_pair[f"{A}:{B}"] = pp
            tabs.append(c)
    allc = pd.concat(tabs, ignore_index=True)
    allc["hard_global"] = S.mark_hard(allc.assign(gene_a=allc.pair + "|" + allc.gene_a, gene_b=allc.pair + "|" + allc.gene_b))
    allc.to_parquet(f"{W}/treefam_pairs.parquet", index=False)
    res = dict(benchmark="TreeFam-A dev (359 families), QfO 2022", n_candidates=int(len(allc)),
               n_orth=int(allc.label.sum()), n_hard_pairs_global=int(allc.hard_global.sum()),
               species_coverage=cov, per_pair=per_pair, coverage={}, evaluation={})
    for ch in CHANNELS:
        ev = allc[f"{ch}_eval"].values; ho = allc.hard_global.values
        res["coverage"][ch] = dict(frac_candidates_evaluable=round(float(ev.mean()), 4),
                                   frac_orth_pairs_evaluable=round(float(ev[allc.label == 1].mean()), 4),
                                   frac_hard_pairs_evaluable=round(float(ev[ho].mean()), 4), n_evaluable=int(ev.sum()))
        res["evaluation"][ch] = S.evaluate(allc, [f"{ch}_count", f"{ch}_cos"], f"{ch}_eval", ch)
        print(ch, json.dumps({k: v for k, v in res["evaluation"][ch].items() if k != "queries"})[:1200], flush=True)
    res["wall_s"] = round(time.time() - t0, 1)
    man = S.manifest(t0)
    man["script"] = os.path.abspath(__file__)
    man["inputs"].update(anchors=f"{RUN}/pairs.tsv", families=f"{RUN}/stages/scores/families_final.npy",
                         labels="splits/treefamA_split.json + QfO 2022 ReconciledTrees_TreeFam-A (make_treefamA_split)")
    res["manifest"] = man
    json.dump(res, open(f"{S.OUT}/standalone_treefam_dev.json", "w"), indent=1, default=str)
    print("wall", res["wall_s"])


if __name__ == "__main__":
    main()
