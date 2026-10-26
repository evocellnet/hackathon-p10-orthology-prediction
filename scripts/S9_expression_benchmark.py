#!/usr/bin/env python
"""S9 exploratory sub-benchmark: does Bgee expression-profile similarity separate orthologs
from paralogs, and does it re-rank the sequence hard cases correctly?

Truth / candidates / hard cases are REUSED from S5 (scripts/S5_refog_labels.py ->
work/S5/refog_pair_labels.parquet; hard queries and decision rules from
scripts/S5_pair_benchmark.py: query_table, query_decisions, pooled_queries, bootstrap_auc).

Pairs
* primary (pre-registered): Homo_sapiens:Mus_musculus, Homo_sapiens:Gallus_gallus,
  Homo_sapiens:Danio_rerio, with the S5 per-pair tables (whole-proteome DIAMOND sequence
  scores, S5 hard flags) work/S5/pair_scores_*.parquet. Tetraodon is not in Bgee, so the
  S5 Danio/Tetraodon pair cannot be used.
* extended: all 28 pairs of the 8 Bgee vertebrates in OrthoBench, sequence score = canonical
  bitscore from the S7 all-vs-all DIAMOND over DEV-RefOG genes (work/S7/hits.tsv).
Expression score: Spearman of Bgee expression scores over Uberon/CL entities shared by the
two species and called in both genes (>= 5); secondary: Jaccard of present-call sets.
Shuffled control: gene <-> profile permuted genome-wide within each species (seeds 11/23/47).
Pre-registered primary endpoint: pooled hard-query top-1 gain of comb@0.9 (Spearman) over
sequence; negative if its 95% CI includes 0 or the shuffled gain is similar.
Outputs: results/S9_expression/expression_benchmark.json (+ per-pair parquet work/S9).
"""
import json, os, pickle, sys, time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, f"{ROOT}/scripts")
from ancestra.expression import load_profiles, pair_similarity, protein_to_gene, shared_columns  # noqa: E402
from S5_pair_benchmark import auc, bootstrap_auc, pooled_queries, query_decisions, query_table  # noqa: E402

BGEE = {"Homo_sapiens": "Homo_sapiens", "Pan_troglodytes": "Pan_troglodytes", "Mus_musculus": "Mus_musculus",
        "Rattus_norvegicus": "Rattus_norvegicus", "Canis_familiaris": "Canis_lupus_familiaris",
        "Monodelphis_domestica": "Monodelphis_domestica", "Gallus_gallus": "Gallus_gallus", "Danio_rerio": "Danio_rerio"}
PRIMARY = [("Homo_sapiens", "Mus_musculus"), ("Homo_sapiens", "Gallus_gallus"), ("Homo_sapiens", "Danio_rerio")]
SEEDS = (11, 23, 47)
MIN_SHARED = 5
HARD = 0.9
OUT = f"{ROOT}/results/S9_expression"
W = f"{ROOT}/work/S9"


def profiles():
    pk = f"{W}/profiles.pkl"
    if os.path.exists(pk):
        return pickle.load(open(pk, "rb"))
    prof = {sp: load_profiles(f"{ROOT}/data/extra/bgee_v15_2/{b}_expr_simple.tsv.gz", sp) for sp, b in BGEE.items()}
    pickle.dump(prof, open(pk, "wb"), protocol=4)
    return prof


def canonical_seq_scores():
    """Canonical-vs-canonical max bits per unordered OrthoBench-ID pair (S7 DIAMOND)."""
    h = pd.read_csv(f"{ROOT}/work/S7/hits.tsv", sep="\t", header=None, usecols=[0, 1, 2], names=["q", "t", "bits"])
    q = h.q.str.split("|", n=1, expand=True); t = h.t.str.split("|", n=1, expand=True)
    m = (q[0] == q[1]) & (t[0] == t[1]) & (q[0] != t[0])
    a, b = q[0][m].values, t[0][m].values
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    return pd.Series(h.bits[m].values).groupby([lo, hi]).max()


def candidates(A, B, lab, canon):
    s5 = f"{ROOT}/work/S5/pair_scores_{A}__{B}.parquet"
    if (A, B) in PRIMARY and os.path.exists(s5):
        c = pd.read_parquet(s5)[["gene_a", "gene_b", "sp_a", "sp_b", "refog_a", "refog_b", "label", "kind", "seq", "hard"]]
        return c.reset_index(drop=True), "S5_whole_proteome"
    c = lab[((lab.sp_a == A) & (lab.sp_b == B)) | ((lab.sp_a == B) & (lab.sp_b == A))].copy()
    sw = (c.sp_a != A).values
    for x, y in (("gene_a", "gene_b"), ("sp_a", "sp_b"), ("refog_a", "refog_b")):
        c.loc[sw, [x, y]] = c.loc[sw, [y, x]].values
    c = c.drop_duplicates(["gene_a", "gene_b"]).reset_index(drop=True)
    lo = np.where(c.gene_a < c.gene_b, c.gene_a, c.gene_b); hi = np.where(c.gene_a < c.gene_b, c.gene_b, c.gene_a)
    c["seq"] = canon.reindex(pd.MultiIndex.from_arrays([lo, hi])).fillna(0.0).values
    c["hard"] = False
    for side in ("gene_a", "gene_b"):
        for g, d in c.groupby(side):
            o, p = d[d.label == 1], d[d.label == 0]
            if len(o) and len(p) and p.seq.max() >= HARD * o.seq.max():
                c.loc[d.index, "hard"] = True
    return c[["gene_a", "gene_b", "sp_a", "sp_b", "refog_a", "refog_b", "label", "kind", "seq", "hard"]], "S7_devrefog_canonical"


def score_pair(A, B, c, prof, p2g, shuf):
    pa, pb = prof[A], prof[B]
    cols = shared_columns(pa, pb)
    ra = pa.rows([p2g[A].get(g) for g in c.gene_a]); rb = pb.rows([p2g[B].get(g) for g in c.gene_b])
    s = pair_similarity(pa, pb, ra, rb, MIN_SHARED, cols)
    c = pd.concat([c, s], axis=1)
    for sd in SEEDS:
        s2 = pair_similarity(shuf[(A, sd)], shuf[(B, sd)], ra, rb, MIN_SHARED, cols)
        c[f"shuf{sd}_spearman"] = s2.spearman.values; c[f"shuf{sd}_jaccard"] = s2.jaccard.values
    # null: random cross-species gene pairs with profiles (for conservation / break calls)
    rng = np.random.default_rng(11)
    xa = rng.choice(len(pa.genes), 4000); xb = rng.choice(len(pb.genes), 4000)
    null = pair_similarity(pa, pb, xa, xb, MIN_SHARED, cols).spearman.dropna().values
    return c, dict(n_shared_entities=int(len(cols[2])), null_n=int(len(null)),
                   null_q50=round(float(np.quantile(null, .5)), 4), null_q95=round(float(np.quantile(null, .95)), 4),
                   null_mean=round(float(null.mean()), 4)), null


def conservation(c, q50, q95):
    out = {}
    for name, d in (("orthologs", c[c.label == 1]), ("paralogs", c[c.label == 0]),
                    ("paralogs_within", c[(c.label == 0) & (c.kind == "within")]),
                    ("paralogs_deep", c[(c.label == 0) & (c.kind == "deep")])):
        x = d.spearman.dropna()
        if len(x) == 0:
            continue
        out[name] = dict(n=int(len(x)), median_spearman=round(float(x.median()), 4),
                         frac_conserved_gt_null95=round(float((x > q95).mean()), 4),
                         frac_break_le_null50=round(float((x <= q50).mean()), 4))
    return out


def query_rows(c, col, near):
    """S5 query_decisions over candidates with expression data; returns (summary, rows)."""
    qs = [x + (False,) for x in query_table(c)]
    return query_decisions(c, qs, col, near=near, return_rows=True)


def boot_diff(rows_a, rows_b, key, B=1000, seed=11):
    """Paired RefOG bootstrap of mean(rows_a[key]) - mean(rows_b[key]) (same query order)."""
    a = np.array([r[key] for r in rows_a], float); b = np.array([r[key] for r in rows_b], float)
    cl = np.array([r["refog"] for r in rows_a])
    u = np.unique(cl); g = {k: np.flatnonzero(cl == k) for k in u}
    rng = np.random.default_rng(seed); ds = []
    for _ in range(B):
        ix = np.concatenate([g[k] for k in rng.choice(u, len(u), replace=True)])
        ds.append(a[ix].mean() - b[ix].mean())
    return round(float(a.mean() - b.mean()), 4), [round(float(np.percentile(ds, 2.5)), 4), round(float(np.percentile(ds, 97.5)), 4)]


def evaluate(cs, label):
    """Pooled AUCs and hard-query decisions over a list of per-pair tables."""
    c = pd.concat(cs, ignore_index=True)
    c["cluster"] = c.refog_a.astype(str)
    has = c.spearman.notna().values
    e = c[has].reset_index(drop=True)
    for col in [x for x in e.columns if x.endswith("jaccard")]:
        e[col] = e[col].fillna(0.0)
    for col in [x for x in e.columns if x.startswith("shuf") and x.endswith("spearman")]:
        e[col] = e[col].fillna(0.0)
    h = e.hard.values
    res = dict(set=label, n_candidates=int(len(c)), n_with_expression=int(has.sum()),
               frac_with_expression=round(float(has.mean()), 4), n_orth=int((e.label == 1).sum()),
               n_par=int((e.label == 0).sum()), n_hard_pairs=int(h.sum()))
    res["auc"] = {}
    for col in ["seq", "spearman", "jaccard"]:
        res["auc"][col] = dict(all=auc(e.label, e[col]), all_ci95=bootstrap_auc(e, col, np.ones(len(e), bool)),
                               hard=auc(e.label[h], e[col][h]), hard_ci95=bootstrap_auc(e, col, h))
        if col != "seq":
            res["auc"][col]["shuf_all"] = [auc(e.label, e[f"shuf{s}_{col}"]) for s in SEEDS]
            res["auc"][col]["shuf_hard"] = [auc(e.label[h], e[f"shuf{s}_{col}"][h]) for s in SEEDS]
    # query-level decisions (hard queries recomputed on the expression-covered candidates;
    # 'all' = every query with >= 1 ortholog and >= 1 paralog candidate)
    res["queries"] = {}
    for qset, near in (("hard", HARD), ("all", 0.0)):
        allrows = {}
        for col in ["spearman", "jaccard"] + [f"shuf{s}_{m}" for s in SEEDS for m in ("spearman", "jaccard")]:
            rows = []
            for pr, d in e.groupby("pair", sort=False):
                _, r = query_rows(d.reset_index(drop=True), col, near)
                for x in r:
                    x.update(pair=pr, col=col)
                rows += r
            allrows[col] = rows
        pq = pooled_queries([{"_qrows": sum(allrows.values(), [])}])
        q = {}
        for m in ("spearman", "jaccard"):
            d = {k: pq[m][k] for k in ("n", "seq_acc", "syn_correct", "syn_wrong", "syn_abstain",
                                       "syn_precision_when_decides", "comb0.9_acc", "comb0.9_minus_seq",
                                       "comb0.9_minus_seq_ci95", "comb0.75_acc", "comb0.75_minus_seq",
                                       "comb0.75_minus_seq_ci95")} if m in pq else {}
            if d:
                d["expr_top1_acc_alone"] = round(d["syn_correct"] / max(1, d["n"]), 4)
                for r in ("0.9", "0.75"):
                    d[f"shuf_comb{r}_minus_seq"] = [pq[f"shuf{s}_{m}"][f"comb{r}_minus_seq"] for s in SEEDS]
                    d[f"shuf_comb{r}_minus_seq_ci95"] = [pq[f"shuf{s}_{m}"][f"comb{r}_minus_seq_ci95"] for s in SEEDS]
                    d[f"shuf_precision_when_decides"] = [pq[f"shuf{s}_{m}"]["syn_precision_when_decides"] for s in SEEDS]
                    # real - shuffled (paired, per seed)
                    d[f"comb{r}_real_minus_shuf"] = [boot_diff(allrows[m], allrows[f"shuf{s}_{m}"], f"comb{r}") for s in SEEDS]
            q[m] = d
        q["per_pair"] = {pr: dict(n=int(len([x for x in allrows["spearman"] if x["pair"] == pr])),
                                  seq_acc=_m([x["seq"] for x in allrows["spearman"] if x["pair"] == pr]),
                                  comb09_spearman=_m([x["comb0.9"] for x in allrows["spearman"] if x["pair"] == pr]),
                                  comb075_spearman=_m([x["comb0.75"] for x in allrows["spearman"] if x["pair"] == pr]),
                                  comb09_shuf11=_m([x["comb0.9"] for x in allrows["shuf11_spearman"] if x["pair"] == pr]))
                         for pr in e.pair.unique()}
        res["queries"][qset] = q
    return res


def _m(x):
    return round(float(np.mean(x)), 4) if len(x) else None


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True); os.makedirs(W, exist_ok=True)
    prof = profiles()
    p2g = {sp: protein_to_gene(f"{ROOT}/data/genomes/orthobench/{sp}.parquet") for sp in BGEE}
    shuf = {(sp, sd): prof[sp].permuted(sd + i) for i, sp in enumerate(BGEE) for sd in SEEDS}
    lab = pd.read_parquet(f"{ROOT}/work/S5/refog_pair_labels.parquet")
    canon = canonical_seq_scores()
    sps = list(BGEE)
    pairs = PRIMARY + [(a, b) for i, a in enumerate(sps) for b in sps[i + 1:] if (a, b) not in PRIMARY]
    per_pair, tabs = {}, {}
    for A, B in pairs:
        c, src = candidates(A, B, lab, canon)
        c, nul, _ = score_pair(A, B, c, prof, p2g, shuf)
        c["pair"] = f"{A}:{B}"
        tabs[(A, B)] = c
        c.to_parquet(f"{W}/pair_expr_{A}__{B}.parquet", index=False)
        e = c[c.spearman.notna()]
        hd = e.hard.values
        per_pair[f"{A}:{B}"] = dict(seq_source=src, n=int(len(c)), n_expr=int(len(e)), n_orth=int((e.label == 1).sum()),
                                    n_par=int((e.label == 0).sum()), n_hard=int(hd.sum()), **nul,
                                    auc_spearman=auc(e.label, e.spearman), auc_jaccard=auc(e.label, e.jaccard.fillna(0.0)),
                                    auc_seq=auc(e.label, e.seq), auc_spearman_hard=auc(e.label[hd], e.spearman[hd]),
                                    auc_seq_hard=auc(e.label[hd], e.seq[hd]),
                                    auc_spearman_shuf=[auc(e.label, e[f"shuf{s}_spearman"].fillna(0)) for s in SEEDS],
                                    conservation=conservation(c, nul["null_q50"], nul["null_q95"]),
                                    ortholog_breaks=c[(c.label == 1) & (c.spearman <= nul["null_q50"])][
                                        ["gene_a", "gene_b", "refog_a", "spearman", "n_shared"]].round(4).values.tolist())
        print(A, B, {k: v for k, v in per_pair[f"{A}:{B}"].items() if k not in ("ortholog_breaks", "conservation")},
              flush=True)
    res = dict(benchmark="orthobench_dev", min_shared=MIN_SHARED, hard_ratio=HARD, shuffle_seeds=SEEDS,
               primary_pairs=[f"{a}:{b}" for a, b in PRIMARY], per_pair=per_pair)
    res["pooled_primary"] = evaluate([tabs[p] for p in PRIMARY], "primary")
    res["pooled_extended"] = evaluate(list(tabs.values()), "extended_28_pairs")
    # conservation pooled over all pairs (ortholog profile conserved / break vs the per-pair null)
    cons = {}
    for k, v in per_pair.items():
        for grp, d in v["conservation"].items():
            x = cons.setdefault(grp, dict(n=0, conserved=0, brk=0))
            x["n"] += d["n"]; x["conserved"] += round(d["n"] * d["frac_conserved_gt_null95"]); x["brk"] += round(d["n"] * d["frac_break_le_null50"])
    res["conservation_pooled"] = {g: dict(n=x["n"], frac_conserved=round(x["conserved"] / x["n"], 4),
                                          frac_break=round(x["brk"] / x["n"], 4)) for g, x in cons.items()}
    res["wall_s"] = round(time.time() - t0, 1)
    json.dump(res, open(f"{OUT}/expression_benchmark.json", "w"), indent=1, default=str)
    for k in ("pooled_primary", "pooled_extended"):
        print(k, json.dumps({kk: vv for kk, vv in res[k].items() if kk != "queries"}))
        print(k, "HARD", json.dumps(res[k]["queries"]["hard"]["spearman"]))
        print(k, "ALL", json.dumps(res[k]["queries"]["all"]["spearman"]))
    print("conservation", res["conservation_pooled"])


if __name__ == "__main__":
    main()
