#!/usr/bin/env python
"""S7 sub-benchmark: best-isoform-pair vs canonical scoring for ortholog/paralog discrimination.

Inputs: work/S7/{isoforms.parquet, hits.tsv} from scripts/S7_isoform_prep.py (all Ensembl 99
isoforms of the DEV-RefOG genes, one all-vs-all DIAMOND run, so canonical and isoform scores
come from the same search). Truth = S5 labels (work/S5/refog_pair_labels.parquet).

Scores per labelled cross-species gene pair:
  canonical = max bits(canonical a, canonical b); best_iso = max over isoform pairs;
  decoy (control, seeds 11/23/47) = max over decoy sets (canonical + the same number of
  non-canonical isoforms drawn from other genes of the same homology component, same species
  first; ancestra.isoforms.decoy_isoform_sets).
Queries (S5 definition): gene a in species A with >= 1 ortholog and >= 1 paralog candidate in B.
Hard queries: best paralog >= 0.9 x best ortholog by the *canonical* score;
"S5 hard" = the S5 hard queries of its 6 species pairs (whole-proteome sequence scores).
Top-1 per query = 1 if max ortholog score > max paralog score, 0.5 on a tie, else 0.
CIs: paired bootstrap over RefOGs (B=2000).
Output: results/S7_isoforms/isoform_benchmark.json (+ work/S7/pair_scores.parquet).
"""
import json, os, sys, time

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, f"{ROOT}/scripts")
from ancestra.isoforms import decoy_isoform_sets, gene_pair_scores, set_pair_scores  # noqa: E402
from S5_pair_benchmark import bootstrap_auc, query_table  # noqa: E402

W = f"{ROOT}/work/S7"
OUT = f"{ROOT}/results/S7_isoforms"
SEEDS = (11, 23, 47)
HARD = 0.9
VERT = {"Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Rattus_norvegicus", "Canis_familiaris",
        "Monodelphis_domestica", "Gallus_gallus", "Danio_rerio", "Tetraodon_nigroviridis"}


def top1(d, col):
    mo = d.loc[d.label == 1, col].max(); mp = d.loc[d.label == 0, col].max()
    return 1.0 if mo > mp else (0.5 if mo == mp else 0.0)


def query_rows(c, cols):
    rows = []
    for pr, d in c.groupby("pair", sort=False):
        for side in ("gene_a", "gene_b"):
            for g, q in d.groupby(side):
                o, p = q[q.label == 1], q[q.label == 0]
                if len(o) == 0 or len(p) == 0:
                    continue
                r = dict(pair=pr, side=side, gene=g, refog=q.refog_a.iloc[0], n_cand=len(q),
                         hard=bool(p.canonical.max() >= HARD * o.canonical.max()),
                         any_alt=bool(((q.best_iso > q.canonical)).any()))
                for col in cols:
                    r[col] = top1(q, col)
                rows.append(r)
    return pd.DataFrame(rows)


def boot(q, a, b, B=2000, seed=11):
    cl = q.refog.values; u = np.unique(cl); g = {k: np.flatnonzero(cl == k) for k in u}
    x, y = q[a].values, q[b].values
    rng = np.random.default_rng(seed); ds = []
    for _ in range(B):
        ix = np.concatenate([g[k] for k in rng.choice(u, len(u), replace=True)])
        ds.append(x[ix].mean() - y[ix].mean())
    return dict(diff=round(float(x.mean() - y.mean()), 4),
                ci95=[round(float(np.percentile(ds, 2.5)), 4), round(float(np.percentile(ds, 97.5)), 4)])


def summarize_queries(q, cols):
    out = dict(n=int(len(q)), n_refogs=int(q.refog.nunique()),
               n_with_any_isoform_gain=int(q.any_alt.sum()))
    for col in cols:
        out[f"acc_{col}"] = round(float(q[col].mean()), 4) if len(q) else None
    if len(q):
        out["best_iso_minus_canonical"] = boot(q, "best_iso", "canonical")
        out["decoy_minus_canonical"] = [boot(q, f"decoy{s}", "canonical") for s in SEEDS]
        out["best_iso_minus_decoy"] = [boot(q, "best_iso", f"decoy{s}") for s in SEEDS]
        out["xfam_minus_canonical"] = [boot(q, f"xfam{s}", "canonical") for s in SEEDS]
        out["best_iso_minus_xfam"] = [boot(q, "best_iso", f"xfam{s}") for s in SEEDS]
        out["changed_decisions"] = int((q.best_iso != q.canonical).sum())
        out["fixed"] = int((q.best_iso > q.canonical).sum()); out["broken"] = int((q.best_iso < q.canonical).sum())
    return out


def auc_block(c, cols, mask):
    out = {}
    for col in cols:
        y, s = c.label.values[mask], c[col].values[mask]
        out[col] = dict(auc=round(float(roc_auc_score(y, s)), 4) if 0 < y.mean() < 1 else None,
                        ci95=bootstrap_auc(c, col, mask, B=1000))
    return out


def within_query_auc(c, col):
    vals = []
    for side in ("gene_a", "gene_b"):
        for (pr, g), d in c.groupby(["pair", side]):
            if 0 < d.label.mean() < 1:
                vals.append(roc_auc_score(d.label, d[col]))
    return round(float(np.mean(vals)), 4), len(vals)


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    iso = pd.read_parquet(f"{W}/isoforms.parquet")
    h = pd.read_csv(f"{W}/hits.tsv", sep="\t", header=None, usecols=[0, 1, 2], names=["q", "t", "bits"])
    prot_gene = dict(zip(iso.uid, iso.ob_id))
    canonical = set(iso.uid[iso.canonical])
    lab = pd.read_parquet(f"{ROOT}/work/S5/refog_pair_labels.parquet")
    # orient: species in alphabetical order; one row per unordered gene pair
    sw = (lab.sp_a > lab.sp_b).values
    c = lab.copy()
    for x, y in (("gene_a", "gene_b"), ("sp_a", "sp_b"), ("refog_a", "refog_b")):
        c.loc[sw, [x, y]] = c.loc[sw, [y, x]].values
    c = c.drop_duplicates(["gene_a", "gene_b"]).reset_index(drop=True)
    c["pair"] = c.sp_a + ":" + c.sp_b
    pairs = list(zip(c.gene_a, c.gene_b))
    s = gene_pair_scores(h.q.values, h.t.values, h.bits.values, prot_gene, canonical, pairs)
    for col in ("canonical", "best_iso", "best_prot_a", "best_prot_b", "best_is_canonical"):
        c[col] = s[col].values
    # control: decoy isoform sets (uids), keyed by OrthoBench gene ID
    comp = pd.read_csv(f"{ROOT}/splits/orthobench_refog_components.tsv", sep="\t")
    r2c = {r: row["#component"] for _, row in comp.iterrows() for r in row["refogs"].split(",")}
    g_ref = dict(zip(lab.gene_a, lab.refog_a)); g_ref.update(dict(zip(lab.gene_b, lab.refog_b)))
    isosets = iso.groupby("ob_id").uid.apply(list).to_dict()
    canon_of = dict(zip(iso.ob_id[iso.canonical], iso.uid[iso.canonical]))
    fam = {g: r2c.get(g_ref.get(g), g_ref.get(g)) for g in isosets}
    spec = dict(zip(iso.ob_id, iso.species))
    hh = h.copy()
    for sd in SEEDS:
        dec = decoy_isoform_sets(isosets, canon_of, fam, spec, seed=sd)
        c[f"decoy{sd}"] = set_pair_scores(hh, dec, pairs)
    # second (non-adversarial) control: decoys from genes of OTHER homology components, same species
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        pool = {}
        for g, ps in isosets.items():
            for p in ps:
                if p != canon_of.get(g):
                    pool.setdefault(spec[g], []).append((fam[g], p))
        xdec = {}
        for g in sorted(isosets):
            k = len(isosets[g]) - 1
            cand = [p for f, p in pool.get(spec[g], []) if f != fam[g]]
            xdec[g] = [canon_of[g]] + (list(rng.choice(cand, size=min(k, len(cand)), replace=False)) if k and cand else [])
        c[f"xfam{sd}"] = set_pair_scores(hh, xdec, pairs)
    c["n_iso_a"] = c.gene_a.map(lambda g: len(isosets.get(g, [])))
    c["n_iso_b"] = c.gene_b.map(lambda g: len(isosets.get(g, [])))
    c["cluster"] = c.refog_a.astype(str)
    c.to_parquet(f"{W}/pair_scores.parquet", index=False)
    cols = ["canonical", "best_iso"] + [f"decoy{s}" for s in SEEDS] + [f"xfam{s}" for s in SEEDS]

    res = dict(benchmark="orthobench_dev", hard_ratio=HARD, seeds=SEEDS, n_pairs=int(len(c)),
               n_orth=int((c.label == 1).sum()), n_par=int((c.label == 0).sum()),
               n_species_pairs=int(c.pair.nunique()))
    # isoform evidence: how often does the best isoform pair differ from the canonical pair?
    ev = {}
    multi = (c.n_iso_a > 1) | (c.n_iso_b > 1)
    for name, m in (("orthologs", c.label == 1), ("paralogs", c.label == 0)):
        d = c[m & (c.canonical > 0)]
        dm = c[m & multi & (c.canonical > 0)]
        ev[name] = dict(n_scored=int(len(d)), frac_either_gene_multi_isoform=round(float(((d.n_iso_a > 1) | (d.n_iso_b > 1)).mean()), 4),
                        frac_best_pair_not_canonical=round(float((d.best_iso > d.canonical).mean()), 4),
                        frac_best_pair_not_canonical_among_multi=round(float((dm.best_iso > dm.canonical).mean()), 4),
                        frac_gain_gt_5pct=round(float((d.best_iso > 1.05 * d.canonical).mean()), 4),
                        median_rel_gain_when_differs=round(float(((d.best_iso / d.canonical) - 1)[d.best_iso > d.canonical].median()), 4)
                        if (d.best_iso > d.canonical).any() else None,
                        frac_canonical_missing_iso_found=round(float(((c.canonical == 0) & (c.best_iso > 0))[m].mean()), 4))
    res["isoform_evidence"] = ev
    # AUC (pooled; raw bits are not comparable across families, see S5 doc) + within-query AUC
    allm = np.ones(len(c), bool)
    res["auc_pooled"] = auc_block(c, cols, allm)
    vm = (c.sp_a.isin(VERT) & c.sp_b.isin(VERT)).values
    res["auc_pooled_vertebrate_pairs"] = auc_block(c, cols, vm)
    res["auc_within_query_mean"] = {col: within_query_auc(c, col) for col in cols}
    # queries
    q = query_rows(c, cols)
    q.to_parquet(f"{W}/query_rows.parquet", index=False)
    res["queries_all"] = summarize_queries(q, cols)
    res["queries_hard"] = summarize_queries(q[q.hard], cols)
    res["queries_hard_vertebrate"] = summarize_queries(
        q[q.hard & q.pair.map(lambda p: all(x in VERT for x in p.split(":")))], cols)
    # S5 hard queries (S5 whole-proteome sequence scores define hardness)
    s5keys = set()
    for f in sorted(os.listdir(f"{ROOT}/work/S5")):
        if not f.startswith("pair_scores_"):
            continue
        d = pd.read_parquet(f"{ROOT}/work/S5/{f}")
        A, B = d.sp_a.iloc[0], d.sp_b.iloc[0]
        for side, qq, ix, no, npar, so, spv in query_table(d):
            if spv >= HARD * so:
                gene = d.loc[ix[0], "gene_a" if side == "a" else "gene_b"]
                s5keys.add((":".join(sorted([A, B])), gene))
    q["s5_hard"] = [(p, g) in s5keys for p, g in zip(q.pair, q.gene)]
    res["queries_s5_hard"] = summarize_queries(q[q.s5_hard], cols)
    res["queries_s5_hard"]["n_s5_hard_keys"] = len(s5keys)
    res["per_pair_hard"] = {pr: dict(n=int(len(d)), canonical=round(float(d.canonical.mean()), 4),
                                     best_iso=round(float(d.best_iso.mean()), 4),
                                     decoy11=round(float(d.decoy11.mean()), 4), xfam11=round(float(d.xfam11.mean()), 4))
                            for pr, d in q[q.hard].groupby("pair") if len(d) >= 5}
    res["wall_s"] = round(time.time() - t0, 1)
    json.dump(res, open(f"{OUT}/isoform_benchmark.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "per_pair_hard"}, indent=1, default=str))


if __name__ == "__main__":
    main()
