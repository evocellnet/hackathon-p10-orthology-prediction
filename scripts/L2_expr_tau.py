#!/usr/bin/env python
"""L2-genome idea 4: tissue-specificity (tau) agreement as a 1:many tie-breaker (docs/L2_genome.md §1.7).

Standalone only, OrthoBench DEV labelled pairs between the 10 Bgee species; sequence = S7 canonical
DIAMOND bits. Writes results/l2_genome/expr_tau.json.
"""
import json
import os
import pickle
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.expression import load_profiles, protein_to_gene  # noqa: E402
from S9_expression_benchmark import BGEE, canonical_seq_scores  # noqa: E402

BG = dict(BGEE, Drosophila_melanogaster="Drosophila_melanogaster", Caenorhabditis_elegans="Caenorhabditis_elegans")
NEAR, MARGIN, MIN_N, SEEDS = 0.9, 0.05, 5, (11, 23, 47)


def tau_of(p):
    s = p.score.astype(float)
    n = np.sum(~np.isnan(s), axis=1)
    mx = np.nanmax(np.where(np.isnan(s), -np.inf, s), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = np.nansum(1 - s / mx[:, None], axis=1) / (n - 1)
    t[(n < MIN_N) | ~(mx > 0)] = np.nan
    return dict(zip(p.genes, t))


def main():
    t0 = time.time()
    pk = os.path.join(ROOT, "work/S9/profiles.pkl")
    prof = pickle.load(open(pk, "rb"))
    for sp in ("Drosophila_melanogaster", "Caenorhabditis_elegans"):
        prof[sp] = load_profiles(os.path.join(ROOT, f"data/extra/bgee_v15_2/{BG[sp]}_expr_simple.tsv.gz"), sp)
    p2g = {sp: protein_to_gene(os.path.join(ROOT, f"data/genomes/orthobench/{sp}.parquet")) for sp in BG}
    tau = {}
    for sp in BG:
        tg = tau_of(prof[sp])
        tau[sp] = {pid: tg.get(g, np.nan) for pid, g in p2g[sp].items()}
        if sp == "Caenorhabditis_elegans":   # OrthoBench worm IDs are WBGene....1
            tau[sp].update({pid + ".1": tg.get(pid, np.nan) for pid in list(tg)})
    lab = pd.read_parquet(os.path.join(ROOT, "work/S5/refog_pair_labels.parquet"))
    lab = lab[lab.sp_a.isin(BG) & lab.sp_b.isin(BG)].reset_index(drop=True)
    canon = canonical_seq_scores()
    lo = np.where(lab.gene_a < lab.gene_b, lab.gene_a, lab.gene_b)
    hi = np.where(lab.gene_a < lab.gene_b, lab.gene_b, lab.gene_a)
    lab["seq"] = canon.reindex(pd.MultiIndex.from_arrays([lo, hi])).fillna(0.0).values
    lab["tau_a"] = [tau[s].get(g, np.nan) for s, g in zip(lab.sp_a, lab.gene_a)]
    lab["tau_b"] = [tau[s].get(g, np.nan) for s, g in zip(lab.sp_b, lab.gene_b)]
    cov = dict(genes_with_tau=round(float(np.mean(~np.isnan(np.r_[lab.tau_a.values, lab.tau_b.values]))), 4))
    # directed candidate table: query q (species A) -> candidates in species B
    d = pd.concat([lab.rename(columns=dict(gene_a="q", gene_b="c", sp_a="sq", sp_b="sc", tau_a="tq", tau_b="tc", refog_a="refog")),
                   lab.rename(columns=dict(gene_b="q", gene_a="c", sp_b="sq", sp_a="sc", tau_b="tq", tau_a="tc", refog_b="refog"))],
                  ignore_index=True)[["q", "c", "sq", "sc", "label", "seq", "tq", "tc", "refog"]]
    rng_perm = {}
    for sd in SEEDS:   # permute tau among genes within species
        rng = np.random.default_rng(sd)
        perm = {}
        for sp in BG:
            ks = sorted(tau[sp])
            vals = np.array([tau[sp][k] for k in ks])
            perm[sp] = dict(zip(ks, rng.permutation(vals)))
        rng_perm[sd] = perm
    arms = {"real": (d.tq.values, d.tc.values)}
    for sd, perm in rng_perm.items():
        arms[f"shuf{sd}"] = (np.array([perm[s].get(g, np.nan) for s, g in zip(d.sq, d.q)]),
                             np.array([perm[s].get(g, np.nan) for s, g in zip(d.sc, d.c)]))
    out = dict(coverage=cov, n_pairs=int(len(lab)), rule=dict(near=NEAR, margin=MARGIN, min_entities=MIN_N))
    groups = d.groupby(["q", "sc"]).indices
    for arm, (tq, tc) in arms.items():
        rows = []
        for (q, sc), ix in groups.items():
            lb = d.label.values[ix]
            if lb.min() == 1 or lb.max() == 0:
                continue
            sq_ = d.seq.values[ix]
            best_o, best_p = sq_[lb == 1].max(), sq_[lb == 0].max()
            hard = best_p >= NEAR * best_o
            seq_ok = int(lb[np.argmax(sq_)] == 1)   # ties -> first (as S5)
            near = sq_ >= NEAR * sq_.max()
            dt = np.abs(tq[ix] - tc[ix])
            ok = near & ~np.isnan(dt)
            dec = 0
            if ok.sum() >= 2 and near.sum() >= 2:
                o = np.argsort(dt[ok])
                if dt[ok][o[1]] - dt[ok][o[0]] >= MARGIN:
                    dec = 1 if lb[ok][o[0]] == 1 else -1
            comb = (1 if dec == 1 else 0) if dec != 0 else seq_ok
            rows.append((hard, seq_ok, dec, comb, d.refog.values[ix[0]]))
        r = pd.DataFrame(rows, columns=["hard", "seq", "dec", "comb", "refog"])
        e = {}
        for name, sub in (("hard", r[r.hard]), ("all", r)):
            dd = sub[sub.dec != 0]
            prec = float((dd.dec == 1).mean()) if len(dd) else np.nan
            rng = np.random.default_rng(11)
            u = sub.refog.unique()
            idx = {k: np.flatnonzero(sub.refog.values == k) for k in u}
            pb, gb = [], []
            for _ in range(1000):
                ix = np.concatenate([idx[k] for k in rng.choice(u, len(u))])
                s2 = sub.iloc[ix]
                d2 = s2[s2.dec != 0]
                if len(d2):
                    pb.append((d2.dec == 1).mean())
                gb.append(s2.comb.mean() - s2.seq.mean())
            e[name] = dict(n=int(len(sub)), seq_acc=round(float(sub.seq.mean()), 4), decides=int(len(dd)),
                           correct=int((dd.dec == 1).sum()), wrong=int((dd.dec == -1).sum()),
                           precision=round(prec, 4), precision_ci=[round(float(np.percentile(pb, 2.5)), 4),
                                                                  round(float(np.percentile(pb, 97.5)), 4)] if pb else None,
                           gain=round(float(sub.comb.mean() - sub.seq.mean()), 4),
                           gain_ci=[round(float(np.percentile(gb, 2.5)), 4), round(float(np.percentile(gb, 97.5)), 4)])
        out[arm] = e
        print(arm, e, flush=True)
    sh = [out[f"shuf{s}"]["hard"]["precision"] for s in SEEDS]
    lo_ci = out["real"]["hard"]["precision_ci"][0] if out["real"]["hard"]["precision_ci"] else np.nan
    out["stop_rule"] = dict(hard_precision_ci_lo=lo_ci, shuffled_precisions=sh,
                            verdict="PASS" if (lo_ci > 0.47 and all(lo_ci > x for x in sh if x == x)) else "STOP")
    out["manifest"] = dict(script="scripts/L2_expr_tau.py", bgee="data/extra/bgee_v15_2 (v15.2 expr_simple, gold+silver)",
                           labels="work/S5/refog_pair_labels.parquet", seq="work/S7/hits.tsv (canonical)",
                           runtime_s=round(time.time() - t0, 1))
    json.dump(out, open(os.path.join(ROOT, "results/l2_genome/expr_tau.json"), "w"), indent=1, default=float)
    print(out["stop_rule"])


if __name__ == "__main__":
    main()
