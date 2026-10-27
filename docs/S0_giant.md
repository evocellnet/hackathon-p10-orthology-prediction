# S0 giant component: Leiden cap and cross-piece edges (G-grid)

Status: in progress. Author: method-engineer agent, from 2026-10-27 03:45 CET (T+13h05). Dev only.
Scripts: `scripts/giant_repartition.py` (S0 stage-2 partition from cached files at another cap),
`scripts/giant_score.py` (incremental S2 scoring of the changed families), `scripts/giant_g2_edges.py`
(G2 cross-piece rep edges), `scripts/giant_merge_eval.py` (pipeline merge + dev evaluation).
Results: `results/giant/` (`runs.tsv` = one JSON line per run). Work: `work/giant/`. No `ancestra/` change.

## 0. Pre-registration (03:55 CET, before any arm result)
Base = current code, cached det S0/S2 (QfO `work/qfo_recall/base` stages; OrthoBench `run_det`). QfO `fast` first;
an arm passing there is confirmed with `accurate` on QfO dev before any recommendation. OrthoBench `accurate`.
Rule (PI): vs G0, TreeFam-A dev TPR up with PPV drop ≤ 0.01, OrthoBench dev OG-F ≥ −1, macro pair F1 ≥ −0.005.

| arm | S0 | S2 | merge |
|---|---|---|---|
| G0 | current (Leiden cap 7,500) | current | current preset |
| G1 | cap 15,000 (same Leiden ladder, seed 11, reassignment) | S2 rescue joins re-applied (`join_families`, max 7,500); new pairs scored with the S2 DIAMOND flags; cap rule for families > 5k (top 2 per other species + 20 same-species) | same |
| G2 | cap 7,500 | + cross-piece rep-vs-rep edges: stage-1 reps of the same split component in different families, from S0's own stage-2 hits (e ≤ 1e-5, span ≥ 0.3 of shorter), ≤ 2 per other species + 20 same-species per rep | same; the merge may join sibling-piece units by mutual best at any level (deviation: not restricted to the root level, see §1) |
| G3 | cap 20,000 | as G1 | same |

**Reproduction checks.** `giant_repartition.py` at cap 7,500 reproduces the cached `families.tsv` exactly on both
OrthoBench and QfO; `join_families` on it reproduces the S2-final families exactly (ARI 1.0, OrthoBench).

**Ceiling diagnosis (dev, before any merge).** QfO has one split component: 626,557 genes (63% of N) in 1,213 pieces
at 7.5k, **857 pieces at 15k, 764 at 20k** (the Leiden ladder picks the coarsest resolution that splits, and that
already yields hundreds of communities; the cap only stops the recursion). Of the 1,286 S0-split TreeFam-A dev FN
records, 665 have both genes in the component (different pieces), 587 have one gene outside it (a separate,
unconnected family: a sequence-sensitivity split, not a Leiden one), 34 neither. Raising the cap puts only
**59 (15k) / 77 (20k)** of the 1,286 records into one family — the maximum G1/G3 can repair (≤ +0.0055 TPR even if
all became TP). G2 can in principle reach all 665 within-component records. OrthoBench: one split component of
118,703 genes, 153 pieces at 7.5k, 71 at 15k, 69 at 20k.
