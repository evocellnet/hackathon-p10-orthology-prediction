# S0 giant component: Leiden cap and cross-piece edges (G-grid)

Status: done (05:42 CET, T+15h01). Verdict: **no arm passes on `fast`; nothing adopted; no `accurate` check needed.** Author: method-engineer agent, from 2026-10-27 03:45 CET (T+13h05). Dev only.
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
unconnected family: a sequence-sensitivity split, not a Leiden one), 34 neither. **Correction (after the scoring
step, before any merge result):** a first count said 59 / 77 records become co-familial at 15k / 20k; that came from
emulating the S2 joins as a transitive union with the old families (it chained pieces into a 74k-gene family on
OrthoBench) and is wrong. With the S2 join rule re-applied properly (`join_families`, reproduces the cached S2
families exactly), **0 of the 1,286 records become co-familial at 15k or 20k** (also 0 in the raw S0 partitions):
the within-component pairs sit in different Leiden communities at every resolution the ladder reaches, not in pieces
that the cap separated. So G1/G3 cannot repair any S0-split record; G2 can in principle reach the 665
within-component ones. OrthoBench: one split component of
118,703 genes, 153 pieces at 7.5k, 71 at 15k, 69 at 20k.

## 1. Results (dev only; `results/giant/summary.{md,json}`, raw `results/giant/runs.tsv`)
QfO: `fast` preset, TreeFam-A dev only (`qfo_treefam_dev_fast.score`). OrthoBench: `accurate` preset. FN classes =
missed TreeFam-A dev speciation pairs (`qfo_merge_iter.fn_breakdown`, with each arm's own families).

| arm | QfO TF-A dev TPR / PPV (fast) | ΔTPR / ΔPPV | FN: diff family / diff root HOG / paralog call | OB OG-F (P/R) | OB macro pair F1 | ΔOG-F / Δmacro | rule |
|---|---|---|---|---|---|---|---|
| G0 | 0.6151 / 0.9224 | – | 1286 / 2306 / 1794 | 79.69 (75.7/84.1) | 0.8064 | – | – |
| G1 | 0.6149 / 0.9214 | -0.0002 / -0.0010 | 1286 / 2408 / 1694 | 81.02 (80.6/81.5) | 0.8024 | +1.33 / -0.0040 | fail |
| G2 | 0.6240 / 0.9120 | +0.0089 / -0.0104 | 1286 / 2072 / 1903 | 79.31 (75.1/84.0) | 0.8023 | -0.38 / -0.0041 | fail |
| G3 | 0.6130 / 0.9209 | -0.0021 / -0.0015 | 1286 / 2284 / 1845 | 81.19 (80.9/81.5) | 0.8024 | +1.50 / -0.0040 | fail |


- **G1 / G3 (cap 15k / 20k)** fail on QfO (TPR −0.0002 / −0.0021). As predicted by §0 the S0-split class does not move
  (1,286). The OrthoBench OG-F gain (+1.3 / +1.5; precision 75.7 → 80.6–80.9, recall 84.1 → 81.5) is a side effect of
  the S2 cap rule, not of joined orthologs: pieces of ≤ 5,000 genes were scored exhaustively, but once inside a
  > 5,000-gene family they fall under the top-2-per-species rule, so the score table shrinks (OrthoBench 13.6 M → 8.7 M
  pairs; QfO 74.2 M → 56.9 M / 55.5 M). Per RefOG: RefOG021 FP 9,314 → 3,241 (but TP 6,793 → 5,609), RefOG004
  FP 12 → 0, RefOG023 splits (1 → 4 groups), RefOG002 FP 48 → 391; micro pair F1 0.804 → 0.786, macro −0.004.
- **G2 (cross-piece rep edges)** raises QfO TPR by +0.0089 but PPV drops by 0.0104, just over the 0.01 limit → **fail**.
  The gain does not come from S0-split records (still 1,286: the rep edges do not reach the *S. pombe* genes'
  orthologs) but from the different-root-HOG class (2,306 → 2,072), while paralog calls rise (1,794 → 1,903) and FP
  rise (722 → 841 records): the new edges link sibling pieces' units and over-merge. OrthoBench: OG-F −0.38, macro
  −0.0041 (RefOG002 FP 48 → 122), inside the OrthoBench limits.
- **Preset change during the grid.** The core agent changed `accurate` at 04:49 (`so_struct` off). OrthoBench G0/G1/G2
  started before, G3 was restarted with `--set so_struct=true` so that all four OrthoBench arms use the same
  (pre-04:49) `accurate` definition (`overrides` field in `runs.tsv`). The `fast` preset did not change. Also: OrthoBench
  G2 ran at 4 threads concurrently with other jobs for ~20 min (11 threads in total for that window, over the 8-thread
  allowance; a waiter-script bug).

## 2. Cost (shared node; RSS = peak of the process)
| step | QfO | OrthoBench |
|---|---|---|
| re-partition from cached S0 files (no search; Leiden on the 98k-rep component) | 5.6–6.0 min, 1 thread | 20 s |
| G1+G3 incremental S2 (new pairs only; DIAMOND on 402k genes in 30 changed groups, 5 blocks) | 78 min wall at 4 threads, **4.7 CPU-h**, **22.8 GB** | 14 min at 6 threads, 1.24 CPU-h, 14.3 GB |
| … a fresh pipeline run instead | S2 DIAMOND cost unchanged (blocks are packed to 100k genes, cost ∝ N·B, not Σg²); only the cap pruning and the pair count change | same |
| G2 edges (reuses S0's own stage-2 hits) | 34 s, 6.3 GB; 113k edges added (+0.15%) | 5 s; 1.8k edges |
| merge, `fast` (1 thread) | G0 356 s / 9.2 GB; G1 286 s / 7.1 GB; G2 349 s / 9.3 GB; G3 280 s / 6.9 GB | – |
| merge, `accurate` (3–4 threads) | – | 1,190–1,640 s wall, 3,585–3,770 CPU-s, ≤ 2.5 GB |
Peak memory of the grid: 22.8 GB, far below the 150 GB budget.

## 3. Verdict and recommendation
- **Nothing passes the rule on `fast`, so nothing goes to the `accurate` check and nothing is adopted.** Presets
  are unchanged.
- The "Leiden cap" hypothesis is refuted: none of the 1,286 S0-split records is separated *by the cap*. The 665
  within-component records sit in different Leiden communities at every resolution, and 587 others have one gene in a
  family that S0 never connected to the component (missing sequence edge, mostly *S. pombe*). Raising the cap only
  changes which pairs the S2 cap rule keeps.
- G2 is the only arm that moves QfO recall (+0.009) but it costs precision by just over the limit; it does not repair
  S0-split records. A guarded variant (cf. R4's RBH check) would be the natural follow-up, but it is outside the
  pre-fixed grid and the code freeze is close, so it is not proposed now.
- What remains for the S0-split class is the *S. pombe* sequence sensitivity (one gene outside the component) and
  community structure inside the component; neither is reachable without a new search or a different clustering
  objective — a post-freeze item.
