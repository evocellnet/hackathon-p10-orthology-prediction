# OrthoBench FINAL evaluation (one-time test; frozen freeze-v1)

Primary endpoint: OrthoBench-**test** F (36 RefOGs; split by homology component). Scores: official Open_Orthobench `benchmark.py` via `bench_repo/scripts/orthobench_eval.py` (subset = sums of per-RefOG normalised TP/FP/FN; all-70 verified identical to benchmark.py). ΔF = paired bootstrap over RefOGs (B=2000, seed 11), percentage points. Non-inferiority (DESIGN §5): *comparable* iff lower 95% bound of Δ(ours − comparator) ≥ −3 pp (−0.03 F). All configs logged in `results/test_log.tsv` before scoring.

## 1. Main table (test)

| method | test F | P | R | exact | dev F | all-70 F | CPU-h | wall-h | peak RSS GB | threads |
|---|---|---|---|---|---|---|---|---|---|---|
| **Ancestra accurate** (default; + gene order) | **76.2** | 69.5 | 84.2 | 10/36 | 79.7 | 77.9 | 3.91 | 0.53 | 14.0 | 10 |
| Ancestra fast (+ gene order)† | **74.8** | 68.7 | 82.1 | 9/36 | 75.7 | 75.2 | 2.88 | 0.40 | 14.0 | 10 |
| OrthoFinder 3.2 | **69.7** | 62.1 | 79.5 | 8/36 | 70.5 | 70.1 | 8.22 | 1.29 | 5.6 | 8 |
| SonicParanoid2 default | **70.4** | 70.8 | 70.1 | 8/36 | 34.2 | 47.4 | 7.85 | 1.31 | 8.1 | 8 |
| ProteinOrtho | **47.5** | 97.7 | 31.4 | 4/36 | 42.6 | 45.1 | 1.59 | 0.23 | 0.8 | 8 |
| FastOMA (LUCA, leaky) ⚑ | **33.5** | 24.3 | 53.7 | 10/36 | 29.6 | 31.7 | 12.74 | 3.50 | 13.6 | 8 |

† fast = S0+S2 stages of the main fresh run + its own merge/output (merge-only rerun on the same cache). ⚑ FastOMA uses the OMA-derived LUCA OMAmer DB (leakage; OMA contains these species) and is not a clean comparator. Competitor runtimes: 8 threads, node load 20–55; Ancestra: 10 threads (load in §3).

### Δ test F vs comparators (paired bootstrap over 36 test RefOGs) and non-inferiority verdict

| Ancestra | comparator | ΔF test (pp) | 95% CI | SE | verdict | ΔF dev | ΔF all-70 [95% CI] |
|---|---|---|---|---|---|---|---|
| accurate_go | OrthoFinder 3.2 | +6.5 | [-6.5, +18.1] | 6.4 | not shown non-inferior | +9.2 | +7.8 [+0.5, +15.1] |
| accurate_go | SonicParanoid2 default | +5.7 | [-4.3, +14.9] | 4.9 | not shown non-inferior | +45.5 | +30.5 [+5.4, +50.4] |
| accurate_go | ProteinOrtho | +28.7 | [+13.1, +41.9] | 7.3 | comparable (non-inferior; CI > 0) | +37.1 | +32.8 [+24.4, +40.7] |
| accurate_go | FastOMA (LUCA, leaky) | +42.7 | [+14.1, +55.2] | 11.2 | comparable (non-inferior; CI > 0) (flagged) | +50.1 | +46.2 [+26.3, +57.7] |
| fast_go | OrthoFinder 3.2 | +5.1 | [-7.5, +16.6] | 6.3 | not shown non-inferior | +5.2 | +5.1 [-1.7, +12.6] |
| fast_go | SonicParanoid2 default | +4.3 | [-6.0, +14.1] | 5.1 | not shown non-inferior | +41.6 | +27.9 [+2.9, +48.2] |
| fast_go | ProteinOrtho | +27.3 | [+12.5, +40.1] | 7.0 | comparable (non-inferior; CI > 0) | +33.1 | +30.2 [+21.6, +37.5] |
| fast_go | FastOMA (LUCA, leaky) | +41.3 | [+9.3, +55.9] | 12.3 | comparable (non-inferior; CI > 0) (flagged) | +46.1 | +43.5 [+22.0, +55.3] |

## 2. Ablations / layer effects (test primary; Δ vs main = accurate + gene order + NCBI tree)

| config | group | test F | P | R | exact | ΔF test vs main [95% CI] | dev F | all-70 F | pair F1 test macro (Δ [CI]) | micro |
|---|---|---|---|---|---|---|---|---|---|---|
| Ancestra accurate + gene order (MAIN, fresh e2e) | main | 76.2 | 69.5 | 84.2 | 10/36 | — | 79.7 | 77.9 | 0.885 | 0.865 |
| Ancestra fast + gene order | main | 74.8 | 68.7 | 82.1 | 9/36 | -1.4 [-7.1, +3.6] | 75.7 | 75.2 | 0.862 (-0.024 [-0.048, +0.001]) | 0.851 |
| accurate, no gene order | layer | 76.1 | 69.3 | 84.2 | 10/36 | -0.1 [-0.4, +0.0] | 78.7 | 77.4 | 0.882 (-0.003 [-0.008, +0.000]) | 0.864 |
| fast, no gene order | layer | 74.8 | 68.6 | 82.3 | 9/36 | -1.4 [-7.1, +3.7] | 71.3 | 73.0 | 0.863 (-0.022 [-0.046, +0.002]) | 0.854 |
| accurate + gene order + family rescue (3Di+profile) + guard | layer | 71.1 | 60.1 | 87.2 | 10/36 | -5.0 [-13.1, +0.0] | 79.9 | 75.2 | 0.885 (+0.000 [+0.000, +0.000]) | 0.865 |
| guide tree ASTRAL-Pro (armB, NCBI-rooted) | tree | 76.2 | 69.4 | 84.5 | 10/36 | +0.0 [-0.3, +0.4] | 76.4 | 76.3 | 0.884 (-0.001 [-0.006, +0.003]) | 0.865 |
| guide tree FastME (armA, NCBI-rooted) | tree | 76.2 | 69.4 | 84.5 | 10/36 | +0.0 [-0.3, +0.4] | 76.4 | 76.3 | 0.884 (-0.001 [-0.006, +0.003]) | 0.865 |
| guide tree NCBI + 1 NNI (seed 11) | tree | 76.2 | 69.5 | 84.2 | 10/36 | -0.0 [-0.0, +0.0] | 79.8 | 77.9 | 0.884 (-0.001 [-0.004, +0.000]) | 0.864 |
| guide tree NCBI + 1 NNI (seed 23) | tree | 72.9 | 67.8 | 78.9 | 10/36 | -3.3 [-8.5, -0.0] | 80.7 | 76.6 | 0.881 (-0.004 [-0.011, +0.001]) | 0.859 |
| guide tree NCBI + 1 NNI (seed 47) | tree | 76.5 | 69.3 | 85.3 | 10/36 | +0.3 [-0.0, +1.0] | 83.9 | 80.0 | 0.887 (+0.001 [-0.000, +0.005]) | 0.863 |
| guide tree random (seed 11; negative control) | tree | 59.8 | 70.3 | 52.0 | 5/36 | -16.4 [-27.1, -5.6] | 60.4 | 60.1 | 0.698 (-0.187 [-0.244, -0.135]) | 0.586 |

Gene order within the fast preset: ΔF test (fast+GO − fast no GO) = -0.0 [-0.1, +0.1].

## 3. Runtime (shared node, 10 threads; `/usr/bin/time -v` + 2-s process-tree sampling)

Fresh end-to-end run of the main config (`work/ob_final/main_accurate_go`):

| stage | wall (s) | CPU (s) | CPU-h | peak RSS tree (GB) | peak RSS own (GB) | node load1 mean |
|---|---|---|---|---|---|---|
| input | 0 | 0 | 0.00 | 0.53 | 0.2 | 10.4 |
| families | 259 | 1798 | 0.50 | 2.1 | 0.97 | 22.8 |
| scores | 1129 | 8516 | 2.37 | 14.28 | 12.81 | 30.8 |
| merge | 491 | 3745 | 1.04 | 6.34 | 2.41 | 26.4 |
| output | 11 | 10 | 0.00 | 2.03 | 1.63 | 30.5 |
| **total** (`time -v`) | 1893 | 14074 | 3.91 | max RSS 14.0 (largest process) | | |

Merge-only variants (cache of the main run; wall/CPU = whole CLI call incl. input, rescue, merge, output):

| variant | wall (s) | CPU-h | max RSS (GB) | merge wall (s) | exit |
|---|---|---|---|---|---|
| fast_go | 69 | 0.02 | 2.2 | 55.8 | 0 |
| accurate_nogo | 507 | 1.08 | 3.2 | 493.5 | 0 |
| fast_nogo | 35 | 0.01 | 2.1 | 21.7 | 0 |
| accurate_go_rescue | 1167 | 2.45 | 3.2 | 1031.6 | 0 |
| tree_astralpro | 518 | 1.04 | 3.2 | 504.3 | 0 |
| tree_fastme | 519 | 1.04 | 3.2 | 505.8 | 0 |
| tree_nni1_s11 | 525 | 1.04 | 3.2 | 511.5 | 0 |
| tree_nni1_s23 | 513 | 1.06 | 3.2 | 498.7 | 0 |
| tree_nni1_s47 | 513 | 1.04 | 3.2 | 500.2 | 0 |
| tree_random_s11 | 696 | 1.52 | 2.4 | 682.3 | 0 |

## 4. Output verification and reproducibility

- `scripts/verify_orthoxml.py` on `accurate_go`: overall **PASS** (xsd, generef_integrity, genes_once, taxrange_lca, pairs_count, flattening_vs_qfo, pyham, orthobench_export).
- `scripts/verify_orthoxml.py` on `fast_go`: overall **PASS** (xsd, generef_integrity, genes_once, taxrange_lca, pairs_count, flattening_vs_qfo, pyham, orthobench_export).
- Fresh S0 families and S2 scores are byte-identical to the pre-freeze deterministic run (`work/core_v0_ob/run_det`): stages/families/families.npy, stages/scores/scores.npz. Main dev F = 79.69, the pre-freeze `accurate` dev value (79.69) (regression check).
- Guide-tree arms reuse the tree-independent S0/S2 checkpoints through `scripts/run_cached.py` (frozen CLI, in-process key adoption, logged in `<run>/reused_stages.json`); same-tree variants hard-link the stages (as `scripts/rescue_ob_arms.sh`).

## 5. Reading (pre-registered rule applied as written)

- **Primary endpoint, OrthoBench-test F:** Ancestra `accurate` 76.2 vs OrthoFinder 3.2 69.7 (Δ +6.5 pp, 95% CI [-6.5, +18.1]) and SonicParanoid2 70.4 (Δ +5.7, [-4.3, +14.9]). The point estimates favour Ancestra, but the lower bounds are below −3 pp, so non-inferiority ("comparable") is **not shown** vs OrthoFinder or SonicParanoid2 on 36 test RefOGs. It is shown (with CI > 0) vs ProteinOrtho. The `fast` preset gives the same verdicts.
- The test-set SE of ΔF (≈5–6 pp = 0.05–0.06 F) is larger than the r2 simulation (0.02–0.04 F). Per-RefOG variance is dominated by a few large/paralog-rich RefOGs. SonicParanoid2 is far better on test (70.4) than on dev (34.2), so the dev-to-test drop of the competitors is not uniform.
- **Gene order (synteny layer):** no measurable test effect (accurate Δ −0.1 [−0.4, +0.0]; fast Δ −0.0 [−0.1, +0.1]). Its dev gain (+1.0 accurate, +4.4 fast) does not replicate on test.
- **Family rescue (3Di + profiles + guard):** test F −5.0 [−13.1, 0.0] (precision 69.5 → 60.1, recall 84.2 → 87.2), against +0.2 on dev; within-RefOG pair F1 is unchanged. So the extra genes attached to RefOG orthogroups are non-members. The decision to keep it off by default is supported.
- **Guide tree:** data-derived trees (ASTRAL-Pro, FastME; one wrong split vs NCBI) and 1-NNI perturbations change test F by ≤0.3, except NNI seed 23 (Gallus+Monodelphis clade; −3.3 [−8.5, −0.0]). The random tree (negative control) costs −16.4 [−27.1, −5.6] F and −0.19 macro pair F1, so the guide tree does real work. ASTRAL-Pro and FastME have the same topology but different Newick child order, and their outputs differ in 18 orthogroups (same scores). That is a small child-order sensitivity of the merge. Dev F moves by −3.3 to +4.2 under the same perturbations, which illustrates the dev volatility (DISCLOSURES 2).
- **Presets:** `fast` is −1.4 [−7.1, +3.6] test F and −0.024 [−0.048, +0.001] macro pair F1 vs `accurate`, for about 1 CPU-h less (merge 56 s vs 491 s).
- **Runtime:** Ancestra `accurate` end to end took 31.5 min wall and 3.9 CPU-h (10 threads, max RSS 14.0 GB, node load ≈ 23–31). OrthoFinder 3.2 took 1.29 h and 8.2 CPU-h (8 threads, load 20–55). The conditions differ (threads, load), so this is indicative only. S2 scoring is 61% of CPU.
- Pair-level test truth: the S5 labelling rule (`scripts/S5_refog_labels.py --subset test`), built for this final evaluation only (`work/ob_final/refog_pair_labels_test.parquet`; 36 RefOGs, 5,545 ortholog / 11,648 within-RefOG paralog pairs). Macro = mean over RefOGs of pair F1; Δ = paired bootstrap over RefOGs (B=2000, seed 11).

## OMA comparators (added at user request, post-freeze, 2026-10-27)

ΔF = F(Ancestra accurate) − F(comparator), paired bootstrap over RefOGs (B=2000, seed 11), NI margin 3 pp. Rows: `fastoma_v2` = FastOMA 0.5.1 + LUCA.h5 + gene-loss patch (pre-declared primary); `fastoma_v2_luca_released` = same run minus rescued HOGs (= unpatched 0.5.1); `oma_browser` = precomputed OMA browser HOGs (All.May2026) cut at Bilateria (primary), `oma_browser_root` = root level; `fastoma_luca` = v1 (unpatched) run.

| comparator | test F | P | R | exact | dev F | all-70 F | ΔF test (Ancestra − comparator) [95% CI] | SE | NI verdict (Ancestra vs comparator) | ΔF all-70 [95% CI] | CPU-h | wall-h | peak RSS GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fastoma_v2 | 33.1 | 23.9 | 53.8 | 10/36 | 26.7 | 29.9 | +43.1 [+14.0, +55.3] | 11.1 | comparable (non-inferior; CI > 0) | +48.0 [+24.3, +60.6] | 19.86 | 2.99 | 12.2 |
| fastoma_v2_luca_released | 33.1 | 23.9 | 53.8 | 10/36 | 25.1 | 29.1 | +43.1 [+14.0, +55.3] | 11.1 | comparable (non-inferior; CI > 0) | +48.8 [+25.3, +61.2] | 19.86 | 2.99 | 12.2 |
| oma_browser | 60.4 | 88.9 | 45.8 | 3/36 | 50.5 | 55.8 | +15.7 [+1.8, +29.2] | 7.2 | comparable (non-inferior; CI > 0) | +22.1 [+13.2, +30.7] | 0.31 | 0.26 | 1.5 |
| oma_browser_root | 7.2 | 3.9 | 49.5 | 2/36 | 7.3 | 7.3 | +69.0 [+22.3, +75.3] | 12.7 | comparable (non-inferior; CI > 0) | +70.6 [+54.6, +77.2] | 0.31 | 0.26 | 1.5 |
| fastoma_luca | 33.5 | 24.3 | 53.7 | 10/36 | 29.6 | 31.7 | +42.7 [+14.1, +55.2] | 11.2 | comparable (non-inferior; CI > 0) | +46.2 [+26.3, +57.7] | 12.74 | 3.50 | 13.6 |

- **LEAKAGE:** OMA (browser HOGs, and the OMAmer LUCA DB used by FastOMA) contains all 12 OrthoBench species/proteins, so these are informed upper-bound-style comparators, not independent ones. Ancestra did not use OMA.
- **Patch verification:** genes in predicted groups rose to 202,691/251,378 (80.6%) vs 198,615 (79.0%) in v1 (unpatched re-run 199,302); 3,795 subHOGs (630 of size>=2, 7,464 genes) were rescued. Big families partly recovered (RefOG005 TP 0.03 -> 3.8, RefOG006 0.06 -> 4.3, RefOG008 0.08 -> 4.3; RefOG013 unchanged). Test F is essentially unchanged (33.1 vs 33.5 v1): the remaining FastOMA gap is not the gene-dropping bug (most absent genes are OMAmer-unplaced / root singletons).
- **Cost (fastoma_v2, 3 segments, 8 threads, shared node):** 19.9 CPU-h, 2.99 h wall, peak RSS 12.2 GB (+~0.75 CPU-h omamer step reused from v1). oma_browser: conversion only, 0.31 CPU-h, 0.26 h, 1.5 GB (OMA inference not measurable). Both are not comparable to Ancestra's 3.9 CPU-h.
- Reading: Ancestra's point estimate is above every OMA row on test; vs oma_browser ΔF +15.7 [+1.8, +29.2] (CI above 0), with the caveat that OMA saw these species.
