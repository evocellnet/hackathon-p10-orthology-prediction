# CLAIMS.md — ledger of quantitative claims

One row per quantitative claim made in `RESEARCH_TREE.md` (plus a few from `docs/` that the tree does not yet repeat, section 9).
Created from the RESEARCH_TREE state at milestone M2 (2026-10-26, T+5h11). Every evidence path was checked for existence
when the ledger was written (`scripts` one-liner at the bottom). **Existence of a file is not verification of the number in it.**

## How to read this

- **dev/test** = which data the number was measured on.
  - `dev` = OrthoBench dev RefOGs (34 RefOGs, `splits/orthobench_split.json`), toy simulator, or a benchmark that has no test split.
  - `test` = OrthoBench test RefOGs or any QfO benchmark. Only one test number exists so far (the S0 completeness, row C-S0-4; logged in `results/test_log.tsv`).
  - `n/a` = timing, data inventory or control; no dev/test split applies.
- **status**
  - `verified` = independently recomputed by an agent or person who did not produce it (recorded in the Reviewer column).
  - `pending` = not yet independently reviewed. The default text is **dev-only, pending review**.
  - `refuted` = falsified, or the pre-registered target was not met.
  - Rows can carry a qualifier in parentheses, e.g. `pending (superseded by C-CORE-12)`.
- **Rule:** a claim moves to `verified` only when a reviewer who did not produce the result reruns the named script (or recomputes at least one number) and writes their name and date in the Reviewer column. Nobody has done this yet, so **no row is `verified`**. (DESIGN.md section 6: "An agent that didn't produce a result reviews it before you rely on it.")
- Numbers are copied from `RESEARCH_TREE.md` as written. Where the evidence file shows something slightly different, the Note column says so.
- Dev numbers were tuned on dev; they are expected to shrink on test. Partition sensitivity of OrthoBench dev OG-F is about +-2 to 4 points (C-CORE-12), and 34 RefOGs give wide CIs.

## 1. Infrastructure and data

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-INF-1 | DIAMOND human-vs-human (20.6k proteins, 24 threads) `--sensitive` takes 12 s; MMseqs2 `-s 7.5` takes 88 s | `ENVIRONMENT.md` (timing table), `scripts/calib_timing.sh` | n/a | pending (dev-only, pending review) | Shared node; single runs. |
| C-INF-2 | The QfO 2022 benchmark runs locally (Darwin binary, no docker). Re-scored OrthoFinder and FastOMA submissions match the published numbers exactly on deterministic tests; STD tests within about 1.4 SE | `bench_repo/results/controls/qfo2022_OrthoFinder_vs_published.tsv`, `bench_repo/results/controls/qfo2022_FastOMA_v0.1.4_nov_vs_published.tsv`, `bench_repo/README.md` (QfO controls table) | n/a | pending (dev-only, pending review) | Positive control of the benchmark, not of Ancestra. |
| C-INF-3 | OrthoBench controls: RefOGs as orthogroups F=100, random groups F=0 | `bench_repo/results/controls/orthobench_controls.tsv`, `bench_repo/README.md` | n/a | pending (dev-only, pending review) | Official `benchmark.py` scoring. |
| C-INF-4 | 21 competitor prediction sets for QfO 2022 downloaded and md5-verified | `bench_repo/reference_results/qfo2022_participant_prediction_files.tsv`, `DATA.md` section 7 | n/a | pending (dev-only, pending review) | Table in `LITERATURE.md` section 7.4. |
| C-INF-5 | Ensembl 99 GFF covers the 12 OrthoBench species; Bgee v15.2 and QfO isoforms/CDS/gcoord (19 species) are available | `DATA.md` sections 2, 3, 9, `data/orthobench/ensembl_99/gff3_coverage.tsv` | n/a | pending (dev-only, pending review) | |
| C-INF-6 | Toy simulator and OrthoXML flattening identical to QfO `map_relations.py`; RefOGs scored through the official `benchmark.py` give F=100 | `tests/test_orthoxml.py`, `tests/test_simulate.py`, `docs/simulator.md`, `bench_repo/results/controls/orthobench_controls.tsv` | n/a | pending (dev-only, pending review) | |
| C-INF-7 | Design v0 (pure progressive pan-proteome merging) rejected: only about 2.5x saving over all-vs-all on QfO | `docs/reviews/r1_algorithms_compute.md`, `docs/reviews/r1_anc_sizes_fastoma_proxy.tsv` | n/a | pending (dev-only, pending review) | Compute-review estimate from a FastOMA proxy, not a measurement of Ancestra. |
| C-INF-8 | ProstT5 on CPU runs at 125-250 residues/s, so full 3Di for QfO would take 140-980 h (infeasible) | `docs/structure_feasibility.md` | n/a | pending (dev-only, pending review) | Measured on a node with load 40-90; CPU-s per residue is the robust quantity (about 0.09-0.10). |
| C-INF-9 | AFDB covers 90% of QfO accessions. Exact-sequence 3Di from Foldseek AFDB/UniProt: QfO 94.5% (928k proteins), OrthoBench 80.5% (202k). Foldcomp 3Di agrees with true 3Di at 97.9% of residues | `docs/structure_feasibility.md`, `data/structure/orthobench_3di_best.coverage.tsv` | n/a | pending (dev-only, pending review) | |

## 2. S0 family partition (docs/S0_families.md)

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-S0-1 | Pre-registered target (>=95% of dev RefOGs in one family) is not attainable by a sequence-similarity partition; even all DIAMOND hits at e<=1e-5 give 27/34; ceiling about 0.82 (28/34) because 6 dev RefOGs have members with no hit to the rest at e<=1e-3 | `results/S0_families/dev_ceiling.json`, `docs/S0_families.md` (Upper bound) | dev | refuted (target not met; deviation recorded) | The target is refuted, the ceiling measurement is pending. |
| C-S0-2 | Chosen config `2sA_dc30_cc_sh30_RA`: dev completeness 0.794 (27/34), fission 1.32, sum(g^2)/N^2 = 0.0075, max family 7,478, 263 s wall, 6.2 GB | `results/S0_families/selection.json`, `results/S0_families/summary.tsv`, `results/S0_families/2sA_dc30_cc_sh30_RA.json`, `results/S0_families/2sA_dc30_cc_sh30_RA_fresh.json` | dev | pending (dev-only, pending review) | Wall time is from the `_fresh` rerun (421 s in the first run; shared node). |
| C-S0-3 | S0 selection rule was fixed before the grid ran | `results/S0_families/configs.json`, `results/S0_families/selection.json` | dev | pending (dev-only, pending review) | Process claim. |
| C-S0-4 | Test completeness 0.722 for the chosen config (single logged test peek) | `results/S0_families/2sA_dc30_cc_sh30_RA.TEST.json`, `results/test_log.tsv` | **test** | pending (dev-only, pending review) | The only test number in the project so far. Not reproduced in REPORT.md. |
| C-S0-5 | Chosen config gave a byte-identical rerun | `results/S0_families/2sA_dc30_cc_sh30_RA_fresh.json` | dev | pending (superseded by C-CORE-10, C-CORE-11) | True for same threads only; later found non-deterministic across thread counts and fixed (C-CORE-11). |
| C-S0-6 | Rep-vs-rep linking (r2 proposal) gave 0.441 completeness; gene-vs-rep (N x R) search was the key | `results/S0_families/summary.tsv` | dev | pending (dev-only, pending review) | 0.441 is the 2s_dc30_cc_min50 row. The 2sR_dc30_cc_nocov row (rep-vs-rep, no guard) gives 0.706. Both are far below 2sA (0.794). |
| C-S0-7 | Root rescue for the 145 Leiden-split families and 15k singletons (6.1%) is still planned | `RESEARCH_TREE.md` (3.4), `results/S0_families/2sA_dc30_cc_sh30_RA.json` | dev | pending (dev-only, pending review) | Counts, not a result. |

## 3. Core engine (docs/core_v0.md)

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-CORE-1 | Toy simulator: pair F1 0.94-0.98 (precision 0.986-0.993); RBH-only arm 0.74 | `results/core_v0/toy_tables.md`, `results/core_v0/toy_grid.jsonl`, `docs/core_v0.md` section 5 | dev (toy) | pending (dev-only, pending review) | Toy results do not replace OrthoBench or QfO. |
| C-CORE-2 | Hybrid scoring (DIAMOND plus MMseqs2 re-align of in-family misses) recovers the 6-10% of pairs DIAMOND alone misses | `results/core_v0/toy_tables.md`, `docs/core_v0.md` section 5.3 | dev (toy) | pending (dev-only, pending review) | `docs/core_v0.md` gives 5-13% on the toy. |
| C-CORE-3 | 33k proteins / 32 species end-to-end in 35 s | `docs/core_v0.md` section 5.5 | dev (toy) | pending (dev-only, pending review) | |
| C-CORE-4 | v0 on OrthoBench dev: OG-F 68.1 (P 71.2, R 65.3; 11/34 exact); 14.7 min wall, 14.9 GB RSS; coverage guard off is the biggest lever (+9 F, 59.0 to 68.1) | `results/core_v0_orthobench/runs.tsv` (rows `B scores.npz default` and `B scores.npz min_cov=0`), `docs/core_v0.md` section 7.2 | dev | pending (dev-only, pending review) | Wall and RSS only in the doc. |
| C-CORE-5 | Dev error classes: merge splits 20/34, over-merges 13/34, S0 splits 7/34 | `docs/core_v0.md` section 7.3, `results/core_v0/diagnose_def32_s23_hybrid_avg.txt` | dev | pending (dev-only, pending review) | Diagnose files are toy-based. The OrthoBench counts are in the doc only. |
| C-CORE-6 | Relative-divergence in-paralog test: better on toy, -1.3 F on dev (not default) | `results/core_v0_orthobench/runs.tsv` (`B scores.npz rd_both min_cov=0` 66.79 vs 68.1), `results/core_v0/sweep_rd.tsv` | dev | pending (dev-only, pending review) | |
| C-CORE-7 | Root parsimony merge (R4): at most 68.3 dev F, off by default; in-family scores sit at the noise floor (median normalised score 0.058 within RefOG vs 0.044 outside) | `results/core_v0_orthobench/runs.tsv` (rows `D ...R4...`), `docs/core_v0.md` section 7.3 | dev | pending (dev-only, pending review) | |
| C-CORE-8 | Synteny R1-only on OrthoBench OG level is null: 68.01 vs 68.10 (shuffled = off) | `results/core_v0_orthobench/runs.tsv` (rows `E syn=...`) | dev | pending (dev-only, pending review) | Expected: OG metric is insensitive to partner choice (C-S5-8 for the pair-level metric). |
| C-CORE-9 | Hybrid local-tree mode: dev OG-F 72.5 (P 63.2, R 85.0), pair F1 0.740 (vs 68.1 / 0.730 seq-only); tree time about 16 s; root tree mode collapses precision (61); toy twoR12 Bilateria HOG pair F1 0.898 vs 0.863 | `docs/core_v0.md` section 7, `results/core_v0/sweep_tree.tsv`, `results/core_v0/sweep_tree2.tsv`, `results/S5_synteny/phase2/eval_base_run_v2.json` | dev | pending (superseded by C-CORE-12 for the partition) | Measured on the old non-deterministic partition `run_v2`. |
| C-CORE-10 | Constrained synteny plus R2/R3 veto plus R3 partner choice: dev OG-F 79.2 (P 73.9), pair F1 0.783 vs off 72.5 / 0.740; shuffled seeds 11/23/47 give 72.5/73.3/72.5 OG-F, 0.740 pair F1 | `docs/core_v0.md` section 7.2 and 7.6 (table), `results/core_v0_orthobench/runs.tsv` | dev | pending (superseded by C-CORE-12) | 79.2 is partition-dependent: S0-agent partition `run_v2`. The headline M2 number. Treat it as optimistic. |
| C-CORE-11 | S0 was not reproducible across thread counts (ARI 0.78 between two runs); that alone moved dev pair F1 0.78 to 0.64 (fresh e2e: OG-F 75.6, pair F1 0.636). Fixed: canonical input order, pinned `-b4 -c1`, sorted edges, Leiden seed 11; OrthoBench `families.tsv` byte-identical across 12/24 threads and input orders | `docs/core_v0.md` section 7.6, `tests/test_s0_determinism.py`, `work/S0det` (not committed) | dev | pending (dev-only, pending review) | The test file is the reproducible evidence. |
| C-CORE-12 | Deterministic partition, dev: default (constrained synteny + veto) OG-F 75.7 (P 72.2, R 79.6), macro pair F1 0.772; off 71.3 / 0.752; shuffled 71.3-72.0. Synteny beats shuffled on all 3 partitions (+3.7...+6.7 OG-F); partition sensitivity +-2-4 OG-F | `results/core_v0_orthobench/det_default_dev_eval.json`, `results/core_v0_orthobench/det_default_manifest.json`, `docs/core_v0.md` section 7.6-7.7 | dev | pending (dev-only, pending review) | Current best dev estimate. Macro pair F1 is the dev pair metric because micro is dominated by RefOG021. |
| C-CORE-13 | Locus-aware fission repair: pair F1 0.7996, S8 fission candidates -45%, but OG-F 75.5; shuffled flat; off by default, kept as a QfO ablation arm | `docs/core_v0.md` section 7.6 (table), `results/core_v0_orthobench/runs.tsv` | dev | pending (dev-only, pending review) | |
| C-CORE-14 | Leak test passes after moving RefOG helpers to `ancestra_eval/` (118 tests at that time; 132 passed on 2026-10-26 after the S0 determinism test, `pytest -m "not slow"`) | `tests/test_no_leak.py`, `ancestra_eval/families_eval.py` | n/a | pending (dev-only, pending review) | Current count: run `pytest`. |
| C-CORE-15 | End-to-end OrthoBench default: 28.5 min wall, 14.2 GB RSS (S0 435 s, scores 1,186 s, merge 75 s); deterministic run 31.5 min, 14.7 GB | `results/core_v0_orthobench/e2e_default_manifest.json`, `results/core_v0_orthobench/det_default_manifest.json` | n/a | pending (dev-only, pending review) | Shared node (load 45-60), 12 threads; upper bounds. det manifest: families 462 s, scores 1,365 s, merge 54 s, output 9 s, total 1,892 s. |

## 4. S1 species tree (docs/S1_speciestree.md)

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-S1-1 | NCBI taxdump tree: 0% split conflict vs the published QfO reference (evaluation only); nRF 0.128 because of polytomies | `results/S1_speciestree/metrics.json`, `results/S1_speciestree/summary.md` | dev (tree benchmark, evaluation reference) | pending (dev-only, pending review) | The published tree is read only by `scripts/S1_evaluate.py`. |
| C-S1-2 | Data-derived trees conflict with the published tree: ASTRAL-Pro 12%, FastME 17%, BIONJ 21% | `results/S1_speciestree/metrics.json`, `results/S1_speciestree/summary.tsv` | dev | pending (dev-only, pending review) | Doc values: 0.120 / 0.173 / 0.213. |
| C-S1-3 | Data-tree errors are classic long-branch artefacts (dog/cow+primates, C. elegans, Nematostella); Bacteria show 25-30% conflict (HGT, poorly resolved phyla) | `results/S1_speciestree/metrics.json` (`conflicting_splits`), `docs/S1_speciestree.md` section 3 | dev | pending (dev-only, pending review) | |
| C-S1-4 | Default guide tree = NCBI; ablation arms = NCBI+ASTRAL-refined, ASTRAL-Pro, FastME, NNI-perturbed (seeds 11/23/47), random (negative control) | `results/S1_speciestree/qfo/`, `results/S1_speciestree/ob/`, `docs/S1_speciestree.md` section 5 | n/a | pending (dev-only, pending review) | Decision, not a measurement. Downstream (full-pipeline) effects of the arms are not yet run. |

## 5. S5 synteny (docs/S5_synteny.md)

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-S5-1 | Gene-order coverage: OrthoBench about 100%; QfO 93% (zebrafish 0.64, amphioxus 0.08) | `docs/S5_synteny.md` section 1, `data/genomes/qfo2022/coverage.tsv`, `data/genomes/orthobench/coverage.json` | n/a | pending (dev-only, pending review) | QfO: 920,254 of 988,778 placed (0.931). |
| C-S5-2 | Synteny discriminates orthologs vs paralogs (OrthoBench dev, AUC at w=20): human/mouse 0.92, human/chicken 0.93, zebrafish/Tetraodon 0.82, human/zebrafish 0.70, human/fly 0.49; shuffled 0.50 | `results/S5_synteny/pair_benchmark_orthobench_dev.json`, `results/S5_synteny/tables.md` | dev | pending (dev-only, pending review) | Doc summary: AUC 0.64-0.93 per pair, 0.75 pooled at w=10 [0.71, 0.80]. |
| C-S5-3 | On 157 hard cases (sequence 50% correct), synteny decides 70 with 84% precision; re-rank gain +0.16 [0.05, 0.26]; shuffled -0.006 | `results/S5_synteny/pair_benchmark_orthobench_dev.json`, `docs/S5_synteny.md` section 0 | dev | pending (dev-only, pending review) | Bootstrap over RefOGs. |
| C-S5-4 | Naive additive bonus in the merge (toy): halves FP paralog pairs but loses 4-8 recall points | `results/S5_synteny/merge_preview_toy.json`, `docs/S5_synteny.md` section 5 | dev (toy) | pending (dev-only, pending review) | Doc: R 0.920 to 0.880. |
| C-S5-5 | Phase 2 pair-level dev metric: 19,146 cross-species pairs in 34 dev RefOGs, labelled from rooted RefOG trees | `results/S5_synteny/phase2/summary.json`, `scripts/orthobench_dev_pairs_eval.py`, `scripts/S5_refog_labels.py` | dev | pending (dev-only, pending review) | |
| C-S5-6 | Synteny R1-only, w10 lambda0.3: +0.032 [+0.009, +0.062] pair F1 over off | `results/S5_synteny/phase2/summary.json`, `results/S5_synteny/phase2/tables.md` | dev | pending (dev-only, pending review) | Config locked at the plateau (less selection bias than the best of 34). |
| C-S5-7 | w10 lambda0.5 core anchors +0.046 [+0.019, +0.072]; vert-teleost +0.092 [+0.029, +0.147] | `results/S5_synteny/phase2/summary.json`, `results/S5_synteny/phase2/tables.md` | dev | pending (dev-only, pending review) | Best of 34 dev configs: expect shrinkage. Doc range for vert-teleost: +0.071 to +0.092. |
| C-S5-8 | 102 shuffled runs stay at 0.716-0.730 pair F1 (seed spread <= 0.012); OG-level F is too noisy on 34 RefOGs for pair-level decisions (shuffled seeds span 66.9-69.8) | `results/S5_synteny/phase2/summary.json`, `docs/S5_synteny.md` section 8.3 | dev | pending (dev-only, pending review) | Keep rule passed. |

## 6. Competitor baselines (OrthoBench)

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-BASE-1 | OrthoFinder default: dev F 70.5 (P 60.8, R 84.0); all-70 F 70.1 | `bench_repo/results/orthobench_competitors/orthofinder/score.json` | dev | pending (dev-only, pending review) | Official `benchmark.py` identical (`verify.identical_to_official`). 8.2 CPU-h, 1.3 h wall at 8 threads on a shared node (`runtime.json`). |
| C-BASE-2 | OrthoFinder -M msa / mcl variant dev 54.2 | `bench_repo/results/orthobench_competitors/orthofinder_mcl/score.json` | dev | pending (dev-only, pending review) | |
| C-BASE-3 | SonicParanoid2 default dev 34.2 (all 47.4) | `bench_repo/results/orthobench_competitors/sonicparanoid_default/score.json` | dev | pending (dev-only, pending review) | Default mode only. |
| C-BASE-4 | T+12h gate passed early: Ancestra dev 75.6-79.2 (two partitions) vs OrthoFinder 70.5; paired bootstrap pending | `results/core_v0_orthobench/det_default_dev_eval.json`, `results/core_v0_orthobench/runs.tsv`, `bench_repo/results/orthobench_competitors/orthofinder/score.json` | dev | pending (dev-only, pending review) | Point estimates only, no CI: "inconclusive" under the DESIGN section 5 rule. The deterministic partition gives 75.7. |

## 7. Extra layers

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-S6-1 | Foldseek 3Di+AA rescue of S0 clusters lifts dev completeness 0.794 to 0.853 (29/34), above the 0.82 sequence ceiling; shuffled 3Di stays at 27/34 | `results/S6_structure/s0_rescue.json`, `results/S6_structure/s0_rescue.tsv`, `results/S6_structure/tables.md` | dev | pending (dev-only, pending review) | sum(g^2)/N^2 rises 0.0075 to 0.0178 and dev fusion rate 0.29 to 0.56 (`docs/S6_structure.md` Verdict). |
| C-S6-2 | Downstream, merging on rescued families drops OG-F 72.5 to 63.8 (RefOG010 absorbs 102 foreign genes) | `results/S6_structure/s0_downstream_scores.json`, `docs/S6_structure.md` | dev | pending (dev-only, pending review) | |
| C-S6-3 | 3Di scores in the merge (5 arms): real > shuffled in all, but no delta over off exceeds seed spread (best +0.015 [-0.004, 0.035]) | `results/S6_structure/merge_runs.jsonl`, `results/S6_structure/tables.md` | dev | pending (dev-only, pending review) | |
| C-S6-4 | 3Di ranks orthologs vs paralogs worse than sequence (deep AUC 0.694 vs 0.751); 654 of 656 structure-only pairs are paralogs | `results/S6_structure/pair_scores_seq.json`, `results/S6_structure/pair_scores_miss.json`, `results/S6_structure/discrimination.json` | dev | pending (dev-only, pending review) | |
| C-S6-5 | S6 cost about 1.35 CPU-h per S0 search | `results/S6_structure/s0_search.json`, `docs/S6_structure.md` | n/a | pending (dev-only, pending review) | |
| C-S8-1 | Of 2,523 apparent single-species mammal losses, 53% are not losses (609 annotated gene split into a target-only HOG plus 199 in another HOG = 808 orthology fission errors; 481 unannotated intact ORFs; 59 relocated); 25% confirmed (351 pseudogenized, 287 deleted with conserved flanks); 21% unresolved | `results/S8_geneloss/final_calls.tsv`, `results/S8_geneloss/summary.json`, `results/S8_geneloss/tables.md` | dev (no gold standard; run on the sequence-only OrthoBench run `run_v2`) | pending (dev-only, pending review) | Computed on the old partition `run_v2`, so the fission share is not re-measured on the deterministic partition. |
| C-S8-2 | Controls: hidden-gene recovery 92.2%, 0/512 random-window hits; loss-call precision 0.80 (human/mouse 0.97, dog 0.67) | `results/S8_geneloss/pos_Homo_sapiens.tsv`, `results/S8_geneloss/neg_Homo_sapiens.tsv`, `results/S8_geneloss/mask_Homo_sapiens.tsv`, `results/S8_geneloss/summary.json` | dev | pending (dev-only, pending review) | Per-species files exist for all 6 mammals. |
| C-S8-3 | Known biology recovered: GULO and UOX pseudogenes in human and chimp | `results/S8_geneloss/known_Homo_sapiens.tsv`, `results/S8_geneloss/known_Pan_troglodytes.tsv`, `results/S8_geneloss/homininae_check.tsv` | dev | pending (dev-only, pending review) | |
| C-S9-1 | Expression (Bgee v15.2), standalone on dev: AUC ortholog vs paralog 0.667 (shuffled about 0.5) | `results/S9_expression/expression_benchmark.json` | dev | pending (dev-only, pending review) | Primary 3 pairs: 0.667 [0.617, 0.721]; extended 28 pairs: 0.657. |
| C-S9-2 | Hard-case gain -0.019 [-0.20, 0.15] (primary) and -0.116 [-0.18, -0.05] (28 vertebrate pairs); when expression decides it is right only 40-47%; KILLED per pre-registered criterion | `results/S9_expression/expression_benchmark.json`, `docs/S9_expression.md` | dev | pending (dev-only, pending review) | Doc extended-set gain: -0.12 [-0.18, -0.05]. |
| C-S9-3 | Expression conservation: 16% of ortholog pairs conserved vs 26% broken (paralogs 4.5% vs 47%); zebrafish 29-45% breaks | `results/S9_expression/expression_benchmark.json`, `docs/S9_expression.md` | dev | pending (dev-only, pending review) | |
| C-S7-1 | Isoforms: best-isoform scoring AUC 0.760 vs 0.759 canonical; hard-case gain +0.013 [-0.018, 0.044]; costs 2-5x more search; KILLED | `results/S7_isoforms/isoform_benchmark.json`, `docs/S7_isoforms.md` | dev | pending (dev-only, pending review) | |
| C-S7-2 | For 21% of ortholog pairs the best-matching isoform pair is not the canonical pair | `results/S7_isoforms/isoform_benchmark.json`, `docs/S7_isoforms.md` | dev | pending (dev-only, pending review) | Orthologous-isoform evidence; not used by the pipeline. |

## 8. Full-pipeline runs in progress

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-QFO-1 | QfO run (arms A gene order, B none) on deterministic S0 with peek 1 on fast tests; **no QfO score exists yet** | `docs/qfo_run.md`, `results/test_log.tsv`, `results/qfo_run/` | test | pending (dev-only, pending review) | No numeric claim in the tree. Do not add one until a peek is logged. |
| C-QFO-2 | Scaling: nested QfO subsets 8/16/32/64/78 vs OrthoFinder 3.2 and SonicParanoid2, CPU-h primary, 6 h cap; **no result claimed yet** | `results/scaling/runs/`, `results/scaling/subsets/`, `scripts/scaling_queue.py`, `scripts/scaling_analyze.py` | n/a | pending (dev-only, pending review) | Partial runs exist (e.g. `results/scaling/runs/anc_8`, `anc_16`); not analysed. |

## 9. Claims in docs/ that RESEARCH_TREE does not yet repeat

| ID | Claim | Evidence file(s) | dev/test | Status | Note |
|---|---|---|---|---|---|
| C-DOC-1 | QfO dry run (old S0) S0->output about 2.1 h at 12 threads: families 36 min (3.2 CPU-h, 6.2 GB), scores 1.5 h (8.8 CPU-h, 29.9 GB), merge 131-348 s | `results/qfo_run/full_nogo_summary.json`, `results/qfo_run/full_go_summary.json`, `docs/qfo_run.md` section 5 | n/a | pending (dev-only, pending review) | Old non-deterministic S0; timing and memory only. Node load 28-65. |
| C-DOC-2 | QfO outputs pass all checks: XSD, geneRef integrity, TaxRange = LCA, flattening equal to QfO `map_relations.py`, pyham loads, OrthoBench reader accepts | `results/qfo_run/full_go_verify_orthoxml.json`, `results/qfo_run/full_nogo_verify_orthoxml.json`, `scripts/verify_orthoxml.py` | n/a | pending (dev-only, pending review) | |
| C-DOC-3 | Input check: 988,778 QfO proteins, 100% mapping coverage, 0 isoform leaks | `results/qfo_run/verify_ids.json`, `scripts/verify_ids.py` | n/a | pending (dev-only, pending review) | |
| C-DOC-4 | S1 runtime: arm A on QfO 757 s wall / 0.84 CPU-h; ASTRAL-Pro3 78 s; concatenation FastTree 926 s | `results/S1_speciestree/qfo/runtime.json`, `docs/S1_speciestree.md` section 4 | n/a | pending (dev-only, pending review) | Upper bounds, shared node. |

## Reviewer log

| Date | Reviewer | Rows verified | Command rerun | Outcome |
|---|---|---|---|---|
| (none yet) | | | | |

## Check that evidence paths exist

```bash
grep -o '`[^`]*`' CLAIMS.md | tr -d '`' | grep -E '^(results|docs|scripts|tests|bench_repo|data|splits|ancestra|ancestra_eval|work|DATA|RESEARCH|ENVIRONMENT)' \
  | sed -E 's/ \(.*//' | sort -u | while read p; do [ -e "$p" ] || echo "MISSING: $p"; done
```
`work/S0det` and some `data/` paths are gitignored and exist only in the working tree. Paths ending with `/` are directories.
