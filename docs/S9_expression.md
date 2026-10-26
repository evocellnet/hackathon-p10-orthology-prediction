# S9: expression-profile evidence (Bgee), exploratory

Status: 2026-10-26, T+2h50 to T+3h15 (`scripts/clock.sh`), well inside the 2 h box. Standalone species-pair sub-benchmark on OrthoBench DEV only. No core file was edited.

All numbers come from `results/S9_expression/expression_benchmark.json`. Checksums are in `results/S9_expression/manifest.json`.

**Verdict: KILL (negative).** The pre-registered negative criterion is met. Expression carries real ortholog signal overall (AUC 0.66 against 0.50 shuffled). On the sequence hard cases, however, it does not re-rank correctly: the gain is −0.02 with a 95% CI of [−0.20, +0.15] on the primary pairs, and −0.12 [−0.18, −0.05] on all 28 pairs.

## 1. Protocol (fixed before the run)
- **Truth, candidates and hard cases** are reused from S5:
  - Labels: `work/S5/refog_pair_labels.parquet` (`scripts/S5_refog_labels.py`): rooted DEV RefOG trees, species overlap, plus deep paralogs from the same homology component.
  - Query and decision code is imported from `scripts/S5_pair_benchmark.py` (`query_table`, `query_decisions`, `pooled_queries`, `bootstrap_auc`).
- **Pairs.**
  - *Primary* (pre-registered): human–mouse, human–chicken, human–zebrafish. These use the S5 per-pair tables (whole-proteome DIAMOND bits and S5 hard flags).
  - Tetraodon is not in Bgee, so zebrafish–Tetraodon is impossible.
  - *Extended*: all 28 pairs of the 8 Bgee vertebrates in OrthoBench (human, chimp, mouse, rat, dog, opossum, chicken, zebrafish). Sequence score = canonical bits from the S7 all-vs-all DIAMOND over the DEV-RefOG genes.
- **Expression score** (`ancestra/expression.py`), from Bgee v15.2 `expr_simple` calls with gold/silver quality, mapped OrthoBench protein → Ensembl gene through the gene-order tables.
  - Gene coverage of DEV genes is 0.89 to 1.00. 78% of primary and 86% of extended candidates have a score.
  - Entities = Uberon/CL IDs present in both species and called in both genes; at least 5 are required.
  - Primary measure: Spearman of the Bgee expression score. Secondary: Jaccard of the present-call sets.
  - Shared entities per pair: human–mouse 210, mouse–zebrafish 402, human–zebrafish 74, human–chicken 13. Most pairs have only 8 to 25, because the chicken, rat, chimp and opossum Bgee sets are small.
- **Decision rule**, the S5 rule: among a query's candidates with seq ≥ r × best seq, pick the highest expression similarity (r = 0.9 primary; 0.75 also reported).
  - Hard query: best paralog bits ≥ 0.9 × best ortholog bits, recomputed on the expression-covered candidates.
- **Shuffled control:** gene ↔ profile permuted genome-wide within each species (seeds 11, 23, 47).
- **CIs:** paired bootstrap over RefOGs.
- **Negative criterion** (pre-registered): the hard-case gain CI includes 0, or the shuffled gain is similar.
- **No circularity.** Bgee is used only as input evidence, never as labels. Labels come from the curated OrthoBench trees. OMA and others use Bgee downstream, but that does not feed OrthoBench.

## 2. Results

### Pooled
| set | candidates (w/ expr) | AUC Spearman [95% CI] | AUC Jaccard | AUC shuffled (3 seeds) | AUC hard: seq / Spearman [CI] |
|---|---|---|---|---|---|
| primary (3 pairs) | 2,123 (1,665) | **0.667** [0.617, 0.721] | 0.612 | 0.52 / 0.47 / 0.48 | 0.620 / 0.578 [0.467, 0.714] |
| extended (28 pairs) | 18,059 (15,562) | **0.657** [0.618, 0.700] | 0.648 | 0.49 / 0.50 / 0.49 | 0.675 / 0.586 [0.526, 0.665] |

### Hard queries, top-1 partner choice (Spearman)
| set | n | seq top-1 | expr alone: correct / wrong (precision) | shuffled precision | comb@0.9 − seq [CI] | comb@0.75 − seq [CI] | shuffled comb@0.9 − seq | real − shuffled (comb@0.9) [CI] |
|---|---|---|---|---|---|---|---|---|
| primary | 53 | 0.509 | 25 / 28 (0.47) | 0.40 / 0.16 / 0.17 | **−0.019 [−0.200, +0.149]** | −0.019 [−0.25, +0.21] | 0.00 / −0.13 / −0.30 | −0.02 [−0.18, +0.12]; +0.11 [−0.13, +0.37]; +0.28 [+0.13, +0.43] |
| extended | 484 | 0.552 | 192 / 287 (0.40) | 0.30 / 0.23 / 0.25 | **−0.116 [−0.183, −0.052]** | −0.134 [−0.224, −0.044] | −0.13 / −0.13 / −0.20 | +0.02 [−0.09, +0.11]; +0.01 [−0.05, +0.09]; +0.09 [+0.01, +0.15] |

All queries with ≥1 ortholog and ≥1 paralog candidate show the same picture:
- Primary (n=326): seq 0.920; comb@0.9 −0.003 [−0.029, +0.029]; comb@0.75 −0.071 [−0.118, −0.007].
- Extended (n=2,990): seq 0.927; comb@0.9 −0.019 [−0.033, −0.005].

Per primary pair (hard queries, seq → comb@0.9):
- human–mouse: 0.60 → 0.80 (n=5)
- human–chicken: 0.45 → 0.45 (n=11)
- human–zebrafish: 0.51 → 0.46 (n=37). This is the ohnolog case expression was supposed to resolve.

**Reading.**
- Expression similarity is real ortholog evidence: AUC 0.66 against a 0.50 shuffled control, and above shuffled in every one of the 28 pairs.
- It is much weaker than sequence, and very weak exactly where it would be needed. On hard pairs the AUC is 0.58 to 0.59. Among competing near-equal paralogs, especially teleost ohnologs and 2R paralogs, the paralog is about as often the closer profile as the ortholog: precision 0.40 to 0.47.
- Real beats shuffled when expression decides (precision 0.40 against 0.23 to 0.30), so the signal is genuine. But it is below the 0.5 to 0.55 that sequence already achieves, so using it as a tie-breaker lowers accuracy.
- **Pre-registered criterion met → negative.** Not kept as a layer.

## 3. Conservation and conservation breaks (project brief)
For each species pair, the null is the Spearman similarity of about 3,000 random cross-species gene pairs (genes with profiles). The null median is about 0 to 0.09 and the 95th percentile is 0.49 to 0.90; the 95th percentile is higher when few entities are shared.
- **Conserved:** ρ above the null 95th percentile.
- **Break:** ρ at or below the null median, i.e. no more similar than a random gene pair.

| group (28 pairs pooled) | n | conserved | break |
|---|---|---|---|
| orthologs | 2,182 | **16.0%** | **26.2%** |
| paralogs (all) | 13,380 | 4.5% | 47.4% |
| within-RefOG paralogs | 9,913 | 5.0% | 44.0% |
| deep (cross-RefOG) paralogs | 3,467 | 3.0% | 57.1% |

**Ortholog profiles by species pair.**
- Mammal–mammal: median ρ 0.20 to 0.51; conserved 6 to 36%; breaks 6 to 39%.
  - Closest pairs: mouse–rat median 0.51, conserved 36%, breaks 8%; human–chimp 0.42, 35%, 6%.
- Chicken pairs: median ρ 0.26 to 0.45; conserved 15 to 28%; breaks 17 to 33%.
- Zebrafish pairs: median ρ 0.09 to 0.30; conserved 4 to 12%; **breaks 29 to 45%**. Profile conservation therefore decays with divergence and with the teleost WGD (sub- and neofunctionalisation of ohnologs).

Every ortholog break (gene pair, RefOG, ρ, number of shared entities) is listed per pair in `per_pair.*.ortholog_breaks` in the JSON. The per-pair table is in §5.

Caveat: when only 8 to 25 entities are shared, ρ is noisy. Many "breaks" in pairs that involve chicken, rat, chimp or opossum may be sampling noise rather than biology. Human–mouse (210 entities) is the most reliable estimate: conserved 6%, break 17%, median ρ 0.26 against 0.11 for paralogs.

## 4. Recommendation
- **Kill S9 as a merge layer.** The hard-case gain CI includes 0 (primary) or is significantly negative (extended), and the gain over shuffled is inconsistent across seeds.
- Keep `ancestra/expression.py` only as an *annotation* of the HOG output, as asked by the brief: per-HOG expression conservation and breaks. It must not feed merge decisions.
- A better design would need many more shared entities, or within-species paralog-vs-paralog profile comparison. That is beyond the time box.

## 5. Per-pair table
| pair | shared entities | cand w/ expr | orth | hard pairs | AUC seq | AUC ρ | AUC ρ hard | AUC Jaccard | AUC shuf | null q95 | orth median ρ | orth conserved | orth break | par median ρ | par conserved | par break |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H. sapiens / M. musculus | 210 | 586 | 81 | 25 | 0.871 | 0.713 | 0.833 | 0.664 | 0.507 | 0.603 | 0.263 | 0.062 | 0.173 | 0.113 | 0.026 | 0.477 |
| H. sapiens / G. gallus | 13 | 413 | 60 | 48 | 0.855 | 0.719 | 0.714 | 0.608 | 0.472 | 0.662 | 0.398 | 0.150 | 0.167 | 0.048 | 0.065 | 0.462 |
| H. sapiens / D. rerio | 74 | 666 | 93 | 263 | 0.785 | 0.594 | 0.540 | 0.573 | 0.493 | 0.653 | 0.168 | 0.118 | 0.376 | 0.042 | 0.038 | 0.442 |
| H. sapiens / P. troglodytes | 25 | 676 | 80 | 14 | 0.910 | 0.813 | 0.546 | 0.692 | 0.466 | 0.592 | 0.425 | 0.350 | 0.063 | −0.010 | 0.039 | 0.552 |
| H. sapiens / R. norvegicus | 22 | 584 | 79 | 21 | 0.864 | 0.718 | 0.982 | 0.633 | 0.471 | 0.648 | 0.370 | 0.139 | 0.203 | 0.096 | 0.030 | 0.485 |
| H. sapiens / C. familiaris | 50 | 514 | 73 | 26 | 0.858 | 0.743 | 0.733 | 0.651 | 0.492 | 0.487 | 0.317 | 0.164 | 0.110 | 0.075 | 0.036 | 0.415 |
| H. sapiens / M. domestica | 17 | 483 | 72 | 36 | 0.826 | 0.672 | 0.657 | 0.639 | 0.477 | 0.613 | 0.259 | 0.069 | 0.264 | 0.018 | 0.051 | 0.506 |
| P. troglodytes / M. musculus | 23 | 732 | 87 | 112 | 0.870 | 0.638 | 0.454 | 0.651 | 0.484 | 0.604 | 0.244 | 0.149 | 0.345 | 0.025 | 0.031 | 0.516 |
| P. troglodytes / R. norvegicus | 11 | 765 | 104 | 107 | 0.844 | 0.591 | 0.376 | 0.574 | 0.484 | 0.672 | 0.227 | 0.164 | 0.394 | 0.050 | 0.047 | 0.496 |
| P. troglodytes / C. familiaris | 17 | 685 | 83 | 123 | 0.854 | 0.688 | 0.643 | 0.656 | 0.513 | 0.588 | 0.284 | 0.133 | 0.217 | 0.018 | 0.040 | 0.495 |
| P. troglodytes / M. domestica | 8 | 610 | 88 | 81 | 0.829 | 0.550 | 0.474 | 0.667 | 0.521 | 0.714 | 0.200 | 0.125 | 0.330 | 0.113 | 0.081 | 0.414 |
| P. troglodytes / G. gallus | 8 | 503 | 69 | 72 | 0.857 | 0.656 | 0.920 | 0.617 | 0.500 | 0.691 | 0.286 | 0.217 | 0.333 | −0.071 | 0.046 | 0.562 |
| P. troglodytes / D. rerio | 9 | 266 | 45 | 144 | 0.735 | 0.550 | 0.527 | 0.599 | 0.522 | 0.900 | 0.200 | 0.044 | 0.400 | 0.100 | 0.014 | 0.498 |
| M. musculus / R. norvegicus | 22 | 633 | 83 | 34 | 0.870 | 0.812 | 0.651 | 0.667 | 0.468 | 0.600 | 0.513 | 0.361 | 0.084 | 0.111 | 0.058 | 0.420 |
| M. musculus / C. familiaris | 55 | 553 | 74 | 29 | 0.861 | 0.745 | 0.550 | 0.701 | 0.516 | 0.536 | 0.292 | 0.203 | 0.095 | 0.062 | 0.031 | 0.422 |
| M. musculus / M. domestica | 22 | 526 | 71 | 82 | 0.817 | 0.611 | 0.600 | 0.746 | 0.502 | 0.713 | 0.289 | 0.141 | 0.296 | 0.132 | 0.024 | 0.402 |
| M. musculus / G. gallus | 15 | 444 | 71 | 56 | 0.837 | 0.670 | 0.757 | 0.628 | 0.520 | 0.664 | 0.261 | 0.169 | 0.282 | 0.054 | 0.035 | 0.469 |
| M. musculus / D. rerio | 402 | 713 | 95 | 324 | 0.771 | 0.635 | 0.580 | 0.648 | 0.512 | 0.625 | 0.215 | 0.095 | 0.326 | 0.008 | 0.036 | 0.505 |
| R. norvegicus / C. familiaris | 17 | 609 | 81 | 95 | 0.833 | 0.775 | 0.551 | 0.640 | 0.506 | 0.573 | 0.380 | 0.309 | 0.062 | 0.098 | 0.059 | 0.403 |
| R. norvegicus / M. domestica | 10 | 529 | 80 | 79 | 0.799 | 0.613 | 0.394 | 0.673 | 0.509 | 0.684 | 0.249 | 0.163 | 0.325 | 0.071 | 0.033 | 0.470 |
| R. norvegicus / G. gallus | 11 | 440 | 66 | 54 | 0.833 | 0.726 | 0.700 | 0.631 | 0.449 | 0.682 | 0.418 | 0.227 | 0.197 | 0.027 | 0.062 | 0.471 |
| R. norvegicus / D. rerio | 11 | 609 | 86 | 288 | 0.747 | 0.614 | 0.563 | 0.649 | 0.505 | 0.800 | 0.300 | 0.047 | 0.291 | 0.029 | 0.027 | 0.526 |
| C. familiaris / M. domestica | 13 | 498 | 75 | 116 | 0.780 | 0.630 | 0.571 | 0.675 | 0.506 | 0.612 | 0.245 | 0.173 | 0.200 | 0.099 | 0.062 | 0.407 |
| C. familiaris / G. gallus | 12 | 407 | 57 | 63 | 0.829 | 0.702 | 0.841 | 0.687 | 0.490 | 0.617 | 0.378 | 0.281 | 0.263 | 0.005 | 0.057 | 0.497 |
| C. familiaris / D. rerio | 25 | 673 | 88 | 312 | 0.780 | 0.546 | 0.579 | 0.645 | 0.509 | 0.624 | 0.089 | 0.057 | 0.398 | 0.033 | 0.055 | 0.489 |
| M. domestica / G. gallus | 11 | 350 | 57 | 53 | 0.839 | 0.699 | 0.778 | 0.706 | 0.495 | 0.714 | 0.446 | 0.228 | 0.175 | 0.055 | 0.072 | 0.468 |
| M. domestica / D. rerio | 12 | 544 | 96 | 209 | 0.750 | 0.561 | 0.533 | 0.601 | 0.501 | 0.771 | 0.093 | 0.094 | 0.448 | 0.061 | 0.058 | 0.489 |
| G. gallus / D. rerio | 11 | 551 | 88 | 256 | 0.807 | 0.605 | 0.618 | 0.623 | 0.507 | 0.691 | 0.155 | 0.102 | 0.364 | 0.036 | 0.050 | 0.488 |

The "AUC seq" column uses raw bits pooled within the pair. Raw bits are not comparable across families (see the S5 doc), so the query-level tables in §2 are the relevant comparison.

## 6. Files
- **Code:** `ancestra/expression.py`, `scripts/S9_expression_benchmark.py` (it imports S5 helpers) and `scripts/S7_S9_manifest.py`.
- **Tests:** `tests/test_expression_isoforms.py` (Bgee parsing and quality filter, shared entities, Spearman/Jaccard, `min_shared`, shuffle).
- **Work:** `work/S9/profiles.pkl` (cached matrices) and `work/S9/pair_expr_*.parquet` (per-pair scores incl. shuffled).
- **Rerun:** `python scripts/S9_expression_benchmark.py`. It needs `work/S7/hits.tsv` and takes about 4 min.
