# S6 — structure (3Di) evidence (DESIGN.md §3 row S6)

Status: in progress. Author: S6 sub-agent (structural bioinformatics). Started T+3h26 (scripts/clock.sh).

## 0. Pre-registration (written at T+3h36, before any S6 result was looked at)
**Hypothesis.** 3Di helps at deep/remote levels: fewer S0 fissions (remote homologs) and better deep-node
decisions; little effect at shallow nodes.

**Data.** `data/structure/orthobench_3di_best.fa`: exact-sequence AFDB-v6 3Di, 202,475 / 251,378 OrthoBench
proteins (80.5 %). DEV RefOG genes: 948, of which 66.7 % have 3Di (lower than overall: chicken, dog, zebrafish,
rat and non-canonical isoforms are the gaps). Proteins without 3Di fall back to sequence everywhere.

**Negative control.** Shuffled 3Di, seeds 11/23/47 (`structure.shuffled_3di`): proteins sorted by length, cut into
bins of 20, 3Di strings permuted within each bin as a derangement (no protein keeps its own string), donor string
resampled to the recipient length (nearest index). The AA sequence is unchanged, so the shuffled arm keeps
Foldseek's AA-substitution component and isolates the 3Di information. Per-residue identity of shuffled to real
3Di: 0.161 (composition background), vs 1.0.

**Exp 1 — S0 rescue.** Units: the selected S0 partition (`work/S0/2sA_dc30_cc_sh30_RA/families.tsv`).
One structure representative per S0 stage-1 cluster (DIAMOND rep if it has 3Di, else longest member with 3Di;
38,802 of 42,732 clusters), searched all-vs-all with `foldseek search --alignment-type 2` (3Di+AA, -s 9.5,
e ≤ 1e-2 stored, max-seqs 1000). Families joined by union-find over hits in increasing e-value.
Grid: e ≤ {1e-3, 1e-5} × small_max {10 (one side ≤ 10 genes), any}; aligned span ≥ 0.3 of the shorter protein
(S0's rule); joins that would exceed 7,500 genes skipped (S0 gmax). Metrics on DEV RefOGs only:
completeness, mean fission, best frac, fusion proxies (fusion rate, foreign-RefOG fraction), Σg²/N², max family.
Selection: max DEV completeness, ties → lower Σg². Keep: completeness(real) − mean(shuffled) > 0 and
completeness(real) − off > seed spread (max − min over the 3 shuffles).

**Exp 2 — structure-aware scores for deep/ambiguous decisions.** Pair set: every in-family pair of the
sequence ScoreTable (`work/core_v0_ob/run_v2`, hybrid DIAMOND+MMseqs2) in the families holding DEV RefOG genes
(evaluation scope; other families cannot change DEV metrics since the merge is family-local), both genes with
3Di, plus self pairs, aligned with `foldseek structurealign --alignment-type 2` on an explicit `tsv2db` pair list.
Score: S_3di = bits/max(self_a, self_b) (same normalisation as the merge's norm_max). Calibration to the
sequence scale: least-squares line d_seq ≈ a + b·d_3di (d = −ln S) on pairs with both; combined
`avg`: S = sqrt(S_seq·S_3di'), `max`: max(S_seq, S_3di'). Written into a ScoreTable (bits' = S·max(self_seq)),
so the unchanged merge (default config = hybrid local-tree mode; its NJ distance is −ln S) consumes it — **no core
change needed**. Arms: scope deep (pairs whose species LCA is Euteleostomi, Chordata, Ecdysozoa or Bilateria: the
tree-mode nodes and the root) or all; mode avg or max; × real / shuf11 / shuf23 / shuf47; off = sequence only.
Metrics: DEV OG-level F (scripts/orthobench_dev_eval.py) and pair P/R/F1 (scripts/orthobench_dev_pairs_eval.py),
pair F1 split deep (vert–teleost, vert–invert, invert–invert) vs shallow (mammal, tetrapod, teleost–teleost).
Primary endpoint: DEV pair F1 (OG-F on 34 RefOGs is too noisy, S5). Keep rule as above with pair F1;
paired bootstrap CI over RefOGs reported.

**Exp 3 — uncovered proteins.** Pair metrics of each Exp-2 run split by how many of the two genes have 3Di
(2 / 1 / 0); effect = Δ vs off per subset.
