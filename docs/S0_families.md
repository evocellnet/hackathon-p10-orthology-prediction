# S0 — family partition (DESIGN.md §1, §3 row S0)

Status: done (dev tuning complete, test scored once). Author: S0 sub-agent, 2026-10-26 15:15–16:32 CET (T+0h35 to T+1h52).
Code: `ancestra/families.py` (module + CLI), `scripts/run_S0_grid.py` (grid runner + dev scoring + one-shot test),
`scripts/summarize_S0.py` (table + pre-stated rule), `scripts/S0_dev_ceiling.py` (upper bound),
`scripts/eval_S0_toy.py` (toy simulation), `scripts/make_orthobench_split.py` (split). Tests: `tests/test_families.py` (14 tests).
Results: `results/S0_families/*.json` (each with a `manifest`: config, config hash, commit, input sha256s, split sha256,
families.tsv sha256, command, time, node load), `summary.tsv`, `selection.json`, `dev_ceiling.json`, `toy/`.
Large intermediates (families.tsv, hit tables) are in `work/S0/<config>/` (gitignored, reproducible).

## 1. Split (DESIGN §4), v2 — homology-aware
v1 (random 35/35 over RefOG ids, seed 11) was replaced on PI instruction (review r2_stats_ai_failure #1: paralogous
RefOGs leak across dev/test). v2 procedure (`scripts/make_orthobench_split.py`, docstring is the spec):
1. All 1,945 RefOG sequences (`Supporting_Data/Data_for_RefOGs/refog_sequences`) all-vs-all with `diamond blastp
   --sensitive -k 0 -e 1e-5`. Two RefOGs are linked by any hit with e < 1e-5 between their genes, or by a shared gene
   (FBpp0309618 is listed in both RefOG021 and RefOG068; the first implementation missed this and found 61
   components; fixed before any use). Result: **60 components**, the same five multi-RefOG components as the reviewer:
   {005,006}, {015,023}, {016,026,027,032,037,039,058}, {021,068}, {038,052}.
2. Components sorted by mean RefOG size, cut into strata of 2; strata visited in `random.Random(11)` order; within a
   stratum a random unit goes to the smaller subset (by #RefOGs; coin flip at equality), the other to the other subset.
3. **dev = 34 RefOGs / 960 genes, test = 36 RefOGs / 985 genes.** Files: `splits/orthobench_split.json`,
   `splits/orthobench_refog_components.tsv`, hashes in `splits/SHA256SUMS`.
All dev numbers below are on the v2 dev set (the five configs scored before the switch were re-scored from their saved
families.tsv with `run_S0_grid.py --rescore`).

## 2. Metrics (computed only on dev RefOGs during tuning)
Per RefOG R (low-certainty genes removed, as in benchmark.py; "keep uncertain" variants are stored too):
- **completeness** (primary): all of R in one family. **fission**: number of families R is spread over.
- **best frac**: largest share of R in one family (also gene-weighted).
- **fusion proxies**: size ratio = size of R's best family / |R|; fusion rate = share of RefOGs whose best family also
  holds genes of another *dev* RefOG (test RefOGs are never looked at during tuning, so this proxy is weak).
- **pair R/P**: OrthoBench-style pair recall/precision of the partition (q_even weighting). Precision is ~0.01 by
  design: a family is a superfamily, not an orthogroup.
- Global: Σ g_f², Σ g_f²/N², top-10 sizes, singletons; wall time and peak RSS from `/usr/bin/time -v` (24 threads).
**Pre-stated selection rule** (written into `configs.json` before the first grid run): among configs with
Σg² ≤ 0.05·N² and max family ≤ 0.03·N (N = 251,378, so ≤ 7,541 genes), maximise dev completeness; ties → higher dev
pair recall → lower wall time.

## 3. Backends and configurations
- Stage 1 (near-linear): DIAMOND 2.2.8 `cluster` (`--approx-id 30 --member-cover 50`, 5-round cascade up to
  more-sensitive) or `deepclust` (`--approx-id 0`, member-cover 50 or 20), or MMseqs2 `cluster` (-s 7.5 / 5.7).
- Stage 2 (two-level, reps ≪ N): DIAMOND `--sensitive -e 1e-5 -k 500` search of either the stage-1 representatives
  against themselves (**2sR**, R×R) or **every gene against the representatives (2sA, N×R)**. Gene-level hits are
  lifted to cluster–cluster edges; optional coverage guard (`min`: both sequences ≥ c; `shorter`: alignment covers ≥ c of
  the shorter sequence); connected components (`cc`) or MCL; components above `s2_gmax` genes are split by Leiden
  (RBConfiguration, coarsest resolution in a fixed ladder that splits, recursive; seed 11) and flagged
  (`families_split_flags.tsv`); optional **RA** = genes of split parts move to the part holding their best
  filter-passing representative hit.
- Name code: `2sA_dc30_cc_sh30_RA` = stage 1 DIAMOND cluster 30% / 50% cover, N×R search, CC, shorter-coverage ≥ 0.3,
  Leiden cap 7,500, reassignment. Configs reusing a stage-1 result (`s1_reuse`) report wall time including the reused
  stage 1 and peak RSS = max of both runs.
- Grid: 23 configs + 1 fresh rerun, larger than the ~12 initially planned because the PI's round-2 instruction
  (two-level scheme, ≥95 % target, Leiden cap) arrived mid-grid. All are listed below; none were dropped.

## 4. Results — OrthoBench 12 proteomes (N = 251,378), DEV RefOGs (n = 34)
Node shared with other agents (load average 27–81 on 24 CPUs at run start, logged per run): times are upper bounds.

| config | dev complete | fission (mean #fam) | best frac | fusion rate | size ratio med/max | pair R/P | Σg²/N² | max fam (%N) | singletons | wall s (RSS GB) | eligible |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2sA_dc30_cc_nocov_RA | 0.794 (27/34) | 1.32 | 0.983 | 0.412 | 9.68/676.1 | 0.970/0.007 | 0.0084 | 7245 (2.9%) | 14547 | 422 (6.0) | yes |
| 2sA_dc30_cc_sh30_RA **(selected)** | 0.794 (27/34) | 1.32 | 0.988 | 0.294 | 10.01/250.4 | 0.973/0.010 | 0.0075 | 7478 (3.0%) | 15344 | 421 (6.0) | yes |
| 2sA_dc30_cc_sh30_RA_fresh | 0.794 (27/34) | 1.32 | 0.988 | 0.294 | 10.01/250.4 | 0.973/0.010 | 0.0075 | 7478 (3.0%) | 15344 | 263 (6.2) | yes |
| 2sA_dc30_cc_nocov | 0.765 (26/34) | 1.32 | 0.983 | 0.412 | 9.65/687.9 | 0.970/0.007 | 0.0084 | 7094 (2.8%) | 14547 | 514 (6.0) | yes |
| 2sA_dc30_cc_nocov_G20k | 0.765 (26/34) | 1.32 | 0.983 | 0.559 | 147.97/1809.6 | 0.970/0.002 | 0.0256 | 18605 (7.4%) | 14547 | 408 (6.0) | no |
| 2sA_dc30_cc_sh30 | 0.765 (26/34) | 1.32 | 0.987 | 0.294 | 10.04/250.4 | 0.973/0.010 | 0.0075 | 7462 (3.0%) | 15344 | 459 (6.0) | yes |
| 2sA_dc30_mcl_I1.2 | 0.765 (26/34) | 1.32 | 0.981 | 0.235 | 12.03/1047.8 | 0.968/0.006 | 0.0082 | 10789 (4.3%) | 14547 | 425 (6.0) | no |
| 2sA_dc30_mcl_I1.4 | 0.765 (26/34) | 1.41 | 0.970 | 0.206 | 3.30/112.2 | 0.958/0.030 | 0.0038 | 9687 (3.9%) | 14547 | 440 (6.0) | no |
| 2sA_deep20_cc_nocov | 0.765 (26/34) | 1.35 | 0.976 | 0.353 | 5.40/555.5 | 0.963/0.009 | 0.0088 | 7375 (2.9%) | 14024 | 303 (6.4) | yes |
| 2sA_deep20_cc_nocov_RA | 0.765 (26/34) | 1.35 | 0.977 | 0.324 | 5.50/546.2 | 0.963/0.009 | 0.0087 | 7101 (2.8%) | 14024 | 313 (6.4) | yes |
| 2sA_deep20_cc_sh30_RA | 0.735 (25/34) | 1.41 | 0.968 | 0.324 | 3.45/433.8 | 0.957/0.010 | 0.0068 | 6586 (2.6%) | 14436 | 267 (6.4) | yes |
| 2sA_deep20_cc_sh50 | 0.735 (25/34) | 1.47 | 0.962 | 0.412 | 2.86/404.4 | 0.945/0.010 | 0.0053 | 6017 (2.4%) | 14885 | 330 (6.4) | yes |
| 2sA_mm75_cc_sh30_RA | 0.735 (25/34) | 1.41 | 0.953 | 0.353 | 13.73/695.9 | 0.940/0.006 | 0.0078 | 7394 (2.9%) | 13655 | 260 (2.8) | yes |
| 2sA_dc30_cc_sh50 | 0.706 (24/34) | 1.47 | 0.975 | 0.412 | 6.45/429.4 | 0.953/0.008 | 0.0064 | 7165 (2.9%) | 16962 | 425 (6.0) | yes |
| 2sR_dc30_cc_nocov | 0.706 (24/34) | 1.56 | 0.960 | 0.412 | 6.11/811.1 | 0.889/0.007 | 0.0065 | 7300 (2.9%) | 16091 | 436 (6.0) | yes |
| 2s_dc30_mcl_I1.4 | 0.706 (24/34) | 1.68 | 0.948 | 0.206 | 2.23/97.2 | 0.878/0.037 | 0.0029 | 8950 (3.6%) | 16091 | 153 (5.5) | no |
| deep_c20 | 0.588 (20/34) | 2.26 | 0.901 | 0.088 | 1.55/41.2 | 0.802/0.065 | 0.0012 | 5456 (2.2%) | 17309 | 224 (6.4) | yes |
| mm_s7.5_c50 | 0.588 (20/34) | 2.09 | 0.901 | 0.088 | 1.50/44.8 | 0.803/0.116 | 0.0010 | 4706 (1.9%) | 16486 | 151 (2.8) | yes |
| mm_s5.7_c30_cc | 0.529 (18/34) | 2.15 | 0.923 | 0.235 | 1.59/4485.3 | 0.853/0.001 | 0.0262 | 40368 (16.1%) | 20025 | 95 (1.5) | no |
| 2s_dc30_cc_min30 | 0.471 (16/34) | 2.68 | 0.873 | 0.147 | 1.00/60.0 | 0.693/0.083 | 0.0035 | 9935 (4.0%) | 19566 | 193 (5.6) | no |
| 2s_dc30_mcl_min30_I1.4 | 0.471 (16/34) | 2.74 | 0.864 | 0.147 | 1.00/60.0 | 0.687/0.104 | 0.0011 | 5446 (2.2%) | 19566 | 218 (6.3) | yes |
| 2s_dc30_cc_min50 | 0.441 (15/34) | 3.00 | 0.844 | 0.000 | 1.00/19.4 | 0.655/0.284 | 0.0014 | 8119 (3.2%) | 21942 | 190 (6.2) | no |
| dc_a30_c50 | 0.441 (15/34) | 3.65 | 0.822 | 0.000 | 1.00/4.5 | 0.634/0.492 | 0.0006 | 4581 (1.8%) | 26032 | 335 (6.0) | yes |
| deep_c50 | 0.412 (14/34) | 2.91 | 0.850 | 0.118 | 1.00/17.9 | 0.684/0.241 | 0.0009 | 5404 (2.1%) | 20867 | 350 (6.2) | yes |


Trade-off (dev): single-level clustering (DIAMOND cluster/deepclust, MMseqs2) tops out at 0.59 completeness with
Σg²/N² ≤ 0.0012. A reps-only second stage (2sR, 2s_*) adds little (0.44–0.71) because centroids of one family are often
too far apart for each other while their members are not. Searching **all genes against the representatives** (2sA)
is the key step (0.71–0.79), at Σg²/N² ≈ 0.004–0.009 (≈ 120–250× fewer pairs than all-vs-all). Coverage guards
(`min` 0.3/0.5) keep families small but lose completeness; `shorter` 0.3 and no guard reach the same completeness.
The 2sA graph has one giant component (102–136k genes, 40–54 % of N, the multi-domain "hairball"); without the
Leiden cap single linkage is unusable (max 136k; PI's 20k cap → max family 18.6k = 7.4 % N, ineligible). MCL on the
same graph keeps max family 9.7–10.8k (ineligible) for the same completeness.

### Upper bound (`results/S0_families/dev_ceiling.json`)
- Single linkage over all 2sA hits (e ≤ 1e-5, no guard, no cap) gives **27/34 = 0.794**; this is the selected
  config's value, so Leiden + reassignment lose nothing on dev.
- Each remaining stray gene was searched against all 251k proteins with `--ultra-sensitive -e 1 -k 2000`: in 6
  RefOGs (002, 021, 054, 055, 057, 066) some member has **no hit with e ≤ 1e-3 to any member of its RefOG's main
  family** (e.g. a 16 kaa fly protein in RefOG002, a 92-aa Ciona fragment in RefOG066 with best e = 0.01, two fly
  paralogs in RefOG054 that hit only each other). So **≤ 28/34 = 0.82 is reachable by any pairwise-hit partition** even
  with ultra-sensitive search at e ≤ 1e-3, and on the toy simulation that setting already creates false fusions between
  unrelated random families (pair precision 0.954 / 0.942). **The PI's pre-registered ≥ 95 % dev completeness is not
  attainable by sequence-similarity S0 on OrthoBench** (the curated RefOGs contain members found by HMM profiles and
  trees). Gene-level the selected config keeps 98.8 % of RefOG genes in the RefOG's family (mean best frac), and 88 % of
  dev RefOGs have at most one stray gene.

### Toy simulation (`ancestra.simulate`, 16 species, 500 root families, seed 11; `results/S0_families/toy/`)
Truth = root families (independent random sequences, so every fusion is false).
| config | height 1.0: complete / pair R | height 2.0: complete / pair R | false fusions |
|---|---|---|---|
| all-vs-all CC, --sensitive e≤1e-5 (ceiling) | 0.876 / 0.953 | 0.600 / 0.778 | 0 |
| all-vs-all CC, --ultra-sensitive e≤1e-3 | 0.922 / 0.972 | 0.692 / 0.834 | yes (P = 0.954 / 0.942) |
| dc_a30_c50 (stage 1 only) | 0.633 / 0.803 | 0.322 / 0.523 | 0 |
| deep_c20 | 0.781 / 0.900 | 0.469 / 0.669 | 0 |
| 2s_dc30_cc_min50 (current pipeline default) | 0.783 / 0.903 | 0.478 / 0.677 | 0 |
| **2sA_dc30_cc_sh30_RA (selected)** | **0.871 / 0.950** | **0.590 / 0.766** | 0 |
The N×R two-level scheme recovers 99 % (h1.0) and 98 % (h2.0) of what a full all-vs-all single linkage at the same
sensitivity achieves. Residual fission on the toy is mostly saturated fast-evolving families (per-family rate
multiplier Gamma(2)).

## 5. Recommendation
**`2sA_dc30_cc_sh30_RA`** (rule-selected; tie at 27/34 with `2sA_dc30_cc_nocov_RA`, broken by pair recall 0.973 vs 0.970):
```json
{"backend": "diamond", "dmnd_cmd": "cluster", "approx_id": 30, "member_cover": 50,
 "stage2": "cc", "s2_query": "all", "s2_sens": "--sensitive", "s2_evalue": 1e-5,
 "s2_cov": 0.3, "s2_cov_mode": "shorter", "s2_gmax": 7500, "s2_reassign": true}
```
- Dev: completeness **0.794 (27/34)** = the single-linkage ceiling; mean fission 1.32; mean best frac 0.988; pair
  recall 0.973; 22,631 families; **15,344 singletons (6.1 % of genes)**; Σg² = 4.77e8 (Σg²/N² = 0.0075); top-10 sizes
  7478, 7380, 7217, 5672, 5537, 5246, 4831, 4688, 4566, 4252 (max = 2.97 % N, just under the 3 % bound); 145 families
  (119,638 genes) come from the Leiden split of one giant component and are flagged for root rescue; 706 genes reassigned.
- Cost (fresh end-to-end rerun, `2sA_dc30_cc_sh30_RA_fresh`): **263 s wall, 6.2 GB peak RSS, 0.55 CPU-h**
  (stage 1 148 s, stage 2 113 s, 42,732 representatives). The rerun's families.tsv is byte-identical
  (sha256 9b63fa5d…), and so is a third run with the final `ancestra/families.py` (sha256 45918626…).
- **Test (scored once, logged in `results/test_log.tsv`): completeness 0.722 (26/36)**; mean fission 1.42; mean best
  frac 0.969; pair recall 0.945; fusion rate 0.25. Incomplete test RefOGs: 011, 014, 039, 047, 053, 058, 060, 063,
  064, 069 (eight of them keep ≥ 90 % of their genes in one family).
- For larger inputs, `s2_gmax` accepts a fraction of N (`0.03`) — untested beyond OrthoBench. The Leiden-split
  giant component is the main precision risk downstream: its parts are flagged so the merge stage can treat
  cross-part matches (root rescue).
- **Action for the pipeline owner**: `ancestra/pipeline.py` `DEFAULT_FAMILY_CONFIG` is currently
  `2s_dc30_cc_min50` (dev 0.441). It should become the config above (I did not edit pipeline.py).

## 6. Caveats
- Timings are on a shared node and partly with reused stage 1; the fresh rerun is the reference number.
- The fusion proxy only sees dev RefOGs; the size ratio (median 10, max 250 on dev) shows families are superfamilies,
  which is intended (S0 should over-merge rather than split).
- The ≥ 95 % target cannot be met by a sequence-similarity partition (§4 upper bound). Structure-based (3Di) or
  profile (HMM) linking of singletons and small families is the only route to the remaining RefOGs, and should be
  checked against the toy negative control (false fusions).

## Determinism fix (core agent, 2026-10-26 ~T+5h; details in docs/core_v0.md §7.6)
- S0 results depended on input order and ID names, because DIAMOND `cluster` greedy ties follow input order (ARI 0.78 between the agent run on `ob_all.fa` and the pipeline's integer-ID FASTA).
- They also depended on `--threads` in the stage-2 `diamond blastp`, through DIAMOND's automatic block and chunk sizing: 0.016% of hits changed between 12 and 24 threads.
- `run_families` now clusters a canonical copy of the input (normalised sequences, sorted by length desc / sequence / ID, renamed by rank). The stage-2 search pins `-b 4 -c 1` (`s2_block`, `s2_chunks`), and edges are sorted before CC and Leiden.
- On OrthoBench, `families.tsv` is byte-identical for 12 vs 24 threads and for agent vs pipeline input order.
- **The selected config's partition changed:** 22,643 families, max 7,269. Families TSVs produced before the fix are not reproducible and should not be reused (`s1_reuse` refuses pre-fix stage-1 outputs).
- Test: `tests/test_s0_determinism.py`.
