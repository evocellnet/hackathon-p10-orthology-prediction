# L2-genome: synteny fission repair and recall at the LCA, intron positions, expression

Agent `genome` (L2 round, post-freeze, dev only). Branch `l2-genome`, worktree `work/l2_genome_wt`.
Protocol: docs/L2_protocol.md. Nothing here changes the frozen method (`freeze-v1`) or any reported number.

## 1. Pre-registration (written and committed 2026-10-27 ~18:15 CET, before any dev integration result)

### 1.0 Context that shaped the design (read before writing this section)
- S5: synteny separates orthologs from paralogs (dev AUC 0.92 human/mouse), but its constrained R1 re-ranking
  did not replicate on test, most likely because local FastTree trees already fix those cases and synteny only
  re-ranks partners inside synteny-off root HOGs.
- core_v0 §7.3 item 8 (prior art in this repo): `locus_repair` (a still-alone single-species unit joins the
  group lacking its child with the highest *unmasked* synteny support >= 2, no sequence gate, applied before
  the local-tree resolver) gave pair F1 +0.016 and S8 fission candidates -45% on the old NJ base, but OG-F
  79.2 -> 75.5 (OG precision -4). It was never run on the frozen `accurate` base.
- S8: 808 orthology fission errors among 2,523 apparent mammal losses; on the accurate dev base the S8-style
  diagnostic ("5/6-mammal Theria HOG whose missing species has a species-only Theria HOG in the same root
  HOG") counts 706 candidates (results/core_v0_orthobench/runs.tsv, row `T[det] tree=on:ft-medoid off:nj`).
- RESEARCH_TREE 8b3/8b4: 78% of missed TreeFam-A dev pairs are merge failures at the pair's LCA; a
  sequence-only "fission repair 0.5" (R3 ratio 0.9 -> 0.5) cost OrthoBench dev OG-F -1.55.

The new arms below differ from `locus_repair` in four pre-stated ways that target its failure (OG-precision
loss through fusions and joins to the wrong group): (i) a **sequence gate** (the partner must be within
rho of the unit's best cross-child score), (ii) a **unique-locus margin** (the chosen group's synteny support
must exceed every other candidate group's by >= 1), (iii) **masked support** by default (joins only inside
the synteny-off root HOG, so no root-HOG fusion is possible), and (iv) it runs **after** the local-tree
resolver, so it only acts on what the frozen preset's trees left alone (the presumed reason S5 did not
replicate).

### 1.1 Rule shared by ideas 1 and 2: gated synteny repair (GSR), `MergeConfig.gsr_*` (default off)
At each internal node X, after R1-R3 and after the local-tree resolver:
1. Eligible units: still-alone units u (group size 1). Idea 1 (F arms): units whose genes all come from one
   species. Idea 2 (R arms): any still-alone unit.
2. Candidate edges (u, v): cross-child, in-family scored pairs (the merge's candidate set), where group
   g = label(v) contains no unit of child(u) (the join is a speciation, never a new duplication).
3. Sequence gate: score(u, v) >= rho * best_any(u), rho = **0.75** (the existing S5 `ctx_window`; not tuned).
4. Synteny gate: support s(u, g) = max over v in g of the S5 context support (anchors from `SyntenyContext`,
   w=10, core anchors; same object as the frozen default). s(u, g) >= **2** (existing `ctx_veto_min` /
   `locus_min_support` convention) and s(u, g) >= s(u, g') + **1** for every other candidate group g' of u.
5. Masking: support is masked to the synteny-off root HOG (`ConstrainedContext.support`) in the `_c` arms;
   unmasked (`locus_support`) in the `_u` arms.
6. Greedy order: support desc, then score desc; each group gains at most one unit per child in one pass.

Shuffled control of the layer: the GSR support (only) comes from a `SyntenyContext` with `shuffle_seed`
11/23/47; the frozen R1 synteny bonus and veto stay real, so the control isolates the new rule.

### 1.2 Idea 1: synteny-driven fission repair (F arms), OrthoBench dev
- Base (`off`) = frozen `accurate` preset with gene order, exactly as the dev reference row
  `T[det] tree=on:ft-medoid off:nj syn constrained+veto+r3` (run `work/core_v0_ob/run_det`, dev OG-F 79.69,
  macro pair F1 0.8026); re-run here to confirm byte-identity of the flags-off path.
- Arms: `F_c` (single-species, masked), `F_u` (single-species, unmasked); each + 3 shuffled seeds.
  Information-only (not eligible for keep): existing `locus_repair` (min support 2) on the same base.
- (a) Standalone signal: for every eligible (u, best-supported g) candidate logged by the rule, label the join
  from the dev RefOG pair labels (work/S5/refog_pair_labels.parquet; join correct if the labelled u x g pairs
  are mostly orthologs, wrong if mostly paralogs). Report AUC of support (real vs shuffled seeds), bootstrap CI
  over RefOGs, coverage (share of candidates with a label); and the S8 fission-candidate count (positive
  control: real should lower it, shuffled should not).
- (b) **Primary dev metric: OrthoBench dev macro pair F1** (34 RefOGs; paired bootstrap over RefOGs, B=2000,
  seed 11). Guards: OrthoBench dev micro pair precision (not worse by > 0.01), OG-F (not worse by > 1.0),
  TreeFam-A dev PPV on the `fast` QfO base (not worse by > 0.01), OrthoXML validity (`verify_orthoxml`).
  Secondary: S8 fission candidates, TreeFam-A dev TPR.

### 1.3 Idea 2: synteny recall anchors at the LCA (R arms), TreeFam-A dev
- Base = QfO `fast` preset merge on `work/qfo_recall/base` (`scripts/qfo_merge_iter.py` path; about 6 min),
  i.e. the frozen fast config with gene order. Confirmation of a passing arm on `accurate` (about 1 h per run)
  only if time allows, reported separately.
- Arms: `R_c` (any still-alone unit, masked), `R_u` (any still-alone unit, unmasked); each + 3 shuffled seeds.
  The F arms are also run on this base (guard for idea 1).
- (a) Standalone: synteny support coverage among missed TreeFam-A dev speciation pairs (share with a QfO gene
  order and support >= 2), real vs shuffled.
- (b) **Primary dev metric: TreeFam-A dev TPR** (359 dev families; paired bootstrap over families, B=2000,
  seed 11, official pseudocount formula). Guards: TreeFam-A dev PPV (not worse by > 0.01); OrthoBench dev
  (accurate base) OG-F (not worse by > 1.0) and macro pair F1 (not worse; CI reported); validity.
- Expectation stated in advance: deep pairs dominate the TreeFam-A misses and gene order is near-random across
  them, so a small effect is the prior.

### 1.4 Idea 3: gene structure (intron positions and phases)
- Data: Ensembl 99 GFF3 CDS segments for all 12 OrthoBench species (data/orthobench/ensembl_99/gff3). Protein
  -> transcript by `protein_id`; C. elegans (gene IDs) and any miss by gene -> transcript whose CDS length
  matches 3 x protein length (+3 with stop). QC: CDS length vs protein length.
- Intron position: each CDS-segment junction -> (residue index, phase 0/1/2) in protein coordinates.
- Pair score: global protein alignment (Biopython PairwiseAligner, BLOSUM62, gap open -10, extend -0.5);
  an intron of a is shared with b if an intron of b with the same phase sits on the aligned column (exact
  column; phase-1/2 introns split the same codon). Score = Dice = 2 x shared / (n_a + n_b) over introns that
  fall inside the aligned region of both proteins; undefined (coverage loss) when both have 0 such introns.
- (a) Standalone: all labelled dev pairs (26,204: within-RefOG + deep paralogs). ROC-AUC of Dice for
  ortholog vs paralog with bootstrap CI over RefOGs (B=1000, seed 11); hard queries as in S5 (best paralog
  sequence score >= 0.9 x best ortholog; sequence = canonical DIAMOND bits, work/S7/hits.tsv) with the S5
  decision rule (top-1 among candidates >= 0.9 x best seq; abstain on ties) and its gain over sequence.
  Shuffled control: intron positions of every gene redrawn uniformly over its CDS (same count, random phase),
  seeds 11/23/47. **Stop rule:** if the hard-case AUC (or the hard-query re-rank gain) is not above the
  shuffled max with a CI excluding it, stop and record a negative result.
- (b) Only if (a) passes: R1 re-ranking bonus inside the frozen constrained window (multiplicative, lambda 0.3,
  window 0.75, masked to synteny-off root HOGs, added to the synteny bonus), real vs 3 shuffled seeds;
  primary OrthoBench dev macro pair F1 with the same guards as 1.2.

### 1.5 Idea 4: expression rethink (only if cheap; otherwise "not pursued")
If pursued: tissue-specificity (tau) agreement as a 1:many tie-breaker, standalone only, on the S9 primary
pairs with the S9 decision rule; stop unless decision accuracy on hard queries beats the S9 40-47% and the
shuffled control.

### 1.6 Keep rule (protocol §6, applied per idea to its primary metric)
Delta(real - off) paired-bootstrap 95% CI lower bound > 0 (B=2000, seed 11); real beats the mean shuffled arm
by more than the shuffled seed spread (max - min); no guard worse than its limit. With two masked/unmasked arms
per idea, an arm is kept only if it passes alone; if both pass, the masked arm (fewer fusions) is the frozen
candidate. Every arm run is reported.

### 1.7 Addendum (committed ~19:20 CET, before running idea 4; ideas 1-3 results were known, idea 4 had not been run)
Idea 4 details fixed now: tau (Yanai et al. 2005 tissue-specificity index) per gene from the Bgee v15.2
`expr_simple` gold/silver calls of its own species (expression score of every called entity; >= 5 entities),
so no shared-entity requirement. Candidates: all labelled DEV pairs between the 10 Bgee species of OrthoBench
(work/S5/refog_pair_labels.parquet; sequence = S7 canonical DIAMOND bits, as S9's extended set). Hard query:
best paralog seq >= 0.9 x best ortholog seq. Rule: among candidates with seq >= 0.9 x best, pick the smallest
|tau_q - tau_c|; decide only if the runner-up's |delta tau| exceeds it by >= 0.05, else abstain. Shuffled: tau
permuted among genes within species (seeds 11/23/47). Stop unless the precision when deciding has a 95% CI
(bootstrap over RefOGs, B=1000, seed 11) whose lower bound is above 0.47 (S9's best) and above every shuffled
seed's precision; report the re-rank gain over sequence on hard queries as well.

## 2. Summary and verdict (dev only; 2026-10-27, 18:05-19:40 CET)
- **No arm passes the keep rule. Nothing is proposed for the test agent; no frozen config is shipped.**
- **Ideas 1 and 2: gated synteny repair (GSR).** The rule clearly works as a mechanism, but the dev benchmarks cannot see its effect.
  - Positive control: it removes 29-40% of the S8-style fission candidates on the frozen `accurate` base (706 -> 426-499), while the shuffled-order control removes none (706-709).
  - Real joins happen at vertebrate and grass nodes (Tetrapoda, Euteleostomi, Theria, Poaceae), where gene order is conserved. Shuffled joins are few (0-65) and fall at deep nodes (Bilateria, Eukaryota).
  - OrthoBench dev macro pair F1 changes by -0.0014 to +0.0059. Every CI includes 0. The best arm is `R_u` at +0.0059 [-0.0015, +0.0159], with OG-F +3.2.
  - TreeFam-A dev TPR changes by -0.0030 to +0.0042. Every CI includes 0, and one shuffled seed alone moves TPR by +0.0097.
  - Why the benchmarks are blind: fewer than 1% of GSR decisions touch a labelled OrthoBench dev pair, and TreeFam-A dev pairs are dominated by deep LCAs.
- **Idea 3: intron positions.** Stopped at the standalone step.
  - There is a real but weak signal: hard-pair AUC 0.578 [0.514, 0.661] against 0.497-0.499 shuffled, and below sequence (0.654) on the same pairs.
  - The S5-style re-rank on hard queries *loses* accuracy: -0.071 [-0.149, +0.007]. Young paralogs (teleost WGD copies, in-paralogs) keep intron positions as well as orthologs do.
- **Idea 4: expression tissue-specificity (tau) tie-breaker.** Stopped at the standalone step.
  - When it decides on a hard query, it is right 0.448 [0.380, 0.554] of the time, inside S9's 40-47% band. The shuffled seeds score 0.36-0.39.
  - The re-rank gain is -0.076 [-0.131, -0.010].
- **Cost.** About 3.5 CPU-h of merges for 34 merge runs. The intron alignments took 0.22 CPU-h; tau and labelling under 0.05 CPU-h. No GPU.
  - GSR itself adds about 10-25 s per QfO `fast` merge (Python, after each node's rules).

## 3. Data
No new downloads. Everything already in the project (DATA.md) was reused:
- Ensembl 99 GFF3 (DATA.md §, `data/orthobench/ensembl_99/gff3`, md5 in `md5sums.txt` there; the chr_patch_hapl_scaff files for human, mouse and zebrafish);
- OrthoBench input proteomes;
- gene-order tables `data/genomes/orthobench` and `data/inputs/qfo2022/gene_order` (S5);
- Bgee v15.2 `expr_simple` (`data/extra/bgee_v15_2`, md5 in `md5sums.txt`; CC0);
- cached run stages `work/core_v0_ob/run_det` (OrthoBench, deterministic S0) and `work/qfo_recall/base` (QfO 2022);
- the S5 dev pair labels `work/S5/refog_pair_labels.parquet`, the S7 canonical DIAMOND hits `work/S7/hits.tsv`, and the TreeFam-A dev reference pairs `work/qfo_recall/treefamA_dev_refpairs.json` (dev families only).

No test RefOG, TreeFam-A test family, SwissTree or other QfO test was read.

## 4. Methods
### 4.1 GSR (`ancestra/merge.py`: `MergeConfig.gsr`, `gsr_mask`, `gsr_rho`, `gsr_min`, `gsr_margin`; `_gsr_repair`)
- The rule is the one in §1.1. It runs per node after R1-R3, after the local-tree resolver and after deep species-overlap absorption, and before R4.
- It is off by default. The new fields are left out of `MergeConfig.key()` at their defaults, so stage hashes of the frozen presets do not change.
- The `off` arm reproduces the frozen dev reference exactly:
  - OrthoBench: OG-F 79.69, macro pair F1 0.8026, S8 candidates 706;
  - QfO `fast`: TreeFam-A dev TPR 0.6151, PPV 0.9224.
- Margin reading, fixed before any result: the runner-up support is taken over *all* groups that u has a candidate edge to. This is the stricter reading.
- Logging: `merge.GSR_LOG`, used only for diagnostics, records every decision with support at least 1.
- Shuffled control: `GSRShuffled` in `scripts/L2_gsr_ob.py`. The R1 bonus, veto and R3 stay real; only the GSR support comes from shuffled gene order.
- Unit tests: `tests/test_gsr.py`. `pytest -q`: 160 passed, 1 skipped.
- Speed-up, opt-in and with no effect on results: `ANCESTRA_FT_CACHE` memoises FAMSA + FastTree Newicks by the sha256 of their input.
  - It is valid because both tools are deterministic on identical input.
  - Evidence: rerunning `fast_R_u` with logging reproduced the TreeFam-A counts exactly; the cache-backed OrthoBench `off` arm reproduced the reference exactly.
- Drivers:
  - `scripts/L2_gsr_ob.py`: OrthoBench, `accurate` base. The synteny-off pass is computed once per process.
  - `scripts/L2_gsr_qfo.py`: QfO `fast`, TreeFam-A dev counts per family.
  - `scripts/L2_compare.py`: keep-rule table.
  - `scripts/L2_gsr_standalone.py` and `scripts/L2_gsr_standalone_qfo.py`: labelling of decisions.
### 4.2 Introns (`scripts/L2_introns_table.py`, `scripts/L2_introns_benchmark.py`)
- The method follows §1.4.
- QC:
  - all 957 dev genes are mapped (C. elegans via gene -> transcript by CDS length);
  - CDS length matches protein length for 953 of 957;
  - genes carry 14.9 introns on average;
  - 99.4% of the 26,049 aligned labelled pairs have a defined score.
- Alignments: 31 ms of CPU per pair.
### 4.3 Tau (`scripts/L2_expr_tau.py`)
- The method follows §1.7.
- The 10 Bgee species include fly and worm, with Bgee gene IDs mapped through the gene-order tables.

## 5. Results
### 5.1 Ideas 1/2 on OrthoBench dev (`accurate` base, 34 RefOGs)
Primary metric: macro pair F1, paired bootstrap over RefOGs (B=2000, seed 11). Source: `results/l2_genome/keep_ob.json` and `ob_arms/*.json`.

| arm | OG-F (Δ) | pair P micro (Δ) | macro pair F1 (Δ [95% CI]) | shuffled macro (s11/23/47) | Δ real − mean shuffled vs seed spread | S8 fission cand. (shuffled) | GSR joins (shuffled) | keep |
|---|---|---|---|---|---|---|---|---|
| off (frozen accurate) | 79.69 | 0.8203 | 0.8026 | - | - | 706 | - | - |
| F_c (single-species, masked) | 77.80 (−1.88) | 0.8221 (+0.0018) | 0.8012 (−0.0014 [−0.0050, +0.0008]) | 0.8024, 0.8026, 0.8026 | −0.0013 vs 0.0002 | 499 (708, 706, 706) | 886 (3, 1, 0) | no |
| F_u (single-species, unmasked) | 81.19 (+1.50) | 0.8227 (+0.0024) | 0.8032 (+0.0006 [−0.0046, +0.0061]) | 0.8024, 0.8026, 0.8026 | +0.0007 vs 0.0002 | 446 (707, 706, 706) | 1300 (22, 18, 21) | no |
| R_c (any unit, masked) | 79.33 (−0.36) | 0.8184 (−0.0019) | 0.8072 (+0.0046 [−0.0026, +0.0144]) | 0.8024, 0.8026, 0.8026 | +0.0046 vs 0.0002 | 476 (708, 706, 706) | 1329 (9, 5, 1) | no |
| R_u (any unit, unmasked) | 82.88 (+3.20) | 0.8188 (−0.0015) | 0.8085 (+0.0059 [−0.0015, +0.0159]) | 0.8024, 0.8026, 0.8026 | +0.0059 vs 0.0002 | 426 (708, 709, 707) | 1873 (65, 54, 48) | no |
| locus_ref (prior-art `locus_repair`, info only) | 81.82 (+2.13) | 0.8190 (−0.0013) | 0.8080 (+0.0054 [−0.0019, +0.0134]) | not run | - | 363 | 2406 locus joins | n/a |

- **OG-F on 34 RefOGs is moved by a handful of RefOGs.**
  - The shuffled seed 11, with 3 joins, already gives OG-F 76.04 (−3.65). A tiny upstream change reorganises RefOG021, the largest RefOG (non-RefOG genes join its OG).
  - Real arms:
    - RefOG007 is fully repaired at the OG level in every real arm (5 groups -> 1, FN 919 -> 0).
    - RefOG023 gains non-RefOG genes and loses 29-35 ortholog pairs, a cascade from the join.
    - RefOG035 and RefOG055 gain pairs (macro F1 +0.04 to +0.13 in R arms).
  - The OG-F differences between arms (F_c −1.9 vs F_u +1.5) sit within this noise.
- **Standalone labelling** (`results/l2_genome/gsr_standalone.json`): only 90-98 of 11.5k-17.8k GSR decisions touch a labelled dev pair (coverage 0.5-0.8%).
  - Support AUC over those: 0.71-0.73, CI [0.54, 0.84].
  - Executed joins with a label: 2 (F_c) and 5 (F_u).
  - OrthoBench dev therefore has almost no power for this rule.
- **Guards:** OG-F (≥ −1.0) passes for R_u, F_u and R_c, and fails for F_c. Pair P (≥ −0.01) passes for all.
  - OrthoXML validity: `verify_orthoxml` PASS for R_u and F_u (`results/l2_genome/verify_ob_*.json`).

### 5.2 Ideas 1/2 on TreeFam-A dev (QfO 2022, `fast` base, 359 dev families)
Primary metric: TPR, paired bootstrap over families (B=2000, seed 11). Source: `results/l2_genome/keep_qfo.json` and `qfo_arms/*.json`.

| arm | TPR (Δ [95% CI]) | PPV (Δ [95% CI]) | shuffled TPR (s11/23/47) | Δ real − mean shuffled vs spread | GSR joins (shuffled) | keep |
|---|---|---|---|---|---|---|
| off (frozen fast) | 0.6151 | 0.9224 | - | - | - | - |
| R_c | 0.6137 (−0.0014 [−0.0141, +0.0098]) | 0.9242 (+0.0017 [−0.0044, +0.0077]) | 0.6151, 0.6151, 0.6151 | −0.0014 vs 0.0000 | 3634 (19, 26, 27) | no |
| R_u | 0.6193 (+0.0042 [−0.0049, +0.0134]) | 0.9205 (−0.0020 [−0.0089, +0.0046]) | 0.6248, 0.6151, 0.6151 | +0.0010 vs 0.0097 | 4277 (46, 63, 64) | no |
| F_c (guard) | 0.6121 (−0.0030 [−0.0120, +0.0041]) | 0.9246 (+0.0021 [−0.0047, +0.0091]) | not run | - | 2516 | no |
| F_u (guard) | 0.6127 (−0.0024 [−0.0132, +0.0064]) | 0.9240 (+0.0016 [−0.0031, +0.0071]) | not run | - | 2930 | no |

- **Where the joins happen.**
  - Real GSR joins (top nodes): Tetrapoda 839-858, Euteleostomi 733-814, Poaceae 471-582, Theria 210-248, Clupeocephala, Amniota.
  - Shuffled joins: Bilateria, Eukaryota, Opisthokonta and the root.
  - The shuffled seed-11 TPR jump (+0.0097) comes from about 46 deep joins. None of them touches a labelled dev pair directly (0 labelled executed joins); the movement is a cascade at deep nodes, where TreeFam-A pairs concentrate.
- **Labelled decisions** (`results/l2_genome/gsr_standalone_qfo.json`):
  - 24 executed joins in R_u and R_c touch TreeFam-A dev pairs. They are speciation joins 79% and 88% of the time, against a 60-62% base rate over all labelled decisions.
  - Over all labelled decisions the support does not rank correct above wrong (AUC 0.46-0.47, CI about [0.37, 0.56]).
  - The shuffled arms' logged decisions include 0 or 1 executed labelled joins (the single one was wrong, at Bilateria).
  - These numbers are small and are not used for a decision.
- **Expectation check.** The pre-registered prior was a small effect on TreeFam-A dev. The prior holds, and the effect is below the shuffled-seed noise.

### 5.3 Idea 3: intron positions, standalone (`results/l2_genome/introns_standalone.json`)
| subset | intron Dice AUC [95% CI] | shuffled AUC (3 seeds) | sequence AUC |
|---|---|---|---|
| all labelled pairs (n=25,885) | 0.655 [0.573, 0.731] | 0.508-0.512 | 0.760 |
| within-RefOG pairs | 0.557 [0.508, 0.648] | 0.512-0.516 | 0.709 |
| hard pairs (n=3,565) | 0.578 [0.514, 0.661] | 0.497-0.499 | 0.654 |

- Mean Dice: orthologs 0.56, within-RefOG paralogs 0.48, deep paralogs 0.01. The signal mostly separates deep paralogs, which sequence already does.
- Hard-pair AUC by class ranges from 0.53 to 0.62. Invert-invert pairs show no signal (0.49).
- **Hard queries (n=695).**
  - Sequence top-1 accuracy is 0.473. Intron re-ranking gives 0.403, a change of −0.071 [−0.149, +0.007]. The intron score is right 167 times, wrong 284 times and abstains 244 times.
  - The shuffled arms change accuracy by −0.012 to −0.026.
- **Stop-rule reading.** The pre-registered stop rule names "hard-case AUC (or the hard-query re-rank gain)".
  - The AUC clears shuffled narrowly. The re-rank test, which is what integration step (b) would do, fails.
  - We read the rule conservatively and **stopped**: integration was not run. This decision was taken before any integration was attempted, and it is recorded here as the interpretation.

### 5.4 Idea 4: tau tie-breaker, standalone (`results/l2_genome/expr_tau.json`)
- Tau is defined for 94.9% of labelled genes.
- Hard queries (n=636): sequence accuracy is 0.550.
- Tau decides 286 of them, with precision 0.448 [0.380, 0.554]. Shuffled precision is 0.394, 0.367 and 0.360.
- Re-rank gain: −0.076 [−0.131, −0.010].
- The CI lower bound of 0.380 is not above 0.47, so the pre-registered rule says **STOP**.
  - Tau agreement carries a little more than shuffled, but less than a coin flip among near-tied candidates, exactly as the S9 profile correlation did.

## 6. Interpretation
1. **Fission repair is real but benchmark-invisible.** GSR fixes many vertebrate-level fissions that the frozen pipeline leaves: S8 candidates fall 29-40% under real gene order and 0% under shuffled. It does so without breaking validity and, in the unmasked arms, with OG-F gains on OrthoBench.
   - The cases are lineage-specific (Tetrapoda, Theria, teleost and grass nodes). OrthoBench dev (34 RefOGs) and TreeFam-A dev labels hardly cover them, and TreeFam-A's deep pairs make it sensitive to tiny deep-node perturbations instead.
   - A test run of these arms would have no power on the current benchmarks. It is not recommended.
   - A fair evaluation needs a benchmark that labels shallow 1:1 orthologs in many vertebrate genomes (for example the QfO VGNC benchmark, or gene-loss/fission counts like S8), pre-registered as a dev metric *before* looking.
2. **Masking vs not.** Masked support (within synteny-off root HOGs) is the safer design in principle, but here it did no better. On OrthoBench it was slightly worse (F_c OG-F −1.9; R_c OG-F −0.4 vs R_u +3.2).
   - Unmasked arms may join across synteny-off root HOGs. The OrthoBench export loses 278 (F_u) and 392 (R_u) top-level groups against `off` (45,237); the masked arms lose 34 (F_c) and 74 (R_c), from singletons joining groups. The OG-F changes do not show these fusions as harmful, but OG-F here is noise-dominated.
3. **The prior-art `locus_repair` no longer hurts on the accurate base.** It gives OG-F +2.1 and macro +0.0054 [−0.0019, +0.0134], against OG-F −3.7 on the old NJ base. The local trees probably absorb its bad joins. This is information only, with no shuffled run, so it is not a candidate.
4. **Gene structure and expression do not help where sequence fails.**
   - Both carry above-shuffled signal overall.
   - Both are worse than a coin among near-tied candidates (intron 37%, tau 45% precision when deciding). The near-tied candidates are mostly young paralogs that keep both gene structure and broad expression patterns.
   - This agrees with S9 and with the expected biology: intron positions are conserved over long timescales in both orthologs and paralogs (Rogozin et al. 2003).

## 7. Cost
| item | wall | CPU |
|---|---|---|
| OrthoBench merges, 18 arms + 1 off pass (accurate; FastTree memo after the first arm) | 0.6 h | 2.7 CPU-h (off arm 1.04 CPU-h without memo) |
| QfO `fast` merges, 11 + 6 logged arms + 2 off passes | 1.1 h | 1.0 CPU-h |
| `verify_orthoxml` (2 arms) | 2 min | 0.03 CPU-h |
| intron table + alignments (3 workers) | 5 min | 0.22 CPU-h |
| tau + labelling + bootstraps | 3 min | < 0.05 CPU-h |
| **total** | about 2 h (in parallel) | **about 4 CPU-h, 0 GPU-h** |
- GSR run-time overhead: QfO `fast` merge 185 s -> 200-223 s. OrthoBench `accurate` merge without the memo: unchanged within noise. Peak RAM < 15 GB.

## 8. Literature credit (all checked on PubMed)
- **Synteny for ortholog assignment.**
  - Jun J, Mandoiu II, Nelson CE (2009) Identification of mammalian orthologs using local synteny. *BMC Genomics* 10:630, doi:10.1186/1471-2164-10-630.
  - Shi G, Zhang L, Jiang T (2010) MSOAR 2.0: incorporating tandem duplications into ortholog assignment based on genome rearrangement. *BMC Bioinformatics* 11:10, doi:10.1186/1471-2105-11-10.
  - Shi G, Peng M-C, Jiang T (2011) MultiMSOAR 2.0. *PLoS ONE* 6:e20892, doi:10.1371/journal.pone.0020892.
  - Synteny-aware refinement in Proteinortho PoFF (Lechner et al. 2014) and OrthoRefine (Ludwig et al. 2024) is credited in LITERATURE.md §5.
  - GSR differs from these: it acts inside a bottom-up HOG merge, only on still-alone units, and is gated by sequence.
- **Intron positions as orthology evidence.**
  - Henricson A, Forslund K, Sonnhammer ELL (2010) Orthology confers intron position conservation. *BMC Genomics* 11:412, doi:10.1186/1471-2164-11-412. This is the closest prior work: an intron-position-conservation score that is higher for orthologs than for the closest non-orthologs. We confirm a weak pairwise signal, and find that it reverses among near-tied hard candidates.
  - Rogozin IB, Wolf YI, Sorokin AV, Mirkin BG, Koonin EV (2003) Remarkable interkingdom conservation of intron positions and massive, lineage-specific intron loss and gain in eukaryotic evolution. *Curr Biol* 13:1512-1517, doi:10.1016/s0960-9822(03)00558-x.
- **Tissue specificity.** Yanai I et al. (2005) Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification. *Bioinformatics* 21:650-659, doi:10.1093/bioinformatics/bti042 (the tau index).

## 9. Failed or not pursued (for S1 Text)
- **Gated synteny repair** (ideas 1/2; joins of still-alone units to their unique best syntenic group, gated by sequence): removes 29-40% of S8 fission candidates (shuffled 0%), but OrthoBench dev macro pair F1 (−0.001 to +0.006) and TreeFam-A dev TPR (−0.003 to +0.004) do not move beyond their CIs or shuffled-seed noise. The benchmarks barely label lineage-specific fissions (<1% of decisions touch a dev label). Not kept.
- **Intron-position conservation** (Dice of shared positions and phases): above-shuffled AUC (hard pairs 0.578 vs 0.50) but below sequence. As a tie-breaker it loses 7 points of hard-query accuracy, because young paralogs keep intron positions. Stopped before integration.
- **Expression tissue-specificity (tau) tie-breaker:** 45% precision when deciding on hard queries (shuffled 36-39%), with re-rank gain −0.076 [−0.131, −0.010]. Stopped, which confirms S9.
- **Not pursued:** QfO `accurate` confirmation runs (no arm passed on `fast`) and intron/tau integration into the merge (stopped by the standalone rule).

## 10. Files
- Code (branch `l2-genome`):
  - `ancestra/merge.py`: GSR, flags off by default.
  - `ancestra/localtree.py`: opt-in FastTree memo.
  - `tests/test_gsr.py`.
  - `scripts/L2_*.py`.
- Results:
  - `results/l2_genome/keep_ob.json` and `keep_qfo.json`: keep-rule tables.
  - `ob_arms/*.json` and `qfo_arms/*.json`: per arm, with a manifest (git commit, config, inputs, runtime) and per-RefOG or per-family counts.
  - `gsr_standalone.json` and `gsr_standalone_qfo.json`.
  - `introns_standalone.json`, `expr_tau.json`.
  - `verify_ob_R_u.json`, `verify_ob_F_u.json`.
- Large outputs (gitignored): `work/l2_genome/` (OrthoXML per arm, GSR decision logs, intron tables, FastTree memo).
- Manifest caveat: the arms ran from the working tree. Its `ancestra/` was committed unchanged in 6745564, so the arm manifests record the parent commit 041e3ee.
