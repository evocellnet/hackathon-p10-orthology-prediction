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
