# Review R1: statistics and benchmarking (DESIGN.md v0)

Clock: 15:06 = T+0h26 by START.txt, not ~1.5 h as the task states.

## Verdict
The test set is not yet protected. Tuning sets overlap the test labels in four places, no uncertainty is quantified, "comparable or better" has no decision rule, and runtime ignores node contention. All are cheap to fix now.

## 1. Dev/test leakage (most important)
1. **VGNC is a QfO test**, yet it is the dev truth for S3, S5 and S9. That is direct leakage. VGNC/HCOP assignments also partly rest on synteny, so synteny gains on VGNC are partly circular.
2. **S2 uses SwissTree and TreeFam-A**, while §4 says they are "not inspected". These contradict each other.
3. **"A subset of QfO species for tuning" protects nothing.** The labels are per pair, and every species subset contains test pairs.
4. **OrthoBench is both the tuning set and a headline result.** After tuning on all 70 RefOGs, its score is a training score.
5. **The S3 truth is wrong, not just leaky.** RefOGs are Bilateria-level orthogroups, so a human–mouse pair inside a RefOG can be a pair of paralogs. Derive speciation/duplication labels from the rooted RefOG trees in `Open_Orthobench/Supporting_Data/Data_for_RefOGs/trees/` using species overlap. Check 5 RefOGs by hand.

**Proposed split** (commit the split file with its hash; pre-register):
- **Dev:**
  - toy simulator;
  - 35 OrthoBench RefOGs (fixed seed, stratified by size);
  - a hash-based 30% of TreeFam-A families. RefPhylo raw output carries the family ID, so per-family scoring works.
- **Test:**
  - the other 35 RefOGs;
  - 70% of TreeFam-A;
  - all of SwissTree, VGNC, GO, EC, FAS, STD and G-STD2.
- Log every test-split evaluation in `results/eval_ledger.tsv` (config hash, split, time, purpose). Allow at most 2 pre-freeze looks, and report the count in REPORT.md.

## 2. Endpoints and decision rule
**Primary endpoints:**
- OrthoBench-test pooled F (with precision and recall);
- SwissTree and TreeFam-A-test TPR/PPV;
- G-STD2-Euk RF together with n samplings (RF without n is meaningless, because n tracks recall);
- EC Schlicker together with n.

**"Comparable"** means that on each primary test, no named comparator dominates us by more than the 95% CI. The comparators are Hieranoid 2, OrthoFinder, FastOMA, SonicParanoid2-sens-g and InParanoid5.

**Caveats to state up front:**
- Species-tree-guided merging makes predictions STD-consistent by construction (Hieranoid has the best STD-Euk). STD gains need this caveat.
- S1 nRF against NCBI is trivially 0 for the NCBI arm, and NCBI has polytomies. Use published phylogenies and RF on resolved splits.

## 3. Uncertainty
- **OrthoBench:** bootstrap over RefOGs (B = 2000). This needs per-RefOG TP/FP/FN from `bench_repo/scripts/orthobench_eval.py`, which currently prints only pooled values rounded to 2 dp. Also report macro-F, because pooled F is dominated by large RefOGs.
- **SwissTree / TreeFam-A:** bootstrap over families. **VGNC / GO / EC:** over genes. **STD:** binomial CI over samplings.
- **Comparisons are always paired.** Use the same resamples for both arms and report the CI of Δ. With 35–70 RefOGs, expect only ~1–3 F points to be resolvable. Say so in advance.
- **Seeds:**
  - Report the seed SD for {11, 23, 47} next to the bootstrap CI.
  - Add a determinism test: same seed and threads must give an identical output hash.
  - The control scripts default to seed 42. Align them with the frozen list.
- **Guide tree and DIAMOND mode** cause more variance than seeds. Treat them as explicit factors.

## 4. Ablations and layer effects
- **Two arms on dev, then one confirmation on test:**
  - add-one: sequence-only + layer;
  - leave-one-out: best − layer.
  - Add a single interaction check (synteny × structure). A full factorial is unaffordable.
- **Mandatory shuffled-layer controls.** Keep each layer's coverage and score distribution but destroy its signal:
  - synteny: permute gene order within each chromosome;
  - 3Di: shuffle strings among proteins of similar length;
  - expression: permute gene labels.
  - Effect = real − shuffled. Otherwise gains may only reflect extra tie-breaking or parameters.
- **Coverage strata.** Synteny exists for only ~15 QfO eukaryotes, and 3Di is used only at deep nodes. Report effects both within the covered stratum and overall.
- **Fix the S5 "ambiguous paralog" stratum before looking**, for example top-2 sequence scores within X%.
- **Multiplicity:** one primary metric per layer, Holm correction across layers. The rest is exploratory.
- **Nulls:** report them with an equivalence bound (e.g. |ΔF| < 1.5).
- **S8:** fix the curated loss set and its source before running, and report an exact binomial CI.

## 5. HOG → pair flattening
- **Why it matters.** QfO pairs genes across the children of every nested orthologGroup. One wrong deep merge adds n₁·n₂ false pairs, which dominate PPV (compare OMA HOGs with OMA Pairs).
- **Unit test.** On toy OrthoXML (nested groups, root paralogGroup, singletons), our flattener must give pair sets identical to `map_relations.py`.
- **One primary export.** It must be a deterministic function of the same HOGs: for example pairs from levels at or below a fixed depth, or pairs restricted by an InParanoid score. Its rule is fixed on dev. Choosing whichever export scores better on test is cherry-picking.
- **Diagnostics:** report pairs per family, and the share of false-positive pairs that come from top-level versus lower HOG levels.
- **OrthoBench:** state the exported HOG level, and pass the official coverage and one-gene-one-OG checks.

## 6. Runtime and scaling
- **Contention.** Parallel agents share the 24 CPUs, so timings need an exclusive window. Log the load before each run, and reject a run that starts with load above 1.
- **Measurement:**
  - `/usr/bin/time -v`: wall time, CPU time, peak RSS;
  - 24 threads for every tool;
  - median of 3 repeats at 16 species or fewer.
- **Subsets:** nested, phylogenetically stratified subsets of 8/16/32/64/78 species, with fixed seeds. Draw 2 independent subsets at 16 and 32 to expose composition variance.
- **Scaling claims:**
  - Plot runtime against total proteins N, not only against S.
  - Fit log–log slopes with CIs, and claim "near-linear" only if the CI upper bound is below ~1.2.
  - Anything beyond 78 species is extrapolation.
- **Merge size is not bounded.** Unmatched genes pass up, so merges near the root approach O(N). Log the time and input size of each merge.
- **Include everything in our runtime:** ProstT5/3Di, database builds and GFF parsing. Report a per-layer breakdown.
- **FastOMA** needs a prebuilt OMAmer database. Report its time without the database build, and flag that.
- **Cap competitors** at 6 h and report capped runs as censored, not dropped.
- **Run competitors on OrthoBench** for same-input accuracy and runtime.

## 7. Additional controls
- Re-scored FastOMA and OrthoFinder must match OpenEBench within 0.005.
- Add these arms:
  - sequence-only Ancestra as the reference for every layer;
  - a random guide tree;
  - RBH-only inside our pipeline, to separate the merge logic from the search.

## 8. Justified claims
- **Justified:**
  - per-test Pareto statements with CIs;
  - measured scaling up to 78 species;
  - layer effects that exceed the shuffled control in the covered stratum.
- **Not justified:**
  - a single "better" verdict across tests;
  - "scales to AFDB";
  - STD gains without the guide-tree caveat.
