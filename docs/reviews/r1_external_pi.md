# Review round 1: external PI (comparative genomics / phylogenomics)

Reviewed: DESIGN.md draft v0. Clock at review: 2026-10-26 15:06 CET, so T+0h26 by START.txt (not 1.5 h). Benchmark data look ready.

## Verdict
The idea is sound and defensible. Its weak point is novelty against Hieranoid 2. The progressive N−1 merge, ancestral representatives and in-paralog rule are all Hieranoid. The claim that the method is "conceptually different" rests on one thing: gene order as a first-class matching signal, reconstructed jointly with orthology. The design should say this plainly and give that component the most rigour. The method is plausibly competitive on QfO precision, because Hieranoid 2 is already on the SwissTree, STD and G-STD2 fronts. Its likely weaknesses are recall and OrthoBench fissions. In 48 h it is achievable only if layers S6 to S9 are cut back hard.

## Scientific soundness: specific concerns

1. **The merge is under-specified where it matters.** A merge between two *ancestral* proteomes is not the 2-species InParanoid problem. Three cases must be handled explicitly:
   - (i) Duplication on a child branch: in-paralogs, which is fine.
   - (ii) A duplication older than the parent node: two ancestral genes on each side must be paired a1–b1 and a2–b2, not grouped. This is a matching problem, and it is exactly where synteny helps.
   - (iii) Differential loss after an old duplication: hidden paralogy produces a false 1:1 match.

   With score rules alone, case (iii) is undetectable. A gene tree with an outgroup is the standard remedy. I recommend a hybrid. Use score-based assignment for unambiguous 1:1 clusters. For ambiguous or many-to-many clusters, build a small tree (FAMSA + FastTree) over the representatives plus the nearest outgroup ancestral gene, and label it with species overlap. The cost stays bounded because only ambiguous clusters need a tree. This makes the method more than "Hieranoid with extra features", and it is where precision is won.

2. **Synteny is circular unless it is staged.** A neighbour's orthology depends on the matching that synteny is meant to inform. Use three passes:
   - (a) Anchors: unambiguous reciprocal-best 1:1 matches.
   - (b) Synteny scored as a GOC-like window (±k ancestral genes, tolerant to gaps and inversions) against the anchors only.
   - (c) Resolution of ambiguous cases.

   Ancestral adjacencies must allow gaps, because unmatched genes break them. Expect the signal to vanish above roughly Vertebrata or Saccharomycetaceae. VGNC is already saturated (TPR/PPV about 0.98/0.99). The measurable gain will therefore come from paralog-resolution precision on SwissTree vertebrate families and on OrthoBench. Pre-register a modest effect size.

3. **"Unmatched passes up" needs a precise semantics.** Under that rule, a gene missed by the search at node P gets its HOG placed too deep or split. On OrthoBench this shows up as fissions. Add a final rescue pass: a profile search of root-level singletons and small HOGs against the HOG representatives, with merges allowed only where species do not overlap.

4. **"Linear in S" is only conditionally true.** Merge cost is roughly |A|·|B|·k² representatives. Ancestral proteome size grows with the clade's pan-genome, and that size does not saturate for prokaryotes. The root merge (Bacteria, Archaea, Eukaryota) will be the largest single job. Log ancestral proteome size and wall time per node. State the claim as "linear in merges, with cost per merge bounded by pan-genome size". At AFDB scale, partition the root merges by family (Linclust on representatives) so that searches stay local.

## Species tree (inadequate as drafted)
S1 benchmarks only the topology (nRF) and the downstream F. The brief asks about *alternatives to a universal tree*, which is a separate question. The design needs:
- **Polytomies.** NCBI taxonomy is multifurcating and merges are binary. The rule for ordering merges inside a polytomy must be stated, for example sketch distance or a multi-way merge.
- **Robustness and negative controls.** Run the full pipeline with the reference tree, the NCBI tree, a tree perturbed by k random NNI/SPR moves, and a random tree. If the random tree barely hurts, the tree is not doing work. If the perturbed tree collapses accuracy, the method is fragile. Both outcomes are publishable findings.
- **Per-family discordance.** At each merge, flag families whose local tree under species overlap contradicts the guide split (HGT or ILS). Then compare (i) strict universal, (ii) universal with delayed merge, and (iii) per-family local labelling for flagged families. Report the fraction flagged, split into prokaryotes and eukaryotes. With 30 of 78 QfO species being prokaryotes, HGT is common, not an edge case. STD_Bacteria and G-STD2-Luca are the relevant tests.
- **HOG levels.** Name levels by taxonomy clade, which is robust to polytomies, and keep the guide tree for merge order only.

## Competitiveness: set concrete bars
Primary comparators are Hieranoid 2 (SwissTree .726/.950, STD-Euk RF .034) and FastOMA/InParanoid5. Success means at least Hieranoid 2 precision with higher recall on SwissTree, TreeFam-A and G-STD2-Euk, with both the HOG and pair exports reported. The OrthoBench bar is OF3's published numbers. Hieranoid scored poorly on OrthoBench in Emms & Kelly 2020, so the gap to OF3 must be shown honestly. FastOMA local runs use the OMA database, which contains the QfO species. Flag this as leakage in every comparison.

## Decomposition and verification
- S3 tests only leaf-level merges. Add **S3b, an internal-node merge benchmark**, using simulated truth. Use ALF (with gene order, duplication, loss and LGT; it needs Darwin, which the QfO pipeline installs anyway) or Zombi (Python, with gene order and DTL). Both are faster and better validated than a home-made simulator. Simulate at two divergence regimes and with LGT so that S1's HGT question can be tested with truth.
- Add a **hierarchy sub-benchmark**: duplication-node placement against reconciled TreeFam-A/SwissTree trees. Run it only at the final freeze, to keep the hold-out clean.
- Tuning: no more than three free weights in the combined score, tuned on simulations plus OrthoBench. Log every tuning run. OrthoBench has only 70 RefOGs, so overfitting is easy.

## Data-layer realism and cuts
- **S7 isoforms are cheap.** The QfO 2022 tarball ships `*_additional.fasta` and `gene2acc`. Keep as *should*.
- **S8 gene-loss validation.** `*_DNA.fasta` holds CDS, not genome assemblies. Genome assemblies for flanking-locus tblastn would have to be downloaded. Cut S8 to a *could* demo on OrthoBench mammals (Ensembl genomes, about 5 species), time-boxed to 3 h.
- **S6 structure.** ProstT5 fp32 weights were not downloaded, and CPU throughput on about 1M proteins is unknown. Use AFDB 3Di where it exists, and apply it only to unmatched or ambiguous representatives at deep merges. Run a 1 h feasibility probe before committing.
- **S9 expression.** Time-box to 2 h with a pre-stated negative-result criterion. Bgee/OMA circularity should be noted.
- **Runtime fairness.** On a shared 24-core node, runtime comparisons are confounded by concurrent jobs. Run competitor timings sequentially at a fixed thread count (for example 8 threads), on 8/16/32-species subsets, starting now in the background. OF3 on all 78 species is optional.

## Schedule
Add a **go/no-go gate at T+12 h**. If v0 is below Hieranoid 2 on the QfO dev subset or OrthoBench, switch to the hybrid tree-at-ambiguous-clusters merge (concern 1) before adding any layer. The synteny layer (S5) should start by T+14 h at the latest, because it carries the novelty claim.
