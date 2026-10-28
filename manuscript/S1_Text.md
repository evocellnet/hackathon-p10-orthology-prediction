---
title: "S1 Text. Methods and design choices that did not work"
subtitle: "Supplement to: Ancestra — orthology inference by level-synchronous ancestral-genome reconstruction inside de novo protein families"
geometry: margin=2.2cm
fontsize: 10pt
---

**Note on authorship.** An AI agent planned, ran and wrote up this project, as an exercise in testing whether AI agents can run a bioinformatics research project (see main text, Methods). No human has verified these results yet.

**How to read this file.** Each entry gives the idea, the test, the result and the decision.

- **DEV** marks results on development data: the 34 OrthoBench dev RefOGs, the toy simulator, or the 30% TreeFam-A dev split.
- **TEST** marks results on held-out data. Test data were scored once with the frozen method, unless stated.
- OG-F is the OrthoBench orthogroup F score, in percentage points. Pair F1 is the macro-averaged F1 over within-RefOG ortholog pairs.
- Intervals are 95% paired-bootstrap CIs (B = 2000, seed 11).

The OrthoBench dev set was reused in more than 240 runs, so all DEV numbers are optimistic. On dev, OG-F changes by ±2–4 points between equally good family partitions. Evidence paths refer to the research repository.

---

# 1. Overall design

## 1.1 Progressive merging of whole ancestral proteomes (design v0)

- **Idea.** Walk the species tree from the leaves. At each node, search the two child "ancestral proteomes" against each other, so that cost grows with the number of merges rather than with the square of the number of proteomes.
- **Test.** Before writing code, a reviewer agent estimated the cost on QfO 2022. It used FastOMA HOG counts as a proxy for the size of each ancestral proteome.
- **Result.** With one representative per ancestral gene, the search needs 2.0 × 10^11^ pairs, against 4.9 × 10^11^ for all-vs-all: only about 2.5× fewer. With three representatives it costs 3.6× more than all-vs-all. Four deep nodes carry 78% of the cost, because about 60% of deep ancestral genes have no partner and are passed upwards.
- **Decision.** Rejected. The final design first partitions proteins into families and then merges inside families. Ancestral genes are never searched again: their scores are aggregated from scores between leaf proteins.
- **Caveat.** This is an estimate from a proxy, not a measurement (`docs/reviews/r1_algorithms_compute.md`).

## 1.2 Novelty claims dropped after review

The first review found close precedents for several planned "novel" parts:

- representative sequences per orthogroup come from OrthoFinder/SHOOT;
- in-paralog absorption comes from InParanoid and Hieranoid;
- bottom-up merging with sequence and synteny comes from SYNERGY.

We therefore describe Ancestra as a scalable, synteny-aware descendant of SYNERGY and Hieranoid, not as a new paradigm (`docs/reviews/r1_prior_work_domain.md`).

# 2. Family partition (S0)

## 2.1 Pre-registered completeness target (failed)

- **Target.** At least 95% of dev RefOGs should fall entirely inside one family.
- **Result.** The best of 23 configurations reached 0.794 (27/34) on DEV and 0.722 (26/36) on TEST.
- **Reason.** Six dev RefOGs have a member with no hit at e $\le$ 10^-3^ to the rest of the RefOG. So no partition built from pairwise sequence hits can exceed 0.82.
- **Decision.** The target was declared unreachable with sequence data, and the near-ceiling configuration was accepted (`docs/S0_families.md`).

## 2.2 Single-level clustering (DEV)

One clustering pass, with DIAMOND cluster/DeepClust or MMseqs2 cluster, reached at most 0.59 completeness (e.g. DeepClust at 20% identity: 0.588). The decisive change was to search *all* proteins against the cluster representatives and then take connected components (0.794). The reason is that the representatives of one family are often too distant from each other to produce hits, while their members are not.

## 2.3 Other partition variants (DEV)

| Variant | Result | Why rejected |
|---|---|---|
| Representative-vs-representative linking | 0.706 completeness | below 0.794 |
| Coverage guard (min. coverage 0.3 / 0.5) | 0.471 / 0.441 | loses completeness |
| MCL instead of connected components | 0.765 | largest family exceeded 3% of proteins |
| Leiden cap of 20,000 genes | — | largest family exceeded 3% of proteins |

## 2.4 Non-deterministic partition (bug, fixed)

Two runs of the same configuration gave different partitions (adjusted Rand index 0.78). Dev OG-F moved from 79.2 to 75.6. Two things caused this:

1. ties in DIAMOND's greedy clustering depend on input order;
2. DIAMOND's automatic block size depends on the thread count.

We fixed the input order, pinned the block and chunk sizes, sorted the edges and fixed the Leiden seed. The partition is now byte-identical across thread counts and input orders (`tests/test_s0_determinism.py`). All merge settings tuned on the old partition were re-checked on three partitions. Some earlier gains turned out to be partition noise.

## 2.5 Larger cap for giant families (DEV; failed pre-registered rule)

- **Hypothesis.** False negatives that are split between families come from the Leiden cap of 7,500 genes, which cuts the giant QfO component (63% of all proteins) into 1,213 pieces.
- **Arms.** Raising the cap to 15,000 or 20,000 genes, or adding edges between the pieces.
- **Result.** All arms failed the rule. None of the 1,286 split records became members of the same family.
- **Decision.** The hypothesis was refuted. Most split pairs lie in different Leiden communities at every resolution, or in families that never connect to the component (often *S. pombe*). The small OrthoBench gain of the larger caps came from a side effect in the scoring stage (`docs/S0_giant.md`).

## 2.6 Rescue of families with 3D structure and profiles (DEV positive; TEST negative)

- **Structure-only rescue.** Foldseek 3Di searches of family representatives against each other raised dev completeness from 0.794 to 0.853, above the ceiling for sequence data. The shuffled-3Di control stayed at 0.794. But false joins doubled, and downstream OG-F fell from 72.5 to 63.8: one 3-gene RefOG absorbed 102 foreign genes.
- **Re-test with the final merge.** Combining 3Di, MMseqs2 profiles and a guard in the merge (only reciprocal best hits may attach) passed the dev rule by a small margin: OG-F +0.24, TreeFam-A TPR +0.002. Both changes are well inside noise.
- **TEST.** OrthoBench F fell from 76.2 to 71.1 (Δ −5.0 [−13.1, 0.0]), with precision falling from 69.5 to 60.1. QfO TreeFam-A TPR changed by −0.004.
- **Decision.** This is an option that is off by default. Structural similarity detects homology but does not separate orthologs from paralogs.

# 3. Scoring inside families (S2)

## 3.1 DIAMOND alone, and other score normalisations (DEV)

- **DIAMOND alone.** DIAMOND (`--sensitive`) missed 6–10% of within-family pairs, mostly the divergent pairs that deep merges need. Re-aligning these pairs with MMseqs2 raised toy pair F1 by 0.01–0.03, so the hybrid is the default.
- **Normalisation by shorter self-score.** Dividing bit scores by the self-score of the shorter protein lowered OrthoBench OG-F to 51.1.
- **OrthoFinder-style length correction.** A per-species-pair length correction broke the comparison between same-child and cross-child scores (−0.04 to −0.07 toy F1).

## 3.2 Coverage guard on scored pairs (DEV; reversed)

Dropping pairs with alignment coverage below 50% was meant to stop multi-domain proteins from chaining unrelated genes. On OrthoBench dev it reduced OG-F from 68.1 to 59.0, because OrthoBench contains many fragments and alternative gene models. The fusion risk did not appear. The guard was removed: this was the largest single improvement during development. One side effect of having no coverage filter is that low-coverage true pairs get low scores, and these account for 743 TreeFam-A dev false negatives.

# 4. Merge rules (S3)

## 4.1 Strict reciprocal best hits (DEV)

Without absorption of in-paralogs, toy pair F1 was 0.74 (precision 0.99, recall 0.60), against at least 0.95 with absorption.

## 4.2 Maximum vs average linkage (DEV)

When scores between ancestral genes were aggregated by the maximum over member pairs, OrthoBench dev OG-F was 6.3 points lower than with the average (61.8 vs 68.1). With the coverage guard on, the maximum was slightly better (+1.4). The verdict therefore depends on other settings. The average is the default.

## 4.3 Variants of the in-paralog rule (DEV)

- **Looser rule (`either`).** Lowered OG-F.
- **Outgroup four-point test.** Doubled paralog false positives on hard toy cases.
- **Relaxed absorption factor (0.8).** Failed the pre-stated precision floor on the toy, but gave +3.9 OG-F on OrthoBench under the coverage guard. It was never re-tested after the guard was removed.

## 4.4 Relative-divergence tests (DEV)

The test divides the distance between two genes by their mean distance to an outgroup gene, so that the evolutionary rate of the family cancels out. It was used to veto hidden paralogs and to absorb fast-evolving young duplicates.

- **Toy.** Precision and recall both improved (F1 0.977 vs 0.974).
- **OrthoBench dev.** OG-F fell by 1.3 points, and by 4.7 points together with local trees.
- **Decision.** Off. Like several other rules, it did not transfer from the toy simulator to real data.

## 4.5 Parsimony re-joining at the root (DEV)

This rule was meant to re-join subfamilies that arose from the two vertebrate whole-genome duplications (2R) and were never united at the root. The unrestricted version reduced OG-F to 60.4. The version restricted to mutual best hits was neutral (68.3 vs 68.1). In-family scores between these subfamilies are at noise level. Local gene trees fixed these splits instead (RefOG021 true positives: 506 → 7,261).

## 4.6 Where to build local gene trees (DEV)

| Placement | Dev OG-F |
|---|---|
| Root only (no outgroup) | 61 |
| Shallow nodes only | 65.7 |
| Every node | OG precision 41% |
| **Divergent nodes with an outgroup (adopted)** | **72.5** (vs 68.1 without trees) |

The toy simulator favoured trees at every node; OrthoBench did not.

## 4.7 Rules to raise QfO recall (DEV; pre-fixed six-arm grid)

QfO peek 1 showed low recall (TreeFam-A TPR 0.55 vs 0.74 for OrthoFinder). On the TreeFam-A dev split, 78% of missed pairs were merge failures at the pair's last common ancestor. Six arms were tested against a rule fixed in advance: TPR must rise, PPV may fall by at most 0.01, and OrthoBench OG-F may fall by at most 1.

| Arm | TreeFam-A TPR / PPV | OrthoBench OG-F | Verdict |
|---|---|---|---|
| Species-overlap structure + global absorption | 0.611 / 0.916 | 79.0 | fail (PPV) |
| + relaxed repair ratio | 0.613 / 0.921 | 78.1 | fail (OG-F) |
| + global tolerance 0.8 | 0.648 / 0.916 | 71.2 | fail |
| + shorter-self-score normalisation | 0.587 / 0.925 | 51.1 | fail |
| Deep-node tolerance 0.8 / 0.85 | 0.649 / 0.890; 0.635 / 0.896 | 79.7 | fail (PPV) |

Two options passed and were adopted: species-overlap structure (in the `fast` preset only, see 4.8) and species-overlap absorption at deep nodes. Looser mutual-best rules recover TreeFam-A recall but fuse OrthoBench orthogroups. This is the central trade-off that we did not resolve.

## 4.8 Species-overlap structure in the `accurate` preset (DEV; failed)

This option passed on the `fast` preset (NJ trees). On the `accurate` preset (FastTree) it cost 0.025 PPV and failed the rule. FastTree already resolves many of these cases, and relabelling duplication-with-loss as speciation adds false positives. The two presets therefore use different options. The same option at first produced malformed OrthoXML nesting, which our validity check caught before the final run.

## 4.9 3D-structure scores in the merge (DEV; failed pre-registered rule)

We tested five ways of adding Foldseek 3Di+AA scores to merge decisions, each with three shuffled-3Di controls. Real 3Di beat shuffled 3Di, but no arm beat the no-structure baseline by more than the spread between seeds. The best arm gained +0.015 [−0.004, 0.035] pair F1. 3Di+AA ranks orthologs above paralogs *worse* than sequence does (deep AUC 0.694 vs 0.751), because paralogs share the fold.

# 5. Local gene-tree method (S4; DEV)

On 32 dev RefOGs, FAMSA + FastTree improved the accuracy of duplication/speciation labels over NJ (0.955 vs 0.928; +0.027 [+0.004, +0.060]). The pre-registered verdict was still "keep NJ", because projected cost exceeded the bar by 1%. A later end-to-end arm (dev OG-F 80.4 vs 75.7) was adopted as the `accurate` preset under a separate rule set by the PI agent. **The pre-registered S4 verdict was therefore overridden.**

| Tree method | Δ label accuracy vs NJ (DEV) | Decision |
|---|---|---|
| IQ-TREE LG+G4 | +0.020 [−0.010, 0.044] | rejected: not significant, more cost |
| Partitioned AA + 3Di IQ-TREE | +0.029 [+0.013, 0.064] | rejected: 3.3 CPU-h on the subset |
| 3Di-only trees | −0.043 | rejected |
| Foldseek-distance NJ (FoldTree-like) | −0.004 (better topology, nRF 0.289 vs 0.339) | rejected |
| FastTree with all genes as leaves | OG-F 79.5 vs 80.4 at 3× cost | rejected |

# 6. Guide species tree (S1)

## 6.1 Trees inferred from the data (DEV benchmark; TEST ablations)

| Tree | Fraction of splits conflicting with the published QfO species tree |
|---|---|
| NCBI taxonomy | 0.000 |
| ASTRAL-Pro | 0.120 |
| Concatenation (FastTree) | 0.133 |
| FastME | 0.173 |
| BIONJ | 0.213 |

The inferred trees repeat known long-branch artefacts (for example *C. elegans* outside Bilateria) and conflict at 25–30% among bacteria. The NCBI taxonomy was chosen as the default before any downstream test. Downstream TEST ablations agreed:

- **OrthoBench.** ΔF +0.0 [−0.3, +0.4] for ASTRAL-Pro and FastME.
- **QfO.** TreeFam-A TPR +0.028 [−0.009, +0.077], but EC function consistency −0.12 to −0.13.

A random tree (negative control) cost −16.4 [−27.1, −5.6] OrthoBench F.

## 6.2 Policy arms that were planned but never run

The design planned three guide-tree policies (strict, delayed merge, per-family) and a check with an independent simulator (Zombi or ALF). Neither was done. The per-family local gene trees in S3 partly cover the per-family policy. All toy results come from our own simulator, which a reviewer agent warned "proves little".

## 6.3 Sensitivity to child order

ASTRAL-Pro and FastME produced the same topology with different child order in the Newick string. Their Ancestra outputs differ in 18 OrthoBench orthogroups. We found this after the freeze and did not fix it.

# 7. Gene order (S5)

Constrained gene order is part of the default when gene positions are given. The variants below were dropped.

## 7.1 Synteny bonus added to all scores (DEV)

Adding or multiplying a synteny bonus into every candidate score cut paralog false positives by half. But recall fell by 4–8 points, because the inflated scores broke the in-paralog comparison (toy F1 0.933 vs 0.953). Gene order is therefore used only to choose between candidate partners. This redesign was made after seeing the result.

## 7.2 Partner choice judged at orthogroup level (DEV)

Gene order changed partner choice but not OG-F (68.0 vs 68.1), as expected for a metric based on orthogroups. On within-RefOG pairs it gave +0.032 [+0.009, +0.062] F1. The best of 34 configurations gave +0.046, which is biased by selection. 102 shuffled runs stayed at 0.716–0.730.

## 7.3 Unconstrained gene order inside the tree-based merge (DEV)

Pair F1 rose by 0.04, but orthogroups fused (OG-F 70.2 vs 72.5), because partner swaps linked two root HOGs. The adopted fix was a two-pass scheme. Gene order may only re-rank partners *inside* the root HOG from the pass without gene order.

## 7.4 Fission repair using gene loci (DEV)

The genome-DNA analysis showed that many apparent single-species losses are HOG fissions. Joining a lone single-species unit to the group with the strongest synteny support improved pair F1 (+0.016) and cut fission candidates by 45%. But OG-F fell from 79.2 to 75.5 (old partition; 74.5 vs 75.7 on the deterministic one). Off.

## 7.5 The adopted variant did not replicate on TEST

- **DEV gains.** +4.4 OG-F (`fast`) and +1.0 (`accurate`). Real gene order beat shuffled gene order on all three partitions.
- **TEST.**
  - OrthoBench: Δ −0.1 [−0.4, 0.0] when gene order is removed.
  - QfO TreeFam-A: −0.012 [−0.028, +0.001] (family bootstrap).

The standalone signal is real: the AUC for separating orthologs from paralogs is 0.92 for human–mouse and 0.50 when shuffled. Inside the merge, however, sequence trees and species-overlap rules already resolve these cases.

# 8. Other data layers

## 8.1 Isoforms (S7; DEV; killed)

Scoring gene pairs by their best isoform pair instead of the canonical one:

- AUC 0.760 vs 0.759;
- hard-case gain +0.013 [−0.018, 0.044];
- 2–5× more search.

Killed. A side result: for 21% of ortholog pairs, the best-matching isoform pair is not the canonical pair.

## 8.2 Expression (S9; DEV; killed)

Bgee expression calls were tested as a tie-breaker between paralogs.

- Standalone AUC: 0.667.
- Gain on hard cases: −0.019 [−0.200, +0.149] for the primary set and −0.116 [−0.183, −0.052] for 28 species pairs.

Expression was worse than sequence on the cases it was meant for.

## 8.3 ProstT5 3Di on CPU (feasibility)

On CPU, ProstT5 predicted 120–160 residues per second, so QfO would need 490–980 CPU-hours. Its 3Di states matched AlphaFold DB 3Di at only 69% per protein. We used precomputed AlphaFold DB 3Di instead, which covers 94.5% of QfO proteins. ProstT5 on GPU is evaluated in the main text.

# 9. Evaluation and process

## 9.1 Leaky dev/test split (fixed before tuning)

The first OrthoBench split was random over RefOGs. A reviewer agent found five groups of paralogous RefOGs that straddled dev and test. We replaced it with a split by homology component (dev 34 / test 36 RefOGs). Three early S0 configurations had already been scored on RefOGs that later fell into the test set.

## 9.2 Other deviations

- **TreeFam-A dev split.** It was created after the first QfO peek. Test-family values were computed once at the second peek, one minute before the freeze.
- **Dev metric.** We switched to macro pair F1 because micro pair F1 was dominated by one RefOG.
- **Design review.** The third round was skipped when the freeze moved earlier.
- **Multiple testing.** The planned Holm correction was not applied; all ablations are reported as exploratory.

Full list: `docs/DISCLOSURES.md`.

## 9.3 Pre-registered endpoints that were not met (TEST)

| Endpoint | Result |
|---|---|
| OrthoBench non-inferiority vs OrthoFinder (margin −3 points) | ΔF +6.5 [−6.5, +18.1]: not shown |
| QfO TreeFam-A TPR vs OrthoFinder | −0.160 [−0.210, −0.106]: fails |
| QfO SwissTree TPR vs OrthoFinder | −0.103 [−0.177, −0.029]: fails |
| QfO G-STD2 Eukaryota RF | 0.257: worst of the main comparators |
| QfO EC | 0.899 vs OrthoFinder 0.931 |

## 9.4 Problems with running FastOMA

- **Gene-dropping bug.** FastOMA 0.5.1 dropped 5,497 OrthoBench genes when splitting HOGs. We patched the bug and re-ran FastOMA: test F was 33.1, against 33.5 unpatched, so the bug does not explain its low score.
- **OMA browser HOGs.** Our first labelling rule for the precomputed HOGs gave test F 45.5. We corrected the rule after seeing this score (60.4), a change in OMA's favour.
- **Scaling runs.** The first FastOMA runs failed silently because Nextflow could not list files in a symbolic-link directory.

# 10. Second round of data layers (post-freeze; DEV only)

These arms ran after the freeze, under one shared protocol (`docs/L2_protocol.md`):
- each layer's arms and keep rule were committed before its integration results;
- shuffled controls used seeds 11/23/47;
- the same keep rule as before applied;
- no test data were used.

No arm passed. Details are in `docs/L2_{network,domains,plm,genome,gpu}.md`.

## 10.1 Interaction networks and co-fractionation

- **Idea.** Orthologs keep their interaction partners. For a candidate pair, the score counts the partners that are linked through the run's own confident 1:1 orthologs.
- **Data.** STRING v12, using only each organism's own experimental and co-expression evidence. Transferred, text-mining and homology channels were excluded, since they are derived from orthology. The human co-fractionation and AP-MS complexes of hu.MAP 3.0 were added.
- **Result.** All-pair AUC was 0.630 [0.532, 0.699] on OrthoBench and 0.622 on TreeFam-A, against 0.50 for shuffled networks. On hard pairs it was 0.501 [0.344, 0.687] and 0.569 [0.501, 0.625]. Only 3.4% of OrthoBench candidates have any network data. On hard queries the layer was right in 4 of 9 (OrthoBench) and 27 of 52 (TreeFam-A) of its decisions. Pairs added as recall rescue were orthologs only 48–54% of the time.
- **Decision.** Stopped at the standalone gate. Most hard cases lie in species with almost no interaction data (fish, chicken, dog, opossum, Ciona). Where data exist, paralogs share partners.

## 10.2 Protein abundance (PaxDb)

- **Result.** All-pair AUC 0.676 [0.604, 0.732] (hard 0.689). On hard queries it was right in 39 of 116 decisions, and it lowered accuracy (−0.131 [−0.246, −0.046]).
- **Decision.** Stopped, as expression was earlier.

## 10.3 Domain architecture (Pfam)

- **Idea and data.** Jaccard similarity of Pfam clan sets, from `hmmsearch --cut_ga` with Pfam-A 38.2. We annotated the 58k proteins in OrthoBench dev families and the 446k proteins in TreeFam-A dev families.
- **Standalone.** Hard-case AUC was 0.637 [0.547, 0.738], against 0.743 for sequence. Architecture separates different orthogroups almost perfectly (0.979) but not orthologs from in-paralogs within one orthogroup (0.54, below its within-family shuffle of 0.57–0.58). Re-ranking hard cases changed accuracy by +0.001.
- **Integration.** The best arm relaxed the mutual-best tolerance to 0.8 only for identical architectures. TreeFam-A dev F1 rose by +0.018 [+0.003, +0.034] and TPR by +0.038, and the arm beat its shuffled control. But PPV fell by 0.025 and OrthoBench OG-F by 2.9, so it broke both guards. An existing sequence-only arm, run on a different base, had higher TPR at similar PPV. A bonus arm cost 3.0 OrthoBench OG-F. The veto arm's lower CI bound was exactly 0.
- **Cost.** 31 CPU-h for annotation alone. Annotating all of QfO would cost about 60 CPU-h.
- **Decision.** Rejected for the merge. The between-orthogroup signal could serve as a family-fusion guard in S0; that is untested.

## 10.4 Protein language models

- **Pretrained models.** We tested ESM-2 (35M, 650M, 3B) and ProtT5, each as mean-pooled cosine and as an embedding-based alignment score. On all triplets they beat shuffled embeddings by 0.23–0.33 in ranking accuracy. On hard triplets the difference was −0.048 to +0.139, and most CIs included 0. All stayed below sequence: 0.61–0.72 on deep pairs, against 0.757. On hard cases where sequence ranks wrongly, six of seven arms scored 0.24–0.41, below chance, and ProtT5 alignment scored 0.53, at chance. Their errors tend to follow the sequence errors. The best combination with sequence gained +0.022 [−0.001, +0.046] on deep hard cases, one of about 40 comparisons.
- **A model trained on our own data.** We trained a contrastive head on frozen ESM-2 3B embeddings. Positives were synteny-supported 1:1 pairs among OrthoBench species, and every family with a dev or test gene was excluded. The head improved on raw embeddings over all triplets (+0.040 [+0.010, +0.056]; shallow triplets +0.051 [+0.026, +0.074]), but not on hard ones (−0.007 [−0.065, +0.056]), with seed 11. The reason is that the rule that makes labels trustworthy (an unambiguous anchor) removes all sequence near-ties from training. A post-hoc variant trained on near-ties decided by synteny did not beat its random-winner control.
- **Nearest-neighbour family rescue.** Mean-pooled cosines cluster near 1, so the pre-registered guard let 90% of links through, and 229k of 251k genes chained into one family. Killed.
- **Not pursued.** ESM-C (licence and package version), LoRA fine-tuning, 3Di tokens as extra input, and pLM distances for NJ trees.

## 10.5 Gene order and gene structure, second attempt

- **Gated synteny repair.** A unit still unmatched at a node may join its uniquely best syntenic partner, if the sequence score is within 0.75 of its best and the join creates no duplication. With real gene order, the number of recomputed S8-style fission candidates fell by 29–40% (706 → 426–499); shuffled gene order left it unchanged. Only 7 executed joins touched a labelled OrthoBench pair, and 5 were correct. On the benchmarks the change was within noise: OrthoBench dev macro pair F1 +0.0059 [−0.0015, +0.0159] and TreeFam-A dev TPR +0.0042 [−0.0049, +0.0134], with PPV within 0.002. Fewer than 1% of its decisions touched a labelled pair. Rejected under the keep rule. It would need a benchmark that labels shallow 1:1 orthologs.
- **Synteny anchors for recall at the LCA.** Same machinery and same outcome (included in the numbers above).
- **Intron positions** (Ensembl 99 GFF3). AUC 0.578 [0.514, 0.661], against 0.50 shuffled and 0.654 for sequence. As a tie-breaker it changed hard-case accuracy from 0.473 to 0.403 (−0.071 [−0.149, +0.007]). Young paralogs keep intron positions as well as orthologs do. The coded stop rule passed on AUC, but integration was skipped on a conservative reading of the rule. That was a judgement call, logged in `docs/L2_genome.md`.
- **Expression, second attempt** (tissue specificity τ as a tie-breaker). Precision 0.448 [0.380, 0.554] against 0.394 shuffled (CIs overlap); accuracy fell by 0.076 [−0.131, −0.010].

## 10.6 GPU ideas left open in the first GPU round

- **Sparsifying the all-GPU score graph.** We ran seven arms, a random-sparsification control and a score-jitter floor. Every all-GPU arm lost OrthoBench dev accuracy. On that benchmark, a 1% random jitter of the CPU scores alone gave −3.2 to −7.3 OG-F (two of three seed CIs include 0), while TreeFam-A barely moved. The best QfO arm (keep pairs scoring at least 0.207 × the query's best hit in that species) gave TreeFam-A dev TPR/PPV 0.637/0.930 at 0.59 CPU-h, against 0.615/0.922 at 9.05 CPU-h, but it lost 8.4 OrthoBench OG-F. The 0.59 CPU-h is an estimate and excludes 0.2–0.3 GPU-h.
- **Batched FastTree.** `FastTree -n` saved 1.6% of merge wall time and changed 4 of 1.39 M pairs. The first GPU report had blamed process start-up for the cost; in fact about 93% of each call is the ML optimisation (measured on 20 small jobs). In-process FAMSA (pyfamsa) was not faster and changed 8 of 998 alignments. Not enabled.
- **Structure rescue at full 3Di coverage.** ProstT5 on the GPU filled the 63k QfO proteins without AlphaFold DB structures in 5 h. Metrics were the same as with AlphaFold DB coverage alone (OrthoBench macro pair F1 +0.0042 [−0.0005, +0.011]; shuffled 3Di kept +0.0021 of it), so the rescue still fails the keep rule. Foldseek's GPU search was slower than its CPU search here (47.6 vs 29.0 s for 3,000 queries) and returned a different hit set.
