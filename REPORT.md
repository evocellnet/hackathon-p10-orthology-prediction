# REPORT.md — Ancestra (skeleton at M2, T+~6h; to be completed)

Status: **no test-set number appears here.** All numbers are dev-only and not independently reviewed. The claim IDs in brackets (e.g. [C-CORE-12]) map to evidence files in `CLAIMS.md`. Do not add a number without an ID or a path. "Dev" means OrthoBench dev RefOGs (34, `splits/orthobench_split.json`), the toy simulator, or a sub-benchmark with no test split.

## 0. Summary (TODO at end)
TODO: a verdict against the non-inferiority rule (DESIGN.md §5), the primary endpoints and runtime. Draft:
- OrthoBench dev OG-F with a deterministic partition is 75.7 with gene order and 71.3 without. OrthoFinder scores 70.5 on the same 34 RefOGs [C-CORE-12, C-BASE-1]. These are point estimates with no CI.
- Only synteny beat its shuffled control. Structure helps S0 but not the merge. Expression and isoforms were rejected by pre-registered rules [C-S5-6, C-S6-1, C-S6-3, C-S9-2, C-S7-1].
- QfO and OrthoBench test performance are unknown.

## 1. Approach
TODO: one paragraph from DESIGN.md §1: family partition, in-family scoring, level-synchronous merge along the guide tree, ancestral gene-order context. Pure progressive pan-proteome merging was rejected at design review (about 2.5x saving over all-vs-all) [C-INF-7].

| Stage | Method | Sub-task benchmark | Doc | Code |
|---|---|---|---|---|
| S0 families | DIAMOND cluster, genes vs representatives, connected components, Leiden cap 7,500 | RefOG completeness | docs/S0_families.md | ancestra/families.py |
| S1 guide tree | NCBI taxdump (default) vs data-derived arms | conflict vs published QfO tree | docs/S1_speciestree.md | ancestra/speciestree.py |
| S2 scoring | in-family DIAMOND plus MMseqs2 rescue | toy truth, runtime | docs/core_v0.md | ancestra/scoring.py |
| S3 merge | k-way level-synchronous merge, in-paralog absorption, fission repair, local NJ trees | toy truth; dev OG-F and pair F1 | docs/core_v0.md | ancestra/merge.py, localtree.py |
| S5 synteny | ancestral adjacency, constrained re-rank, veto | species-pair AUC; pair-level dev F1 | docs/S5_synteny.md | ancestra/synteny.py |
| S6 structure | 3Di rescue of S0 families (optional, off) | completeness; merge arms | docs/S6_structure.md | ancestra/structure.py |
| S7 / S9 | isoforms / Bgee expression (both rejected) | hard-case gain | docs/S7_isoforms.md, docs/S9_expression.md | ancestra/isoforms.py, expression.py |
| S8 gene loss | locus search with miniprot (add-on) | control recovery | docs/S8_geneloss.md | ancestra/geneloss.py |

Not benchmarked on their own: gene-tree inference with AA vs AA+3Di, and a reconciliation step. The only tree step is the local NJ inside S3 [C-CORE-9]. TODO: run an AA vs 3Di local-tree arm before the freeze, or list it as a follow-up.
TODO: frozen config hash, git commit, command line and seeds (11, 23, 47) from `manifest.json`. Freeze is planned for T+36h.

## 2. Best methods
TODO: name the frozen default and the reported ablation arms. Current dev-selected default (not frozen):
- S0 `2sA_dc30_cc_sh30_RA`, hybrid scoring, average-link merge with local trees, NCBI guide tree.
- With `--gene-order`, constrained synteny (w=10, λ=0.3, core anchors, R2/R3 veto) [C-CORE-10, C-CORE-12].
- Off or rejected: coverage guard, root parsimony merge, relative-divergence test, structure in the merge, expression, isoforms. Locus-aware fission repair is kept only as an ablation arm [C-CORE-4, -6, -7, -13, C-S6-3, C-S9-2, C-S7-1].

## 3. Test-set results (TODO: do not fill from dev)
**3.1 OrthoBench test (36 RefOGs)**

| Method | F | P | R | exact | paired-bootstrap ΔF vs Ancestra [95% CI] |
|---|---|---|---|---|---|
| Ancestra, gene order | TODO | TODO | TODO | TODO | n/a |
| Ancestra, no gene order | TODO | TODO | TODO | TODO | TODO |
| OrthoFinder | TODO | TODO | TODO | TODO | TODO |
| SonicParanoid2 | TODO | TODO | TODO | TODO | TODO |
| FastOMA | TODO (if run) | | | | |

Rule: comparable if the lower 95% bound of Δ is ≥ −0.03. Without a CI, write "inconclusive".

**3.2 QfO 2022 (78 species).** Local scoring was validated against published submissions [C-INF-2].

| Endpoint | Ancestra +GO | Ancestra −GO | OrthoFinder | FastOMA | Hieranoid 2 |
|---|---|---|---|---|---|
| SwissTree TPR/PPV | TODO | TODO | TODO | TODO | TODO |
| TreeFam-A TPR/PPV | TODO | TODO | TODO | TODO | TODO |
| G-STD2 Eukaryota (RF, n) | TODO | TODO | TODO | TODO | TODO |
| EC (Schlicker, n) | TODO | TODO | TODO | TODO | TODO |
| Exploratory: GO, VGNC, SP, STD_* | TODO | TODO | TODO | TODO | TODO |

Peeks used: TODO, from `results/test_log.tsv`. Competitor rows come from `bench_repo/reference_results/openebench_QfO_2022.tsv`. VGNC is partly synteny-derived, so it is not evidence for the synteny layer.

## 4. Sub-task benchmarks (all dev)
**S0.** Selected config: completeness 0.794 (27/34), fission 1.32, Σg²/N² 0.0075, max family 7,478, 263 s, 6.2 GB [C-S0-2]. The pre-registered target of ≥95% is refuted: the ceiling is about 0.82 [C-S0-1, `results/S0_families/dev_ceiling.json`]. The key step is searching all genes against the representatives; single-level clustering reached ≤0.59 [C-S0-6]. S0 was made deterministic after an ARI of 0.78 was found across thread counts [C-CORE-11]. Test completeness is a single logged peek (C-S0-4); TODO copy it here only with the final test tables.

**S1.** Conflict vs the published QfO tree (`results/S1_speciestree/metrics.json`):

| Tree | Conflict | ID |
|---|---|---|
| NCBI | 0.000 (nRF 0.128) | C-S1-1 |
| ASTRAL-Pro | 0.120 | C-S1-2 |
| FastME | 0.173 | C-S1-2 |
| BIONJ | 0.213 | C-S1-2 |
| random | 0.987 | add to ledger |

**S2/S3.** On the toy, pair F1 is 0.94–0.98 and the RBH-only arm is 0.74 [C-CORE-1]. Hybrid scoring recovers 6–10% of the pairs DIAMOND misses [C-CORE-2].

| OrthoBench dev step | OG-F | Source |
|---|---|---|
| coverage guard on | 59.0 | `results/core_v0_orthobench/runs.tsv` |
| guard off | 68.1 (P 71.2, R 65.3) | C-CORE-4 |
| + local trees, old partition | 72.5 | C-CORE-9 |
| + synteny + veto, old partition | 79.2 (superseded) | C-CORE-10 |
| deterministic, no gene order | 71.3 | C-CORE-12 |
| deterministic, gene order | 75.7; macro pair F1 0.772 | C-CORE-12 |

- Partition sensitivity is ±2–4 OG-F, so 79.2 vs 75.7 is noise [C-CORE-12].
- Dev error classes: merge splits 20/34, over-merges 13/34, S0 splits 7/34 [C-CORE-5].
- Average link beat max link, 68.1 vs 61.8 (`runs.tsv`; add to ledger).
- TODO: paired bootstrap and replication on ≥2 partitions for every arm in the final table.

**S5.**
- Coverage: OrthoBench about 100%, QfO 0.931 [C-S5-1].
- Ortholog-vs-paralog AUC at w=20: human/mouse 0.92, human/chicken 0.93, zebrafish/Tetraodon 0.82, human/zebrafish 0.70, human/fly 0.49. Shuffled order gives 0.50 [C-S5-2].
- Hard cases (n=157): synteny decides 70 with precision 0.84; re-rank gain +0.16 [0.05, 0.26] [C-S5-3].
- Pair F1 in the merge: +0.032 [+0.009, +0.062] at w10 λ0.3. The best arm, +0.046 [+0.019, +0.072], was chosen among 34 configs. 102 shuffled runs stay at 0.716–0.730 [C-S5-6, -7, -8].
- A naive additive bonus costs 4–8 recall points on the toy [C-S5-4].

**S6.**
- S0 rescue: completeness 0.794 → 0.853, shuffled 27/34 [C-S6-1].
- Downstream the merge falls from 72.5 to 63.8 OG-F [C-S6-2].
- In merge scores, the best gain is +0.015 [−0.004, 0.035], within seed spread [C-S6-3].
- Deep AUC is 0.694 for 3Di vs 0.751 for sequence [C-S6-4].
- CPU ProstT5 is infeasible at QfO scale; AFDB 3Di covers 94.5% of QfO and 80.5% of OrthoBench [C-INF-8, -9].

**S7.** AUC 0.760 vs 0.759; hard-case gain +0.013 [−0.018, 0.044]; 2–5x search cost; rejected [C-S7-1]. For 21% of ortholog pairs the best isoform pair is not the canonical pair [C-S7-2].

**S8.**
- Controls: 92.2% recovery, 0/512 random-window hits, loss-call precision 0.80 [C-S8-2].
- Of 2,523 apparent losses, 53% are not losses (808 are fission errors); 25% are confirmed; 21% are unresolved [C-S8-1].
- GULO and UOX are recovered [C-S8-3].
- This was measured on the old-partition HOGs. TODO: rerun on the final default.

**S9.** AUC 0.667 [C-S9-1]. The hard-case gain is −0.019 [−0.20, 0.15] and −0.116 [−0.18, −0.05]; rejected [C-S9-2]. Conservation: 16% of ortholog pairs conserved vs 26% broken [C-S9-3].

**Toy controls.** TODO: positive and negative controls, and the Zombi/ALF independent-simulator check if done.

## 5. Runtime
Walls come from a shared node (load 20–90) and are upper bounds.

| Run | Proteins | Threads | Wall | CPU-h | RSS | Source |
|---|---|---|---|---|---|---|
| Ancestra OrthoBench, deterministic, gene order | 251k | 12 | 31.5 min | TODO | 14.7 GB | C-CORE-15 |
| OrthoFinder 3.2 | 251k | 8 | 1.30 h | 8.2 | 5.6 GB | C-BASE-1 |
| SonicParanoid2 | 251k | TODO | TODO | 7.9 | TODO | `bench_repo/results/orthobench_competitors/sonicparanoid_default/runtime.json` (add to ledger) |
| Ancestra QfO dry run (old S0) | 989k | 12 | about 2.1 h | S0 3.2 + scores 8.8 only | 29.9 GB | C-DOC-1 |
| Ancestra QfO frozen | 989k | TODO | TODO | TODO | TODO | TODO |
| Nested subsets 8/16/32/64/78 | | 8 | TODO | TODO | TODO | `results/scaling/runs/` (partial, unreviewed) |

TODO: log-log slope with CI, stage breakdown, preprocessing counted, GPU outlook.

## 6. Effect of each extra data layer

| Layer | Sub-task benchmark effect (dev) | Full pipeline: OrthoBench test | Full pipeline: QfO | Decision |
|---|---|---|---|---|
| Synteny | pair F1 +0.032 [+0.009, +0.062]; hard-case +0.16 [0.05, 0.26]; OG-F 71.3 → 75.7, shuffled 71.3–72.0 [C-S5-3, -6, C-CORE-12] | TODO | TODO | kept |
| Structure (3Di) | S0 0.794 → 0.853; no merge gain [C-S6-1, -3] | TODO | TODO | optional S0 rescue |
| Expression | hard-case −0.019 [−0.20, 0.15] [C-S9-2] | not run | not run | rejected |
| Isoforms | +0.013 [−0.018, 0.044] [C-S7-1] | not run | not run | rejected |
| Genome DNA | recovery 92.2%; 53% of losses are not losses [C-S8-1, -2] | n/a | TODO (fission-repair arm) | add-on |

TODO: add-one and leave-one-out against shuffled controls, with Holm correction.

## 7. Species-tree alternatives
- Measured: NCBI conflicts 0% with the published tree, data trees 12–21% [C-S1-1, -2]. The data-tree errors are long-branch artefacts, and bacteria conflict 25–30% [C-S1-3]. On OrthoBench all data trees share one error. The NCBI tree costs 6 s; data trees cost about 0.5–0.84 CPU-h [C-DOC-4].
- Decision: NCBI is the default; data trees are ablation arms [C-S1-4].
- TODO, with the full pipeline:
  - **Policy arms:** strict, delayed merge, and per-family local labelling. The S3 local-tree mode is a partial version of the last [C-CORE-9].
  - **Downstream effect of tree choice:** dev and QfO results for the NNI-perturbed and random-tree arms.
  - **Discordant families:** the count, split into prokaryote and eukaryote.
- TODO: discuss deep bacteria and the four eukaryote nodes where discordance should matter.

## 8. Follow-ups
1. Root-level (Bilateria) decisions have no outgroup; add structure evidence or a second outgroup.
2. Over-merging in large families (OG-P 63–75%): add an outgroup witness test.
3. Root rescue for the 145 Leiden-split families and 15k singletons [C-S0-7].
4. Species-tree policy arms.
5. Locus-driven fission repair as a candidate default [C-CORE-13].
6. 3Di gene trees and a reconciliation benchmark.
7. An orthologous-isoform output.
8. Structure at AFDB scale (GPU or precomputed 3Di).
9. A per-family comparison with FastOMA and OrthoFinder.

## 9. Limitations
- The dev set has 34 RefOGs, and OG-F partition noise is ±2–4 points [C-CORE-12, C-S5-8]. Dev-selected configs will shrink on test [C-S5-7].
- Competitor baselines show large dev-vs-test gaps (`bench_repo/results/orthobench_competitors/*/score.json`, fields `dev` and `test`), so a dev margin does not predict a test margin. TODO: quantify.
- No claim is independently reviewed (CLAIMS.md).
- Timings are single runs on a shared node.
- The S0 ceiling is about 0.82 and the pre-registered target was missed [C-S0-1].
- Deep decisions are score-based, and synteny gives nothing at deep levels (human vs fly AUC 0.49) [C-CORE-5, C-S5-2].
- Synteny coverage has gaps (zebrafish 0.64, amphioxus 0.08) [C-S5-1].
- The toy simulator is simplified (`docs/simulator.md`).
- The NCBI tree is an external input.
- Prokaryotic HGT, thousands of proteomes, structure at scale and GPU are not covered.

## 10. Remaining ideas
- Ancestral gene-order graphs shared across families (AGORA-style).
- Profile/HMM rescue for S0 members without pairwise hits.
- Calibrated distances and outgroup quartet tests for fast-evolving young duplicates.
- Structure as a deep-node tie-breaker only.
- Language-model embeddings for family partitioning.
- Epigenetic layers.
- Delayed merge for HGT-prone families.
- An independent-simulator check of multi-level HOGs.

## 11. Reproduction
TODO: one command per table, generated from `results/**/*.json`, plus the env lock (`envs/orth.lock.txt`), seeds, frozen config hash, and the `bench_repo/` runners. Until then, use the paths in `CLAIMS.md`.
