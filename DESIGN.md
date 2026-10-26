# DESIGN.md — v1 (PI, 2026-10-26 15:20 CET, T+0h40; after review round 1)

The project goal and its parts are fixed (PROMPT.md). Round-1 reviews are in docs/reviews/r1_*.md, and the change log is in §10. Time is always read with `scripts/clock.sh`, never estimated (the v0 draft mis-stamped T+1h15 when it was T+0h26).

## 1. Pitch
**Ancestra**: orthology as *level-synchronous ancestral genome reconstruction inside de novo protein families*.

1. **S0, family partition.** Near-linear clustering (DIAMOND DeepClust / MMseqs2 linclust cascade) splits all proteins into homologous families. There is no reference orthology database and no all-vs-all search.
2. **Within-family similarity.** We compute similarity only for pairs inside a family (Σ_f g_f², with large families split or subsampled). Ancestral genes never need to be searched again: their similarity is an aggregate of their members' leaf-pair scores.
3. **Level-synchronous merging.** We walk a guide species tree bottom-up. At each node X, *all families at once* merge their child ancestral genes into ancestral genes of X, and each one is a HOG at X.
   - The merge decision combines three things: (a) sequence evidence (aggregated scores, mutual-best plus in-paralog absorption, with fission repair); (b) a small local gene tree with species-overlap labelling when the cluster is ambiguous or many-to-many; (c) *ancestral gene-order context*. Ancestral genes inherit adjacencies from their members, and a candidate match gets support when its neighbours are matched to each other at X (anchor-based, three-pass, not circular).
   - Structure (3Di) evidence is used only for deep and ambiguous decisions.
4. **Output.** The HOG hierarchy goes to OrthoXML (primary output). A pair export derived from the same HOGs is a labelled secondary output.

Cost model (replaces the v0 "linear in S" claim): C = C_cluster(N) + Σ_f min(g_f², cap·g_f) + Σ_X Σ_f (merge work). This is near-linear in N when families are bounded. We report measured cost per stage, RSS, time per node and empirical log-log slopes, not asymptotic claims.

## 2. Novelty (honest; see r1_prior_work_domain)
Precedents:
- **SYNERGY** (Wapinski 2007): bottom-up along the species tree with sequence plus synteny, all-vs-all, fungi.
- **Hieranoid 1/2**: progressive InParanoid merges with representative profiles.
- **AGORA / edgeHOG / DeCoSTAR**: ancestral adjacencies *after* orthology is known.
- **OrthoFFGC**: a synteny-based QfO entry.
- **OF3 / SHOOT**: medoid representatives.
- **FastOMA**: k-mer placement into an existing database, then bottom-up processing.

Our claimed differences:
1. Reference-free near-linear family partitioning combined with bottom-up ancestral reconstruction at scale.
2. Ancestral gene-order context used *inside* the merge decision at every level, jointly with orthology, for all families synchronously.
3. A hybrid merge: score rules where unambiguous, local trees plus species overlap where ambiguous.
4. Structure (3Di) as a deep-node tie-breaker and rescue.
5. Explicit treatment of guide-tree uncertainty (polytomies, HGT, delayed merging).

We describe the method as "a scalable, synteny-aware descendant of SYNERGY and Hieranoid", not as unprecedented.

## 3. Sub-tasks and their benchmarks
| ID | Sub-task | Arms | Benchmark (dev split) | Metric |
|---|---|---|---|---|
| S0 | Family partition | DeepClust · linclust/cluster cascades · thresholds | OrthoBench dev RefOGs: fusion/fission, RefOG completeness; toy sim truth; runtime | % RefOGs fully within one family, fusion rate, CPU-h |
| S1 | Guide species tree | NCBI taxdump-derived · sketch distance + NJ · single-copy HOGs → ASTRAL-Pro · *perturbed (NNI/SPR)* · *random (negative control)* | nRF vs NCBI and published trees; downstream dev F | nRF; downstream |
| S2 | Within-family scoring | DIAMOND (mode) · MMseqs2 align on in-family pairs · SIMD SW (pyopal/parasail) | toy sim; OrthoBench dev pairs recall; CPU-h | recall, CPU-h |
| S3 | Merge rule at a node | RBH-only · mutual-best + in-paralogs · + fission repair · hybrid with local trees | toy sim (exact); OrthoBench dev rooted RefOG trees (speciation vs duplication labels from Data_for_RefOGs trees) | pair P/R/F1 |
| S3b | Multi-level HOG correctness | same | Zombi or ALF simulation with gene order, DTL, LGT (if installable within 1 h), else in-house simulator | per-level HOG Jaccard, pair F1 |
| S4 | Ancestral representation for scoring/rescue | max-link · avg-link · medoid · k-diverse rescue (capped at 1.2·n) | dev deep-level recall/fission | Δ |
| S5 | Synteny evidence | off · anchors-only tie-break · additive score · **shuffled-order control** | OrthoBench dev + toy sim; also reported on QfO test at freeze only | Δ(real − shuffled) |
| S6 | Structure evidence | off · 3Di (ProstT5/AFDB) for deep/ambiguous · **shuffled-3Di control** | dev deep fission/recall | Δ(real − shuffled) |
| S7 | Isoforms | canonical · best-isoform-pair scoring (QfO `_additional.fasta`, Ensembl isoforms) | dev pairs | ΔF1 |
| S8 | Gene-loss validation | demo: search for the missing HOG member in the syntenic locus of the genome (Ensembl mammals) | small set of known losses | precision of "lost" vs "not annotated" |
| S9 | Expression (exploratory, 2 h box) | Bgee correlation as paralog tie-breaker, vertebrates | dev ambiguous cases | negative if Δ < CI |
| F | Full pipeline | frozen config ± each layer | OrthoBench **test** RefOGs; QfO 2022 all tests (local service) | Pareto vs public methods; runtime |

## 4. Data, splits, leakage control (must)
- Prediction input is QfO 2022_02_with_updated_UP000000437 only (bare accessions, canonical set). `scripts/verify_ids.py` runs before every QfO run (≥99.5% coverage, no isoforms leaking).
- Synteny coordinates come from the 2022 `_gcoord.xml` first; 2026 coordinates are used only when the sequence MD5 matches; Ensembl GFF3 for OrthoBench uses the primary assembly only. Coverage is reported per species.
- **Split** (fixed seed 11, hashes committed in `splits/`):
  - Dev: toy simulations; 35 OrthoBench RefOGs; S0/S2 tuning on unlabeled data.
  - Test: the other 35 RefOGs; every QfO test (SwissTree, TreeFam-A, VGNC, GO, EC, FAS, STD, G-STD2).
- At most 3 logged QfO "peeks" before the freeze (`results/test_log.tsv`: config hash, time, reason).
- VGNC is not used as development evidence for synteny, because VGNC is partly synteny-derived.
- `tests/test_no_leak.py`: method code may read only declared inputs (proteomes, GFF/gcoord, guide tree, optional structure/expression). It fails on references to reference_data, RefOG, SwissTree, TreeFam, lineage_tree or species_tree_samples. The NCBI guide tree comes from the NCBI taxdump.

## 5. Evaluation protocol (pre-registered)
- **Primary comparators:** Hieranoid 2, OrthoFinder, FastOMA, InParanoid 5, SonicParanoid2-sens-g (public QfO rows). Their public predictions are re-scored locally and must match OpenEBench within 0.005, otherwise the local QfO pipeline is not trusted. FastOMA's advantage from the OMA database is flagged.
- **OrthoBench:** all measures, paired bootstrap over RefOGs (B=2000) for the CI of the difference vs comparators run locally (OrthoFinder; FastOMA and SonicParanoid if feasible).
- **QfO:** report each test as (recall proxy, precision proxy) with n. "Comparable" means no comparator dominates us on both axes beyond the 95% CI (bootstrap over families or pairs where the test allows; otherwise point estimates, labelled as such).
- **Layers:** add-one and leave-one-out ablations against shuffled-layer controls, with Holm correction across layers. Null effects are reported with an equivalence bound.
- **Runtime:** `/usr/bin/time -v`, 24 threads, node load logged; nested species subsets (8/16/32/64/78) with fixed seeds; log-log slope with CI. All preprocessing counts (ProstT5, DB builds). Censored runs are reported as "> T".
- **Primary output** is the HOG OrthoXML. The secondary pair export rule is fixed on dev; we never pick the better export per test.
- **Seeds** {11, 23, 47} for every stochastic step. A determinism test hashes the output under a fixed seed.

## 6. Verification
- `scripts/verify_orthoxml.py`: XSD validity, geneRef integrity, TaxRange = LCA of the species, our flattening = QfO `map_relations.py` flattening, pyHam loads, OrthoBench accepts the export.
- Toy simulator (in-house, already being built) with exact truth: positive control (truth → 1.0), negative controls (random sequences, shuffled gene order). Zombi or ALF as an independent simulator for S3b.
- Every reported number traces to `results/**.json` with a `manifest.json` (config hash, git commit, inputs, runtime). CLAIMS.md ledgers each claim with its evidence file. Report tables are generated from JSON.
- Each result is reviewed by an agent that did not produce it, and the reviewer reruns the verify scripts and recomputes at least one number.
- Spot-check 20 random HOGs at each milestone.

## 7. Plan, gates and kill criteria (clock from START.txt 14:40 CET 26 Oct)
| When | Milestone / gate |
|---|---|
| T+3h | Design frozen (≤3 review rounds). Data ready. OrthoBench and QfO local pipelines validated on public submissions |
| T+8h | v0 end-to-end (S0 → S2 → S3 → OrthoXML) on toy and OrthoBench dev. Competitor timings started (sequential, fixed threads) |
| **T+12h gate** | If OrthoBench-dev F < OrthoFinder −0.10, switch the core merge to the hybrid local-tree mode before any layer work |
| T+16h | First QfO full run (peek 1). S1 arms |
| **T+20h gate** | If v0 is dominated by Hieranoid 2 on both SwissTree and TreeFam-A by more than 2 SE, fix the core before layers |
| T+20–34h | Layers S5/S6/S7 (time-box 4 h each), S8 (3 h), S9 (2 h). Each is killed if Δ(real − shuffled) ≤ 0 on dev at its time-box |
| T+24h | REPORT.md skeleton |
| T+36h | Code freeze. T+38h final runs started. T+44h independent review of results. T+47h REPORT.md final |

Checkpointing: atomic, hash-keyed per-stage outputs (families, scores, per-node merge state). Nodes with more than 200k queries are chunked. RESUME.md explains how to restart.

## 8. Engineering choices for v0
Python 3.11 with numpy, numba and pyarrow. Integer gene IDs and a species ID mapping table with tests. Per-stage parquet/npy files. lxml streaming OrthoXML writer. A single search-backend interface (DIAMOND, MMseqs2, Foldseek; GPU-ready later). Polytomies are merged k-way inside families, and HOG levels are named only after real taxonomy clades (synthetic binarisation nodes are dissolved on output).

## 9. Risks
- Recall loss at deep nodes, the known weakness of Hieranoid and SYNERGY → fission repair, root-level rescue search among root HOGs and singletons, 3Di rescue.
- Over-merging at deep levels lowers QfO precision → conservative deep merges, and bidirectional ≥60% coverage guard for multi-domain proteins and fusions (no transitive chaining through non-overlapping hit regions).
- Synteny circularity → anchor-based three passes; abstain for tandem duplicates; never prune co-orthologs.
- HGT in prokaryotes → "unmatched passes up" plus delayed merge; report STD-Bacteria and G_STD2_Luca separately.
- QfO local run without Docker fails → fallback to reimplementing the pairwise tests, validated against public submissions.

## 10. Change log
- v1, after round 1:
  - Added S0 family partition, which fixes the cost blow-up of v0's pan-proteome searches (r1_algorithms_compute: v0 saving only ~2.5× over all-vs-all).
  - Replaced the "linear" claim with a cost model.
  - Added a hybrid merge (r1_external_pi), fission repair and root rescue (r1_prior_work).
  - Added a split and leakage tests; removed VGNC as dev evidence (stats, ai_failure).
  - Fixed S3 truth to come from RefOG trees.
  - Pre-registered comparators and a success rule.
  - Added shuffled controls for each layer.
  - Specified the runtime protocol.
  - Added gates and kill criteria.
  - Rewrote the novelty section (SYNERGY, AGORA, OrthoFFGC).
  - k-way polytomy merge with taxonomy-only levels.
