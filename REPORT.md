# Ancestra: orthology inference as level-synchronous ancestral-genome reconstruction inside de novo protein families, with an evaluation of extra data layers

*AI-generated report (Claude Opus 5.5 as principal investigator, delegating to Opus/Sonnet subagents). It is written in the style of a journal article but is not a peer-reviewed publication; see the human-written header in `README.md`. Every quantitative statement maps to an evidence file through `CLAIMS.md` (claim IDs in brackets) and to `docs/VERIFICATION.md`. All protocol deviations are listed in `docs/DISCLOSURES.md`.*

**Status of this draft:** sections 3.6 (large-scale runtime) and 3.8 (GPU) are filled when the corresponding runs finish. All other numbers are final, frozen-method (`freeze-v1`) results.

---

## Abstract

**Background.** State-of-the-art orthology methods infer gene trees per orthogroup (OrthoFinder), build hierarchical orthologous groups (HOGs) from placements into a reference database (FastOMA), or combine pairwise best hits (SonicParanoid2, InParanoid, Hieranoid). We asked whether a conceptually different formulation — reconstructing ancestral gene sets level by level along the species tree, entirely inside de novo protein families — can match their accuracy at lower cost. We also asked whether data layers beyond protein sequence (gene order, 3D structure, expression, isoforms, genomic DNA) improve orthology inference.

**Method.** *Ancestra* (i) partitions all proteins into families in near-linear time with a two-level DIAMOND clustering (no all-vs-all, no reference database), (ii) scores protein pairs only inside families, and (iii) walks the NCBI taxonomy guide tree bottom-up. At each ancestral node, the ancestral genes of all families are merged with mutual-best seeds, in-paralog absorption and small local gene trees (FastTree or NJ) labelled by species overlap. Gene order can re-rank partners inside an existing root HOG. Output: HOGs in OrthoXML, ortholog pairs and orthogroups.

**Results.** On the held-out OrthoBench test RefOGs (36, split by homology component, scored once), the default preset reaches F = 76.2 vs 69.7 for OrthoFinder 3.2, 70.4 for SonicParanoid2, 60.4 for precomputed OMA browser HOGs and 33.1 for FastOMA 0.5.1. It uses 3.9 CPU-h vs 8.2 for OrthoFinder. The differences against OMA (+15.7 F [+1.8, +29.2]) and FastOMA (+43.1 [+14.0, +55.3]) are significant. The pre-registered non-inferiority criterion against OrthoFinder is not met (Δ +6.5 [−6.5, +18.1]): the confidence interval is too wide for 36 RefOGs. On the Quest for Orthologs 2022 benchmark (78 species, ~1 M proteins; full official suite run locally), Ancestra is non-inferior to FastOMA and OMA Pairs on TreeFam-A recall, but clearly below OrthoFinder on recall (TreeFam-A −0.16, SwissTree −0.10). Its species-tree discordance (G-STD2 Eukaryota) is the worst of the main comparators. CPU time grows as proteins^1.16 [1.10, 1.21], against ~1.4–1.5 for OrthoFinder and SonicParanoid2 (3-point fits). Of the extra layers, gene order separates orthologs from paralogs very well on its own (human–mouse AUC 0.92) but gives no test-measurable gain inside the merge. 3D structure (3Di) raises family completeness beyond the sequence ceiling but does not improve final orthologs. Expression and isoforms were rejected by pre-registered rules. Genomic DNA shows that 53% of apparent gene losses are not real losses. The NCBI taxonomy is the best guide tree; a random tree costs 16 F points.

**Conclusion.** Level-synchronous ancestral reconstruction inside de novo families is a viable, cheaper and more scalable route to HOGs. It is competitive with or better than OMA-based resources and pairwise methods on OrthoBench, but it does not yet match OrthoFinder's recall at deep taxonomic levels.

---

## 1. Introduction

Orthology inference underpins function transfer, comparative genomics and phylogenomics. Current methods fall into three families:
1. **Graph/pairwise methods** compute (reciprocal) best hits between every pair of proteomes and cluster them: InParanoid, OrthoMCL, SonicParanoid2, ProteinOrtho, OMA pairs. Their all-vs-all search costs O(n²) in the number of proteomes.
2. **Tree-based methods** build a gene tree per orthogroup and reconcile it with a species tree: OrthoFinder, PhylomeDB, Ensembl Compara, PANTHER.
3. **Placement methods** map new proteomes into a precomputed reference: FastOMA via OMAmer, eggNOG-mapper. They are fast but depend on the reference database.

Hierarchical orthologous groups (HOGs) — sets of genes descended from one ancestral gene at a given taxonomic level — are the natural output of all three. Most methods derive HOGs post hoc from pairs or trees (Altenhoff et al. 2013; Majidian et al. 2025).

We explored a fourth formulation: **build the HOGs directly, by reconstructing the ancestral gene complement at each internal node of the species tree, one level at a time**. This is how a phylogeneticist reasons about gene families (what did the ancestor of mammals have?). It also has two computational advantages:
- every decision at a node is local to a family and a level;
- extra evidence such as ancestral gene order can be attached to ancestral genes as they are reconstructed.

The brief (`PROMPT.md`) required:
- a conceptually new method;
- scalability to thousands of proteomes;
- OrthoXML HOG output;
- per-sub-task benchmarks;
- exploration of alternatives to a universal species tree;
- tests of extra data layers;
- benchmarking on QfO, OrthoBench and SwissTree against competitors, with at least one OMA resource in every comparison.

---

## 2. Approach and research path

### 2.1 Design history
The first design (v0) merged complete ancestral pan-proteomes progressively. The design review estimated only ~2.5× saving over all-vs-all on QfO, because ancestral pan-proteomes stay large. Design v1/v2 therefore added a **near-linear family partition first** and restricted all subsequent work to within-family comparisons [C-INF-7] (`DESIGN.md`; reviews in `docs/reviews/`).

The protocol was pre-registered before any results:
- dev/test splits: OrthoBench split by homology component; a TreeFam-A 30% dev split was added later (disclosed);
- at most three QfO peeks before the freeze;
- a non-inferiority success rule: the lower 95% bound of Δ ≥ −0.03;
- shuffled controls for every data layer.

The research record with every decision is in `RESEARCH_TREE.md`; process lessons are in `LESSONS.md`.

### 2.2 Final method (frozen as `freeze-v1`)
![Ancestra method overview](docs/figures/fig1_method.png)

*Figure 1. Ancestra pipeline: S0 family partition → S2 in-family scoring → S3 level-synchronous merge along the guide tree, with optional gene-order, 3Di and gene-loss layers.*

**S0 — Families** (`ancestra/families.py`).
- DIAMOND cluster at 30% identity.
- All genes are searched against the cluster representatives, and connected components are taken over hits with ≥30% coverage.
- Components above 7,500 genes are split with Leiden.
- The partition is deterministic: canonical input order, pinned DIAMOND block/chunk parameters and sorted edges. A non-deterministic first version (ARI 0.78 across thread counts) was found and fixed [C-CORE-11].

**S2 — In-family scoring** (`ancestra/scoring.py`). DIAMOND blocks inside each family, plus MMseqs2 re-alignment of pairs missed by DIAMOND. The hybrid recovers 6–10% of missed pairs [C-CORE-2].

**S3 — Level-synchronous merge** (`ancestra/merge.py`, `localtree.py`). The NCBI taxdump tree restricted to the input species is traversed bottom-up. At each internal node, the ancestral genes of the child clades in every family are merged k-way:
1. mutual-best seeds;
2. absorption of in-paralogs;
3. a small local gene tree for ambiguous clusters (FastTree on a FAMSA alignment of cluster medoids in the `accurate` preset, NJ in `fast`), labelled by species overlap;
4. at deep nodes (median seed −ln S ≥ 1.5), species-overlap absorption to limit splitting at the LCA.

**Gene order** (`ancestra/synteny.py`, `--gene-order`). Ancestral adjacency is propagated with the ancestral genes. Gene order may only *re-rank* partners within the HOG that the synteny-off pass placed them in (constrained mode; window 10, λ 0.3), and it can veto merges. This constraint was added after an unconstrained additive bonus cost 4–8 recall points [C-S5-4].

**Outputs.** OrthoXML HOGs, validated by `scripts/verify_orthoxml.py` (XSD, gene uniqueness, LCA consistency, pyham); ortholog pairs; OrthoBench orthogroups; and a manifest with the code hash and seeds.

**Presets.**
- `accurate` (default): FastTree local trees on the gene-order pass, NJ otherwise, plus deep-node absorption.
- `fast`: NJ plus species-overlap structure.

### 2.3 Sub-task benchmarks (dev data only)
Each stage had its own benchmark before end-to-end evaluation (details in `docs/S*_*.md`):

| Stage | Benchmark | Key dev result |
|---|---|---|
| S0 families | Completeness of OrthoBench dev RefOGs in one family | 0.794 (ceiling of any sequence-similarity partition ≈ 0.82; the pre-registered 95% target is refuted) [C-S0-1, C-S0-2] |
| S1 guide tree | Conflict with the published QfO species tree | NCBI 0.000; ASTRAL-Pro 0.120; FastME 0.173; BIONJ 0.213 [C-S1-1, -2] |
| S2/S3 core | Toy simulator truth; OrthoBench dev OG-F | Toy pair F1 0.94–0.98 (RBH-only 0.74); dev OG-F 79.7 (`accurate`) vs OrthoFinder 70.5 [C-CORE-1, C-CORE-12] |
| S4 local trees | Duplication/speciation label accuracy on dev families | FastTree > NJ at deep nodes → `accurate` preset (`docs/S4_trees.md`) |
| S5 gene order | Ortholog-vs-paralog AUC; hard-case gain | AUC 0.92 human–mouse, 0.49 human–fly, shuffled 0.50 [C-S5-2] |
| S6 structure | S0 completeness; merge pair F1 | 0.794 → 0.853 completeness; merge +0.015 [−0.004, 0.035] [C-S6-1, -3] |
| S7 isoforms | Hard-case gain | +0.013 [−0.018, 0.044] → rejected [C-S7-1] |
| S8 genomic DNA | Recovery of planted controls; loss-call precision | 92.2% recovery, 0/512 false hits [C-S8-2] |
| S9 expression | Hard-case gain | −0.019 [−0.20, 0.15] → rejected [C-S9-2] |

### 2.4 Main course corrections
These are recorded in `RESEARCH_TREE.md`.
1. Pure progressive merging was replaced by family-first design (M1).
2. A leaky random OrthoBench split was replaced by a homology-component split before any tuning.
3. The non-deterministic S0 was fixed. All merge arms were then replicated on the deterministic partition, and some earlier gains turned out to be partition noise (±2–4 OG-F) [C-CORE-12].
4. QfO peek 1 showed low recall. A dev-only diagnosis on a TreeFam-A 30% split found that 78% of misses were merge failures at the pair's LCA. This led to the FastTree local trees and deep-node absorption (M3).
5. A malformed-nesting bug in the species-overlap structure option was caught by the OrthoXML validity gate before the final run.

---

## 3. Results

All numbers in this section come from the frozen code (`freeze-v1`, sha256 of every module in `results/freeze.json`). Test sets were scored once, after the freeze.

### 3.1 OrthoBench (held-out test RefOGs)

| Method | Test F | P | R | Exact | CPU-h | Peak RSS (GB) | ΔF (Ancestra − method) [95% CI] |
|---|---|---|---|---|---|---|---|
| **Ancestra `accurate`** | **76.2** | 69.5 | 84.2 | 10/36 | 3.9 | 14.0 | — |
| Ancestra `fast` | 74.8 | 68.7 | 82.1 | 9/36 | 2.9 | 14.0 | +1.4 [−3.6, +7.1] |
| OrthoFinder 3.2 | 69.7 | 62.1 | 79.5 | 8/36 | 8.2 | 5.6 | +6.5 [−6.5, +18.1] |
| SonicParanoid2 | 70.4 | 70.8 | 70.1 | 8/36 | 7.9 | 8.1 | +5.7 [−4.3, +14.9] |
| OMA browser HOGs (Bilateria) ⚑ | 60.4 | 88.9 | 45.8 | 3/36 | n/a (precomputed) | — | **+15.7 [+1.8, +29.2]** |
| ProteinOrtho | 47.5 | 97.7 | 31.4 | 4/36 | 1.6 | 0.8 | **+28.7 [+13.1, +41.9]** |
| FastOMA 0.5.1, patched ⚑ | 33.1 | 23.9 | 53.8 | 10/36 | 19.9 | 12.2 | **+43.1 [+14.0, +55.3]** |

*Table 1. OrthoBench test (36 RefOGs). F, P and R are the official Open Orthobench metrics. ΔF is a paired bootstrap over RefOGs (B = 2000, seed 11), in percentage points. ⚑ OMA contains all 12 OrthoBench species, so OMA-based results are informed rather than independent comparators. Runtimes come from a shared node (8–10 threads, load 20–55) and are indicative. Sources: `results/orthobench_final/summary.md`; bootstrap re-run independently [C-FIN-*, C-OMA-*].*

Under the pre-registered rule (comparable if the lower 95% bound of Δ ≥ −3 pp):
- Ancestra is **comparable to or better than** ProteinOrtho, the OMA browser HOGs and FastOMA, with CIs entirely above 0.
- Non-inferiority to OrthoFinder and SonicParanoid2 is **not shown**. The point estimates favour Ancestra by 5–7 F, but with 36 RefOGs the standard error of ΔF is ~5–6 pp.
- On all 70 RefOGs, which include the dev half used for tuning (optimistic), the difference vs OrthoFinder is +7.8 [+0.5, +15.1].

**FastOMA.** FastOMA 0.5.1 performed poorly on OrthoBench. While diagnosing this we found a bug that drops ~5.5 k genes when groups are split. A patched re-run raised gene coverage from 79.0% to 80.6% but did not change test F (33.5 → 33.1). The gap is therefore methodological: many proteins are not placed by OMAmer, and RefOGs are fused at the root. This is consistent with the large FastOMA fusion rate on OrthoBench reported by Emms et al. (2026). The precomputed OMA browser HOGs, cut at Bilateria, are a fairer OMA representative (60.4): precise (P 88.9) but low-recall (R 45.8).

### 3.2 Quest for Orthologs 2022 (78 species, 988,778 proteins)
The full QfO benchmark suite (SwissTree, TreeFam-A, VGNC, GO/EC function, FAS, STD and G-STD2) was run locally from the official code. Our local scoring reproduces 60/60 published values of competitor submissions [C-INF-2].

| Method | SwissTree TPR / PPV | TreeFam-A TPR / PPV | EC Schlicker (n) | G-STD2 Euk RF (n trees) |
|---|---|---|---|---|
| **Ancestra `accurate`** | 0.721 / 0.926 | 0.583 / 0.921 | 0.899 (181k) | 0.257 (3,379) |
| Ancestra `fast` | 0.731 / 0.923 | 0.608 / 0.916 | 0.896 (197k) | 0.266 (4,753) |
| OrthoFinder | 0.825 / 0.917 | 0.737 / 0.918 | 0.931 (204k) | 0.249 (6,629) |
| SonicParanoid2-sens-g | 0.717 / 0.938 | 0.667 / 0.936 | 0.906 (222k) | 0.232 (4,878) |
| Hieranoid 2 | 0.726 / 0.950 | 0.621 / 0.930 | 0.948 (169k) | 0.222 (3,643) |
| FastOMA v0.1.4 | 0.699 / 0.953 | 0.546 / 0.919 | 0.928 (157k) | 0.239 (3,092) |
| OMA HOGs | 0.696 / 0.924 | 0.580 / 0.923 | 0.829 (308k) | 0.237 (2,981) |
| OMA Pairs | 0.677 / 0.929 | 0.565 / 0.920 | 0.952 (156k) | 0.226 (2,369) |

*Table 2. QfO 2022 primary endpoints for the main comparators. The full 21-method table, VGNC, GO, FAS and STD results are in `results/qfo_final/summary.md`. Competitor values are published OpenEBench submissions.*

The formal non-inferiority analysis used a paired bootstrap over the 749 TreeFam-A test families and the 18 SwissTree cases [`results/verification/qfo_primary_ni.json`]. `accurate` on TreeFam-A test recall:

| Comparator | Δ TPR [95% CI] | Verdict |
|---|---|---|
| FastOMA | +0.043 [−0.005, +0.093] | non-inferior |
| OMA Pairs | +0.018 [−0.028, +0.069] | non-inferior |
| OMA HOGs | +0.002 [−0.058, +0.067] | inconclusive |
| OrthoFinder | −0.160 [−0.210, −0.106] | fails |

On SwissTree recall it fails against OrthoFinder (−0.103 [−0.177, −0.029]) and is inconclusive against the rest.

Ancestra therefore sits with the OMA family of methods on QfO, precise but recall-limited, while OrthoFinder and the tree-based methods (PANTHER-all, pbDOOD) reach higher recall. G-STD2 Eukaryota is the worst of the main comparators. EC precision is also slightly lower, consistent with — but not proven to be caused by — over-merging in large vertebrate families.

**Contrast with OrthoBench.** Ancestra leads on OrthoBench but trails OrthoFinder on QfO. The two benchmarks reward different things:
- OrthoBench scores whole orthogroups in 12 mostly bilaterian species, where Ancestra's family-first partition and in-family merge do well.
- QfO's TreeFam-A and SwissTree reward pairwise recall across deep, eukaryote-wide divergences. There Ancestra's score-based merge decisions at deep nodes (Bilateria, Opisthokonta, Eukaryota), which have no outgroup, split true orthologs at the LCA.

### 3.3 Effect of extra data layers
Every layer was tested at two levels:
1. standalone, on its own sub-task benchmark with a shuffled control;
2. inside the full pipeline, on held-out test data.

| Layer | Standalone (dev) | OrthoBench test (Δ F vs main) | QfO (Δ TreeFam-A TPR, family bootstrap) | Decision |
|---|---|---|---|---|
| Gene order (synteny) | AUC 0.92 human–mouse, 0.82 fish, 0.49 human–fly (shuffled 0.50); hard cases +0.16 [0.05, 0.26] | −0.1 [−0.4, 0.0] (removing it) | −0.012 [−0.028, +0.001] (removing it) | option, default on; no measurable test gain |
| 3D structure (3Di, AFDB/Foldcomp, 94.5% QfO coverage) | S0 completeness 0.794 → 0.853 | −5.0 [−13.1, 0.0] as family rescue | −0.004 | option, off |
| Expression (Bgee) | AUC 0.667; hard-case −0.019 [−0.20, 0.15] | not run | not run | rejected (pre-registered) |
| Isoforms | hard-case +0.013 [−0.018, 0.044]; 2–5× search cost | not run | not run | rejected |
| Genomic DNA (miniprot at syntenic loci) | 92.2% control recovery; 0/512 false hits | n/a (validation layer) | n/a | add-on |

*Table 3. Data-layer effects. Sources: `docs/S5_synteny.md`, `docs/S6_structure.md`, `docs/S7_isoforms.md`, `docs/S8_geneloss.md`, `docs/S9_expression.md`, `results/orthobench_final/summary.md` §2, `results/verification/qfo_ablation_familyboot.txt`.*

**Gene order.** Synteny strongly separates orthologs from paralogs between closely related vertebrates, and its standalone effect is real (shuffled controls are flat). Inside the merge, however, sequence evidence already resolves almost every case that synteny can resolve. The dev gains (+1.0 to +4.4 F) did not replicate on test. Synteny also carries no signal at the deep levels where Ancestra loses recall (human–fly AUC 0.49).

**Structure.** 3Di profile search re-attaches fragments that sequence similarity misses, raising S0 completeness above the sequence ceiling. In the final pipeline, though, the extra genes lowered OrthoBench precision (69.5 → 60.1). Within-RefOG pair F1 was unchanged, so the rescued genes were mostly non-members. Structural similarity at the 3Di level is a homology signal, not an orthology signal.

**Genomic DNA.** Of 2,523 apparent gene losses in mammals:
- 53% were not losses;
- 32% were our own fission errors (a gene split across HOGs);
- 25% were confirmed pseudogenes or deletions, e.g. *GULO* and *UOX* in primates;
- the rest were unresolved [C-S8-1, -3].

This makes genome-level validation a useful quality-control layer for HOG gene-loss claims, rather than an input to inference.

### 3.4 Alternatives to a universal species tree
Four arms were tested:
1. the NCBI taxonomy (default);
2. data-derived trees (ASTRAL-Pro from gene trees; FastME/BIONJ from distances), rooted on NCBI;
3. NNI-perturbed NCBI trees;
4. a random tree as a negative control.

| Guide tree | Conflict vs published QfO tree | OrthoBench test ΔF | QfO TreeFam-A ΔTPR (fast) | QfO EC Δ |
|---|---|---|---|---|
| NCBI (default) | 0.000 | — | — | — |
| ASTRAL-Pro | 0.120 | +0.0 [−0.3, +0.4] | +0.028 [−0.009, +0.077] | −0.125 |
| FastME | 0.173 | +0.0 [−0.3, +0.4] | +0.022 [−0.001, +0.045] | −0.134 |
| NCBI + 1 NNI (3 seeds) | — | −3.3 to +0.3 | — | — |
| NCBI + 5 NNI | — | — | +0.003 [−0.020, +0.026] | +0.003 |
| Random (negative control) | 0.987 | **−16.4 [−27.1, −5.6]** | **−0.243** | −0.060 |

*Table 4. Guide-tree arms. Sources: `docs/S1_speciestree.md`, `results/orthobench_final/summary.md`, `results/qfo_final/ablation.tsv`.*

The guide tree matters: a random tree collapses performance. Among reasonable trees, the curated NCBI taxonomy is best overall:
- Data-derived trees err on long branches, especially among bacteria (25–30% conflict).
- They buy a small, non-significant QfO recall gain at a large cost in EC function consistency (−0.12 to −0.13), probably by merging across wrongly placed prokaryote clades.
- A per-family local labelling — the local-tree step in S3 — is the partial alternative to a single universal tree that we retained.
- ASTRAL-Pro and FastME have the same topology but differ in Newick child order, and their outputs differ in 18 orthogroups. This reveals a small child-order sensitivity of the merge (DISCLOSURES #11).

### 3.5 Runtime and scaling (≤ 78 species)
![Benchmarks and scaling](docs/figures/fig2_benchmarks.png)

*Figure 2. OrthoBench test F with bootstrap CIs, QfO TreeFam-A precision–recall, and CPU-h scaling on nested QfO subsets.*

On nested, domain-stratified QfO subsets (8 threads per run):

| Method | 8 sp (62k) | 16 sp (144k) | 32 sp (355k) | 64 sp (816k) | 78 sp (989k) | Slope vs proteins |
|---|---|---|---|---|---|---|
| Ancestra `accurate` | 0.78 | 2.16 | 6.02 | 16.30 | 18.54 | 1.16 [1.10, 1.21] |
| Ancestra `fast` | 0.47 | 0.99 | 3.20 | 10.63 | 11.63 | 1.21 [1.02, 1.41] |
| OrthoFinder 3.2 | 1.22 | 3.54 | 14.43 | *[T2 pending]* | — | 1.42 (3 pts) |
| SonicParanoid2 | 0.77 | 2.34 | 11.11 | *[T2 pending]* | — | 1.53 (3 pts) |
| FastOMA 0.5.1 | 1.83 | *[T2 pending]* | *[T2 pending]* | *[T2 pending]* | — | — |

*Table 5. CPU-h per run (`results/scaling/summary.md`). Peak RSS for Ancestra is 2–29 GB. FastOMA excludes its precomputed OMAmer database.*

**Where Ancestra spends its time** (QfO): S2 in-family scoring 8.6 CPU-h > `accurate` merge 7.0 CPU-h (61k small FAMSA+FastTree jobs) > S0 3.0 CPU-h. The `fast` merge takes 2.5 CPU-min.

### 3.6 Large-scale demonstration (125–1,000 proteomes)
*[Pending T3: UniProt reference proteomes, 0.5–4.75 M proteins, `fast` preset; `results/scaling_large/`.]*

### 3.7 Output validity and reproducibility
- Every final OrthoXML passes `scripts/verify_orthoxml.py`: XSD, gene-reference integrity, each gene exactly once, LCA consistency, pair counts, flattening equivalence and pyham parsing.
- A fresh end-to-end run reproduces the pre-freeze S0 and S2 checkpoints byte-for-byte.
- 156 unit and integration tests pass, including a leak test that keeps reference orthogroups out of `ancestra/`.

### 3.8 GPU acceleration (post-freeze, exploratory)
*[Pending GPU v2 track on branch `v2-gpu`: MMseqs2-GPU, batched GPU alignment, GPU trees, ProstT5 3Di, ESM-2 + FAISS candidate generation; dev-only evaluation; `docs/GPU_v2.md`.]*

---

## 4. Discussion

**What is new.** Ancestra treats HOG inference as ancestral-genome reconstruction rather than as clustering of pairwise relations or decomposition of large gene trees. Each decision is small: one family, one ancestral node, a handful of ancestral genes. Four consequences follow:
1. the cost is dominated by within-family alignment, which scales close to linearly in the number of proteins;
2. the expensive phylogenetic step is reserved for ambiguous clusters;
3. ancestral gene order can be carried along with the ancestral genes;
4. HOGs at every level are produced directly, not extracted after the fact.

It needs no reference database, unlike placement methods.

**Where it falls short.** Recall at deep levels is the main weakness. At the root of the input (e.g. Bilateria for OrthoBench, Eukaryota for QfO) there is no outgroup to decide whether two ancestral genes are orthologous or derived from an older duplication. Score thresholds there either split true orthologs (QfO recall) or over-merge (large vertebrate families, G-STD2 discordance). OrthoFinder's whole-orthogroup gene trees, rooted by duplication parsimony (STRIDE), handle this better. The most promising remedies are:
- outgroup "witness" tests using a third clade;
- structure as a deep-node tie-breaker;
- a profile-based S0 that is not capped by the sequence-similarity ceiling.

**Extra data layers.** Our layer results are largely negative for orthology *assignment*, despite strong standalone signals. This is an important and somewhat counter-intuitive finding: a layer can be highly informative in isolation yet redundant given good sequence evidence. The cases where the extra signal is unique — distant homologs (structure) and deep levels (synteny) — are exactly where that signal is weakest or blurs orthology and paralogy. Genomic DNA was most useful as validation: it exposed our own fission errors and confirmed known pseudogenisations.

**Benchmarks disagree.** The 36-RefOG OrthoBench test gives ΔF CIs of ±12 pp, so a 6.5-point lead over OrthoFinder is not significant. QfO's pair-level tests are much tighter but measure different things. We report both, and every pre-registered verdict as written, including the failures.

## 5. Limitations
- Dev data were reused heavily. A dev split for TreeFam-A was created after the first QfO peek, and test TreeFam-A families were scored once at peek 2, one minute before the freeze (DISCLOSURES #1). Dev numbers are therefore optimistic. Test numbers come from one evaluation of a frozen method.
- Runtimes were measured once per point on a shared node (load 20–70), so walls are upper bounds. Competitor fits use 3 points.
- OMA-based comparators saw the benchmark species (leakage in their favour). QfO competitor rows are published submissions produced under other conditions.
- The OMA browser labelling rule was corrected after a first score (45.5 → 60.4), in OMA's favour (DISCLOSURES #19).
- The merge is slightly sensitive to Newick child order (18/thousands of orthogroups).
- Prokaryotic HGT is handled only implicitly. The NCBI tree is an external input.
- 57 of 81 claims in `CLAIMS.md` were matched to their evidence files but not independently recomputed.

## 6. Methods (summary)
- **Software:** DIAMOND 2.2.8, MMseqs2 18-8cc5c, FAMSA, FastTree, Foldseek, miniprot, igraph/leidenalg, Python 3.11. Full lock in `envs/orth.lock.txt` and `ENVIRONMENT.md`.
- **Data:** OrthoBench (Open Orthobench), QfO 2022_02 reference proteomes and benchmark service, Ensembl 99 annotations for gene order, AFDB/Foldcomp 3Di, Bgee, UniProt reference proteomes for the scaling demo. Versions, checksums and licences are in `DATA.md`.
- **Splits:** OrthoBench by homology component (dev 34 / test 36, `splits/orthobench_split.json`); TreeFam-A 30% dev (`splits/treefamA_split.json`). Hashes are in `splits/SHA256SUMS`.
- **Statistics:** paired bootstrap over RefOGs or families, B = 2000, seed 11. Non-inferiority margin −0.03, as in `DESIGN.md` §5.
- **Verification:** an agent that did not produce a result recomputed it (`docs/VERIFICATION.md`). Discrepancies were added to `docs/DISCLOSURES.md` and are not hidden.

## 7. Data and code availability
- Method code, research record and results are in this repository (tag `freeze-v1`).
- Benchmark runners and competitor outputs are in the separate repository `bench_repo/`.
- Restricted or large data are not committed; `DATA.md` describes how to obtain them.

## 8. Follow-ups
1. Outgroup witness tests for deep-node merge decisions (QfO recall, G-STD2).
2. Outgroup-aware splitting of giant vertebrate HOGs (EC/VGNC precision).
3. A profile-based S0 rescue that keeps precision (the guarded 3Di rescue lost precision on test).
4. Removing the Newick child-order sensitivity.
5. A clean, unseen benchmark for v2 (QfO and OrthoBench test are now consumed).

## References
Key references (full list with verification status in `LITERATURE.md`):
- Altenhoff et al. 2013 *PLoS ONE* 8:e53786 (GETHOGs)
- Altenhoff et al. 2016 *Nat Methods* 13:425 (QfO benchmarking)
- Altenhoff et al. 2024 *NAR Genom Bioinform* 6:lqae167 (QfO 2022)
- Emms & Kelly 2019 *Genome Biol* 20:238 (OrthoFinder 2)
- Emms et al. 2026 *Nat Methods* doi:10.1038/s41592-026-03126-6 (OrthoFinder 3)
- Majidian et al. 2025 *Nat Methods* 22:269 (FastOMA)
- Cosentino et al. 2024 *Genome Biol* 25:195 (SonicParanoid2)
- Kaduk & Sonnhammer 2017 *Bioinformatics* 33:1154 (Hieranoid 2)
- Buchfink et al. 2021 *Nat Methods* 18:366 (DIAMOND)
- Steinegger & Söding 2017 *Nat Biotechnol* 35:1026 (MMseqs2)
- Kallenborn et al. 2025 *Nat Methods* (MMseqs2-GPU)
- van Kempen et al. 2024 *Nat Biotechnol* 42:243 (Foldseek)
- Heinzinger et al. 2024 *NAR Genom Bioinform* 6:lqae150 (ProstT5)
- Price et al. 2010 *PLoS ONE* 5:e9490 (FastTree 2)
- Zhang et al. 2020 *MBE* 37:3292 (ASTRAL-Pro)
- Huerta-Cepas et al. 2007 *Genome Biol* 8:R109 (species overlap)
- Bastian et al. 2021 *NAR* 49:D831 (Bgee)
- Trachana et al. 2011 *BioEssays* and Emms & Kelly 2020 (OrthoBench / Open Orthobench; citation details not re-verified)
