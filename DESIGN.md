# DESIGN.md — draft v0 (PI, 2026-10-26 T+1h15)

Status: DRAFT for review round 1. The project goal and its parts are fixed (PROMPT.md). Everything below is a proposal.

## 1. One-paragraph pitch
**Ancestra**: orthology inference by *progressive reconstruction of ancestral proteomes*. Instead of comparing every proteome with every other (O(S²·G²) for S species of G genes), we walk a guide species tree bottom-up. At each internal node we merge the two child (ancestral) proteomes into the parent ancestral proteome by matching their *ancestral genes*. Each ancestral gene is a HOG at that taxonomic level, by construction. The matching score combines sequence similarity of multiple ancestral representatives with optional evidence: conserved gene neighbourhood (synteny, from GFF), structural similarity (3Di), and isoform-aware representatives. The cost is S−1 merges of bounded size, so it is linear in S. Independent subtrees merge in parallel. No pre-built orthology database is used (unlike FastOMA/OMAmer, PANTHER, eggNOG), and no all-vs-all search (unlike OrthoFinder ≤2, OMA, SonicParanoid). The output is OrthoXML HOGs plus a pair export.

## 2. Positioning and novelty (honest)
- The closest precedent is Hieranoid 1/2 (Kaduk & Sonnhammer 2013/2017): progressive InParanoid merges with profile or consensus representatives. Our differences:
  1. Gene order is a *first-class* matching signal. Each ancestral gene carries ancestral adjacencies (reconstructed from children). Merging becomes a context-aware assignment, close to ancestral genome reconstruction (AGORA-like) done *jointly* with orthology, not as post-hoc QC.
  2. Multiple diverse representatives per ancestral gene (k-medoids) instead of one profile, plus optional 3Di/AA joint search for deep nodes.
  3. Explicit duplication handling at merges: in-paralog co-orthology, and "unmatched passes up" so that guide-tree errors and HGT delay a merge rather than corrupt it.
  4. Per-node guide tree from proteome sketches when no trusted tree exists. The species tree is benchmarked separately against alternatives.
- What is not novel: DIAMOND/MMseqs2 as the search engine, and the InParanoid-style in-paralog rule. We say so.

## 3. Sub-tasks, each with its own benchmark (pre-registered)
| ID | Sub-task | Variants to compare | Benchmark / QC | Primary metric |
|---|---|---|---|---|
| S1 | Guide species tree | NCBI taxonomy (as provided by QfO) · protein k-mer/MinHash distance + NJ · first-pass single-copy HOGs → ASTRAL-Pro/concat FastTree | QfO 78 and OrthoBench 12 vs reference (NCBI and published phylogenies); then downstream full-pipeline effect | nRF; downstream OrthoBench F |
| S2 | Candidate homolog search per merge | DIAMOND fast/sensitive/ultra · MMseqs2 · Foldseek AA+3Di | recall of reference ortholog pairs (OrthoBench RefOG-derived, SwissTree/TreeFam-A) among hits; runtime | recall@hits vs CPU-h |
| S3 | Orthology call at one merge (pairwise level) | RBH · InParanoid-like (seed + in-paralogs) · local-tree species overlap | 2-species truth: OrthoBench RefOG pairs for human–mouse, human–fly; VGNC | pair P/R/F1 |
| S4 | Ancestral representatives | 1 medoid · k=3 diverse · consensus · HMM profile | recall at deep merges (e.g. Bilateria) on OrthoBench | RefOG recall/fission |
| S5 | Synteny evidence | off · tie-break only · additive score | OrthoBench (Ensembl GFF), VGNC; paralog-resolution precision | ΔF, Δprecision |
| S6 | Structure evidence | off · 3Di joint search only at deep nodes (ProstT5 or AFDB) | S2 recall at deep nodes; OrthoBench fissions | Δrecall |
| S7 | Isoforms | canonical only · best-isoform-pair representatives | S3 pairs | ΔF1 |
| S8 | Gene-loss validation (DNA) | none · tblastn of the HOG representative in the genomic locus flanked by syntenic neighbours | small curated set (e.g. known losses in mammals) | precision of "lost" calls |
| S9 | Omics/expression | exploratory: Bgee expression correlation as paralog tie-breaker, for vertebrates only | VGNC pairs where synteny is ambiguous | report effect or a negative result |
| F | Full pipeline | best config ± each layer (ablation) | OrthoBench (all 7 measures) · QfO 2022 local (SwissTree, TreeFam-A, STD/G-STD2, GO, EC, VGNC) | Pareto position vs 21 public methods + runtime |

Pre-registered claims: (a) runtime grows ~linearly with S (8/16/32/64/78 species) versus OrthoFinder and SonicParanoid run on the same node; (b) on QfO we target the precision–recall Pareto front for SwissTree, TreeFam-A and G-STD2-Euk, with no single scalar; (c) on OrthoBench we report all measures. Metrics, test sets and controls are frozen once round 3 of review completes. Any later change goes to RESEARCH_TREE.md with a reason.

## 4. Data
- QfO 2022_02 (with updated UP000000437) proteomes, which is what the benchmark expects. Gene coordinates come from QfO 2026_02 gcoord/gff mapped by accession, with coverage reported.
- OrthoBench 12 Ensembl proteomes plus the matching Ensembl GFF3 (synteny development set).
- Structure: ProstT5 3Di (cost from feasibility probe) or AFDB where available.
- Development vs test discipline: OrthoBench and a *subset* of QfO species are used for tuning. The final QfO run uses frozen parameters. SwissTree/TreeFam-A families are not inspected during development.

## 5. Verification
- Toy simulator: simulate a species tree and gene families with duplications, losses and gene order (adjacencies, inversions), then sequences (e.g. a simple WAG/LG + gamma evolver). This gives exact HOG truth. Unit tests check that the merge reproduces the truth on easy cases (positive control) and degrades gracefully on shuffled gene order (negative control for synteny) and random sequences (negative control).
- Benchmark controls: a perfect-truth input scores 1 on OrthoBench, and random grouping scores ~0. For QfO, a public submission is re-scored locally and must match its published numbers.
- Every result is reviewed by an agent that did not produce it before the PI relies on it.
- Seeds: fixed list {11, 23, 47} recorded in configs.

## 6. Compute plan (48 h, 24 CPU, 150 GB)
- T+1–3 h: design reviews; data and benchmark readiness.
- T+3–10 h: core v0 (sequence only) on OrthoBench; S1, S2, S3 sub-benchmarks; toy simulator.
- T+10–20 h: QfO full run and local QfO benchmark; competitor runtime runs on species subsets.
- T+20–36 h: layers S4–S9; ablations.
- T+36–44 h: frozen final runs; independent result review.
- T+44–48 h: REPORT.md.
- Checkpointing: each merge writes its ancestral proteome to disk (resume from the last completed node).

## 7. Risks
- Recall loss at deep nodes (Hieranoid's known weakness) → S4 and S6.
- Over-merging at deep levels hurts QfO precision (OMA HOGs vs OMA Pairs) → conservative deep merges, export of both HOGs and pairs.
- QfO local run may not work without Docker → fallback is to reimplement the pairwise tests from reference files, validated on public submissions.
- Synteny coverage on QfO may be patchy → report coverage; synteny is optional evidence.
