# RESEARCH TREE — p10 orthology prediction

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding
Clock: START.txt = 2026-10-26T14:40 CET; hard deadline 2026-10-28T14:40 CET (48 h). Check `date` vs START.txt at each milestone.

## Status summary (latest first)

### M1 — 2026-10-26 16:05 (T+1h25) — design v2
Design went through two review rounds (8 reviews in docs/reviews/). v0 (pure progressive ancestral merging) was rejected on cost grounds: the compute review estimated only a ~2.5× saving over all-vs-all on QfO because ancestral pan-proteomes are large. v1/v2 restructure it as **family partition (near-linear, two-level) → within-family scores → level-synchronous ancestral merging across all families, with ancestral gene-order context inside the merge**. Leakage fixes: the OrthoBench split is by homology component; at most 3 QfO peeks before freeze; VGNC is not dev evidence. The success rule is non-inferiority (lower 95% bound of Δ ≥ −0.03).

Running: S0 family-partition benchmark, core engine v0 (sequence only, toy sim), data/QfO local pipeline, structure feasibility probe. Done: env, literature, simulator, OrthoXML tools (48 tests). Decision: review round 3 is deferred to the v0-results checkpoint (~T+8h), so the last design round critiques data rather than plans.


### M0 — 2026-10-26 14:45 (T+0h05)
Project initialised: git repo (root = method + research record), separate benchmark repo in `bench_repo/` (nested, gitignored by root). PROMPT.md saved verbatim, START.txt written. Three parallel agents launched: (1) literature review → LITERATURE.md, (2) software env → conda_envs/orth + ENVIRONMENT.md + timing calibration, (3) benchmark data acquisition (OrthoBench, QfO 2022 proteomes + local QfO benchmark, competitor public results) → DATA.md + bench_repo/.

Next: when literature arrives, draft DESIGN.md and run up to 3 review rounds (external PI + specialist reviewers). Implementation of the family-clustering sub-task can start in parallel as soon as the env exists, because every plausible design needs a sub-quadratic family step.

Decision note M0: benchmark lives in its own repo (bench_repo/) as the prompt requires; raw data never committed (data/ gitignored, checksums in DATA.md).

## Tree
- [~] 0. Infrastructure
  - [x] 0.1 Conda env (conda_envs/orth, all tools incl. orthofinder/fastoma/sonicparanoid) — diamond human-vs-human sensitive 12 s, mmseqs s7.5 88 s
  - [x] 0.2 Benchmark data + local runners (bench_repo/, DATA.md). [!] The QfO 2022 benchmark runs locally (Darwin binary, no docker). Re-scored OrthoFinder/FastOMA/OMA submissions match published numbers exactly on deterministic tests; STD within ~1.4 SE. OrthoBench controls: RefOGs F=100, random F=0. 21 competitor prediction sets downloaded. Ensembl 99 GFF (12 OrthoBench species), Bgee v15.2, QfO isoforms/CDS/gcoord (19 spp) available.
- [x] 1. Literature review → LITERATURE.md (agent). [!] Open space: no method uses structure/synteny/expression inside inference; competitor QfO predictions downloadable (21 sets, table in LIT §7.4); predict on QfO 2022_02 (not 2026) proteomes.
- [~] 2. Design → DESIGN.md v2 after 2 review rounds; round 3 at v0-results checkpoint.
  - [-] 2.0 v0 'pure progressive pan-proteome merging' — rejected: ~2.5x saving only vs all-vs-all (r1_algorithms_compute)
  - [x] 2.3 v1/v2 family-partition + level-synchronous merge; split by homology component (random split leaked: 5 paralogous RefOG groups)
  - [x] 2.1 Structure feasibility (docs/structure_feasibility.md): ProstT5 CPU 125–250 res/s → 140–980 h for full sets (infeasible). AFDB covers 90% of QfO accessions; the 1.8 GB Foldseek AFDB/Proteome DB already gives 3Di for 263k QfO / 98k OrthoBench proteins. Decision: S6 uses AFDB 3Di for deep/ambiguous representatives only; decide on streaming the 491 GB AFDB/UniProt DB when S6 starts
  - [~] 2.1b User suggestion (T+1h20): Foldcomp prebuilt AFDB databases → foldseek createdb for true 3Di of matched sequences (agent evaluating coverage/download)
  - [x] 2.2 Toy simulator + OrthoXML utils — flattening identical to QfO map_relations.py; RefOGs→F=100 in official benchmark.py
- [x] 3. S0 family partition — own benchmark (docs/S0_families.md)
  - [!] 3.1 Ceiling finding: no sequence-similarity partition can reach the pre-registered ≥95% dev RefOG completeness; even all DIAMOND hits at e≤1e-5 give 27/34, ceiling ≈0.82 (6 dev RefOGs have members with no hit to the rest at e≤1e-3 even with ultra-sensitive). Deviation recorded; target not met, method accepted as near-ceiling.
  - [x] 3.2 Chosen: 2sA_dc30_cc_sh30_RA (diamond cluster 30%/cov50 → all genes vs reps → CC with ≥30% cov → Leiden >7.5k → reassignment). Dev completeness 0.794, fission 1.32, Σg²/N² 0.0075, max family 7,478, 263 s wall, 6.2 GB. **Test 0.722** (logged in results/test_log.tsv). Byte-identical rerun.
  - [-] 3.3 Rep-vs-rep linking (r2 proposal) gave 0.441; gene-vs-rep N×R search was the key.
  - [ ] 3.4 Root rescue for the 145 Leiden-split families + 15k singletons (6.1%)
- [~] 3b. Core engine v0 (S2 scoring + S3 merge, sequence-only) on toy sim (agent)
- [ ] 4. Sub-task: gene-tree inference (AA vs AA+3Di) — own benchmark
- [ ] 5. Sub-task: species tree (universal vs inferred vs family-local) — own benchmark
- [ ] 6. Sub-task: reconciliation → HOGs → OrthoXML — own benchmark
- [ ] 7. Extra layers: structure (3Di), synteny (GFF), isoforms, omics, DNA gene-loss validation
- [ ] 8. Full pipeline benchmark vs competitors (accuracy + runtime + scaling)
- [ ] 9. REPORT.md
