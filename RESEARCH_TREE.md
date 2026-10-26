# RESEARCH TREE — p10 orthology prediction

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding
Clock: START.txt = 2026-10-26T14:40 CET; hard deadline 2026-10-28T14:40 CET (48 h). Check `date` vs START.txt at each milestone.

## Status summary (latest first)

### M0 — 2026-10-26 14:45 (T+0h05)
Project initialised: git repo (root = method + research record), separate benchmark repo in `bench_repo/` (nested, gitignored by root). PROMPT.md saved verbatim, START.txt written. Three parallel agents launched: (1) literature review → LITERATURE.md, (2) software env → conda_envs/orth + ENVIRONMENT.md + timing calibration, (3) benchmark data acquisition (OrthoBench, QfO 2022 proteomes + local QfO benchmark, competitor public results) → DATA.md + bench_repo/.

Next: when literature arrives, draft DESIGN.md and run up to 3 review rounds (external PI + specialist reviewers). Implementation of the family-clustering sub-task can start in parallel as soon as the env exists, because every plausible design needs a sub-quadratic family step.

Decision note M0: benchmark lives in its own repo (bench_repo/) as the prompt requires; raw data never committed (data/ gitignored, checksums in DATA.md).

## Tree
- [~] 0. Infrastructure
  - [x] 0.1 Conda env (conda_envs/orth, all tools incl. orthofinder/fastoma/sonicparanoid) — diamond human-vs-human sensitive 12 s, mmseqs s7.5 88 s
  - [~] 0.2 Benchmark data + local runners: OrthoBench, QfO (SwissTree, TreeFam-A, STD, GO/EC), competitor results (agent)
- [x] 1. Literature review → LITERATURE.md (agent). [!] Open space: no method uses structure/synteny/expression inside inference; competitor QfO predictions downloadable (21 sets, table in LIT §7.4); predict on QfO 2022_02 (not 2026) proteomes.
- [~] 2. Design → DESIGN.md, ≤3 review rounds. v0 drafted: 'Ancestra' progressive ancestral-proteome merging along guide tree with synteny/3Di evidence. Review round 1 running (5 reviewers).
  - [~] 2.1 Structure feasibility probe (ProstT5 CPU throughput, AFDB coverage)
  - [~] 2.2 Toy simulator + OrthoXML utils (design-independent)
- [ ] 3. Sub-task: family placement / clustering (sub-quadratic) — own benchmark
- [ ] 4. Sub-task: gene-tree inference (AA vs AA+3Di) — own benchmark
- [ ] 5. Sub-task: species tree (universal vs inferred vs family-local) — own benchmark
- [ ] 6. Sub-task: reconciliation → HOGs → OrthoXML — own benchmark
- [ ] 7. Extra layers: structure (3Di), synteny (GFF), isoforms, omics, DNA gene-loss validation
- [ ] 8. Full pipeline benchmark vs competitors (accuracy + runtime + scaling)
- [ ] 9. REPORT.md
