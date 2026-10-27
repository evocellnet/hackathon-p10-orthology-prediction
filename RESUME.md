# RESUME.md — how to pick up this project (paused 2026-10-27 ~10:20 CET, T+19h40 of the original 48 h budget)

The user paused the project to move it to a **GPU node** and explore GPU speed-ups of the method. Everything is committed in git (root repo; the benchmark code is in the separate repo `bench_repo/`). No processes were left running.

## 1. Read first (in this order)
1. `docs/ONE_PAGE_SUMMARY.md`: main findings and figures (`docs/figures/`).
2. `RESEARCH_TREE.md`: full research record. The top summaries M3/M2 hold the status; the tree holds every experiment with its decision.
3. `docs/DISCLOSURES.md` (18 protocol deviations) and `docs/VERIFICATION.md` (independent recomputation).
4. `LESSONS.md` (14 process rules), `DESIGN.md` (pre-registered protocol v2), `PROMPT.md` (original brief, verbatim).
5. `README.md` (CLI and quick start), `docs/core_v0.md` (algorithm and presets §7.9–7.14), `docs/qfo_run.md` (QfO runs §1–10).
6. `CLAIMS.md` lists every quantitative claim with its evidence path. Statuses still need updating from `docs/VERIFICATION.md`; that is open task T5.

## 2. Environment and data
- Conda env: `source /opt/conda/etc/profile.d/conda.sh && conda activate <root>/conda_envs/orth` (specs in `envs/orth.yml`, `envs/orth.lock.txt`; see `ENVIRONMENT.md`). The env and the data are NOT in git. On a new node, rebuild the env from `envs/orth.lock.txt` and re-download data per `DATA.md` (every dataset lists source, version, md5 and license). If the new node mounts the same /scratch path, nothing needs re-downloading.
- Large data under `data/`:
  - OrthoBench, QfO 2022_02 (+gcoord), Ensembl 99 GFF/genomes, Bgee;
  - Foldseek AFDB/UniProt DB (491 GB);
  - Foldcomp afdb_uniprot_v6 (1.2 TB; only a cross-check, so it can be deleted);
  - 3Di tables `data/structure/*_3di_best.fa`;
  - the 1000-proteome scalability set `data/scale1000/`.
- Cached pipeline work dirs live under `work/` (gitignored). Most important: `work/qfo2022_det/` holds the deterministic QfO S0/S2 caches plus the `cand2_accurate` and `cand2_fast` outputs that were evaluated as final.

## 3. Frozen method
- Git tag **`freeze-v1`**. `results/freeze.json` holds sha256 of every `ancestra/*.py`; `ancestra/` was not changed after the freeze.
- Presets (`--preset`):
  - `accurate` (default) = FastTree local trees on the synteny-on pass, NJ on the off pass, plus deep-node species-overlap absorption.
  - `fast` = NJ + so_struct + deep absorption.
- Options: `--gene-order DIR` turns on constrained synteny; `--family-rescue both --family-rescue-guard` turns on the 3Di + profile S0 rescue (off by default).
- Final results:
  - `results/orthobench_final/` (OrthoBench test, one-time);
  - `results/qfo_final/` (full QfO incl. STD/G-STD2/FAS + 7 ablations);
  - `results/scaling/` (nested QfO subsets vs OrthoFinder/SonicParanoid2).
- **Any change to `ancestra/` from now on is post-freeze "v2" work.** Give it a new tag, and evaluate only on dev data unless a new, clean benchmark is used. QfO and OrthoBench test have been consumed (see DISCLOSURES).

## 4. What was in progress at pause (each has its own resume note)
| Task | State | Resume note |
|---|---|---|
| T1 **OMA comparators on OrthoBench** (user requirement: ≥1 OMA resource in every comparison) — corrected FastOMA run + OMA-browser HOGs mapped to OrthoBench IDs | interrupted | `bench_repo/results/orthobench_competitors/OMA_RESUME_NOTE.md` |
| T2 **Runtime scaling** incl. FastOMA (user requirement) on nested QfO subsets 8/16/32/64 | Ancestra 8→78 and OrthoFinder/SonicParanoid2 8→32 done; FastOMA and the remaining competitor runs pending | `results/scaling/RESUME_NOTE.md` |
| T3 **Large-scale demo** (125/250/500/1000 UniProt proteomes, 0.5–4.75 M proteins; data + trees ready) | n=125 started and killed in S0; nothing measured | `results/scaling_large/RESUME_NOTE.md` (projected RSS for n=1000 scoring ≈142 GB on CPU → the stage needs streaming/chunking or more RAM) |

## 5. Next plans (priority order)
1. **T1–T3 above.** On the GPU node, re-measure competitor and Ancestra runtimes on an *idle* node. Current timings come from a shared node with load 20–65, so they are upper bounds.
2. **GPU speed-up exploration** (user goal). Profile first: on QfO the CPU time is S2 scoring 8.9 CPU-h > S0 3.0 CPU-h > `accurate` merge 7 CPU-h (61k FAMSA+FastTree jobs) ≫ `fast` merge 6 min. Candidates:
   - (a) GPU homology search for S0/S2: MMseqs2-GPU (`mmseqs ... --gpu 1`, ungapped prefilter + GPU SW, available in MMseqs2 ≥16) as a drop-in backend behind `ancestra/scoring.py`'s backend interface; or DIAMOND GPU builds if available.
   - (b) Batched GPU Smith–Waterman for the explicit within-family pair lists (S2 is "align these pairs", which suits GPU batching: libraries like CUDASW++4 / pyopal-GPU / custom CUDA via cupy/numba-cuda).
   - (c) GPU ProstT5 → 3Di for proteins without AFDB structures. CPU was infeasible at 125–250 res/s; GPU makes the structure layer complete, see `docs/structure_feasibility.md`.
   - (d) Replace the 61k small FAMSA+FastTree jobs in the `accurate` merge with a batched GPU distance/NJ (or a pLM-embedding distance) and compare with S4's label-accuracy benchmark (`docs/S4_trees.md`).
   - (e) pLM embeddings (ESM-2/ProstT5) + FAISS-GPU kNN as an alternative S0 candidate generator, which targets the sequence-similarity ceiling (0.82). Benchmark with the S0 protocol (`docs/S0_families.md`, dev RefOGs only).

   Keep outputs byte-identical to the CPU path where the algorithm is unchanged, and test determinism (LESSONS #8).
3. **Known method weaknesses for v2** (dev-only work):
   - vertebrate over-merging (giant root HOGs; VGNC/EC precision; 6k-gene Chordata HOG);
   - recall vs OrthoFinder at Bilateria/Eukaryota levels (merge splits at the LCA);
   - S0 splits involving *S. pombe* (needs a different clustering objective);
   - Newick child-order sensitivity of the merge (DISCLOSURES #11).
4. **REPORT.md** (the skeleton exists): fill it from `docs/ONE_PAGE_SUMMARY.md`, `results/*_final/summary.md`, `docs/DISCLOSURES.md` and `docs/VERIFICATION.md`; update `CLAIMS.md` statuses. The subagent Write tool refused to write report files, so the PI must assemble REPORT.md.
5. Milestone hygiene: update the RESEARCH_TREE status summary, prune LESSONS, and get a fresh-agent retrospective at each milestone (see `PROMPT.md`).

## 6. Clock
`START.txt` holds the original start. If the user grants a new budget, rewrite START.txt or add a line, and use `scripts/clock.sh`. The original 48 h budget had 28h35m left at pause.
