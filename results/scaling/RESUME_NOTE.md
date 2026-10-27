# RESUME NOTE: runtime scaling benchmark (paused 2026-10-27 10:10 CET, T+19h30)

All my processes are stopped (queue runners killed, interrupted run's whole session killed; checked with ps). Nothing else of mine runs.

## Completed (kept; per-run files in `results/scaling/runs/<label>/`, raw table `scaling_raw_runs.tsv`, aggregate `scaling.tsv`, `summary.md`, `scaling.png`)
| label | what | wall h | CPU-h | peak RSS GB |
|---|---|---|---|---|
| anc_8, anc_16, anc_32, anc_32_go | fresh Ancestra runs (S0+S2 are the real costs; merge was the pre-preset NJ default) | 0.07 / 0.17 / 0.53 / 0.52 | 0.47 / 0.99 / 3.20 / 3.25 | 0.8 / 2.1 / 5.7 / 5.6 |
| anc_8_acc, anc_16_acc, anc_32_acc, anc_32_go_acc | merge+output re-run with `--preset accurate` on cached S0+S2 of the base run (merge-only cost) | - | 0.31 / 1.17 / 2.84 / 2.92 (merge only) | - |
| anc_*_fast (8,16,32,64,78) | same, `--preset fast` | - | merge only, 0.002-0.04 | - |
| anc_64, anc_78 | fresh, accurate (default) | 2.65 / 2.83 | 16.30 / 18.54 | 25.1 / 28.9 |
| of_8, of_16, of_32 | OrthoFinder 3.2 | 0.19 / 0.57 / 2.27 | 1.22 / 3.54 / 14.43 | 2.4 / 2.9 / 4.8 |
| sp_8, sp_16, sp_32 | SonicParanoid2 default | 0.15 / 0.38 / 1.80 | 0.77 / 2.34 / 11.11 | 8.0 / 8.1 / 8.1 |
| fo_8 | FastOMA 0.5.1 (LUCA.h5) | 0.37 | 1.83 | 13.6 |

## Interrupted
- **fo_16** (FastOMA, 16 species), killed after ~20 min at the PI's pause request. Must be rerun from scratch. Its partial dirs were moved to `work/scaling/interrupted/` (not needed).

## Not run
fo_16 (interrupted), fo_32, fo_64, of_64, sp_64 (dropped for budget; extrapolated 38-45 CPU-h), all of them at 8 threads. Ancestra and the other competitors are complete for the planned sizes (Ancestra also 78).

## Resume (exact commands)
All from `/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction`, env `source /opt/conda/etc/profile.d/conda.sh && conda activate $PWD/conda_envs/orth`. Inputs already exist: `work/scaling/inputs/sub{8,16,32,64,78}` (symlinked FASTA), `fo_{8,16,32,64}/` (FastOMA: real `proteome/` dir of file symlinks + `species_tree.nwk`), trees `results/scaling/subsets/sub{n}.nwk` and `.fastoma.nwk`. `LUCA.h5` at `data/orthobench/runs/fastoma/db/LUCA.h5` (9.9 GB, must be reachable from the new node; nextflow home `data/orthobench/runs/fastoma/nxf_home`).
1. Edit `work/scaling/queue3.txt` (one job per line `label|kind|n_species|cap_h|extras`); currently: `fo_8` (done), `fo_16`, `fo_32`, `anc_32_go_acc` (done), `fo_64` commented out. Remove done lines or leave them (done = `results/scaling/runs/<label>/runtime.json` exists).
2. Start the runner (ONE run at a time, 8 threads, 60 s uptime log, /usr/bin/time -v, timeout = cap, fresh work dir required, so `rm -rf work/scaling/runs/fo_16*` first if present):
   `nohup python scripts/scaling_queue3.py > work/scaling/queue3.log 2>&1 &`
   Single FastOMA run by hand: `/usr/bin/time -v -o t.txt bash scripts/scaling_fastoma.sh 16 work/scaling/runs/fo_16` (nextflow, `data/orthobench/runs/fastoma/local.config`, 8 cpus).
   OrthoFinder: `orthofinder -f work/scaling/inputs/sub64 -t 8 -a 8 -mog -o work/scaling/runs/of_64` (kind `of`); SonicParanoid: `sonicparanoid -i work/scaling/inputs/sub64 -o work/scaling/runs/sp_64 -p scal64 -t 8 -m default` (kind `sp`). Both only read their queue lines `of_64|of|64|<cap_h>` / `sp_64|sp|64|<cap_h>`; the runner writes `runtime.json` with `censored: true` on timeout.
3. After every batch: `python scripts/scaling_analyze.py` regenerates `scaling.tsv`, `summary.md` (+ `status_note.md`, edit it), `scaling.png`, `fits.json`.
Notes: node loads differ (new node: wall-time not comparable with the shared-node runs; rerun at least one anchor, e.g. of_8/anc_8, there if timing across nodes is to be mixed). Code for Ancestra runs is frozen in `work/scaling/code_snapshot` (taken 2026-10-27 05:00).
