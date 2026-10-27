# OMA comparators on OrthoBench — PAUSED 2026-10-27T10:09:12+01:00 (user request: move to GPU node)

Task (user request, post-freeze): include >=1 OMA resource in all comparisons. Two resources: (1) FastOMA re-run fixed, (2) OMA browser HOGs.

## Done
- **OMA browser HOGs (All.May2026) — scored, final.** Files: `oma_browser/` (primary, Bilateria-level HOGs, pre-declared) and
  `oma_browser_root/` (root-level HOGs). score.json, per_refog.tsv, runtime.json (conversion only), version.txt, mapping_stats.txt,
  mapping_final.tsv.gz, PREREGISTRATION.txt (+ addendum: a first, buggy Bilateria-cut rule scored F all 36.1 / test 45.5; fixed rule reported).
  Provenance in root `DATA.md` §14 (md5s match Zenodo; CC BY 4.0). LEAKAGE: OMA contains all 12 species.

| comparator | test F | P | R | exact | dev F | all-70 F | ΔF test (Ancestra − comparator) [95% CI] | SE | NI verdict (Ancestra vs comparator) | ΔF all-70 [95% CI] | CPU-h | wall-h | peak RSS GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| oma_browser | 60.4 | 88.9 | 45.8 | 3/36 | 50.5 | 55.8 | +15.7 [+1.8, +29.2] | 7.2 | comparable (non-inferior; CI > 0) | +22.1 [+13.2, +30.7] | 0.31 | 0.26 | 1.5 |
| oma_browser_root | 7.2 | 3.9 | 49.5 | 2/36 | 7.3 | 7.3 | +69.0 [+22.3, +75.3] | 12.7 | comparable (non-inferior; CI > 0) | +70.6 [+54.6, +77.2] | 0.31 | 0.26 | 1.5 |
| fastoma_luca | 33.5 | 24.3 | 53.7 | 10/36 | 29.6 | 31.7 | +42.7 [+14.1, +55.2] | 11.2 | comparable (non-inferior; CI > 0) | +46.2 [+26.3, +57.7] | 12.74 | 3.50 | 13.6 |

  (ΔF = F(Ancestra accurate, results/orthobench_final/scores/accurate_go.json) − F(comparator), paired bootstrap B=2000 seed 11, NI margin 3 pp.)
- **FastOMA v1 diagnosis (complete):** inputs, species tree (12 leaves = file names) and export level were correct (top-level OrthoXML HOGs ==
  RootHOGs.tsv, root = Bilateria). OMAmer DB is not the issue: unplaced 4.10% with LUCA.h5 vs 4.17% with Metazoa.h5 (`fastoma_v2/metazoa_arm.txt`).
  The 21% absent genes = 23.0k never in a FastOMA family (9.1%) + 24.3k root-level singletons (9.7%) + **5,497 genes (2.2%) silently LOST by a
  FastOMA 0.5.1 bug**: `FastOMA/_hog_class.py::split_hog` rebuilds a split HOG only from subHOGs containing a representative and drops the rest
  (same code on upstream `dev`). Losses concentrate in big families (RefOG005/006/008/013 almost entirely missing in v1). Reproduced on one
  family (F0800523: 238/356 genes lost). Patch (`data/orthobench/runs/fastoma_v2/patched_pkg/FastOMA/_infer_subhog.py`, marked
  "ANCESTRA-BENCH PATCH", used via PYTHONPATH; site-packages untouched) keeps such subHOGs as frozen HOGs and logs them to
  `rescued_subhogs.tsv`; frozen HOGs never merge, so "released 0.5.1" output = patched output minus rescued HOGs (verified on F0800523).
- Pre-declared primary FastOMA row: `fastoma_v2` = FastOMA 0.5.1 + LUCA.h5 + patch; `fastoma_v2_luca_released` = same run minus rescued HOGs.
  Metazoa arm stopped for time (no score).

## In progress / remaining
1. **FastOMA LUCA+patch run, ~90% done** (59/66 hog_rest batches + all upstream steps cached; 0/6 hog_big; 8 batches were killed by the pause).
   Run dir `/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/data/orthobench/runs/fastoma_v2/luca_patched/` (logs_seg1/, logs_seg2/ = /usr/bin/time + load of the two segments so far;
   seg1 ran at 4 cpus in parallel with the Metazoa arm and was killed/resumed at 8 cpus). The nextflow work dir lived on node-local /tmp;
   a full copy is in `luca_patched/nxwork_backup/` (921 MB). rescued_subhogs.tsv is already cleaned of the killed batches
   (raw copies *.seg1_raw.tsv / *.seg2_raw.tsv; killed rhog lists killed_rhogs*.txt).
2. Then: convert + score + runtime, update SUMMARY.md, append the OMA section to results/orthobench_final/summary.md, commit.

## Exact commands to resume (<=8 threads)
```bash
R=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction; D=$R/data/orthobench/runs/fastoma_v2
# nextflow -resume needs the cached task dirs at their ORIGINAL path -> restore it (any node):
OLD=/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-p10-orthology-prediction/d3baa247-6354-5748-9b99-4ccb63a14cf9/scratchpad/nxwork_luca_patched
mkdir -p $(dirname $OLD); [ -d $OLD ] || ln -s $D/luca_patched/nxwork_backup $OLD
mkdir -p $D/luca_patched/logs
RESUME=-resume nohup $D/launch.sh luca_patched $R/data/orthobench/runs/fastoma/db/LUCA.h5 $OLD > $D/luca_patched/launch3.out 2>&1 &
# (if resume finds no cache: rm the symlink and rerun without RESUME, ~2-3 h at 8 cpus on local disk; delete rescued_subhogs.tsv first)
# when logs/exit_code.txt == 0:
cd $R/bench_repo && source /opt/conda/etc/profile.d/conda.sh && conda activate $R/conda_envs/orth
python scripts/fastoma_v2_to_og.py $D/luca_patched/run/Output/FastOMA_HOGs.orthoxml $D/luca_patched/rescued_subhogs.tsv $D/luca_patched/fastoma_v2 $R/data/orthobench/Open_Orthobench/BENCHMARKS/Input
B=$R/data/orthobench/Open_Orthobench/BENCHMARKS
for v in patched:fastoma_v2 released:fastoma_v2_luca_released; do mkdir -p results/orthobench_competitors/${v#*:}
  python -W ignore scripts/orthobench_eval.py $B $D/luca_patched/fastoma_v2_${v%%:*}.og --label ${v#*:} --json results/orthobench_competitors/${v#*:}/score.json --per-refog-tsv results/orthobench_competitors/${v#*:}/per_refog.tsv
  python scripts/make_runtime_json.py $D/luca_patched/logs_seg1 results/orthobench_competitors/${v#*:}/runtime.json --threads 8 --include $D/luca_patched/logs_seg2 $D/luca_patched/logs --note "FastOMA 0.5.1 + LUCA.h5 + gene-loss patch; 3 segments (pause/resume); omamer step reused v1 hogmaps (add ~0.75 CPU-h, 13.6 GB RSS from v1 trace for a full run); local NVMe work dir; LEAKAGE"
done
python scripts/summarize_competitors.py
python scripts/oma_comparators_vs_ancestra.py $R/results/orthobench_final/scores/accurate_go.json $R/splits/orthobench_split.json results/orthobench_competitors/{fastoma_v2,fastoma_v2_luca_released,oma_browser,oma_browser_root,fastoma_luca}
# -> paste under "## OMA comparators (added at user request, post-freeze)" in $R/results/orthobench_final/summary.md; write version.txt; commit.
```
