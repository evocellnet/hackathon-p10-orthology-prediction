# Scaling-large: paused at 2026-10-27T10:06+01:00 (T+19h25m; 28h34m of the 48 h budget remained)

## Finished
Nothing. No scaling result exists yet; results/scaling_large/ holds only the pre-existing subsets/ and trees/.
Prepared and kept: `scripts/scale_large_run.sh` (new; ancestra/ untouched, tag freeze-v1 intact).

## Interrupted
- Subset n=125 (523,225 proteins), `--preset fast`, 12 threads, work dir `work/scale_large/n125/`.
- Started 10:04:40; killed (process group 241773) at ~10:05:55, about 70 s in. Stage `input` completed (0.3 s, 302 MB RSS).
  Stage `families` (S0) was in its first DIAMOND cluster step and was lost; only the partial `out/stages/families.tmp247688/` remains.
- No processes of mine remain (ps check for ancestra|diamond|mmseqs|FastTree|famsa|scale_large was empty).

## What scripts/scale_large_run.sh does
`scripts/scale_large_run.sh N [PRESET] [TAG]`: wipes and recreates `work/scale_large/<TAG>/` (default n<N>), builds `input/` of symlinks from
`results/scaling_large/subsets/n<N>.txt` to `data/scale1000/proteomes/`, runs `python -m ancestra run --preset fast --threads 12`
under `/usr/bin/time -v` (`logs/time_v.txt`, `logs/run.log`), logs `uptime` + `free -g` + my process-tree RSS every 60 s (`logs/monitor.log`),
and kills the run if my tree exceeds 100 GB RSS (`LIMIT_GB`). Pipeline output goes to `work/scale_large/<TAG>/out/` (manifest.json there).

## Resume (detached)
```
cd /scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
setsid nohup scripts/scale_large_run.sh 125 fast > /dev/null 2>&1 &     # then 250, 500, 1000 in order
# optional last: setsid nohup scripts/scale_large_run.sh 250 accurate n250_accurate > /dev/null 2>&1 &
```
The wrapper deletes the work dir of the tag it starts, so n=125 simply restarts from S0 (only ~70 s lost).
The time budget must be recomputed on resume (the 48 h clock in START.txt keeps running): the T+36h launch stop and T+40h finish are
defined on that clock (`scripts/clock.sh`); if the pause ends later, say so in the report.

## Remaining plan after each run
`python scripts/verify_orthoxml.py work/scale_large/n<N>/out --json results/scaling_large/verify_n<N>.json` (structural subset if too slow),
collect stage wall/CPU/RSS from `out/manifest.json` + `logs/time_v.txt`, families, root HOGs, HOGs per level.
Then write scaling.tsv, summary.md, scaling_large.png (load the dataviz skill first).

## Extrapolation so far (no scaling runs yet; from docs/qfo_run.md only, 988,778 proteins, 12 threads, loaded node)
QfO: S0 1,534 s / 3.0 CPU-h / 6.4 GB; scores 3,886 s / 8.9 CPU-h / 29.5 GB; merge (fast) ~6 min / 9.5 GB; accurate merge 55 min / 7 CPU-h.
Subset totals: 125 -> 0.53x QfO, 250 -> 1.10x, 500 -> 2.28x, 1000 -> 4.80x proteins (4.75 M).
Linear scaling from QfO would give n=1000: S0 ~2 h, scores ~5.2 h wall, ~43 CPU-h, scores RSS ~142 GB (above the 100 GB cap),
so the n=1000 memory projection must be re-checked from the n=125/250/500 fits before launching. Not a measurement.

## RESUME 2026-10-27 11:10: lane B running (LIMIT_GB raised to 300 in scale_large_run.sh)
Chain `work/scale_large_chain.sh` runs n=125,250,500 (fast, 12 threads). n=125 done 11:11-11:35: 23.9 min wall, 3.26 CPU-h, 7.9 GB RSS. n=1000 launched manually only after RSS/time projection check. Collector: scripts/scaling_large_collect.py (to write).

DEADLINE UPDATE (PI, 11:45): soft deadline 2026-10-28T14:40; n=1000 must finish before 10-28 06:40; final report by 10-28 08:40.

## 16:00 status: n125/250/500 done (0.40/1.12/3.32 h wall; 7.6/23.4/68.8 GB RSS; collector run). Projection for n1000 from fits: ~9.7 h wall, ~212 GB RSS (<300, ends ~02:00 < 06:40 limit) -> launched 16:00 via scripts/scale_large_run.sh 1000 fast (LIMIT_GB 300). Afterwards: n250 accurate (`scripts/scale_large_run.sh 250 accurate n250_accurate`), verify_orthoxml on each, rerun `python scripts/scaling_large_collect.py`.
PI UPDATE (01:25 10-28): n1000 may run until 11:30; if unfinished record censored '> T' with stage times. Report n500->n1000 per-stage slopes honestly (S0 4922s -> 17614s).
02:36 10-28: n1000 DONE (10.6 h wall, 64.3 CPU-h, 221 GB RSS). n250_accurate launched 02:35 (work/scale_large/n250_accurate). verify_n125/250/500 PASS; verify_n1000 running. Collector reruns: python scripts/scaling_large_collect.py (collects accurate tag too).
04:35: lane B complete (n125/250/500/1000 fast + n250_accurate, all verify PASS); summary/png/fits written. Only sp_64 pending in lane A.
