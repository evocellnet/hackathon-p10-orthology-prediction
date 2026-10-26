# Runtime scaling (DESIGN §5)

8 threads per run, one run at a time, shared 24-core node (load logged per run in `runs/*/load.log`; wall times from a busy node are upper bounds, CPU-h = user+sys from `/usr/bin/time -v` is the primary metric). Fresh work dir per run. Censored runs: `> cap`.

Subsets: nested, domain-stratified (random.Random(11)), `subsets/sub{n}.{txt,nwk,json}`; guide tree = NCBI tree pruned to the subset.

| method | 8 sp (62k) | 16 sp (144k) |
|---|---|---|
| Ancestra (sequence-only) | 0.07 h / 0.47 CPU-h / 0.8 GB | 0.17 h / 0.99 CPU-h / 2.1 GB |
| Ancestra (+gene order) | not run | not run |
| OrthoFinder 3.2 | 0.19 h / 1.22 CPU-h / 2.4 GB | 0.57 h / 3.54 CPU-h / 2.9 GB |
| SonicParanoid2 (default) | 0.15 h / 0.77 CPU-h / 8.0 GB | 0.38 h / 2.34 CPU-h / 8.1 GB |

Cells: wall / CPU-h / peak RSS (max over processes as reported by `time -v`).

## Log-log slope of CPU-h (uncensored points)

| method | x | n pts | slope | OLS 95% CI | bootstrap 95% CI | R² |
|---|---|---|---|---|---|---|
| Ancestra (sequence-only) | proteins | 2 | 0.88 | n<3 | n<3 | 1.000 |
| Ancestra (sequence-only) | species | 2 | 1.06 | n<3 | n<3 | 1.000 |
| OrthoFinder 3.2 | proteins | 2 | 1.26 | n<3 | n<3 | 1.000 |
| OrthoFinder 3.2 | species | 2 | 1.53 | n<3 | n<3 | 1.000 |
| SonicParanoid2 (default) | proteins | 2 | 1.32 | n<3 | n<3 | 1.000 |
| SonicParanoid2 (default) | species | 2 | 1.59 | n<3 | n<3 | 1.000 |

Bootstrap: percentile, B=5000, seed 11, resampling the points (resamples with <3 distinct x skipped). With 3-5 points both intervals are rough.

## Per-stage cost, Ancestra (CPU s / wall s / peak RSS MB)

| run | input | families | scores | merge | output |
|---|---|---|---|---|---|
| anc_8 | 0 / 0 / 80 | 989 / 154 / 268 | 702 / 100 / 629 | 3 / 3 / 453 | 1 / 1 / 420 |
| anc_16 | 0 / 0 / 129 | 1599 / 313 / 578 | 1931 / 276 / 2108 | 11 / 12 / 838 | 3 / 4 / 769 |
