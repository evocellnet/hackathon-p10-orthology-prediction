# Runtime scaling (DESIGN §5)

8 threads per run, one run at a time, shared 24-core node (load logged per run in `runs/*/load.log`; wall times from a busy node are upper bounds, CPU-h = user+sys from `/usr/bin/time -v` is the primary metric). Fresh work dir per run. Censored runs: `> cap`.

Subsets: nested, domain-stratified (random.Random(11)), `subsets/sub{n}.{txt,nwk,json}`; guide tree = NCBI tree pruned to the subset.

| method | 8 sp (62k) | 16 sp (144k) | 32 sp (355k) | 64 sp (816k) |
|---|---|---|---|---|
| Ancestra default (accurate preset, seq-only) | 0.12 h / 0.78 CPU-h / 2.2 GB | 0.35 h / 2.16 CPU-h / 2.3 GB | 0.92 h / 6.02 CPU-h / 5.7 GB | 2.65 h / 16.30 CPU-h / 25.1 GB |
| Ancestra fast preset (NJ merge, seq-only) | 0.07 h / 0.47 CPU-h / 0.8 GB | 0.17 h / 0.99 CPU-h / 2.1 GB | not run | not run |
| Ancestra default + gene order | not run | not run | not run | not run |
| OrthoFinder 3.2 | 0.19 h / 1.22 CPU-h / 2.4 GB | 0.57 h / 3.54 CPU-h / 2.9 GB | 2.27 h / 14.43 CPU-h / 4.8 GB | not run |
| SonicParanoid2 (default) | 0.15 h / 0.77 CPU-h / 8.0 GB | 0.38 h / 2.34 CPU-h / 8.1 GB | 1.80 h / 11.11 CPU-h / 8.1 GB | not run |

Cells: wall / CPU-h / peak RSS (max over processes as reported by `time -v`).

## Log-log slope of CPU-h (uncensored points)

| method | x | n pts | slope | OLS 95% CI | bootstrap 95% CI | R² |
|---|---|---|---|---|---|---|
| Ancestra default (accurate preset, seq-only) | proteins | 4 | 1.18 | [1.14, 1.21] | [1.17, 1.18] | 1.000 |
| Ancestra default (accurate preset, seq-only) | species | 4 | 1.46 | [1.44, 1.49] | [1.46, 1.47] | 1.000 |
| Ancestra fast preset (NJ merge, seq-only) | proteins | 2 | 0.88 | n<3 | n<3 | 1.000 |
| Ancestra fast preset (NJ merge, seq-only) | species | 2 | 1.06 | n<3 | n<3 | 1.000 |
| Ancestra S0+S2 only (families+scoring) | proteins | 4 | 1.22 | [0.84, 1.61] | [1.08, 1.38] | 0.989 |
| Ancestra S0+S2 only (families+scoring) | species | 4 | 1.52 | [1.03, 2.01] | [1.35, 1.72] | 0.989 |
| OrthoFinder 3.2 | proteins | 3 | 1.42 | [0.33, 2.51] | [1.42, 1.42] | 0.996 |
| OrthoFinder 3.2 | species | 3 | 1.78 | [-0.03, 3.59] | [1.78, 1.78] | 0.994 |
| SonicParanoid2 (default) | proteins | 3 | 1.53 | [-0.00, 3.07] | [1.53, 1.53] | 0.994 |
| SonicParanoid2 (default) | species | 3 | 1.92 | [-0.48, 4.33] | [1.92, 1.92] | 0.990 |

Bootstrap: percentile, B=5000, seed 11, resampling the points (resamples with <3 distinct x skipped). With 3-5 points both intervals are rough.

## Per-stage cost, Ancestra (CPU s / wall s / peak RSS MB)

| run | input | families | scores | merge | output |
|---|---|---|---|---|---|
| anc_8_acc (anc) | 0 / 0 / 80 | 989 / 154 / 268 | 702 / 100 / 629 | 1126 / 186 / 464 | 1 / 1 / 402 |
| anc_16_acc (anc) | 0 / 0 / 129 | 1599 / 313 / 578 | 1931 / 276 / 2108 | 4220 / 664 / 755 | 4 / 4 / 703 |
| anc_32_acc (anc) | 0 / 0 / 256 | 4126 / 738 / 1733 | 7345 / 1106 / 5807 | 10194 / 1455 / 1810 | 11 / 11 / 1746 |
| anc_64 (anc) | 0 / 1 / 522 | 8914 / 1533 / 5282 | 29198 / 4632 / 25676 | 20539 / 3338 / 6954 | 32 / 33 / 4823 |
| anc_8_fast (anc_fast) | 0 / 0 / 80 | 989 / 154 / 268 | 702 / 100 / 629 | 3 / 3 / 364 | 1 / 1 / 310 |
| anc_16_fast (anc_fast) | 0 / 0 / 129 | 1599 / 313 / 578 | 1931 / 276 / 2108 | 11 / 11 / 755 | 4 / 4 / 696 |
