# Runtime scaling (DESIGN §5)

8 threads per run, one run at a time, shared node (24 cores until 2026-10-27 11:10, 48 cores after; load logged per run in `runs/*/load.log`; wall times from a busy node are upper bounds, CPU-h = user+sys from `/usr/bin/time -v` is the primary metric). Fresh work dir per run. Censored runs: `> cap`.

**FastOMA** relies on a precomputed OMA-derived OMAmer database (LUCA.h5, 9.9 GB); its build cost is not counted (download excluded); the timed workflow is the full nextflow pipeline incl. OMAmer placement, roothog inference, HOG inference and the report notebook. OMA contains these species (leakage flag for accuracy, irrelevant for timing). FastOMA places only the proteins OMAmer can assign.

Subsets: nested, domain-stratified (random.Random(11)), `subsets/sub{n}.{txt,nwk,json}`; guide tree = NCBI tree pruned to the subset.

| method | 8 sp (62k) | 16 sp (144k) | 32 sp (355k) | 64 sp (816k) | 78 sp (989k) |
|---|---|---|---|---|---|
| Ancestra default (accurate preset, seq-only) | 0.12 h / 0.78 CPU-h / 2.2 GB | 0.35 h / 2.16 CPU-h / 2.3 GB | 0.92 h / 6.02 CPU-h / 5.7 GB | 2.65 h / 16.30 CPU-h / 25.1 GB | 2.83 h / 18.54 CPU-h / 28.9 GB |
| Ancestra fast preset (NJ merge, seq-only) | 0.07 h / 0.47 CPU-h / 0.8 GB | 0.17 h / 0.99 CPU-h / 2.1 GB | 0.53 h / 3.20 CPU-h / 5.7 GB | 1.75 h / 10.63 CPU-h / 25.1 GB | 1.88 h / 11.63 CPU-h / 28.9 GB |
| Ancestra default + gene order | not run | not run | 0.91 h / 6.14 CPU-h / 5.6 GB | not run | not run |
| OrthoFinder 3.2 | 0.19 h / 1.22 CPU-h / 2.4 GB | 0.57 h / 3.54 CPU-h / 2.9 GB | 2.27 h / 14.43 CPU-h / 4.8 GB | > 6 h wall (45.9 CPU-h used, 5.3 GB) | not run |
| SonicParanoid2 (default) | 0.15 h / 0.77 CPU-h / 8.0 GB | 0.38 h / 2.34 CPU-h / 8.1 GB | 1.80 h / 11.11 CPU-h / 8.1 GB | FAILED rc=251 | not run |
| FastOMA 0.5.1 (OMAmer LUCA, DB build not counted) | 0.37 h / 1.83 CPU-h / 13.6 GB | 0.76 h / 5.27 CPU-h / 14.9 GB | 2.70 h / 15.73 CPU-h / 20.8 GB | 5.45 h / 36.64 CPU-h / 24.6 GB | not run |

Cells: wall / CPU-h / peak RSS (max over processes as reported by `time -v`).

## Log-log slope of CPU-h (uncensored points)

| method | x | n pts | slope | OLS 95% CI | bootstrap 95% CI | R² |
|---|---|---|---|---|---|---|
| Ancestra default (accurate preset, seq-only) | proteins | 5 | 1.16 | [1.10, 1.21] | [1.11, 1.18] | 0.999 |
| Ancestra default (accurate preset, seq-only) | species | 5 | 1.42 | [1.30, 1.53] | [1.29, 1.47] | 0.998 |
| Ancestra fast preset (NJ merge, seq-only) | proteins | 5 | 1.21 | [1.02, 1.41] | [1.10, 1.37] | 0.993 |
| Ancestra fast preset (NJ merge, seq-only) | species | 5 | 1.49 | [1.23, 1.74] | [1.38, 1.72] | 0.991 |
| Ancestra S0+S2 only (families+scoring) | proteins | 5 | 1.21 | [1.02, 1.41] | [1.10, 1.37] | 0.992 |
| Ancestra S0+S2 only (families+scoring) | species | 5 | 1.49 | [1.23, 1.74] | [1.38, 1.72] | 0.991 |
| OrthoFinder 3.2 | proteins | 3 | 1.42 | [0.33, 2.51] | n<4 | 0.996 |
| OrthoFinder 3.2 | species | 3 | 1.78 | [-0.03, 3.59] | n<4 | 0.994 |
| SonicParanoid2 (default) | proteins | 3 | 1.53 | [-0.00, 3.07] | n<4 | 0.994 |
| SonicParanoid2 (default) | species | 3 | 1.92 | [-0.48, 4.33] | n<4 | 0.990 |
| FastOMA 0.5.1 (OMAmer LUCA, DB build not counted) | proteins | 4 | 1.17 | [1.00, 1.34] | [1.11, 1.24] | 0.998 |
| FastOMA 0.5.1 (OMAmer LUCA, DB build not counted) | species | 4 | 1.45 | [1.21, 1.70] | [1.38, 1.55] | 0.997 |

OLS CI: t-interval, df = n-2 (primary; very wide for 3 points). Bootstrap: percentile, B=5000, seed 11, resampling the points, only for n>=4 (resamples with <3 distinct x skipped); with 4-5 points it is optimistic, read the OLS CI first. The CIs capture only fit scatter around a power law, not run-to-run or node-load variation (one run per point).

## Extrapolated CPU-h where no completed run exists (power-law point estimate, NOT a measurement)

| method | size | predicted CPU-h |
|---|---|---|
| OrthoFinder 3.2 | 64 sp (816k) | 45 |
| OrthoFinder 3.2 | 78 sp (989k) | 59 |
| SonicParanoid2 (default) | 78 sp (989k) | 51 |
| FastOMA 0.5.1 (OMAmer LUCA, DB build not counted) | 78 sp (989k) | 49 |

## Per-stage cost, Ancestra (CPU s / wall s / peak RSS MB)

| run | input | families | scores | merge | output |
|---|---|---|---|---|---|
| anc_8_acc (anc) | 0 / 0 / 80 | 989 / 154 / 268 | 702 / 100 / 629 | 1126 / 186 / 464 | 1 / 1 / 402 |
| anc_16_acc (anc) | 0 / 0 / 129 | 1599 / 313 / 578 | 1931 / 276 / 2108 | 4220 / 664 / 755 | 4 / 4 / 703 |
| anc_32_acc (anc) | 0 / 0 / 256 | 4126 / 738 / 1733 | 7345 / 1106 / 5807 | 10194 / 1455 / 1810 | 11 / 11 / 1746 |
| anc_64 (anc) | 0 / 1 / 522 | 8914 / 1533 / 5282 | 29198 / 4632 / 25676 | 20539 / 3338 / 6954 | 32 / 33 / 4823 |
| anc_78 (anc) | 1 / 1 / 608 | 10707 / 1807 / 6410 | 30928 / 4760 / 29635 | 25040 / 3558 / 9093 | 41 / 42 / 5753 |
| anc_8_fast (anc_fast) | 0 / 0 / 80 | 989 / 154 / 268 | 702 / 100 / 629 | 3 / 3 / 364 | 1 / 1 / 310 |
| anc_16_fast (anc_fast) | 0 / 0 / 129 | 1599 / 313 / 578 | 1931 / 276 / 2108 | 11 / 11 / 755 | 4 / 4 / 696 |
| anc_32_fast (anc_fast) | 0 / 0 / 256 | 4126 / 738 / 1733 | 7345 / 1106 / 5807 | 33 / 33 / 1835 | 11 / 11 / 1748 |
| anc_64_fast (anc_fast) | 0 / 1 / 522 | 8914 / 1533 / 5282 | 29198 / 4632 / 25676 | 103 / 104 / 6282 | 30 / 31 / 4564 |
| anc_78_fast (anc_fast) | 1 / 1 / 608 | 10707 / 1807 / 6410 | 30928 / 4760 / 29635 | 150 / 151 / 8524 | 42 / 42 / 5687 |
| anc_32_go_acc (anc_go) | 0 / 0 / 256 | 4210 / 658 / 1728 | 7384 / 1125 / 5778 | 10480 / 1482 / 2222 | 12 / 12 / 1858 |

## Status and deviations (resumed 2026-10-27 11:10 CET on a 48-core node; lane A finished 2026-10-28 06:12)

- **Run status (resumed):** FastOMA 16/32/64 completed (5.3 / 15.7 / 36.6 CPU-h). OrthoFinder 64 hit the 6 h cap: censored, > 6 h wall, 45.9 CPU-h used (lower bound). SonicParanoid2 64 FAILED (rc=251) after 4.12 h wall / 36.5 CPU-h: after the 874 s MMseqs profile search it aborted with `The directory with the architectures contains 37 files when 64 are expected` (pfam profile-DB errors in the log, `runs/sp_64/cmd_tail.txt`); the 8/16/32 runs worked, cause not diagnosed, no time to rerun; its CPU-h is a lower bound for the cost, not a completed measurement. Neither point enters the fits. Rows below with 8 threads, 48-core node shared with other agents (load 20-45).
- **Fits:** OrthoFinder and SonicParanoid2 have 3 points (very wide OLS CIs); FastOMA now has 4 points and a fit; Ancestra 5.
- **Merge preset drift.** `ancestra/pipeline.py` and `merge.py` were edited during the benchmark (presets `accurate` = default, FastTree local trees; `fast` = NJ). anc_8/16/32 (and 32_go) ran first with the pre-preset NJ merge; anc_64 and anc_78 ran with `accurate`. For a consistent series, merge+output were re-measured with `--preset accurate|fast` on the cached S0+S2 of the fresh base run (hard links, `anc_*_acc`, `anc_*_fast`); an entry = base input/families/scores (fresh) + re-measured merge/output (`composite` column in scaling.tsv). anc_64/anc_78 base runs used the code of 02:45 (anc_64) / the frozen snapshot `work/scaling/code_snapshot` (anc_78 and all remerges); S0/S2 code was identical (same stage keys). Gene-order tables were prebuilt (build cost not counted).
- **Accurate-preset merge is expensive:** merge = 1126 / 4220 / 10194 / 20539 / 25040 CPU-s at 8/16/32/64/78 species vs 3 / 11 / 33 / 103 / 150 CPU-s for the fast preset (see per-stage table). Gene order adds about 0.1 CPU-h at 32 species with the accurate preset (6.14 vs 6.02 CPU-h) and about 50 s of merge wall with the fast/NJ merge (86 vs 34 s).
- **Shared node:** load average 25-70 on 24 cores (another agent's QfO pipeline + competitor runs); walls are upper bounds, CPU-h is the primary metric but is also inflated by contention (cache, hyperthreads). One run per point.
- FastOMA: first launches (8/16/32) ended after check_input in seconds because `proteome/` was a symlinked directory (nextflow listFiles() sees nothing); fixed with a real directory of file symlinks. Those bogus results are in `work/scaling/failed_attempts/`. FastOMA species tree = pruned NCBI tree with spaces replaced by `_` (`subsets/sub{n}.fastoma.nwk`), polytomies kept.
