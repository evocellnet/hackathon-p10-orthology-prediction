# Large-scale demo (T3): Ancestra on UniProt proteome subsets

12 threads, node shared with other agents (load logged per run). Measured values only in the first table.

| run | preset | species | proteins | wall h | CPU-h | peak RSS GB | S0 h | S2 h | merge h | mean/max load1 | families | root groups |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n250_accurate | accurate | 250 | 1,082,777 | 1.91 | 17.7 | 23.4 | 0.42 | 0.61 | 0.87 | 19/22 | 190,737 | 49,879 |
| n125 | fast | 125 | 523,225 | 0.40 | 3.3 | 7.6 | 0.17 | 0.20 | 0.02 | 29/46 | 111,690 | 30,552 |
| n250 | fast | 250 | 1,082,777 | 1.12 | 8.6 | 23.4 | 0.42 | 0.63 | 0.05 | 32/45 | 190,737 | 50,362 |
| n500 | fast | 500 | 2,257,658 | 3.32 | 23.0 | 68.8 | 1.37 | 1.77 | 0.14 | 26/50 | 355,370 | 84,865 |
| n1000 | fast | 1000 | 4,751,551 | 10.60 | 64.3 | 221.2 | 4.89 | 5.17 | 0.44 | 24/57 | 619,799 | 146,079 |

## Log-log slopes vs number of proteins (fast preset, OLS, 95% t-interval)

| quantity | slope | 95% CI | points | R2 |
|---|---|---|---|---|
| wall_h | 1.49 | [1.39, 1.58] | 4 | 1.000 |
| cpu_h | 1.35 | [1.31, 1.39] | 4 | 1.000 |
| max_rss_gb | 1.52 | [1.46, 1.58] | 4 | 1.000 |
| families_wall_h | 1.54 | [1.23, 1.85] | 4 | 0.996 |
| scores_wall_h | 1.46 | [1.37, 1.55] | 4 | 1.000 |
| merge_wall_h | 1.40 | [1.30, 1.50] | 4 | 0.999 |

## EXTRAPOLATION (not measured) to 200 M proteins / ~10,000 proteomes

Power-law extrapolation from the fits above over a ~40-400x range beyond the largest measured point; the interval reflects the slope CI only and ignores model misspecification (e.g. memory is not a power law once chunked). Treat as order of magnitude.

| quantity | extrapolated | slope-CI range |
|---|---|---|
| wall_h | 2,679 h | 439 - 16,355 h |
| cpu_h | 9,804 h | 4,619 - 20,809 h |
| max_rss_gb | 65,347 GB | 20,911 - 204,207 GB |
| families_wall_h | 1,435 h | 4 - 532,457 h |
| scores_wall_h | 1,221 h | 209 - 7,142 h |
| merge_wall_h | 80 h | 11 - 574 h |

## Local slopes between the two largest measured points (n500 -> n1000, 2.10x proteins)

Report these alongside the global fit: the cost is superlinear; the local exponent for S0 (families) is higher than the global one.

| quantity | n_prev | n_last | ratio | local slope |
|---|---|---|---|---|
| wall_h | 3.32 | 10.6 | 3.19 | 1.56 |
| cpu_h | 23 | 64.3 | 2.80 | 1.38 |
| max_rss_gb | 68.8 | 221 | 3.22 | 1.57 |
| families_wall_h | 1.37 | 4.89 | 3.58 | 1.71 |
| scores_wall_h | 1.77 | 5.17 | 2.92 | 1.44 |
| merge_wall_h | 0.144 | 0.437 | 3.03 | 1.49 |
| output_wall_h | 0.0275 | 0.07 | 2.55 | 1.26 |

Caveats: timings come from a shared node (load1 mean/max in the table); S0 (families) at n1000 was 3.6x the n500 time for 2.1x proteins, i.e. the global power-law slope understates the growth at the top end, so the extrapolation below is more likely an under- than an over-estimate for time. Peak RSS is the S2 (scores) stage; a 65 TB extrapolation is not meaningful, it only says S2 would need chunking/streaming at that scale.
