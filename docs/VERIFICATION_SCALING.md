# Independent verification: runtime and scaling claims (commit b51138b)

Verifier: an agent that did not produce the scaling runs or the text. Date 2026-10-28. No metric, test or result file was changed.

**Scope.** `manuscript/body.tex` (Runtime and scaling: Table `tab:scaling`, the paragraph after it, Table `tab:large`, Fig `fig:scaling` caption; abstract; author summary; Discussion paragraph 1; Limitations; Methods "Scaling"), `REPORT.md` (Abstract scaling sentences, §3.5, §3.6, Discussion item 1, Limitations runtime bullet), `README.md` "Status (final…)", `docs/ONE_PAGE_SUMMARY.md` "Scaling (final)", `docs/figures/fig3_scaling.png` / `scripts/fig3_scaling.py`, and a search for remaining "near-linear" claims.

**Method.** I recomputed everything from the raw files, not from `scaling.tsv` or `summary.md`. CPU-h = (User + System) / 3600 and peak RSS = "Maximum resident set size" from each `/usr/bin/time -v` file (`results/scaling/runs/*/time.txt`, `work/scale_large/<tag>/logs/time_v.txt`). Per-stage times came from `stage_accounting.json` (QfO) and `logs/run.log` plus `out/stages/*/stage.json` (large runs). Ancestra QfO composites were rebuilt my own way: base run CPU − base merge/output CPU + remerge run CPU. A pure sum of stage CPU agrees to within 0.003 CPU-h. Fits are OLS of log10(y) on log10(proteins) with a 95% t-interval (df = n − 2), using conda_envs/orth/bin/python (numpy/scipy). The scripts are in the session scratchpad and were not committed.

**Result.** I checked 52 claims. Every number agrees with the raw logs to the precision printed. Problems are in interpretation and wording: 1 mismatch (the extrapolation's "~10,000 proteomes"), 1 misleading statement ("we did not profile why", although the sub-step timings were already logged), 6 overstatements or comparisons that are not like-for-like, 3 figure defects and 4 minor wording or unit issues. They are listed in the fixes below.

## Claim table

| # | Location | Stated | Recomputed (raw) | Verdict |
|---|---|---|---|---|
| 1 | body Tab scaling / REPORT T5, Ancestra accurate 8/16/32/64/78 sp | 0.78 / 2.16 / 6.02 / 16.3 / 18.5 | 0.784 / 2.155 / 6.023 / 16.304 / 18.536 | OK |
| 2 | same, Ancestra fast | 0.47 / 0.99 / 3.20 / 10.6 / 11.6 | 0.472 / 0.986 / 3.201 / 10.629 / 11.625 | OK |
| 3 | same, OrthoFinder 8/16/32 | 1.22 / 3.54 / 14.4 | 1.225 / 3.543 / 14.434 | OK |
| 4 | same, OrthoFinder 64 censored | > 45.9, stopped at 6 h cap | exit 124 (`timeout`), wall 6.0017 h (cap 21,600 s), 45.865 CPU-h; cmd_tail "stopped after 21604 s" | OK |
| 5 | same, SonicParanoid2 8/16/32 | 0.77 / 2.34 / 11.1 | 0.774 / 2.336 / 11.112 | OK |
| 6 | same, SonicParanoid2 64 failed | rc=251, 4.1 h, 36.5 CPU-h, "architectures contains 37 files when 64 are expected" | exit 251, wall 4.119 h, 36.521 CPU-h; message present in `cmd_tail.txt`, after the 874 s MMseqs profile search | OK |
| 7 | same, FastOMA 8/16/32/64 | 1.83 / 5.27 / 15.7 / 36.6 | 1.830 / 5.273 / 15.728 / 36.640 | OK |
| 8 | Slope Ancestra accurate | 1.16 [1.10, 1.21] | 1.155 [1.098, 1.213], n = 5, R² 0.999 | OK |
| 9 | Slope Ancestra fast | 1.21 [1.02, 1.41] | 1.212 [1.019, 1.406], n = 5 | OK |
| 10 | Slope OrthoFinder (3 pts) | 1.42 | 1.420 [0.33, 2.51] | OK |
| 11 | Slope SonicParanoid2 (3 pts) | 1.53 | 1.535 [−0.00, 3.07] | OK |
| 12 | Slope FastOMA | 1.17 [1.00, 1.34] | 1.171 [1.003, 1.339], n = 4 | OK (but mixes two nodes, see #47) |
| 13 | body ¶ / REPORT §3.5: at 64 sp, Ancestra acc 16.3, FastOMA 36.6, OF > 45.9 | as stated | 16.30 / 36.64 / > 45.87 | OK (cross-node, see #47) |
| 14 | body ¶ / REPORT §3.5: QfO 78 accurate 18.5 CPU-h, 29 GB; S2 8.6, merge 7.0, S0 3.0 | as stated | 18.54 CPU-h, 28.94 GiB; scores 8.59, merge 6.96, families 2.97 CPU-h | OK |
| 15 | fast merge 2.5 CPU-min | 2.5 | 2.50 min (anc_78_fast merge 150.2 CPU-s) | OK |
| 16 | REPORT T5 caption: peak RSS Ancestra 2–29 GB | 2–29 | accurate series 2.24–28.94 GiB (fast composites 0.8–28.9) | OK |
| 17 | Tab large / REPORT T6, proteins 125/250/500/1000 | 523,225 / 1,082,777 / 2,257,658 / 4,751,551 | identical (`subset_totals.tsv`; re-summed from `selection.tsv`) | OK |
| 18 | wall h 125/250/500/1000 | 0.40 / 1.12 / 3.32 / 10.6 | 0.398 / 1.120 / 3.319 / 10.596 | OK |
| 19 | CPU-h 125/250/500/1000 | 3.3 / 8.6 / 23.0 / 64.3 | 3.281 / 8.624 / 22.979 / 64.305 | OK |
| 20 | peak memory 125/250/500/1000 | 7.6 / 23.4 / 68.8 / 221 | 7.561 / 23.426 / 68.782 / 221.155 GiB (= 237.5 GB decimal) | OK (unit, see fix 13) |
| 21 | Root HOGs | 30,552 / 50,362 / 84,865 / 146,079 | `merge/stage.json` n_root_groups identical | OK |
| 22 | REPORT T6 S0 h | 0.17 / 0.42 / 1.37 / 4.89 | 0.167 / 0.420 / 1.367 / 4.893 | OK |
| 23 | REPORT T6 S2 h | 0.20 / 0.63 / 1.77 / 5.17 | 0.204 / 0.631 / 1.768 / 5.165 | OK |
| 24 | REPORT T6 merge h | 0.02 / 0.05 / 0.14 / 0.44 | 0.020 / 0.052 / 0.144 / 0.437 | OK |
| 25 | REPORT T6 n250 accurate row | 1.91 h, 17.7 CPU-h, 23.4 GB, S0 0.42, S2 0.61, merge 0.87, 49,879 | 1.914, 17.657, 23.416, 0.417, 0.614, 0.866, 49,879 | OK |
| 26 | Large CPU-h exponent | 1.35 [1.31, 1.39] | 1.347 [1.308, 1.387], R² 0.9999 | OK |
| 27 | Large wall exponent | 1.49 [1.39, 1.58] | 1.486 [1.392, 1.581] | OK |
| 28 | Large memory exponent | 1.52 [1.46, 1.58] | 1.524 [1.464, 1.583] | OK |
| 29 | Families wall global | 1.54 [1.23, 1.85] | 1.539 [1.229, 1.848] | OK |
| 30 | Families n500→n1000 | 3.6×, 2.1× proteins, local 1.71 | 3.58×, 2.105×, 1.71 | OK |
| 31 | "family step grew fastest" | fastest | wall-time local slopes n500→n1000: families 1.71 > merge 1.49 > scores 1.44 > CPU 1.38. Global: families 1.54 [1.23, 1.85] vs scores 1.46 [1.37, 1.55] (CIs overlap). Wall only, since per-stage CPU was not logged for the large runs | OK as a point estimate; say "wall time" (fix 11) |
| 32 | "We did not profile why" / "did not profile which sub-step" | not known | The S0 sub-step timings are already logged (`out/stages/families/stage.json`, `s0_timings`) and identify the cause (fix 2) | MISLEADING |
| 33 | Peak memory in the scoring stage, which holds all scores at once | S2 | n1000 S2 VmHWM of the Python process itself was 226,463 MB = 221 GiB (other stages ≤ 114 GiB). The process-tree monitor peaks at 218 GB at 00:43–01:02, inside S2 | OK (location verified; "holds all scores" is consistent with the in-memory `ScoreTable` but was not measured directly) |
| 34 | accurate vs fast n250 | 17.7 vs 8.6 CPU-h, "twice", mostly merge (0.87 vs 0.05 h) | 17.66 / 8.62 = 2.05×; merge wall 0.866 vs 0.052 h; S0 and S2 walls equal (0.42/0.42, 0.61/0.63); total wall difference 0.79 h ≈ merge wall difference 0.81 h | OK (the hours are wall time, not CPU; fix 12) |
| 35 | "All five outputs passed OrthoXML validation" | 5/5 | `verify_n{125,250,250_accurate,500,1000}.json`: `ok: true`, all 8 checks true; logs "overall: PASS" (no separate log for n250_accurate, JSON only) | OK |
| 36 | Selection 800/100/100, one per genus, UniProt 2026_03, seed 11, nested | as stated | `selection.tsv`: 800/100/100. No genus repeated within a domain (the only shared label is the placeholder "candidate", once in Bacteria and once in Archaea). Subsets nested, with 100/12/13, 200/25/25, 400/50/50, 800/100/100. Release 2026_03 per DATA.md §15 | OK |
| 37 | Methods: guide trees from NCBI taxonomy; competitors capped at 6 h | as stated | `trees/ncbi_n*.nwk`; competitor `cap_s` 21,600 | OK |
| 38 | Abstract (body): 1.16 up to 1 M, then 1.35 from 0.5 to 4.75 M; "64 CPU-hours but 221 GB" | read as "exponent rose" | 1.16 is the **accurate** preset on QfO; 1.35 is the **fast** preset on UniProt. Like-for-like (fast): QfO 1.21 [1.02, 1.41] vs 1.35 [1.31, 1.39]. The CIs overlap, and data, threads (8 vs 12) and node all differ | OVERSTATED (fix 1) |
| 39 | body Results: "At this scale the cost grew faster than on QfO" | faster | same issue as #38 | OVERSTATED (fix 1) |
| 40 | body Discussion: "exponent rose from 1.16 on QfO to 1.35" | rose | same issue | OVERSTATED (fix 1) |
| 41 | body Discussion: "memory … is the first limit to larger inputs" | first limit | Inference, not a measurement. On the 500 GB node the memory fit (1.52) predicts exhaustion at about 8 M proteins (about 1,700 proteomes of this composition), where the wall fit predicts about 24 h, so the inference is plausible | OVERSTATED (minor, fix 6) |
| 42 | body author summary: "ran on 1,000 proteomes on one computer" | one computer | one node, 12 threads, but 221 GiB of RAM. "One computer" suggests ordinary hardware | OVERSTATED (minor, fix 7) |
| 43 | body Limitations: no competitor beyond 64 species; OF cap, SP2 failed | as stated | true for runtime. The QfO-78 competitor results are published OpenEBench submissions, not local runs | OK |
| 44 | REPORT abstract: "On 125–1,000 … the exponent rises to 1.35" | rises | same issue as #38 | OVERSTATED (fix 1) |
| 45 | REPORT §3.6 bullet: "These are higher than the 1.16 measured on QfO subsets" | higher | same issue as #38 | OVERSTATED (fix 1) |
| 46 | REPORT Discussion 1: "dominated by the family step and within-family alignment … clearly faster … memory the first limit" | — | For `accurate` on QfO the merge is 7.0 of 18.5 CPU-h (38%), so S0+S2 dominate only for `fast`. "Clearly faster" fails like-for-like (#38) | OVERSTATED (fix 8) |
| 47 | body Tab scaling footnote and Fig caption: "Single runs on a shared node" | one node | Two machines, according to the uptime counters in `load.log`. Machine 1 (uptime 60–61 d) ran all Ancestra runs, OF and SP2 8–32, and FastOMA 8. Machine 2 (66 d) ran FastOMA 16/32/64, OF 64, SP2 64 and all large-scale runs. `os.cpu_count()` reported 48 on both | MISMATCH (wording, fix 3) |
| 48 | REPORT §3.6 extrapolation: 200 M proteins "(~10,000 proteomes)", 9,800 CPU-h, "40× beyond" | ~10,000 proteomes | 9,804 CPU-h and 42× check out. But 200 M proteins at this set's mean of 4,752 proteins per proteome is **about 42,000 proteomes**. "~10,000" only holds for eukaryote-sized proteomes (~20,000 proteins) | MISMATCH (fix 5) |
| 49 | REPORT Limitations: one run per point, shared node (load 20–70); OF/SP2 3 points; no competitor on 125–1,000 | as stated | load1 means 18.5–32 and maxima up to 73.6 over all runs. Two nodes, not one | OK, but add the two nodes (fix 3) |
| 50 | README Status scaling bullet | 1.16 up to 1 M; 10.6 h / 64 CPU-h / 221 GB; "steeper exponent (1.35)" | numbers OK; "steeper" compares accurate with fast | OVERSTATED (fix 1) |
| 51 | ONE_PAGE "Scaling (final)" | FastOMA 1.17, Ancestra 1.16; 10.6 h, 64 CPU-h, 221 GB, 1.35 | numbers OK; the preset is not named for 1.35 | minor (fix 1) |
| 52 | Fig 3 (png and script) | caption: open triangle (lower bound), cross for failed; B memory on same axis | Data points match the TSVs, which match the raw logs. Defects: (a) the upward arrow for the OF censored point is **not drawn**: its head (y = 73.4) is above the axis limit (57.7), so matplotlib clips it (`_check_xy` = False); (b) the "stopped at 6 h cap" text extends above the panel frame; (c) the SP2 "failed" cross at 36.5 sits on top of the FastOMA 64-species marker at 36.6 and cannot be told apart; (d) panel B has no y-axis label although it shows CPU-h and GB; (e) the "1000 proteomes" label runs past the right frame | Figure defects (fix 4) |

Remaining "near-linear" or "linear" wording (search over body.tex, REPORT.md, README.md, ONE_PAGE_SUMMARY.md, CLAIMS.md):
- body.tex: none left. The abstract's "first splits all proteins into families with a clustering step" and the intro's "fast clustering step" are fine.
- REPORT.md:48 (design history): "Design v1/v2 therefore added a **near-linear family partition first**". This now contradicts §3.6 (fix 9).
- REPORT.md:264: "S0 was designed to be near-linear; at 4.75 M proteins it is not". Correct as written.
- docs/ONE_PAGE_SUMMARY.md:20: "On nested QfO subsets, CPU grows about as proteins^1.16". Not wrong, but needs the range (fix 10).
- README.md and CLAIMS.md: none. CLAIMS C-SCALE-1 and C-SCALE-2 are still `PENDING-VERIFICATION` (fix 14).
- Out of scope but stale: DESIGN.md lines 8, 27 and 110 ("Near-linear clustering", "Reference-free near-linear family partitioning"). These are design-time statements and could carry a pointer to REPORT §3.6.

## Required fixes

1. **Compare like with like (preset).** 1.16 is the `accurate` preset on QfO and 1.35 the `fast` preset on UniProt. The fast QfO slope 1.21 [1.02, 1.41] overlaps 1.35 [1.31, 1.39]. Do not say that the exponent "rose".
   - body.tex abstract: replace "CPU time grew with the number of proteins to the power 1.16 \CI{1.10}{1.21} up to one million proteins, and to the power 1.35 \CI{1.31}{1.39} from 0.5 to 4.75 million; 1,000 proteomes took 64 CPU-hours but 221 GB of memory." with
     "Up to one million proteins, CPU time grew with the number of proteins to the power 1.16 \CI{1.10}{1.21}. On 0.5 to 4.75 million UniProt proteins, the \texttt{fast} preset grew to the power 1.35 \CI{1.31}{1.39}, and 1,000 proteomes took 64 CPU-hours but 221 GB of memory."
   - body.tex Results: replace "At this scale the cost grew faster than on QfO: the CPU-time exponent was 1.35 \CI{1.31}{1.39} and the memory exponent 1.52." with
     "The CPU-time exponent was 1.35 \CI{1.31}{1.39}, compared with 1.21 \CI{1.02}{1.41} for the same preset on QfO subsets (8 threads, different data and node), so we cannot say that the exponent increased. Wall time (1.49) and memory (1.52) grew faster than CPU time."
   - body.tex Discussion: replace "The cost is not linear: its exponent rose from 1.16 on QfO to 1.35 between 0.5 and 4.75 million proteins, the family step grew fastest, and memory (221 GB at 1,000 proteomes) is the first limit to larger inputs." with
     "The cost is not linear. The CPU-time exponent was 1.16 on QfO (\texttt{accurate}) and 1.35 between 0.5 and 4.75 million proteins (\texttt{fast}; 1.21 on QfO). The graph step of the family partition grew almost quadratically. Memory (221 GB at 1,000 proteomes) is the likely first limit to larger inputs."
   - REPORT.md abstract: replace "On 125–1,000 UniProt reference proteomes (0.5–4.75 M proteins, `fast` preset) the exponent rises to 1.35 [1.31, 1.39];" with "On 125–1,000 UniProt reference proteomes (0.5–4.75 M proteins) the `fast` preset's exponent is 1.35 [1.31, 1.39] (1.21 [1.02, 1.41] for `fast` on QfO subsets);".
   - REPORT.md §3.6 bullet 1: replace "These are higher than the 1.16 measured on QfO subsets." with "For the same preset on QfO subsets the CPU-h slope was 1.21 [1.02, 1.41]. The intervals overlap, and data, thread count and node differ, so a rise in the exponent is not established. The cost is clearly superlinear at this scale."
   - README.md Status: replace "with a steeper exponent (1.35) at that scale" with "CPU-time exponent 1.35 [1.31, 1.39] for the `fast` preset at that scale (1.21 for `fast` on QfO subsets)".
   - ONE_PAGE_SUMMARY.md "Scaling (final)": replace "exponent 1.35 at that scale" with "`fast`-preset exponent 1.35 at that scale (1.21 on QfO subsets)".

2. **Remove "We did not profile why"; the logged S0 sub-step timings already show the cause.** From `work/scale_large/n*/out/stages/families/stage.json` (`s0_timings`), wall seconds at n125/n250/n500/n1000:
   - DIAMOND clustering (stage1): 292 / 568 / 1,134 / 2,546. Local exponents 0.91, 0.94, 1.09, so linear.
   - Representative search (stage2_search): 175 / 411 / 1,129 / 3,159. Local exponents 1.17–1.38.
   - Post-search graph step (stage2_s − stage2_search_s: hit filtering, lifting edges to clusters, connected components, Leiden split of oversized components, reassignment): 125 / 510 / 2,606 / 11,784. Local exponents 1.94, 2.22, 2.03. Its share of S0 rose from 21% to 67%.
   - On QfO, S0 wall scaled as proteins^0.90 [0.86, 0.93].
   - body.tex: replace "The family step, which we designed to be near-linear, grew fastest. Its wall time rose 3.6-fold from 500 to 1,000 proteomes, for 2.1 times more proteins (local exponent 1.71). We did not profile why." with
     "The wall time of the family step grew fastest. It rose 3.6-fold from 500 to 1,000 proteomes for 2.1 times more proteins (local exponent 1.71), although on QfO it had grown sublinearly (0.90). The DIAMOND clustering inside it stayed linear. The cause was the graph step that follows the second search (edge filtering, connected components and Leiden splitting of oversized components), whose time grew roughly with the square of the number of proteins (local exponents 1.9--2.2). It took 67\% of the family step at 1,000 proteomes."
   - REPORT.md §3.6 bullet 2: replace "We did not profile which of its sub-steps is responsible." with "The logged sub-step timings show the cause: DIAMOND clustering stays linear (local exponents 0.9–1.1), while the post-search graph step (edge filtering, components, Leiden split, reassignment) grows about quadratically (local exponents 1.9–2.2: 2,606 → 11,784 s from n500 to n1000, 67% of S0 at n1000)."
   - Optionally, add a sentence to Discussion item 1 or Follow-ups: "rewrite the S0 graph step so that it scales linearly".

3. **Two machines, not one shared node** (QfO scaling).
   - body.tex Table `tab:scaling` footnote: replace "Single runs on a shared node; CPU time is the primary measure." with "Single runs on two shared nodes: all Ancestra runs, OrthoFinder and SonicParanoid2 at 8--32 species and FastOMA at 8 species on the first; FastOMA at 16--64 species and OrthoFinder and SonicParanoid2 at 64 species on the second. CPU time is the primary measure."
   - Fig `fig:scaling` caption: replace "Measured points only; single runs on a shared node." with "Measured points only; single runs on shared nodes (panel A spans two nodes; see Table~\ref{tab:scaling})."
   - REPORT.md Limitations: replace "Runtimes were measured once per point on a shared node (load 20–70)" with "Runtimes were measured once per point on shared nodes (load 20–70). The QfO scaling series spans two nodes: FastOMA 16–64, OrthoFinder 64 and SonicParanoid2 64 ran on the second node, together with all 125–1,000-proteome runs."
   - REPORT.md T5 caption: add "Two nodes (see Limitations)."

4. **Fig 3 (`scripts/fig3_scaling.py`).**
   - (a) Draw the lower-bound arrow: add `a.set_ylim(0.35, 150)` before plotting (or pass `annotation_clip=False` to `a.annotate`).
   - (b) Keep the label inside the axes: `a.text(x * 0.55, y * 1.9, "stopped at 6 h cap", fontsize=7, color=col)`.
   - (c) Separate the SonicParanoid2 failure from the FastOMA marker, for example by plotting the cross at `x * 1.08` with a short note, or by stating in the caption "the SonicParanoid2 cross coincides with the FastOMA point".
   - (d) `b.set_ylabel("CPU-hours / peak memory (GB)")`.
   - (e) Raise `b.set_xlim(4e5, 9e6)` or put the n1000 label on the left (`xytext=(-60, 6)`).
   - Regenerate the PNG and the PDF.

5. **Extrapolation size (REPORT §3.6 last bullet).** Replace "to 200 M proteins (~10,000 proteomes)" with "to 200 M proteins (about 42,000 proteomes of this set's average size, or about 10,000 eukaryote-sized proteomes)". The same edit applies to the `results/scaling_large/summary.md` heading if regenerated. I did not change that file.

6. **Memory as "first limit" (body Discussion; REPORT §3.6 bullet 3; REPORT Discussion 1).** Mark it as an inference. In REPORT §3.6 replace "**Memory is the first limit.**" with "**Memory is the likely first limit.** On our 500 GB node, the memory fit (1.52) predicts that RAM runs out at about 8 M proteins (about 1,700 proteomes like these), when wall time would still be about a day (extrapolated, not measured)." In body.tex use "the likely first limit" (already included in the fix 1 Discussion text).

7. **Author summary "one computer".** Replace "ran on 1,000 proteomes on one computer" with "ran on 1,000 proteomes on a single server node (12 threads, 221 GB of memory)".

8. **REPORT Discussion item 1.** Replace the whole item with "1. with the `fast` preset the cost is dominated by the family step and within-family alignment; with `accurate` the merge adds a large share (7.0 of 18.5 CPU-h on QfO). CPU time grows slightly faster than linearly up to 1 M proteins (exponent 1.16 [1.10, 1.21], `accurate`; 1.21 [1.02, 1.41], `fast`) and as 1.35 [1.31, 1.39] (`fast`) between 0.5 and 4.75 M proteins. The family step's graph sub-step grows about quadratically at the top end, and memory is the likely first limit;".

9. **REPORT.md:48 (design history).** Replace "added a **near-linear family partition first**" with "added a **family partition first, designed to be near-linear** (sublinear on QfO subsets; superlinear above about 2 M proteins, §3.6)".

10. **ONE_PAGE_SUMMARY.md:20.** Replace "On nested QfO subsets, CPU grows about as proteins^1.16" with "On nested QfO subsets (up to 1 M proteins), CPU grows about as proteins^1.16".

11. **"Grew fastest" refers to wall time.** No per-stage CPU was logged for the large runs. The fix 2 text already says "The wall time of the family step grew fastest". In REPORT §3.6 bullet 2, the heading "The family step (S0) grows fastest at the top" can stay, since the next sentence names wall time.

12. **accurate vs fast at n250.** body.tex: replace "mostly in the merge (0.87 vs 0.05 hours)" with "mostly in the merge (merge wall time 0.87 vs 0.05 hours; the other stages took the same time)". REPORT: replace "almost all in the merge (0.87 vs 0.05 h)" with "almost all in the merge (merge wall 0.87 vs 0.05 h; S0 and S2 walls equal)".

13. **Memory units (optional).** All "GB" values are GiB (2^30 bytes) from `/usr/bin/time` kB / 1024². 221 GiB = 237 GB. Add to the `tab:large` footnote: "Memory in GiB (2$^{30}$ bytes)", or convert the values.

14. **CLAIMS.md.** After the fixes above, set C-SCALE-1 and C-SCALE-2 to VERIFIED-WITH-CAVEAT, citing this file. Caveats: two nodes (C-SCALE-1), and the preset mismatch in the 1.16 vs 1.35 comparison (C-SCALE-2). Note that `CLAIMS.md` currently has uncommitted edits by another agent.

15. Mirror every body.tex fix in `manuscript/ancestra_manuscript.tex`, which contains a copy of the same text (for example "We did not profile why" appears in both), and rebuild the PDF.

Minor, outside the strict scope: the REPORT abstract Conclusion says that Ancestra "used less CPU time than OrthoFinder and SonicParanoid2". At 8 species, SonicParanoid2 used 0.774 CPU-h and Ancestra `accurate` 0.784 (equal within noise); `fast` was cheaper at every size.
