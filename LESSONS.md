# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at M3 by a fresh-agent retrospective)

1. Launch independent setup streams (literature, env, data) as parallel agents at T0. Evidence: M0, three streams within 5 min.
2. Cost the design with real proxies before building it. Evidence: r1 rejected v0 progressive merging at ~2.5× saving.
3. Split benchmarks by homology component, never at random. Evidence: r2 found 5 leaked paralogous RefOG groups.
4. Read time only from scripts/clock.sh; the node is shared, so log CPU time and load and rerun the final runtime table in a quiet window. Evidence: v0 mis-stamped T+1h15 when it was T+0h26; ProstT5 identical wall at 12/24 threads.
5. Iterate on the fast QfO tests only; run STD tests once at the end. Evidence: STD 15 min–1.5 h each.
6. Re-check toy-tuned defaults on real dev data. Evidence: coverage guard 0.5 best on toy, −9 F on OrthoBench dev.
7. Use a metric sensitive to the decision; prefer macro/pair-level metrics when n < 50. Evidence: synteny null at OG level, +0.03–0.05 pair F1; RefOG021 dominated micro F1.
8. Check a stage is deterministic before tuning on its cached output. Evidence: S0 thread-dependence moved pair F1 0.78→0.64.
9. Give every layer a shuffled control and a pre-registered keep rule. Evidence: expression, isoforms and 3Di-in-merge were killed cleanly.
10. Fix grids in advance with CIs; treat differences inside partition noise (±2–4 OG-F) as null. Evidence: 243 OrthoBench dev runs.
11. Have comparator numbers on dev before any gate. Evidence: OrthoFinder finished 17:58 but was unrecorded until T+5h.
12. Re-check adoption on every preset, and run verify_orthoxml before every peek. Evidence: so_struct passed on fast, failed on accurate; the gate caught malformed nesting.
13. Create the dev split before the first look at a test benchmark. Evidence: the TreeFam-A dev split came after peek 1 (disclosure #1).
14. Keep long runs as detached, logged jobs, and never edit a running script. Evidence: the API limit at T+11h40 killed 4 agents but not their jobs; 2 QfO runs were lost to script edits.
15. A strong standalone signal does not predict a gain inside the merge; test the hard cases and integrate before investing in data. Evidence: L2 round, no layer passed the keep rule; domain architecture separated orthogroups with AUC 0.98 but orthologs from in-group paralogs only at 0.54, and most layers were near chance on hard cases.
16. Measure scaling at the target scale, not only on the benchmark set. Evidence: exponent 1.16 on QfO subsets became 1.35 at 0.5–4.75 M proteins, the "near-linear" S0 had local exponent 1.71, and memory (221 GB) became the limit.
17. Check that a competitor run completes at every size before the deadline night. Evidence: SonicParanoid2 at 64 species failed after 4.1 h with no time to rerun.
