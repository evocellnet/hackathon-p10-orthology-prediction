# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at M2 by a fresh-agent retrospective)

1. Launch independent setup work (literature, env, data) as parallel agents at T0. Evidence: M0, three streams started within 5 min.
2. Have a compute reviewer cost a design with real proxies before building. Evidence: r1 rejected v0 progressive merging at ~2.5× saving.
3. Split benchmarks by homology component, never at random. Evidence: r2 found 5 leaked paralogous RefOG groups.
4. Read time only from scripts/clock.sh. Evidence: DESIGN v0 stamped T+1h15 when it was T+0h26.
5. The node is shared (load 40–90). Log CPU time and load beside wall time; walls from busy periods are upper bounds; rerun the final runtime table in a quiet window. Evidence: ProstT5 identical wall at 12 and 24 threads.
6. Never edit a running script; copy it to a new version. Evidence: 2 QfO runs lost.
7. Iterate only on the fast QfO tests (~5 min); STD tests take 15 min–1.5 h each. Evidence: data-agent timings.
8. Re-check toy-tuned defaults on real dev data. Evidence: coverage guard 0.5 best on toy, −9 F on OrthoBench dev.
9. Pick a metric sensitive to the decision under test. Evidence: synteny null at OG level, +0.03–0.05 at pair level.
10. Check a stage is deterministic before tuning downstream on its cached output. Evidence: S0 thread-dependence moved pair F1 0.78→0.64.
11. Every layer needs a shuffled or negative control and a pre-registered keep rule. Evidence: expression and isoforms killed; 3Di merge gains within seed spread.
12. With <50 dev items, fix the config grid in advance and report CIs and macro metrics. Evidence: shuffled seeds span 66.9–69.8 OG-F on 34 RefOGs; RefOG021 dominates micro pair F1.
13. Get comparator numbers on the dev set as soon as they exist; a gate without a comparator is no gate. Evidence: OrthoFinder finished 17:58 but its dev F was not recorded by T+5h.
