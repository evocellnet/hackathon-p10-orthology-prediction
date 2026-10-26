# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at each milestone)

- Launch independent setup work (literature, env, data) as parallel background agents at T0. Evidence: M0, three streams started within 5 min.
- Get a compute reviewer to estimate cost with real proxies (e.g. competitor HOG sizes) before building. Evidence: r1 caught that v0 progressive merging saved only ~2.5x, which saved a full implementation cycle.
- Never split benchmark families at random; split by homology component. Evidence: r2 found 5 paralogous RefOG groups that leaked across a 35/35 split.
- Read the time from scripts/clock.sh; never estimate it. Evidence: DESIGN v0 was stamped T+1h15 when it was T+0h26.
- The node is shared (load 40–90 seen against our 24-CPU allocation). Always record CPU time (user+sys) and load next to wall time; wall times from busy periods are upper bounds. Evidence: ProstT5 probe saw identical wall time at 12 and 24 threads.
- Never edit a running bash script; copy it to a new version instead. Evidence: the data agent lost 2 QfO runs this way.
- For iteration, use only the fast QfO tests (GO EC VGNC SwissTrees TreeFam-A, ~5 min); STD tests take 15 min–1.5 h each. Evidence: data agent timings.
- The toy simulator lacks fragments/bad gene models and 2R-like WGD, so check toy-tuned defaults on real dev data before trusting them. Evidence: coverage guard 0.5 was best on toy and −9 F on OrthoBench dev.
- Pick metrics that are sensitive to the decision being tested. Evidence: synteny changes which paralog is chosen as the ortholog, so it was null at OrthoBench OG level but +0.03–0.05 pair F1 on the dev pair metric.
