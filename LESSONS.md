# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at each milestone)

- Launch independent setup work (literature, env, data) as parallel background agents at T0. Evidence: M0, three streams started within 5 min.
- Get a compute reviewer to estimate cost with real proxies (e.g. competitor HOG sizes) before building. Evidence: r1 caught that v0 progressive merging saved only ~2.5x, which saved a full implementation cycle.
- Never split benchmark families at random; split by homology component. Evidence: r2 found 5 paralogous RefOG groups that leaked across a 35/35 split.
- Read the time from scripts/clock.sh; never estimate it. Evidence: DESIGN v0 was stamped T+1h15 when it was T+0h26.
- The node is shared (load 40–90 seen against our 24-CPU allocation). Always record CPU time (user+sys) and load next to wall time; wall times from busy periods are upper bounds. Evidence: ProstT5 probe saw identical wall time at 12 and 24 threads.
