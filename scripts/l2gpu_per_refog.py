#!/usr/bin/env python
"""Per-dev-RefOG OG-level TP/FP/FN (normalised by n-1, as benchmark.py) for several runs.
    python scripts/l2gpu_per_refog.py RUN [RUN ...]"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from orthobench_dev_eval import dev_eval  # noqa: E402
runs = sys.argv[1:]
R = [dev_eval(os.path.join(r, "orthobench_ogs.txt"))["per_refog"] for r in runs]
print("refog n " + " ".join(os.path.basename(r.rstrip("/")) for r in runs))
for i, row in enumerate(R[0]):
    cells = []
    for x in (Rk[i] for Rk in R):
        d = max(x["n"] - 1, 1)
        cells.append(f"{x['tp']/d:.0f}/{x['fp']/d:.0f}/{x['fn']/d:.0f}")
    if len(set(cells)) > 1:
        print(row["refog"], row["n"], *cells)
