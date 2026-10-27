#!/usr/bin/env python
"""Diagnostic (post hoc): precision of 'rescue' = add TreeFam-A dev reference pairs NOT predicted by the frozen
QfO accurate run whose interolog support (own) >= k; real vs shuffled networks.
Output results/l2_network/treefam_rescue_precision.json. Needs work/l2_network/treefam/treefam_pairs.parquet."""
import json, os
import pandas as pd
MAIN = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results/l2_network")
c = pd.read_parquet(f"{MAIN}/work/l2_network/treefam/treefam_pairs.parquet")
genes = set(c.gene_a) | set(c.gene_b); P = set()
for l in open(f"{MAIN}/work/qfo2022_det/cand2_accurate/pairs.tsv"):
    if l[0] == "#":
        continue
    a, b = l.split("\t", 2)[:2]
    if a in genes and b in genes:
        P.add((a, b)); P.add((b, a))
c["pred"] = [(a, b) in P for a, b in zip(c.gene_a, c.gene_b)]
out = {}
for k in (1, 2, 3, 5):
    for col in ("own_count", "shuf11_own_count", "shuf23_own_count", "shuf47_own_count"):
        m = (~c.pred) & c.own_eval & (c[col] >= k)
        out[f"{col}>={k}"] = dict(n_rescued=int(m.sum()), FN_fixed=int((m & (c.label == 1)).sum()),
                                  FP_added=int((m & (c.label == 0)).sum()),
                                  precision=round(float(c.label[m].mean()), 3) if m.any() else None)
json.dump(out, open(f"{OUT}/treefam_rescue_precision.json", "w"), indent=1)
print(json.dumps(out, indent=1))
