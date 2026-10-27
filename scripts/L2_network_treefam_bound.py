#!/usr/bin/env python
"""Oracle upper bound on TreeFam-A dev (frozen QfO accurate run): if every network-evaluable error with the
'right' interolog sign were fixed (FN with count>0 -> TP; FP evaluable with count==0 -> TN), how far could
TPR/PPV move? Also per-species-pair AUC of own_count. Output results/l2_network/treefam_bound.json."""
import json, os
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
MAIN = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results/l2_network")
c = pd.read_parquet(f"{MAIN}/work/l2_network/treefam/treefam_pairs.parquet")
genes = set(c.gene_a) | set(c.gene_b)
P = set()
for line in open(f"{MAIN}/work/qfo2022_det/cand2_accurate/pairs.tsv"):
    if line[0] == "#":
        continue
    a, b = line.split("\t", 2)[:2]
    if a in genes and b in genes:
        P.add((a, b)); P.add((b, a))
c["pred"] = [(a, b) in P for a, b in zip(c.gene_a, c.gene_b)]
def rates(tp, fn, fp):   # official formula (uniform prior +1, half counts)
    tp_, fn_, fp_ = 1 + tp / 2, 1 + fn / 2, 1 + fp / 2
    return round(tp_ / (tp_ + fn_), 4), round(tp_ / (tp_ + fp_), 4)
out = {}
for ch in ("own", "exp150"):
    ev, cnt = c[f"{ch}_eval"].values, c[f"{ch}_count"].values
    y, p = c.label.values == 1, c.pred.values
    tp, fn, fp = int((y & p).sum()), int((y & ~p).sum()), int((~y & p).sum())
    fixfn = int((y & ~p & ev & (cnt > 0)).sum()); fixfp = int((~y & p & ev & (cnt == 0)).sum())
    out[ch] = dict(TP=tp, FN=fn, FP=fp, FN_evaluable=int((y & ~p & ev).sum()), FP_evaluable=int((~y & p & ev).sum()),
                   fixable_FN=fixfn, fixable_FP=fixfp, TPR_PPV_off=rates(tp, fn, fp),
                   TPR_PPV_oracle=rates(tp + fixfn, fn - fixfn, fp - fixfp))
    pp = {}
    e = c[ev]
    for pr, d in e.groupby("pair"):
        if len(d) >= 50 and d.label.nunique() == 2:
            pp[pr] = dict(n=int(len(d)), orth=int(d.label.sum()), auc=round(roc_auc_score(d.label, d[f"{ch}_count"]), 3),
                          auc_shuf11=round(roc_auc_score(d.label, d[f"shuf11_{ch}_count"]), 3),
                          auc_seq=round(roc_auc_score(d.label, d.seq), 3))
    out[ch]["per_pair_auc"] = pp
json.dump(out, open(f"{OUT}/treefam_bound.json", "w"), indent=1)
print(json.dumps(out, indent=1))
