#!/usr/bin/env python
"""L2-genome step (a) for GSR: label every logged GSR decision with the dev RefOG pair labels.

    python scripts/L2_gsr_standalone.py ARM [ARM ...]   (reads work/l2_genome/ob/<arm>/gsr_log.json)

A decision = (still-alone unit u, its best-supported group g). Label: ortholog share among the
labelled u x g gene pairs (work/S5/refog_pair_labels.parquet, DEV RefOGs; within-RefOG and deep
paralog pairs); correct if >= 0.5. Reports, per arm: AUC of support s1 over labelled decisions
(bootstrap over RefOGs, B=1000, seed 11), and precision / count of the joins the rule executed.
Writes results/l2_genome/gsr_standalone.json.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"


def load_labels():
    lab = pd.read_parquet(os.path.join(MAIN, "work/S5/refog_pair_labels.parquet"))
    d = {}
    for a, b, y, r in zip(lab.gene_a, lab.gene_b, lab.label, lab.refog_a):
        d[(a, b)] = d[(b, a)] = (int(y), r)
    return d


def label_rows(rows, L):
    out = []
    for r in rows:
        ys, rs = [], []
        for x in r["u"]:
            for y in r["g"]:
                v = L.get((x, y))
                if v is not None:
                    ys.append(v[0])
                    rs.append(v[1])
        if not ys:
            continue
        out.append(dict(s1=r["s1"], s2=r["s2"], margin=r["s1"] - r["s2"], ratio=r["ratio"], done=r["done"],
                        ok=r["ok_c"] and r["ok_s"] and r["ok_m"], y=int(np.mean(ys) >= 0.5),
                        refog=pd.Series(rs).mode()[0], node=r["node"], n_u=len(r["u"])))
    return pd.DataFrame(out)


def boot_auc(d, col, B=1000, seed=11):
    rng = np.random.default_rng(seed)
    cl = d.refog.values
    u = np.unique(cl)
    idx = {k: np.flatnonzero(cl == k) for k in u}
    vals = []
    for _ in range(B):
        ix = np.concatenate([idx[k] for k in rng.choice(u, len(u))])
        y = d.y.values[ix]
        if 0 < y.sum() < len(y):
            vals.append(roc_auc_score(y, d[col].values[ix]))
    return [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)] if vals else None


def main():
    L = load_labels()
    res = {}
    for arm in sys.argv[1:]:
        rows = json.load(open(os.path.join(ROOT, "work/l2_genome/ob", arm, "gsr_log.json")))
        d = label_rows(rows, L)
        e = dict(n_decisions=len(rows), n_labelled=len(d), coverage=round(len(d) / max(1, len(rows)), 4))
        if len(d) and 0 < d.y.sum() < len(d):
            e["auc_s1"] = round(float(roc_auc_score(d.y, d.s1)), 4)
            e["auc_s1_ci"] = boot_auc(d, "s1")
            e["auc_margin"] = round(float(roc_auc_score(d.y, d.margin)), 4)
        ex = d[d.done] if len(d) else d
        e["executed_labelled"] = int(len(ex))
        e["executed_correct"] = int(ex.y.sum()) if len(ex) else 0
        e["executed_precision"] = round(float(ex.y.mean()), 4) if len(ex) else None
        e["executed_total"] = int(sum(r["done"] for r in rows))
        e["gate_pass_labelled"] = int(d.ok.sum()) if len(d) else 0
        e["gate_pass_precision"] = round(float(d[d.ok].y.mean()), 4) if len(d) and d.ok.any() else None
        e["all_decisions_precision"] = round(float(d.y.mean()), 4) if len(d) else None
        res[arm] = e
        print(arm, e, flush=True)
    out = os.path.join(ROOT, "results/l2_genome/gsr_standalone.json")
    old = json.load(open(out)) if os.path.exists(out) else {}
    old.update(res)
    json.dump(old, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
