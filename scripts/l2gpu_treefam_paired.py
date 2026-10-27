#!/usr/bin/env python
"""TreeFam-A DEV: TPR/PPV of runs vs a reference run, with a paired bootstrap over dev families
(B = 2000, seed 11) for delta TPR / delta PPV. Official pair formula of qfo_treefam_dev_fast.py
(counts / 2 + 1 pseudocount). Test families are never read (dev_refpairs filters at parse time).

    python scripts/l2gpu_treefam_paired.py REF_RUN RUN [RUN ...] --json out.json
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from qfo_treefam_dev_fast import dev_refpairs, predicted_orthologs  # noqa: E402


def per_family(run, refs, fams, genes):
    src = os.path.join(run, "pairs.tsv")
    if not os.path.exists(src):
        src = os.path.join(run, "hogs.orthoxml")
    pred = predicted_orthologs(src, genes)
    fi = {f: i for i, f in enumerate(fams)}
    C = np.zeros((len(fams), 4))       # TP FP FN TN
    for a, b, t, f in refs:
        p = ((a, b) if a < b else (b, a)) in pred
        C[fi[f], (0 if p else 2) if t == "S" else (1 if p else 3)] += 1
    return C


def rates(C):
    tp, fp, fn = 1 + C[..., 0] / 2, 1 + C[..., 1] / 2, 1 + C[..., 2] / 2
    return tp / (tp + fn), tp / (tp + fp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--json")
    a = ap.parse_args()
    refs = dev_refpairs()
    fams = sorted({r[3] for r in refs})
    genes = {x for r in refs for x in r[:2]}
    R = per_family(a.ref, refs, fams, genes)
    rng = np.random.default_rng(11)
    W = rng.multinomial(len(fams), np.full(len(fams), 1 / len(fams)), size=2000).astype(float)
    tpr0, ppv0 = rates(R.sum(0))
    out = [dict(run=a.ref, TPR=round(tpr0, 4), PPV=round(ppv0, 4))]
    for run in a.runs:
        X = per_family(run, refs, fams, genes)
        tpr, ppv = rates(X.sum(0))
        bt, bp = rates(W @ X)
        rt, rp = rates(W @ R)
        dt, dp = bt - rt, bp - rp
        out.append(dict(run=run, TPR=round(tpr, 4), PPV=round(ppv, 4), dTPR=round(tpr - tpr0, 4),
                        dTPR_ci=[round(float(np.percentile(dt, 2.5)), 4), round(float(np.percentile(dt, 97.5)), 4)],
                        dPPV=round(ppv - ppv0, 4),
                        dPPV_ci=[round(float(np.percentile(dp, 2.5)), 4), round(float(np.percentile(dp, 97.5)), 4)],
                        n_families=len(fams)))
        print(json.dumps(out[-1]))
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
