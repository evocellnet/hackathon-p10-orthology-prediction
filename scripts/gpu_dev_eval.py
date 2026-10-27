#!/usr/bin/env python
"""v2-gpu: OrthoBench DEV comparison of runs against a reference run (docs/GPU_v2.md).

    python scripts/gpu_dev_eval.py REF_RUN RUN [RUN ...] [--json out.json]

Per run: OG-level F/P/R (official benchmark.py on the 34 DEV RefOGs, scripts/orthobench_dev_eval.py),
ortholog-pair P/R/F1 micro and macro (mean over DEV RefOGs; scripts/orthobench_dev_pairs_eval.py on
hogs.orthoxml). Against REF_RUN: paired bootstrap over DEV RefOGs (B = 2000, seed 11) of
delta OG-F (per-RefOG TP/FP/FN pooled) and delta macro pair F1. Test RefOGs are never read.
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from orthobench_dev_eval import dev_eval  # noqa: E402
from orthobench_dev_pairs_eval import evaluate as pair_eval, load_truth  # noqa: E402


def one(run, truth):
    ev = dev_eval(os.path.join(run, "orthobench_ogs.txt"))
    pm, parr, prefogs = pair_eval(os.path.join(run, "hogs.orthoxml"), truth)
    per = parr.sum(axis=1)
    tp, fn, fp = (per[:, k].astype(float) for k in (0, 1, 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        f1_r = 2 * tp / (2 * tp + fp + fn)
    og = np.array([[r["tp"], r["fp"], r["fn"]] for r in ev["per_refog"]], float)
    return dict(run=run, F=round(ev["F"], 2), P=round(ev["precision"], 2), R=round(ev["recall"], 2),
                exact=ev["exact"], pairF1_micro=pm["all"]["F1"], pairF1_macro=round(float(np.nanmean(f1_r)), 4),
                _og=og, _f1r=f1_r, _refogs=[r["refog"] for r in ev["per_refog"]], _prefogs=list(prefogs))


def boot(a, b, B=2000, seed=11):
    rng = np.random.default_rng(seed)
    n = len(a["_og"])
    assert a["_refogs"] == b["_refogs"] and a["_prefogs"] == b["_prefogs"]
    W = rng.multinomial(n, np.full(n, 1.0 / n), size=B).astype(float)

    def ogf(x, w):
        s = w @ x
        return 100 * 2 * s[:, 0] / (2 * s[:, 0] + s[:, 1] + s[:, 2])
    dF = ogf(b["_og"], W) - ogf(a["_og"], W)
    pa, pb = a["_f1r"], b["_f1r"]
    ok = np.isfinite(pa) & np.isfinite(pb)
    m = int(ok.sum())
    W2 = rng.multinomial(m, np.full(m, 1.0 / m), size=B).astype(float) / m
    dM = W2 @ (pb[ok] - pa[ok])
    pt = lambda x: 100 * 2 * x[:, 0].sum() / (2 * x[:, 0].sum() + x[:, 1].sum() + x[:, 2].sum())  # noqa: E731
    return dict(dOG_F=round(pt(b["_og"]) - pt(a["_og"]), 2),
                dOG_F_ci=[round(float(np.percentile(dF, 2.5)), 2), round(float(np.percentile(dF, 97.5)), 2)],
                dmacro=round(float((pb[ok] - pa[ok]).mean()), 4),
                dmacro_ci=[round(float(np.percentile(dM, 2.5)), 4), round(float(np.percentile(dM, 97.5)), 4)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--json")
    a = ap.parse_args()
    truth = load_truth()
    ref = one(a.ref, truth)
    out = []
    for r in [a.ref] + a.runs:
        x = ref if r == a.ref else one(r, truth)
        row = {k: v for k, v in x.items() if not k.startswith("_")}
        if r != a.ref:
            row.update(boot(ref, x))
        out.append(row)
        print(json.dumps(row))
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
