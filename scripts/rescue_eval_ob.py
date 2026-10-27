"""OrthoBench DEV eval of one rescue arm run dir: OG-level (official benchmark.py, dev RefOGs)
+ pair P/R/F1 micro and macro (orthobench_dev_pairs_eval). Writes RUN/dev_eval_rescue.json and
appends to results/rescue/ob_runs.tsv. Optional --base RUN0: paired bootstrap of pair F1.

    python scripts/rescue_eval_ob.py work/rescue/ob/R1 --label R1 [--base work/rescue/ob/R0]
"""
import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from orthobench_dev_eval import dev_eval  # noqa: E402
from orthobench_dev_pairs_eval import evaluate as pair_eval, load_truth  # noqa: E402

LOG = os.path.join(ROOT, "results/rescue/ob_runs.tsv")


def run_eval(run):
    ev = dev_eval(os.path.join(run, "orthobench_ogs.txt"))
    truth = load_truth()
    pm, arr, refogs = pair_eval(os.path.join(run, "hogs.orthoxml"), truth)
    per = arr.sum(axis=1)
    tp, fn, fp = per[:, 0].astype(float), per[:, 1].astype(float), per[:, 2].astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        f1 = 2 * tp / (2 * tp + fp + fn)
    return ev, pm, f1, refogs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    ev, pm, f1, refogs = run_eval(a.run)
    man = json.load(open(os.path.join(a.run, "manifest.json")))
    st = man["stages"]
    res = dict(label=a.label, F=round(ev["F"], 2), P=round(ev["precision"], 2), R=round(ev["recall"], 2),
               exact=ev["exact"], n_split=ev["n_refogs_split"], n_with_fp=ev["n_refogs_with_fp"],
               pairP=pm["all"]["P"], pairR=pm["all"]["R"], pairF1=pm["all"]["F1"],
               pairF1_macro=round(float(np.nanmean(f1)), 4),
               rescue_s=st.get("rescue", {}).get("wall_s"), merge_s=st["merge"]["wall_s"],
               rescue_stats=st.get("rescue", {}).get("stats", {}).get("attach"),
               guard=st["merge"].get("stats", {}).get("rescue_guard"))
    out = dict(res, per_refog_pairF1={r: (None if not np.isfinite(x) else round(float(x), 4))
                                      for r, x in zip(refogs, f1)},
               per_refog_og={r["refog"] if "refog" in r else i: r for i, r in enumerate(ev.get("per_refog", []))})
    json.dump(out, open(os.path.join(a.run, "dev_eval_rescue.json"), "w"), indent=1, default=str)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    cols = ["time", "label", "F", "P", "R", "exact", "n_split", "n_with_fp", "pairP", "pairR", "pairF1",
            "pairF1_macro", "rescue_s", "merge_s", "rescue_stats", "guard"]
    row = dict(time=datetime.now().isoformat(timespec="seconds"), **res)
    new = not os.path.exists(LOG)
    with open(LOG, "a") as fh:
        if new:
            fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(json.dumps(row[c]) if isinstance(row[c], dict) else str(row[c]) for c in cols) + "\n")
    print(json.dumps(res, default=str))


if __name__ == "__main__":
    main()
