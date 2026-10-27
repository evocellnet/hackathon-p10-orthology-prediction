#!/usr/bin/env python
"""L2-genome step (a) for GSR on QfO: label logged GSR decisions with TreeFam-A DEV reference pairs.

    python scripts/L2_gsr_standalone_qfo.py fast_R_u_log fast_R_u_s11_log ...

A decision (still-alone unit u, best-supported group g) is labelled by the TreeFam-A DEV reference
pairs between u's and g's genes: correct if >= 50% are speciation (S) pairs. Reports AUC of the
support s1 over labelled decisions (bootstrap over families, B=1000, seed 11), and the precision
of the joins the rule executed (+ of all gate-passing decisions). Appends to
results/l2_genome/gsr_standalone_qfo.json.
"""
import json
import os
import sys
from collections import Counter

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from qfo_treefam_dev_fast import dev_refpairs  # noqa: E402


def main():
    L = {}
    for a, b, t, f in dev_refpairs():
        L[(a, b)] = L[(b, a)] = (t == "S", f)
    res = {}
    for arm in sys.argv[1:]:
        d = json.load(open(os.path.join(ROOT, "work/l2_genome/qfo", arm, "gsr_log.json")))
        rows = []
        for r in d["rows"]:
            v = [L[(x, y)] for x in r["u"] for y in r["g"] if (x, y) in L]
            if not v:
                continue
            fam = Counter(f for _, f in v).most_common(1)[0][0]
            rows.append((r["s1"], r["s1"] - r["s2"], int(np.mean([s for s, _ in v]) >= 0.5), r["done"],
                         r["ok_c"] and r["ok_s"] and r["ok_m"], fam, len(v), r["node"]))
        e = dict(n_decisions=d["n_decisions"], n_labelled=len(rows))
        if rows:
            s1 = np.array([x[0] for x in rows]); y = np.array([x[2] for x in rows])
            done = np.array([x[3] for x in rows]); ok = np.array([x[4] for x in rows])
            fam = np.array([x[5] for x in rows]); npairs = np.array([x[6] for x in rows])
            if 0 < y.sum() < len(y):
                e["auc_s1"] = round(float(roc_auc_score(y, s1)), 4)
                rng = np.random.default_rng(11)
                u = np.unique(fam)
                idx = {k: np.flatnonzero(fam == k) for k in u}
                vals = []
                for _ in range(1000):
                    ix = np.concatenate([idx[k] for k in rng.choice(u, len(u))])
                    if 0 < y[ix].sum() < len(ix):
                        vals.append(roc_auc_score(y[ix], s1[ix]))
                e["auc_s1_ci"] = [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)]
            e.update(all_precision=round(float(y.mean()), 4),
                     executed_labelled=int(done.sum()), executed_precision=round(float(y[done].mean()), 4) if done.any() else None,
                     executed_ref_pairs=int(npairs[done].sum()),
                     executed_ref_pairs_S=int(sum(npairs[i] for i in np.flatnonzero(done & (y == 1)))),
                     gate_pass_labelled=int(ok.sum()), gate_pass_precision=round(float(y[ok].mean()), 4) if ok.any() else None,
                     executed_by_node=dict(Counter(rows[i][7] for i in np.flatnonzero(done)).most_common(10)))
        res[arm] = e
        print(arm, e, flush=True)
    out = os.path.join(ROOT, "results/l2_genome/gsr_standalone_qfo.json")
    old = json.load(open(out)) if os.path.exists(out) else {}
    old.update(res)
    json.dump(old, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
