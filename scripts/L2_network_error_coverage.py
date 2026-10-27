#!/usr/bin/env python
"""Upper bound for step (b): which of the frozen accurate run's OrthoBench-dev pair errors (FN, FP)
are network-evaluable at all, and does the interolog count point the right way on them?
(Diagnostic only; no merge is run. Output results/l2_network/error_coverage.json.)"""
import json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import L2_network_standalone as S

t0 = time.time()
sp_of, fam = S.load_species()
anc = S.anchors(sp_of)
nets = S.networks()
lab = pd.read_parquet(f"{S.MAIN}/work/S5/refog_pair_labels.parquet")
lab = lab[lab.kind == "within"].reset_index(drop=True)          # the 19,146 dev truth pairs
pred = pd.read_csv(f"{S.RUN}/pairs.tsv", sep="\t", comment="#", header=None, usecols=[0, 1], names=["a", "b"])
P = set(zip(pred.a, pred.b)) | set(zip(pred.b, pred.a))
lab["pred"] = [(a, b) in P for a, b in zip(lab.gene_a, lab.gene_b)]
lab["err"] = np.where((lab.label == 1) & ~lab.pred, "FN", np.where((lab.label == 0) & lab.pred, "FP",
                      np.where(lab.label == 1, "TP", "TN")))
out = dict(n_truth=int(len(lab)), counts=lab.err.value_counts().to_dict())
for ch in ("own", "db400_CIRCULAR"):
    adj = {sp: S.adj_from(a, b) for (sp, c), (a, b, nodes) in nets.items() if c == ch and len(nodes) >= S.MIN_NODES}
    ev = np.zeros(len(lab), bool); cnt = np.zeros(len(lab))
    for (A, B), d in lab.groupby(["sp_a", "sp_b"]):
        if A not in adj or B not in adj:
            continue
        k = tuple(sorted((A, B)))
        dd = d if (A, B) == k else d.rename(columns={"gene_a": "gene_b", "gene_b": "gene_a"})
        M = anc.get(k, {})
        c, _, na_, nb_ = S.interolog(dd, adj[k[0]], adj[k[1]], M, fam)
        ev[d.index] = (na_ > 0) & (nb_ > 0); cnt[d.index] = c
    r = {}
    for e in ("TP", "FN", "FP", "TN"):
        m = (lab.err == e).values
        r[e] = dict(n=int(m.sum()), evaluable=int(ev[m].sum()), count_gt0=int((cnt[m] > 0).sum()))
    # best case: every evaluable FN with count>0 becomes TP, every evaluable FP with count==0 becomes TN
    tp, fn, fp = r["TP"]["n"], r["FN"]["n"], r["FP"]["n"]
    fix_fn = int(((lab.err == "FN").values & ev & (cnt > 0)).sum())
    fix_fp = int(((lab.err == "FP").values & ev & (cnt == 0)).sum())
    f0 = 2 * tp / (2 * tp + fn + fp); f1 = 2 * (tp + fix_fn) / (2 * (tp + fix_fn) + (fn - fix_fn) + (fp - fix_fp))
    r["oracle_upper_bound"] = dict(fixable_FN=fix_fn, fixable_FP=fix_fp, micro_F1_off=round(f0, 4),
                                   micro_F1_if_all_fixed=round(f1, 4), delta=round(f1 - f0, 4))
    by = lab.assign(ev=ev).groupby("err").ev.mean().round(4).to_dict()
    r["frac_evaluable_by_class"] = by
    out[ch] = r
out["wall_s"] = round(time.time() - t0, 1)
json.dump(out, open(f"{S.OUT}/error_coverage.json", "w"), indent=1)
print(json.dumps(out, indent=1))
