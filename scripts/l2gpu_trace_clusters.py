#!/usr/bin/env python
"""L2-gpu mechanism trace (diagnostic): wraps localtree.resolve_ambiguous during the one-family merge of
l2gpu_trace_family.py and records, at one node, the ambiguous clusters that contain genes of the RefOG
(units, RefOG genes per unit, labels before/after the local tree).
    python scripts/l2gpu_trace_clusters.py NODE --family 77 --refog RefOG021 --tables A B --out X.json"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ancestra.localtree as LT  # noqa: E402
import l2gpu_trace_family as TF  # noqa: E402

NODE = sys.argv.pop(1)
REC = []
_orig = LT.resolve_ambiguous


def wrapped(n, child, lu, lv, score, labels, out_edges=None, tau=0.9, *a, **k):
    new = _orig(n, child, lu, lv, score, labels, out_edges, tau, *a, **k)
    members = a[5] if len(a) > 5 else k.get("members")
    if members is None or WANT["ref"] is None:
        return new
    ref = WANT["ref"]
    s = np.asarray(score, float)
    best = np.full(n, -np.inf)
    np.maximum.at(best, lu, s)
    np.maximum.at(best, lv, s)
    amb = (labels[lu] != labels[lv]) & (s >= tau * np.minimum(best[lu], best[lv]))
    par = np.arange(n)

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for x, y in zip(lu[amb].tolist(), lv[amb].tolist()):
        rx, ry = f(x), f(y)
        if rx != ry:
            par[rx] = ry
    nref = np.array([len(ref.intersection(members(i).tolist())) for i in range(n)])
    roots = {}
    for i in np.flatnonzero(nref).tolist():
        roots.setdefault(f(i), []).append(i)
    isref = nref > 0
    bestref = np.full(n, -np.inf)
    m = isref[lu] & isref[lv]
    np.maximum.at(bestref, lu[m], s[m])
    np.maximum.at(bestref, lv[m], s[m])
    bestpart = np.full(n, -1)
    for e in np.argsort(s, kind="stable").tolist():   # ascending: last write = best partner
        bestpart[lu[e]] = lv[e]
        bestpart[lv[e]] = lu[e]
    deg = np.bincount(np.concatenate([lu, lv]), minlength=n)
    unit_info = {int(i): dict(nref=int(nref[i]), size=int(len(members(i))), child=int(child[i]), best=round(float(best[i]), 4),
                              best_to_ref=round(float(bestref[i]), 4), best_partner_is_ref=bool(isref[bestpart[i]]) if bestpart[i] >= 0 else None,
                              best_partner_size=int(len(members(int(bestpart[i])))) if bestpart[i] >= 0 else None,
                              degree=int(deg[i]))
                 for i in np.flatnonzero(isref).tolist()}
    rec = [{"units": unit_info}]
    for r, us in roots.items():
        comp = [i for i in range(n) if f(i) == r]
        rec.append(dict(comp_units=len(comp), ref_units=[(int(nref[i]), int(labels[i]), int(new[i]), int(child[i])) for i in us],
                        n_unit_labels_before=len({int(labels[i]) for i in comp}),
                        n_unit_labels_after=len({int(new[i]) for i in comp})))
    REC.append(rec)
    return new


LT.resolve_ambiguous = wrapped
WANT = {"ref": None}
_orig_run = TF.run_merge_constrained


def run_mc(sgt, tree, st, merge, ctx, off_cfg=None, trace=None):
    return _orig_run(sgt, tree, st, merge, ctx, off_cfg=off_cfg, trace=trace)


# hook: set the RefOG gene set once the family is loaded (TF.main computes it locally) -> patch np.array use
_orig_main = TF.main


def main():
    import ancestra.merge as M
    orig_rm = M.run_merge
    out_all = {}

    def rm(gt, tree, st, cfg, context=None, **kw):
        names = gt.names
        base = {nm.split(".")[0]: i for i, nm in enumerate(names)}
        ref = [base.get(l.strip().split(".")[0]) for l in open(os.path.join(TF.OB, "RefOGs", ARGS_REFOG + ".txt")) if l.strip()]
        WANT["ref"] = {r for r in ref if r is not None}
        REC.clear()
        node_names = [nd.name for nd in tree.postorder() if nd.children]
        res = orig_rm(gt, tree, st, cfg, context=context, **kw)
        # REC has one entry per tree-mode node, in postorder of nodes where tree mode ran
        out_all.setdefault("passes", []).append([{"rec": r} for r in REC])
        out_all.setdefault("tree_nodes", []).append([s["node"] for s in res.node_stats if "tree_s" in s])
        return res
    M.run_merge = rm
    TF.run_merge_constrained.__globals__["run_merge"] = rm
    _orig_main()
    json.dump(out_all, open(OUT2, "w"), indent=0)


ARGS_REFOG = sys.argv[sys.argv.index("--refog") + 1]
OUT2 = sys.argv[sys.argv.index("--out") + 1].replace(".json", "_clusters.json")
main()
