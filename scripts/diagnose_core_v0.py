"""Classify core-v0 errors on a toy simulation (guide tree = true species tree).

    python scripts/diagnose_core_v0.py SIMDIR [--arm max+repair] [--score norm_max] [--backend mmseqs]

For each missed true ortholog pair (a, b) with true level L (speciation node),
look at the merge decision at L: the ancestral genes ua, ub of a and b at L's
children ("units") were not merged. Unit pairs (L, ua, ub) are counted once,
classified by what happened to each unit at L:

* ``both_unmatched``   both units passed up unchanged
* ``one_unmatched``    one unit joined a group, the other passed up
* ``split_groups``     both joined different groups at L
* ``no_edge``          (sub-flag) no scored leaf pair between ua and ub

For false ortholog pairs that are true paralogs, the predicted level is reported.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ancestra.io import load_families, load_guide_tree, load_proteomes  # noqa: E402
from ancestra.merge import ARMS, leaf_scores, run_merge, to_orthoxml  # noqa: E402
from ancestra.orthoxml import hogs_to_pairs, read_orthoxml  # noqa: E402
from ancestra.scoring import ScoringConfig, score_families  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("sim")
    ap.add_argument("--arm", default="max+repair")
    ap.add_argument("--score", default="norm_max")
    ap.add_argument("--backend", default="mmseqs")
    ap.add_argument("--evalue", type=float, default=10.0)
    ap.add_argument("--rule", default="ortholog")
    ap.add_argument("--cov", type=float, default=0.5)
    ap.add_argument("--work", default="/tmp/ancestra_diag")
    ap.add_argument("--examples", type=int, default=6)
    a = ap.parse_args(argv)
    gt = load_proteomes(os.path.join(a.sim, "proteomes"))
    fam = load_families(os.path.join(a.sim, "truth/gene_families.tsv"), gt, 0, 2)
    tree = load_guide_tree(os.path.join(a.sim, "species_tree.nwk"), gt.species)
    st, _ = score_families(gt, fam, ScoringConfig(backend=a.backend, evalue=a.evalue), a.work)
    cfg = replace(ARMS[a.arm], score=a.score, inpar_rule=a.rule, min_cov=a.cov)
    trace: dict = {}
    res = run_merge(gt, tree, st, cfg, trace=trace)
    ox = to_orthoxml(res, gt, tree)
    truth = read_orthoxml(os.path.join(a.sim, "truth/hogs.orthoxml"))
    tp_lvl = hogs_to_pairs(truth, with_level=True)
    _, tpar = hogs_to_pairs(truth, with_paralogs=True)
    pred = hogs_to_pairs(ox, with_level=True)
    idx = gt.index()
    s = leaf_scores(gt, st, cfg.score)
    keep = np.minimum(st.qcov, st.tcov) >= cfg.min_cov
    sc = {(int(x), int(y)): float(v) for x, y, v, k in zip(st.q, st.t, s, keep) if k}

    def S(x, y):
        return sc.get((x, y) if x < y else (y, x))
    node_maps = {}
    for nd, (genes, before, after) in trace.items():
        node_maps[nd] = (dict(zip(genes.tolist(), before.tolist())),
                         dict(zip(genes.tolist(), after.tolist())))
    fn = [(p, l) for p, l in tp_lvl.items() if p not in pred]
    cls = Counter()
    cls_pairs = Counter()
    seen = {}
    by_level = Counter()
    for (x, y), L in fn:
        gx, gy = idx[x], idx[y]
        bef, aft = node_maps[L]
        ux, uy = bef[gx], bef[gy]
        key = (L, min(ux, uy), max(ux, uy))
        mx, my = aft[gx] != ux, aft[gy] != uy
        c = "split_groups" if (mx and my) else ("one_unmatched" if (mx or my) else "both_unmatched")
        if key not in seen:
            membx = [g for g, u in bef.items() if u == ux]
            memby = [g for g, u in bef.items() if u == uy]
            vals = [S(p, q) for p in membx for q in memby]
            vals = [v for v in vals if v is not None]
            ne = not vals
            seen[key] = (c, ne, len(membx), len(memby), max(vals) if vals else None)
            cls[c + ("+no_edge" if ne else "")] += 1
            by_level[L] += 1
        cls_pairs[c + ("+no_edge" if seen[key][1] else "")] += 1
    print(f"{a.sim} arm={a.arm} score={a.score} backend={a.backend}")
    print(f"FN pairs={len(fn)} of {len(tp_lvl)} true; distinct unit pairs={len(seen)}")
    for k, v in cls.most_common():
        print(f"  {k:28s} unit-pairs={v:5d}  pairs={cls_pairs[k]}")
    depth = {n.name: 0 for n in tree.preorder()}
    for n in tree.preorder():
        for ch in n.children:
            depth[ch.name] = depth[n.name] + 1
    print("  FN unit-pairs by level (depth):",
          ", ".join(f"{L}(d{depth[L]})={v}" for L, v in sorted(by_level.items(), key=lambda kv: depth[kv[0]])))
    fp = [p for p in pred if p not in tp_lvl]
    fpl = Counter(pred[p] for p in fp if p in tpar)
    print(f"FP pairs={len(fp)} (true paralogs {sum(fpl.values())}); by predicted level: {dict(fpl)}")
    # FP examples: predicted ortholog pairs that are true paralogs -> unit-level view
    fp_units = Counter()
    gene_name = gt.names
    for (x, y) in fp:
        if (x, y) not in tpar:
            continue
        L = pred[(x, y)]
        bef, aft = node_maps[L]
        fp_units[(L, bef[idx[x]], bef[idx[y]])] += 1
    print(f"  FP distinct unit pairs={len(fp_units)}; top: {fp_units.most_common(5)}")
    # examples
    ex = [(k, v) for k, v in seen.items()][: a.examples]
    for (L, ux, uy), (c, ne, nx, ny, smax) in ex:
        print(f"  ex L={L} units {ux}({nx} genes) {uy}({ny} genes) class={c} maxS={smax}")


if __name__ == "__main__":
    main()
