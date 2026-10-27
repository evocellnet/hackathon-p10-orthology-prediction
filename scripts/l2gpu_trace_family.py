#!/usr/bin/env python
"""L2-gpu task 1 mechanism trace: rerun the frozen merge (accurate preset, synteny two-pass) on ONE S0
family with two score tables and report, per species-tree node, how the genes of one dev RefOG are
partitioned into ancestral genes, plus the avg-link scores between the RefOG's units at the first node
where the tables disagree.

    python scripts/l2gpu_trace_family.py --family 77 --refog RefOG021 --tables A.npz B.npz --out out.json
"""
import argparse
import json
import os
import sys
from dataclasses import replace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import GeneTable, load_proteomes  # noqa: E402
from ancestra.merge import MergeConfig, run_merge_constrained, leaf_scores, root_labels  # noqa: E402
from ancestra.pipeline import SYNTENY_DEFAULT, apply_preset, PRESETS, load_guide_tree  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from ancestra.synteny import SyntenyContext, gene_positions  # noqa: E402

OB = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS")


def sub_table(st, keep_idx, remap):
    m = np.isin(st.q, keep_idx) & np.isin(st.t, keep_idx)
    return ScoreTable(remap[st.q[m]].astype(np.int32), remap[st.t[m]].astype(np.int32), st.bits[m],
                      st.qcov[m], st.tcov[m], st.pident[m], st.evalue[m], st.self_bits[keep_idx])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", type=int, required=True)
    ap.add_argument("--refog", required=True)
    ap.add_argument("--tables", nargs="+", required=True)
    ap.add_argument("--families", default=os.path.join(ROOT, "..", "gpuv2/ob_gpusw_acc/stages/families/families.npy"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--link", default="avg")
    ap.add_argument("--merge-json", default="{}", help="diagnostic MergeConfig overrides")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    gt = load_proteomes(os.path.join(OB, "Input"))
    tree = load_guide_tree(os.path.join(ROOT, "results/S1_speciestree/ob/ncbi.nwk"), gt.species)
    fam = np.load(a.families)
    keep = np.flatnonzero(fam == a.family)
    remap = np.full(gt.n, -1, np.int64)
    remap[keep] = np.arange(len(keep))
    sgt = GeneTable([gt.names[i] for i in keep], gt.species, gt.species_of[keep], gt.lengths[keep],
                    [gt.seqs[i] for i in keep])
    base = {n.split(".")[0]: i for i, n in enumerate(sgt.names)}
    ref = [base.get(l.strip().split(".")[0]) for l in open(os.path.join(OB, "RefOGs", a.refog + ".txt")) if l.strip()]
    ref = np.array([r for r in ref if r is not None])
    names = [sp if os.path.exists(os.path.join(ROOT, "data/genomes/orthobench", sp + ".parquet")) else sp.split(".")[0]
             for sp in sgt.species]
    cc, oo, _ = gene_positions(sgt.names, sgt.species_of, names, os.path.join(ROOT, "data/genomes/orthobench"))
    merge = apply_preset(MergeConfig(), "accurate", a.threads)
    merge = replace(merge, ctx_veto=True, ctx_r3=True, link=a.link, **json.loads(a.merge_json))
    out = {"merge_overrides": json.loads(a.merge_json), "link": a.link, "family": a.family, "family_size": int(len(keep)), "refog": a.refog, "n_ref_in_family": int(len(ref)),
           "runs": {}}
    traces = {}
    for tp in a.tables:
        st = sub_table(ScoreTable.load(tp), keep, remap)
        ctx = SyntenyContext(cc, oo, family=np.zeros(sgt.n, np.int64), **SYNTENY_DEFAULT)
        tr: dict = {}
        res = run_merge_constrained(sgt, tree, st, merge, ctx,
                                    off_cfg=replace(merge, tree_method=PRESETS["accurate"]["off_tree_method"]),
                                    trace=tr)
        lab = root_labels(res, sgt.n)
        per_node = {}
        for node, val in tr.items():
            genes, at_children, at_node = val
            pos = {int(g): i for i, g in enumerate(genes)}
            idx = [pos[g] for g in ref.tolist() if g in pos]
            if not idx:
                continue
            u_child = np.asarray(at_children)[idx]
            u_node = np.asarray(at_node)[idx]
            per_node[node] = dict(n_ref=len(idx), units_in=int(len(np.unique(u_child))),
                                  units_out=int(len(np.unique(u_node))))
        root_parts = np.unique(lab[ref], return_counts=True)[1]
        out["runs"][tp] = dict(n_pairs=int(len(st)), per_node=per_node,
                               root_parts=sorted(root_parts.tolist(), reverse=True),
                               stats={s["node"]: {k: s.get(k) for k in ("n_units", "n_mb", "r1_unions", "r1_blocked",
                                                                         "r2_absorbed", "r3_repaired", "tree_clusters",
                                                                         "tree_changed", "tree_too_big")}
                                      for s in res.node_stats})
        traces[tp] = tr
        print(tp, out["runs"][tp]["root_parts"], json.dumps(per_node))
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
