"""Merge-only rerun on cached QfO stages + fast TreeFam-A DEV scoring (no Darwin, no test).

    python scripts/qfo_merge_iter.py [--run work/qfo_recall/base] --preset fast --set k=v ... --label X

Replicates the pipeline's merge with gene order (S5 constrained two-pass + veto + R3; the
off pass uses the preset's off_tree_method), then scores TreeFam-A DEV families with
scripts/qfo_treefam_dev_fast.py and appends to results/core_v0_qfo/runs.tsv.
"""
import argparse
import json
import os
import sys
import time
from dataclasses import asdict, replace
from datetime import datetime

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import load_guide_tree, read_gene_table  # noqa: E402
from ancestra.merge import MergeConfig, run_merge_constrained, to_orthoxml  # noqa: E402
from ancestra.pipeline import PRESETS, SYNTENY_DEFAULT, apply_preset  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from qfo_treefam_dev_fast import score as tf_score  # noqa: E402

LOG = os.path.join(ROOT, "results/core_v0_qfo/runs.tsv")
GO = os.path.join(ROOT, "data/inputs/qfo2022/gene_order")
PROT = os.path.join(ROOT, "data/inputs/qfo2022/proteomes")


def fn_breakdown(ox, gt, fam, tree):
    """Missed reference speciation pairs (TreeFam-A DEV) by category x LCA level."""
    from collections import Counter
    from ancestra.newick import TreeIndex
    from ancestra.orthoxml import Group, hogs_to_pairs
    from qfo_treefam_dev_fast import dev_refpairs
    refs = dev_refpairs()
    genes = {x for r in refs for x in r[:2]}
    top = {}
    keep = []
    for k, g in enumerate(ox.groups):
        if isinstance(g, Group):
            gg = g.genes()
            if genes & set(gg):
                keep.append(g)
                for x in gg:
                    top[x] = k
    from ancestra.orthoxml import OrthoXML
    pred = hogs_to_pairs(OrthoXML(genes=ox.genes, groups=keep, species_tree=ox.species_tree))
    idx = gt.index()
    tix = TreeIndex(tree)
    out = {}
    for a, b, t, _ in refs:
        if t != "S" or ((a, b) if a < b else (b, a)) in pred:
            continue
        if a not in idx or b not in idx:
            cat = "missing_gene"
        elif fam[idx[a]] != fam[idx[b]]:
            cat = "i_diff_family"
        elif top.get(a, -1) != top.get(b, -2) or a not in top:
            cat = "ii_diff_root_hog"
        else:
            cat = "iii_paralog_call"
        lvl = tix.mrca([ox.genes[a].species, ox.genes[b].species]) if a in ox.genes and b in ox.genes else "?"
        out.setdefault(cat, Counter())[lvl] += 1
    return {k: dict(v.most_common()) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=os.path.join(ROOT, "work/qfo_recall/base"))
    ap.add_argument("--preset", default="fast", choices=["fast", "accurate"])
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--label", required=True)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--save", help="write hogs.orthoxml here")
    ap.add_argument("--outdir", help="write a run-shaped dir (hogs.orthoxml, pairs.tsv, orthobench_ogs.txt, "
                    "stages/input hard-linked) for scripts/verify_orthoxml.py")
    a = ap.parse_args()
    cfg = apply_preset(MergeConfig(), a.preset, a.threads)
    cfg = replace(cfg, **{k: json.loads(v) for k, v in (kv.split("=", 1) for kv in a.set)})
    mcfg = replace(cfg, ctx_veto=True, ctx_r3=True)
    off_cfg = replace(mcfg, tree_method=PRESETS[a.preset]["off_tree_method"])
    gt = read_gene_table(os.path.join(a.run, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(a.run, "stages/input/guide_tree.nwk"), gt.species)
    st = ScoreTable.load(os.path.join(a.run, "stages/scores/scores.npz"))
    fam = np.load(os.path.join(a.run, "stages/scores/families_final.npy"))
    if mcfg.tree_method == "fasttree":
        from ancestra.io import load_proteomes
        full = load_proteomes(PROT)
        assert full.names == gt.names
        gt.seqs = full.seqs
    from ancestra.synteny import SyntenyContext, gene_positions
    names = [sp if os.path.exists(os.path.join(GO, sp + ".parquet")) else sp.split(".")[0] for sp in gt.species]
    cc, oo, _ = gene_positions(gt.names, gt.species_of, names, GO)
    ctx = SyntenyContext(cc, oo, family=fam, **SYNTENY_DEFAULT)
    t0 = time.time()
    c0 = os.times()
    res = run_merge_constrained(gt, tree, st, mcfg, ctx, off_cfg=off_cfg)
    c1 = os.times()
    ms = time.time() - t0
    cpu = (c1.user - c0.user) + (c1.system - c0.system) + (c1.children_user - c0.children_user) + \
        (c1.children_system - c0.children_system)
    ox = to_orthoxml(res, gt, tree)
    if a.save:
        from ancestra.orthoxml import write_orthoxml
        write_orthoxml(ox, a.save)
    if a.outdir:
        import shutil
        from ancestra.orthoxml import to_orthobench_ogs, write_orthoxml
        from ancestra.pipeline import write_pairs_streaming
        os.makedirs(os.path.join(a.outdir, "stages"), exist_ok=True)
        dst = os.path.join(a.outdir, "stages", "input")
        if not os.path.exists(dst):
            shutil.copytree(os.path.join(a.run, "stages/input"), dst, copy_function=os.link)
        write_orthoxml(ox, os.path.join(a.outdir, "hogs.orthoxml"))
        write_pairs_streaming(ox, os.path.join(a.outdir, "pairs.tsv"))
        to_orthobench_ogs(ox, None, os.path.join(a.outdir, "orthobench_ogs.txt"))
    tf = tf_score(ox)
    fnb = fn_breakdown(ox, gt, fam, tree)
    json.dump(fnb, open(os.path.join(ROOT, "results/core_v0_qfo", f"fn_breakdown_{a.label.split()[0]}.json"), "w"),
              indent=1)
    stats = {k: int(sum(s.get(k, 0) for s in res.node_stats)) for k in
             ("so_absorbed", "so_struct_merged", "r3_repaired", "tree_clusters", "tree_fasttree_jobs")}
    row = dict(time=datetime.now().isoformat(timespec="seconds"), label=a.label, preset=a.preset,
               TF_dev_TPR=round(tf["TPR"], 4), TF_dev_PPV=round(tf["PPV"], 4), **{f"TF_{k}": tf[k] for k in
               ("n_TP", "n_FP", "n_FN", "n_TN")}, merge_wall_s=round(ms, 1), merge_cpu_s=round(cpu, 1),
               stats=json.dumps(stats), config=json.dumps(asdict(cfg), sort_keys=True))
    new = not os.path.exists(LOG)
    with open(LOG, "a") as fh:
        if new:
            fh.write("\t".join(row) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")
    print(f"{a.label}: TreeFam-A DEV TPR={tf['TPR']:.4f} PPV={tf['PPV']:.4f} (TP {tf['n_TP']} FN {tf['n_FN']} "
          f"FP {tf['n_FP']}) | merge {ms:.0f}s wall, {cpu / 3600:.2f} CPU-h | {stats}")


if __name__ == "__main__":
    main()
