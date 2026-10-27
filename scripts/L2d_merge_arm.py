"""L2-domains step (b): one merge arm on cached stages (QfO fast / OrthoBench accurate), DEV scoring only.

    python scripts/L2d_merge_arm.py --bench qfo --arm gate --shuffle 11 --domtbl D [--domtbl D2 ...] --out X.json

Replicates the pipeline merge with gene order (constrained two-pass synteny + veto + R3 choice, as
ancestra/pipeline.py s_merge and scripts/qfo_merge_iter.py) and passes an ancestra.domains.DomainArch
(`arch=`) with MergeConfig.arch_mode = arm. QfO: TreeFam-A DEV counts per family (paired bootstrap later);
OrthoBench: DEV OG-F (official benchmark on dev RefOGs) + pair counts per RefOG.
Pre-registration: docs/L2_domains.md §0.4.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from collections import Counter
from dataclasses import asdict, replace

import numpy as np

ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
WT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, WT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.domains import DomainArch, load_clans, parse_domtbl  # noqa: E402
from ancestra.io import load_guide_tree, load_proteomes, read_gene_table  # noqa: E402
from ancestra.merge import MergeConfig, run_merge_constrained, to_orthoxml  # noqa: E402
from ancestra.pipeline import PRESETS, SYNTENY_DEFAULT, apply_preset  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from ancestra.synteny import SyntenyContext, gene_positions  # noqa: E402

CLANS = os.path.join(ROOT, "data/extra/l2_domains/Pfam-A.clans.tsv.gz")
BENCH = {
    "qfo": dict(run=os.path.join(ROOT, "work/qfo_recall/base"), preset="fast",
                go=os.path.join(ROOT, "data/inputs/qfo2022/gene_order"),
                prot=os.path.join(ROOT, "data/inputs/qfo2022/proteomes")),
    "ob": dict(run=os.path.join(ROOT, "work/core_v0_ob/run_det"), preset="accurate",
               go=os.path.join(ROOT, "data/genomes/orthobench"),
               prot=os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/Input")),
}


def tf_counts(ox):
    """TreeFam-A DEV counts per family: {family: [TP, FN, FP, TN]} (raw pair records)."""
    from qfo_treefam_dev_fast import dev_refpairs, predicted_orthologs
    refs = dev_refpairs()
    genes = {x for r in refs for x in r[:2]}
    pred = predicted_orthologs(ox, genes)
    c: dict = {}
    for a, b, t, f in refs:
        p = ((a, b) if a < b else (b, a)) in pred
        k = (0 if p else 1) if t == "S" else (2 if p else 3)
        c.setdefault(f, [0, 0, 0, 0])[k] += 1
    return c


def ob_eval(ox, outdir):
    from ancestra.orthoxml import to_orthobench_ogs
    from orthobench_dev_eval import dev_eval
    from orthobench_dev_pairs_eval import evaluate as pair_eval
    ogs = os.path.join(outdir, "orthobench_ogs.txt")
    to_orthobench_ogs(ox, None, ogs)
    ev = dev_eval(ogs)
    pm, parr, refogs = pair_eval(ox)
    per = parr.sum(axis=1)
    return dict(OG_F=ev["F"], OG_P=ev["precision"], OG_R=ev["recall"], exact=ev["exact"],
                pairP=pm["all"]["P"], pairR=pm["all"]["R"], pairF1=pm["all"]["F1"],
                per_refog={r: [int(x) for x in per[i]] for i, r in enumerate(refogs)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", choices=list(BENCH), required=True)
    ap.add_argument("--arm", choices=["off", "bonus", "gate", "gate_deep", "veto"], required=True)
    ap.add_argument("--shuffle", default="none")
    ap.add_argument("--domtbl", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--save-orthoxml")
    a = ap.parse_args()
    B = BENCH[a.bench]
    t0 = time.time()
    c0 = os.times()
    cfg = apply_preset(MergeConfig(), B["preset"], a.threads)
    cfg = replace(cfg, **{k: json.loads(v) for k, v in (kv.split("=", 1) for kv in a.set)})
    mcfg = replace(cfg, ctx_veto=True, ctx_r3=True, arch_mode=a.arm)
    off_cfg = replace(mcfg, tree_method=PRESETS[B["preset"]]["off_tree_method"])
    run = B["run"]
    gt = read_gene_table(os.path.join(run, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(run, "stages/input/guide_tree.nwk"), gt.species)
    st = ScoreTable.load(os.path.join(run, "stages/scores/scores.npz"))
    fam = np.load(os.path.join(run, "stages/scores/families_final.npy"))
    if mcfg.tree_method == "fasttree":
        full = load_proteomes(B["prot"])
        assert full.names == gt.names
        gt.seqs = full.seqs
    names = [sp if os.path.exists(os.path.join(B["go"], sp + ".parquet")) else sp.split(".")[0]
             for sp in gt.species]
    cc, oo, _ = gene_positions(gt.names, gt.species_of, names, B["go"])
    ctx = SyntenyContext(cc, oo, family=fam, **SYNTENY_DEFAULT)
    arch = None
    ann = {}
    if a.arm != "off":
        clans = load_clans(CLANS)
        for d in a.domtbl:
            ann.update(parse_domtbl(d, clans))
        seed = None if a.shuffle == "none" else int(a.shuffle)
        arch = DomainArch(gt.names, gt.species_of, ann, shuffle_seed=seed)
    t1 = time.time()
    res = run_merge_constrained(gt, tree, st, mcfg, ctx, off_cfg=off_cfg, arch=arch)
    merge_s = time.time() - t1
    ox = to_orthoxml(res, gt, tree)
    outdir = os.path.splitext(a.out)[0] + "_files"
    os.makedirs(outdir, exist_ok=True)
    if a.save_orthoxml:
        from ancestra.orthoxml import write_orthoxml
        write_orthoxml(ox, a.save_orthoxml)
    stats = Counter()
    for s in res.node_stats:
        for k in ("arch_bonus", "arch_gate_added", "arch_veto", "n_mb", "r3_repaired", "so_absorbed"):
            stats[k] += int(s.get(k, 0))
    out = dict(bench=a.bench, arm=a.arm, shuffle=a.shuffle, stats=dict(stats))
    if a.bench == "qfo":
        c = tf_counts(ox)
        tot = np.array(list(c.values())).sum(axis=0)
        tp, fn, fp, tn = (int(x) for x in tot)
        out.update(TF_TP=tp, TF_FN=fn, TF_FP=fp, TF_TN=tn, TF_TPR=tp / (tp + fn), TF_PPV=tp / (tp + fp),
                   TF_F1=2 * tp / (2 * tp + fp + fn), per_family=c)
        # official-formula values (pseudocount, /2) for comparison with results/core_v0_qfo/runs.tsv
        TP, FP, FN = 1 + tp / 2, 1 + fp / 2, 1 + fn / 2
        out.update(TF_TPR_official=TP / (TP + FN), TF_PPV_official=TP / (TP + FP))
    else:
        out.update(ob_eval(ox, outdir))
    c1 = os.times()
    git = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain", "ancestra"], capture_output=True,
                           text=True).stdout.strip()
    out["manifest"] = dict(git_commit=git, ancestra_dirty=bool(dirty), argv=sys.argv, run=run,
                           preset=B["preset"], config=mcfg.key(), off_tree_method=off_cfg.tree_method,
                           synteny=SYNTENY_DEFAULT, domtbl=a.domtbl, n_annotated=int(arch.n_annotated) if arch
                           else 0, merge_s=round(merge_s, 1), wall_s=round(time.time() - t0, 1),
                           cpu_s=round((c1.user - c0.user) + (c1.system - c0.system) + (c1.children_user -
                                       c0.children_user) + (c1.children_system - c0.children_system), 1))
    json.dump(out, open(a.out, "w"), indent=1)
    keys = [k for k in out if k.startswith(("TF_", "OG_", "pair"))]
    print(a.bench, a.arm, a.shuffle, {k: (round(out[k], 4) if isinstance(out[k], float) else out[k]) for k in keys},
          dict(stats), f"{merge_s:.0f}s")


if __name__ == "__main__":
    main()
