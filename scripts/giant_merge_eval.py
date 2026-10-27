"""Merge (exactly the pipeline's merge stage with gene order: constrained synteny two-pass + veto + R3,
preset applied) on given scores/families, then DEV evaluation: QfO -> TreeFam-A dev (fast scorer, dev
families only) + FN breakdown; OrthoBench -> dev OG-level + pair P/R/F1 micro/macro.
Appends to results/giant/runs.tsv.

    python scripts/giant_merge_eval.py {qfo,ob} --preset fast --label G1 --scores S.npz --families F.npy
"""
import argparse
import json
import os
import resource
import sys
import time
from dataclasses import replace
from datetime import datetime

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import load_guide_tree, load_proteomes, read_gene_table  # noqa
from ancestra.merge import MergeConfig, run_merge_constrained, to_orthoxml  # noqa
from ancestra.orthoxml import to_orthobench_ogs, write_orthoxml  # noqa
from ancestra.pipeline import PRESETS, SYNTENY_DEFAULT, apply_preset  # noqa
from ancestra.scoring import ScoreTable  # noqa
from ancestra.synteny import SyntenyContext, gene_positions  # noqa

B = {"qfo": dict(run="work/rescue/qfo/base", prot="data/inputs/qfo2022/proteomes", go="data/inputs/qfo2022/gene_order"),
     "ob": dict(run="work/core_v0_ob/run_det", prot="data/orthobench/Open_Orthobench/BENCHMARKS/Input",
                go="data/genomes/orthobench")}
LOG = os.path.join(ROOT, "results/giant/runs.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bench", choices=["qfo", "ob"])
    ap.add_argument("--preset", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--scores")
    ap.add_argument("--families")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    C = B[a.bench]
    run = os.path.join(ROOT, C["run"])
    out = os.path.join(ROOT, "work/giant", a.bench, "runs", f"{a.label}_{a.preset}")
    os.makedirs(out, exist_ok=True)
    gt = read_gene_table(os.path.join(run, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(run, "stages/input/guide_tree.nwk"), gt.species)
    st = ScoreTable.load(a.scores or os.path.join(run, "stages/scores/scores.npz"))
    fam = np.load(a.families or os.path.join(run, "stages/scores/families_final.npy"))
    merge = apply_preset(MergeConfig(), a.preset, a.threads)
    if merge.tree_method == "fasttree":
        full = load_proteomes(os.path.join(ROOT, C["prot"]))
        assert full.names == gt.names
        gt.seqs = full.seqs
    go = os.path.join(ROOT, C["go"])
    names = [sp if os.path.exists(os.path.join(go, sp + ".parquet")) else sp.split(".")[0] for sp in gt.species]
    cc, oo, _ = gene_positions(gt.names, gt.species_of, names, go)
    ctx = SyntenyContext(cc, oo, family=fam, **SYNTENY_DEFAULT)
    mcfg = replace(merge, ctx_veto=True, ctx_r3=True)
    t0, c0 = time.time(), os.times()
    res = run_merge_constrained(gt, tree, st, mcfg, ctx, off_cfg=replace(mcfg, tree_method=PRESETS[a.preset]["off_tree_method"]))
    c1 = os.times()
    wall = time.time() - t0
    cpu = sum(getattr(c1, k) - getattr(c0, k) for k in ("user", "system", "children_user", "children_system"))
    ox = to_orthoxml(res, gt, tree)
    xml = os.path.join(out, "hogs.orthoxml")
    write_orthoxml(ox, xml)
    to_orthobench_ogs(ox, None, os.path.join(out, "orthobench_ogs.txt"))
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    row = dict(time=datetime.now().isoformat(timespec="seconds"), bench=a.bench, label=a.label, preset=a.preset,
               merge_wall_s=round(wall, 1), merge_cpu_s=round(cpu, 1), peak_rss_gb=round(rss, 1),
               scores=a.scores, families=a.families)
    if a.bench == "qfo":
        from qfo_merge_iter import fn_breakdown
        from qfo_treefam_dev_fast import score as tf_score
        tf = tf_score(ox)
        fnb = fn_breakdown(ox, gt, fam, tree)
        row.update(TPR=round(tf["TPR"], 4), PPV=round(tf["PPV"], 4), n_TP=tf["n_TP"], n_FP=tf["n_FP"],
                   n_FN=tf["n_FN"], fn_classes={k: sum(v.values()) for k, v in fnb.items()},
                   fn_breakdown=fnb)
    else:
        from orthobench_dev_eval import dev_eval
        from orthobench_dev_pairs_eval import evaluate as pair_eval
        ev = dev_eval(os.path.join(out, "orthobench_ogs.txt"))
        pm, arr, refogs = pair_eval(xml)
        per = arr.sum(axis=1)
        tp, fn, fp = per[:, 0].astype(float), per[:, 1].astype(float), per[:, 2].astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            f1 = 2 * tp / (2 * tp + fp + fn)
        row.update(OG_F=round(ev["F"], 2), OG_P=round(ev["precision"], 2), OG_R=round(ev["recall"], 2),
                   exact=ev["exact"], pairF1=pm["all"]["F1"], pairF1_macro=round(float(np.nanmean(f1)), 4))
        json.dump(ev, open(os.path.join(out, "dev_eval.json"), "w"), indent=1)
    json.dump(row, open(os.path.join(out, "result.json"), "w"), indent=1, default=str)
    short = {k: v for k, v in row.items() if k not in ("fn_breakdown", "scores", "families")}
    with open(LOG, "a") as fh:
        fh.write(json.dumps(short, default=str) + "\n")
    print(json.dumps(short, default=str))


if __name__ == "__main__":
    main()
