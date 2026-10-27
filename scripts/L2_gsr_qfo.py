#!/usr/bin/env python
"""L2-genome: GSR arms on QfO 2022 cached stages, scored on TreeFam-A DEV only (docs/L2_genome.md).

    python scripts/L2_gsr_qfo.py --preset fast --arms off R_c R_c_s11 ... [--threads 8]

Replicates scripts/qfo_merge_iter.py (preset merge with gene order: constrained two-pass synteny +
veto + R3, off pass = the preset's off_tree_method). Synteny-off pass computed once per process.
Arm suffix _sNN: GSR support from shuffled gene order (seed NN), everything else real.
Outputs results/l2_genome/qfo_arms/<preset>_<arm>.json with per-family TreeFam-A DEV counts
(for the paired family bootstrap) and work/l2_genome/qfo/<preset>_<arm>/hogs.orthoxml (with --save).
"""
import argparse
import json
import os
import sys
import time
from collections import Counter
from dataclasses import asdict, replace
from datetime import datetime

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import ancestra.merge as M  # noqa: E402
from ancestra.io import load_guide_tree, read_gene_table  # noqa: E402
from ancestra.merge import ConstrainedContext, MergeConfig, root_labels, run_merge, to_orthoxml  # noqa: E402
from ancestra.pipeline import PRESETS, SYNTENY_DEFAULT, apply_preset  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from ancestra.synteny import SyntenyContext, gene_positions  # noqa: E402
from L2_gsr_ob import ARMS, GSRShuffled, git_commit  # noqa: E402
from qfo_treefam_dev_fast import dev_refpairs, predicted_orthologs  # noqa: E402

MAIN = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
RUN = os.path.join(MAIN, "work/qfo_recall/base")
GO = os.path.join(MAIN, "data/inputs/qfo2022/gene_order")
PROT = os.path.join(MAIN, "data/inputs/qfo2022/proteomes")
OUTJ = os.path.join(ROOT, "results/l2_genome/qfo_arms")
OUTW = os.path.join(ROOT, "work/l2_genome/qfo")


def tf_counts(ox, refs):
    genes = {x for r in refs for x in r[:2]}
    pred = predicted_orthologs(ox, genes)
    per = {}
    for a, b, t, f in refs:
        p = ((a, b) if a < b else (b, a)) in pred
        k = ("TP" if p else "FN") if t == "S" else ("FP" if p else "TN")
        per.setdefault(f, Counter())[k] += 1
    tot = Counter()
    for c in per.values():
        tot.update(c)
    tp, fp, fn, tn = (1 + tot[k] / 2 for k in ("TP", "FP", "FN", "TN"))
    return dict(TPR=tp / (tp + fn), PPV=tp / (tp + fp), n=dict(tot)), \
        {f: [c["TP"], c["FN"], c["FP"], c["TN"]] for f, c in per.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="fast", choices=["fast", "accurate"])
    ap.add_argument("--arms", nargs="+", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--save", action="store_true")
    ap.add_argument("--log", action="store_true", help="write GSR decisions touching TreeFam-A DEV genes")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(OUTJ, exist_ok=True)
    cfg0 = apply_preset(MergeConfig(), a.preset, a.threads)
    base = replace(cfg0, ctx_veto=True, ctx_r3=True)
    off_cfg = replace(base, tree_method=PRESETS[a.preset]["off_tree_method"])
    gt = read_gene_table(os.path.join(RUN, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(RUN, "stages/input/guide_tree.nwk"), gt.species)
    st = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
    fam = np.load(os.path.join(RUN, "stages/scores/families_final.npy"))
    if base.tree_method == "fasttree":
        from ancestra.io import load_proteomes
        full = load_proteomes(PROT)
        assert full.names == gt.names
        gt.seqs = full.seqs
    names = [sp if os.path.exists(os.path.join(GO, sp + ".parquet")) else sp.split(".")[0] for sp in gt.species]
    cc, oo, _ = gene_positions(gt.names, gt.species_of, names, GO)
    refs = dev_refpairs()
    t0 = time.time()
    off = run_merge(gt, tree, st, off_cfg)
    root_lab = root_labels(off, gt.n)
    t_off = time.time() - t0
    del off
    print(f"off pass {t_off:.0f}s", flush=True)
    real = SyntenyContext(cc, oo, family=fam, **SYNTENY_DEFAULT)
    for arm in a.arms:
        name, _, seed = arm.partition("_s")
        seed = int(seed) if seed else None
        cfg = replace(base, **ARMS[name])
        cons = ConstrainedContext(real, root_lab)
        ctx = cons
        if seed is not None:
            ctx = GSRShuffled(cons, SyntenyContext(cc, oo, family=fam, shuffle_seed=seed, **SYNTENY_DEFAULT),
                              root_lab)
        M.GSR_LOG = [] if (a.log and cfg.gsr != "off") else None
        c0, t1 = os.times(), time.time()
        res = run_merge(gt, tree, st, cfg, context=ctx)
        c1 = os.times()
        ms = time.time() - t1
        cpu = sum(getattr(c1, k) - getattr(c0, k) for k in ("user", "system", "children_user", "children_system"))
        ox = to_orthoxml(res, gt, tree)
        lbl = f"{a.preset}_{arm}{a.tag}"
        if M.GSR_LOG is not None:
            refg = {x for r in refs for x in r[:2]}
            nm = gt.names
            keep = []
            for r in M.GSR_LOG:
                u = [nm[i] for i in r["u"]]
                g = [nm[i] for i in r["g"]]
                if refg.intersection(u) and refg.intersection(g):
                    keep.append(dict(r, u=[x for x in u if x in refg], g=[x for x in g if x in refg]))
            os.makedirs(os.path.join(OUTW, lbl), exist_ok=True)
            json.dump(dict(n_decisions=len(M.GSR_LOG), rows=keep), open(os.path.join(OUTW, lbl, "gsr_log.json"), "w"))
            M.GSR_LOG = None
        if a.save:
            from ancestra.orthoxml import write_orthoxml
            os.makedirs(os.path.join(OUTW, lbl), exist_ok=True)
            write_orthoxml(ox, os.path.join(OUTW, lbl, "hogs.orthoxml"))
        tf, per = tf_counts(ox, refs)
        gs = {k: int(sum(s.get(k, 0) for s in res.node_stats)) for k in
              ("gsr_repaired", "gsr_eligible", "gsr_plans", "locus_repaired", "so_absorbed", "r3_repaired",
               "tree_fasttree_jobs")}
        gsr_by_node = {s["node"]: s.get("gsr_repaired", 0) for s in res.node_stats if s.get("gsr_repaired")}
        rec = dict(arm=lbl, time=datetime.now().isoformat(timespec="seconds"), TF_dev=tf, per_family=per,
                   stats=gs, gsr_by_node=gsr_by_node, merge_wall_s=round(ms, 1), merge_cpu_s=round(cpu, 1),
                   off_pass_wall_s=round(t_off, 1),
                   manifest=dict(git_commit=git_commit(), run=RUN, config=asdict(cfg), off_config=asdict(off_cfg),
                                 synteny=SYNTENY_DEFAULT, gsr_shuffle_seed=seed, gene_order=GO,
                                 ft_cache=os.environ.get("ANCESTRA_FT_CACHE"), threads=a.threads))
        json.dump(rec, open(os.path.join(OUTJ, lbl + ".json"), "w"), indent=1)
        print(f"{lbl}: TreeFam-A DEV TPR {tf['TPR']:.4f} PPV {tf['PPV']:.4f} {tf['n']} | {gs} | {ms:.0f}s "
              f"{cpu / 3600:.2f} CPU-h", flush=True)


if __name__ == "__main__":
    main()
