#!/usr/bin/env python
"""L2-genome: gated synteny repair (GSR) arms on OrthoBench DEV (accurate base), docs/L2_genome.md.

    python scripts/L2_gsr_ob.py --arms off F_c F_c_s11 ... [--threads 8]

Base = frozen `accurate` preset with gene order, as the dev reference row
`T[det] tree=on:ft-medoid off:nj syn constrained+veto+r3` (run_det; OG-F 79.69, macro 0.8026):
two-pass constrained synteny (off pass NJ), ctx_veto + ctx_r3, FastTree medoid local trees.
The synteny-off pass is computed once per process and shared by all arms.
Arm suffix _sNN: the GSR support (only) comes from shuffled gene order (seed NN); the frozen R1
bonus / veto stay real. Outputs: results/l2_genome/ob_arms/<arm>.json (+ hogs.orthoxml,
orthobench_ogs.txt, gsr_log.json under work/l2_genome/ob/<arm>/). Dev RefOGs only.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from dataclasses import asdict, replace
from datetime import datetime

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import ancestra.merge as M  # noqa: E402
from ancestra.io import load_guide_tree, load_proteomes, read_gene_table  # noqa: E402
from ancestra.merge import ConstrainedContext, MergeConfig, root_labels, run_merge, to_orthoxml  # noqa: E402
from ancestra.orthoxml import to_orthobench_ogs, write_orthoxml  # noqa: E402
from ancestra.pipeline import SYNTENY_DEFAULT, apply_preset  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from ancestra.synteny import SyntenyContext, gene_positions  # noqa: E402

RUN = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/work/core_v0_ob/run_det"
OUTJ = os.path.join(ROOT, "results/l2_genome/ob_arms")
OUTW = os.path.join(ROOT, "work/l2_genome/ob")

ARMS = {
    "off": {},
    "F_c": dict(gsr="single", gsr_mask=True),
    "F_u": dict(gsr="single", gsr_mask=False),
    "R_c": dict(gsr="any", gsr_mask=True),
    "R_u": dict(gsr="any", gsr_mask=False),
    "locus_ref": dict(locus_repair=True, locus_min_support=2.0),   # prior art, information only
}


class GSRShuffled:
    """Real constrained synteny for R1 / veto / R3; GSR support from a shuffled-order context."""

    def __init__(self, cons, shuf, root_lab):
        self.cons, self.shuf, self.root_lab = cons, shuf, root_lab

    def __call__(self, ctx):
        return self.cons(ctx)

    def support(self, ctx):
        return self.cons.support(ctx)

    def locus_support(self, ctx):
        return self.cons.locus_support(ctx)

    def gsr_support(self, ctx, masked):
        s = np.asarray(self.shuf.support(ctx), dtype=np.float64)
        return s * self.cons._same(ctx) if masked else s


def git_commit():
    try:
        return subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True,
                              text=True).stdout.strip()
    except Exception:
        return "?"


def evaluate(ox, gt, out):
    from ob_merge_iter import s8_fission_candidates
    from orthobench_dev_eval import dev_eval
    from orthobench_dev_pairs_eval import evaluate as pair_eval
    ogs = os.path.join(out, "orthobench_ogs.txt")
    to_orthobench_ogs(ox, None, ogs)
    ev = dev_eval(ogs)
    s8 = s8_fission_candidates(ox, gt)
    pm, parr, prefogs = pair_eval(ox)
    per = parr.sum(axis=1)
    tp, fn, fp = per[:, 0].astype(float), per[:, 1].astype(float), per[:, 2].astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        f1_r = 2 * tp / (2 * tp + fp + fn)
    return dict(OG_F=ev["F"], OG_P=ev["precision"], OG_R=ev["recall"], exact=ev["exact"],
                n_refogs_split=ev["n_refogs_split"], pair_micro=pm["all"],
                pair_macro_F1=float(np.nanmean(f1_r)), s8_5of6=s8[0], s8_fission_cand=s8[1],
                pair_classes={k: v for k, v in pm.items() if k != "all"},
                per_refog={r: [int(x) for x in per[i]] for i, r in enumerate(prefogs)},
                per_refog_class_counts=parr.tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="+", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--set", nargs="*", default=[], help="extra MergeConfig overrides k=v (all arms)")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(OUTJ, exist_ok=True)
    gt = read_gene_table(os.path.join(RUN, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(RUN, "stages/input/guide_tree.nwk"), gt.species)
    st = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
    man = json.load(open(os.path.join(RUN, "manifest.json")))
    pdir = man.get("argv_proteomes") or os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/Input")
    full = load_proteomes(pdir)
    assert full.names == gt.names
    gt.seqs = full.seqs
    fam = np.load(os.path.join(RUN, "stages/scores/families_final.npy"))
    stems = [sp.split(".")[0] for sp in gt.species]
    cc, oo, _ = gene_positions(gt.names, gt.species_of, stems, os.path.join(ROOT, "data/genomes/orthobench"))
    base = apply_preset(MergeConfig(), "accurate", a.threads)
    base = replace(base, ctx_veto=True, ctx_r3=True,
                   **{k: json.loads(v) for k, v in (kv.split("=", 1) for kv in a.set)})
    off_cfg = replace(base, tree_method="nj")
    t0 = time.time()
    off = run_merge(gt, tree, st, off_cfg)
    root_lab = root_labels(off, gt.n)
    t_off = time.time() - t0
    print(f"off pass {t_off:.0f}s", flush=True)
    real = SyntenyContext(cc, oo, family=fam, **SYNTENY_DEFAULT)
    for arm in a.arms:
        name, _, seed = arm.partition("_s")
        seed = int(seed) if seed else None
        cfg = replace(base, **ARMS[name])
        cons = ConstrainedContext(real, root_lab)
        ctx = cons
        if seed is not None:
            shuf = SyntenyContext(cc, oo, family=fam, shuffle_seed=seed, **SYNTENY_DEFAULT)
            ctx = GSRShuffled(cons, shuf, root_lab)
        M.GSR_LOG = [] if cfg.gsr != "off" else None
        c0, t1 = os.times(), time.time()
        res = run_merge(gt, tree, st, cfg, context=ctx)
        c1 = os.times()
        ms = time.time() - t1
        cpu = sum(getattr(c1, k) - getattr(c0, k) for k in ("user", "system", "children_user", "children_system"))
        ox = to_orthoxml(res, gt, tree)
        lbl = arm + a.tag
        out = os.path.join(OUTW, lbl)
        os.makedirs(out, exist_ok=True)
        write_orthoxml(ox, os.path.join(out, "hogs.orthoxml"))
        ev = evaluate(ox, gt, out)
        gs = {k: int(sum(s.get(k, 0) for s in res.node_stats)) for k in
              ("gsr_repaired", "gsr_eligible", "gsr_plans", "locus_repaired", "tree_fasttree_jobs", "r3_repaired")}
        if M.GSR_LOG is not None:
            names = gt.names
            log = [dict(r, u=[names[i] for i in r["u"]], g=[names[i] for i in r["g"]]) for r in M.GSR_LOG]
            json.dump(log, open(os.path.join(out, "gsr_log.json"), "w"))
        rec = dict(arm=lbl, time=datetime.now().isoformat(timespec="seconds"), **ev, stats=gs,
                   merge_wall_s=round(ms, 1), merge_cpu_s=round(cpu, 1), off_pass_wall_s=round(t_off, 1),
                   manifest=dict(git_commit=git_commit(), run=RUN, config=asdict(cfg),
                                 off_config=asdict(off_cfg), synteny=SYNTENY_DEFAULT, gsr_shuffle_seed=seed,
                                 inputs=dict(scores_md5=hashlib.md5(open(os.path.join(RUN, "stages/scores/scores.npz"), "rb").read()).hexdigest()),
                                 ft_cache=os.environ.get("ANCESTRA_FT_CACHE"), threads=a.threads,
                                 orthoxml=os.path.join(out, "hogs.orthoxml")))
        json.dump(rec, open(os.path.join(OUTJ, lbl + ".json"), "w"), indent=1)
        print(f"{lbl}: OG-F {ev['OG_F']:.2f} P {ev['OG_P']:.2f} R {ev['OG_R']:.2f} | pair micro P "
              f"{ev['pair_micro']['P']} R {ev['pair_micro']['R']} F1 {ev['pair_micro']['F1']} macro "
              f"{ev['pair_macro_F1']:.4f} | S8 cand {ev['s8_fission_cand']} | {gs} | {ms:.0f}s", flush=True)


if __name__ == "__main__":
    main()
