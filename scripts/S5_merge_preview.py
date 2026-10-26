#!/usr/bin/env python
"""S5 integration preview on toy sims: unmodified core merge (run_merge) +/- synteny context.

Arms: off; multiplicative (lam 0.05/0.1/0.2, cap 3, w 5); additive (lam 0.05); shuffled-order
control for the main arm. Scores: cached core_v0 DIAMOND e10 tables (work/core_v0/runs/<sim>/truth),
truth families (isolates the merge). Metric: ortholog-pair P/R/F1 (evaluate_toy).
Output: results/S5_synteny/merge_preview_toy.json
"""
import argparse, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.evaluate_toy import evaluate  # noqa: E402
from ancestra.genomes import assign_order, parse_gff3_proteins, write_gene_order  # noqa: E402
from ancestra.io import load_guide_tree, load_proteomes, load_families  # noqa: E402
from ancestra.merge import ARMS, run_merge, to_orthoxml  # noqa: E402
from ancestra.orthoxml import read_orthoxml  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from ancestra.genomes import to_gene_order  # noqa: E402
from ancestra.synteny import SyntenyContext, gene_positions  # noqa: E402

ARMS_S5 = {"off": None,
           "ovr_m1": dict(mode="override", margin=1),
           "ovr_m2": dict(mode="override", margin=2),
           "ovr_m1_near0.75": dict(mode="override", margin=1, near=0.75),
           "ovr_m1_shuf11": dict(mode="override", margin=1, shuffle_seed=11),
           "ovr_m1_shuf23": dict(mode="override", margin=1, shuffle_seed=23),
           "mult0.05": dict(mode="multiplicative", lam=0.05),
           "mult0.1": dict(mode="multiplicative", lam=0.1),
           "mult0.2": dict(mode="multiplicative", lam=0.2),
           "add0.05": dict(mode="additive", lam=0.05),
           "mult0.1_shuf11": dict(mode="multiplicative", lam=0.1, shuffle_seed=11),
           "mult0.1_shuf23": dict(mode="multiplicative", lam=0.1, shuffle_seed=23),
           "mult0.1_coreadj": dict(mode="multiplicative", lam=0.1, core_adjacency=True),
           "mult0.1_coreanch": dict(mode="multiplicative", lam=0.1, anchor="core")}


def one(sim):
    sd = f"{ROOT}/work/core_v0/sims/{sim}"
    gt = load_proteomes(f"{sd}/proteomes")
    tree = load_guide_tree(f"{sd}/species_tree.nwk", gt.species)
    truth = read_orthoxml(f"{sd}/truth/hogs.orthoxml")
    fam = load_families(f"{sd}/truth/gene_families.tsv", gt, 0, 2)
    st = ScoreTable.load(f"{ROOT}/work/core_v0/runs/{sim}/truth/scores_diamond_e10.npz")
    tdir = f"{ROOT}/work/S5/integ/{sim}/genomes"
    for s in gt.species:
        if not os.path.exists(f"{tdir}/{s}.parquet"):
            df = parse_gff3_proteins(f"{sd}/gff/{s}.gff3").assign(primary=True, source="toy")
            write_gene_order(assign_order(df), f"{tdir}/{s}.parquet")
    cc, oo, cov = gene_positions(gt.names, gt.species_of, gt.species, tdir)
    out = []
    for arm, kw in ARMS_S5.items():
        t0 = time.time()
        if kw is None:
            res = run_merge(gt, tree, st, ARMS["default"])
            stats = None
        else:
            kw = dict(kw)
            if kw.pop("core_adjacency", False):   # adjacency computed by the core from GeneOrder
                ctx = SyntenyContext(None, None, family=fam, **kw)
                res = run_merge(gt, tree, st, ARMS["default"], context=ctx,
                                gene_order=to_gene_order(gt, tdir), adjacency_w=5)
            else:
                ctx = SyntenyContext(cc, oo, family=fam, w=5, **kw)
                res = run_merge(gt, tree, st, ARMS["default"], context=ctx)
            stats = dict(n_anchors=int(sum(s["n_anchors"] for s in ctx.stats)),
                         n_cand=int(sum(s["n_cand"] for s in ctx.stats)),
                         n_supported=int(sum(s["n_supported"] for s in ctx.stats)))
        r = evaluate(to_orthoxml(res, gt, tree), truth)
        p = r["pairs"]
        rec = dict(sim=sim, arm=arm, P=p["precision"], R=p["recall"], F1=p["f1"],
                   fp_paralogs=p["fp_true_paralogs"], wall_s=round(time.time() - t0, 2), ctx=stats)
        print(json.dumps(rec), flush=True)
        out.append(rec)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", nargs="+", default=["def8_s11", "def8_s23", "def8_s47", "hard8_s11", "hard8_s23",
                                                  "hard8_s47", "def32_s11", "hard32_s11"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        rows = [r for rs in ex.map(one, a.sims) for r in rs]
    summ = {}
    for arm in ARMS_S5:
        rr = [r for r in rows if r["arm"] == arm]
        off = {r["sim"]: r for r in rows if r["arm"] == "off"}
        summ[arm] = dict(F1=round(float(np.mean([r["F1"] for r in rr])), 4),
                         P=round(float(np.mean([r["P"] for r in rr])), 4),
                         R=round(float(np.mean([r["R"] for r in rr])), 4),
                         dF1_vs_off=[round(r["F1"] - off[r["sim"]]["F1"], 4) for r in rr])
    json.dump(dict(rows=rows, summary=summ, wall_s=round(time.time() - t0, 1)),
              open(f"{ROOT}/results/S5_synteny/merge_preview_toy.json", "w"), indent=1)
    for k, v in summ.items():
        print(k, v)


if __name__ == "__main__":
    main()
