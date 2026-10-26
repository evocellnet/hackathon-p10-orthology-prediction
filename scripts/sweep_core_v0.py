"""Sweep merge parameters on cached hybrid scores (truth families), pair metrics only.

Pre-stated rule (2026-10-26 ~16:15 CET, before running): for each parameter,
choose the value maximising mean pair F1 over the 12 truth-family toy
conditions subject to mean precision >= (precision of the current default) - 0.003.

    python scripts/sweep_core_v0.py > results/core_v0/sweep.tsv
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.evaluate_toy import pair_metrics  # noqa: E402
from ancestra.io import load_guide_tree, load_proteomes  # noqa: E402
from ancestra.merge import MergeConfig, run_merge, to_orthoxml  # noqa: E402
from ancestra.orthoxml import hogs_to_pairs, read_orthoxml  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402

SIMS = os.path.join(ROOT, "work/core_v0/sims")
RUNS = os.path.join(ROOT, "work/core_v0/runs")
D = MergeConfig()
GRID_NAME = os.environ.get("SWEEP_GRID", "base")
SCORES = os.environ.get("SWEEP_SCORES", "scores_hybrid_e10.npz")
GRID = [("default", D)]
for v in (0.8, 0.9, 0.95, 1.05):
    GRID.append((f"inpar_factor={v}", replace(D, inpar_factor=v)))
for v in (0.9, 0.98, 1.0):
    GRID.append((f"tol={v}", replace(D, tol=v)))
for v in (0.3, 0.5, 0.7, 0.9):
    GRID.append((f"repair_ratio={v}", replace(D, repair=True, repair_ratio=v)))
for v in (0.3, 0.0):
    GRID.append((f"min_cov={v}", replace(D, min_cov=v)))
if GRID_NAME == "rd":   # pre-stated rule as above (P >= default P - 0.003, max F1)
    GRID = [("default", D)]
    for z in (2.0, 3.0, 4.0, 6.0):
        GRID.append((f"rd_hidden z={z}", replace(D, rd_hidden=True, rd_z=z)))
    for f in (0.8, 1.0, 1.2):
        GRID.append((f"rd_inpar f={f}", replace(D, rd_inpar=True, rd_inpar_factor=f)))
    for z, f in ((3.0, 1.0), (4.0, 1.0), (3.0, 0.8)):
        GRID.append((f"rd_both z={z} f={f}", replace(D, rd_hidden=True, rd_z=z, rd_inpar=True,
                                                      rd_inpar_factor=f)))


def one(sim):
    d = os.path.join(SIMS, sim)
    gt = load_proteomes(os.path.join(d, "proteomes"))
    tree = load_guide_tree(os.path.join(d, "species_tree.nwk"), gt.species)
    st = ScoreTable.load(os.path.join(RUNS, sim, "truth", SCORES))
    truth = set(hogs_to_pairs(read_orthoxml(os.path.join(d, "truth", "hogs.orthoxml"))))
    out = []
    for name, cfg in GRID:
        res = run_merge(gt, tree, st, cfg)
        pm = pair_metrics(hogs_to_pairs(to_orthoxml(res, gt, tree)), truth)
        out.append((sim, name, pm["precision"], pm["recall"], pm["f1"]))
    return out


def main():
    sims = sorted(os.listdir(SIMS))
    rows = []
    with ProcessPoolExecutor(12) as ex:
        for r in ex.map(one, sims):
            rows.extend(r)
    print("sim\tparam\tP\tR\tF1")
    for r in rows:
        print("\t".join([r[0], r[1]] + [f"{x:.5f}" for x in r[2:]]))
    print("\n#param\tmeanP\tmeanR\tmeanF1\t(n sims)")
    names = [g[0] for g in GRID]
    for nm in names:
        rs = [r for r in rows if r[1] == nm]
        print(f"#{nm}\t{np.mean([r[2] for r in rs]):.5f}\t{np.mean([r[3] for r in rs]):.5f}\t"
              f"{np.mean([r[4] for r in rs]):.5f}\t{len(rs)}")


if __name__ == "__main__":
    main()
