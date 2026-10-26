#!/usr/bin/env python
"""S4 cost at localtree cluster size: CPU per tree for random 5-protein sets (the mean
localtree cluster has 79,782 / 15,456 = 5.2 units in the OrthoBench e2e default run).
100 sets: a DEV RefOG drawn uniformly, then 5 of its 3Di-covered genes (np.random.default_rng(11)).
Same arm code as scripts/S4_trees.py (do_task); a_nj gets the pipeline scores.
    python scripts/S4_cost.py --workers 6 -> results/S4_trees/cost_5units.json
"""
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import S4_trees as T  # noqa: E402
from ancestra.io import load_proteomes  # noqa: E402
from ancestra import structure as STR  # noqa: E402

N_CLUSTERS = 15456


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--n", type=int, default=100)
    a = ap.parse_args()
    T.W = os.path.join(ROOT, "work/S4/cost")
    split = json.load(open(f"{ROOT}/splits/orthobench_split.json"))
    gp = load_proteomes(f"{T.OB}/BENCHMARKS/Input")
    idx = gp.index()
    ss = STR.load_3di(f"{ROOT}/data/structure/orthobench_3di_best.fa", gp)
    pools = {}
    for r in split["dev"]:
        g = [l.strip() for l in open(f"{T.OB}/BENCHMARKS/RefOGs/{r}.txt") if l.strip()]
        g = sorted(x for x in g if x in ss and x in idx)
        if len(g) >= 5:
            pools[r] = g
    rng = np.random.default_rng(11)
    keys = sorted(pools)
    tasks = []
    for i in range(a.n):
        r = keys[rng.integers(len(keys))]
        gs = sorted(rng.choice(pools[r], 5, replace=False).tolist())
        S = np.full((5, 5), 0.5)
        np.fill_diagonal(S, 1.0)
        tasks.append(dict(key=f"c{i:03d}_{r}", names=gs, seqs=[gp.seqs[idx[x]] for x in gs],
                          ss=[ss[x] for x in gs], S=S, arms=T.ARMS_ALL))
    with Pool(a.workers) as pool:
        res = dict(pool.imap_unordered(T.do_task, tasks))
    arms = [x for x in T.ARMS_ALL]
    out = {}
    for arm in arms:
        c = np.array([res[k][arm]["cpu"] for k in res])
        w = np.array([res[k][arm]["wall"] for k in res])
        out[arm] = dict(cpu_mean_s=float(c.mean()), cpu_median_s=float(np.median(c)), wall_mean_s=float(w.mean()),
                        projected_cpu_h_orthobench=float(c.mean() * N_CLUSTERS / 3600),
                        acceptable=bool(c.mean() * N_CLUSTERS / 3600 <= 1.0))
    out["_note"] = ("a_nj CPU here is the NJ only (scores are precomputed by the pipeline); "
                    "MSA arms include FAMSA; e_* include Foldseek DB build + all-pair alignment")
    json.dump(out, open(f"{ROOT}/results/S4_trees/cost_5units.json", "w"), indent=1)
    for k, v in out.items():
        print(k, v)


if __name__ == "__main__":
    main()
