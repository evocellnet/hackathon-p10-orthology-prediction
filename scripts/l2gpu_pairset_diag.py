#!/usr/bin/env python
"""L2-gpu diagnostic (label-free): how does DIAMOND's sparse in-family pair set relate to the
exhaustive GPU SW pair set in large families, and which *rank-based / self-normalised* rule
reproduces it best?  No benchmark labels are read.

    python scripts/l2gpu_pairset_diag.py CPU_SCORES.npz GPU_SCORES.npz FAMILIES.npy SPECIES_OF.npy OUT.json

For pairs in families with > MINF genes (default 1000):
- pair-set overlap |C|, |G|, |C & G|;
- Jaccard(C, rule(G)) for rules: top-k per (query, target species) either direction; ratio
  bits >= a * best(query, target species) either direction; norm_max >= x; bits >= b.
"""
import json
import sys

import numpy as np


def directed(q, t, b):
    return np.concatenate([q, t]), np.concatenate([t, q]), np.concatenate([b, b])


def rank_and_ratio(q, t, b, sp):
    """Per undirected pair: min over the two directions of the rank of the partner among the
    query's hits to the partner's species (0 = best), and max of bits / best(query, species)."""
    m = len(q)
    dq, dt, db = directed(q, t, b)
    tsp = sp[dt].astype(np.int64)
    order = np.lexsort((-db, tsp, dq))
    k1, k2 = dq[order], tsp[order]
    start = np.r_[True, (k1[1:] != k1[:-1]) | (k2[1:] != k2[:-1])]
    gid = np.cumsum(start) - 1
    first = np.flatnonzero(start)
    rank_s = np.arange(len(order)) - first[gid]
    best_s = db[order][first][gid]
    rank = np.empty(2 * m, np.int64)
    ratio = np.empty(2 * m, np.float64)
    rank[order] = rank_s
    ratio[order] = db[order] / best_s
    return np.minimum(rank[:m], rank[m:]), np.maximum(ratio[:m], ratio[m:])


def main():
    cpu_p, gpu_p, fam_p, sp_p, out = sys.argv[1:6]
    minf = 1000
    fam = np.load(fam_p)
    sp = np.load(sp_p)
    sizes = np.bincount(fam)
    res = {"min_family": minf}
    tabs = {}
    for name, p in (("cpu", cpu_p), ("gpu", gpu_p)):
        z = np.load(p)
        q, t = z["q"].astype(np.int64), z["t"].astype(np.int64)
        big = (sizes[fam[q]] > minf) & (fam[q] == fam[t])
        tabs[name] = dict(q=q[big], t=t[big], bits=z["bits"][big].astype(np.float64),
                          self=z["self_bits"].astype(np.float64))
        res[f"n_{name}"] = int(big.sum())
    N = len(fam)
    kc = tabs["cpu"]["q"] * N + tabs["cpu"]["t"]
    G = tabs["gpu"]
    kg = G["q"] * N + G["t"]
    inC = np.isin(kg, kc)
    res["n_both"] = int(inC.sum())
    res["jaccard_all"] = res["n_both"] / (res["n_cpu"] + res["n_gpu"] - res["n_both"])
    rank, ratio = rank_and_ratio(G["q"], G["t"], G["bits"], sp)
    sb = G["self"]
    norm = G["bits"] / np.maximum(sb[G["q"]], sb[G["t"]])
    nc = res["n_cpu"]

    def jac(mask):
        a = int((mask & inC).sum())
        return dict(kept=int(mask.sum()), in_cpu=a, jaccard=round(a / (nc + int(mask.sum()) - a), 4),
                    recall_of_cpu=round(a / nc, 4), precision=round(a / max(int(mask.sum()), 1), 4))
    res["topk_per_species"] = {k: jac(rank < k) for k in (1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 100)}
    res["ratio_to_best"] = {a: jac(ratio >= a) for a in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)}
    res["norm_max"] = {x: jac(norm >= x) for x in (0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15, 0.2)}
    res["bits"] = {b: jac(G["bits"] >= b) for b in (30, 35, 40, 45, 50, 60, 80)}
    # inclusion probability of a GPU pair in DIAMOND's set by rank bin and by norm bin
    rb = np.minimum(rank, 200)
    edges = [0, 1, 2, 3, 5, 10, 20, 50, 100, 201]
    res["p_in_cpu_by_rank"] = {f"{edges[i]}-{edges[i+1]-1}": round(float(inC[(rb >= edges[i]) & (rb < edges[i + 1])].mean()), 4)
                               for i in range(len(edges) - 1) if ((rb >= edges[i]) & (rb < edges[i + 1])).any()}
    nb = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1.01]
    res["p_in_cpu_by_norm"] = {f"{nb[i]}-{nb[i+1]}": round(float(inC[(norm >= nb[i]) & (norm < nb[i + 1])].mean()), 4)
                               for i in range(len(nb) - 1) if ((norm >= nb[i]) & (norm < nb[i + 1])).any()}
    # self-score agreement (CPU DIAMOND self hit vs GPU SW) in big families
    gb = np.flatnonzero(sizes[fam] > minf)
    a, b = tabs["cpu"]["self"][gb], tabs["gpu"]["self"][gb]
    ok = np.isfinite(a) & np.isfinite(b)
    res["self_ratio_gpu_over_cpu_median"] = round(float(np.median(b[ok] / a[ok])), 4)
    res["self_cpu_missing"] = int((~np.isfinite(a)).sum())
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
