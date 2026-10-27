#!/usr/bin/env python
"""L2-gpu task 1: apply a `gpu_sparsify` rule post hoc to an all-GPU score table (the same code the
pipeline option runs: ancestra.scoring.sparsify_table), or compute density-matched parameters.

    # density matching (label-free): k / alpha / x whose kept count in families > 1000 genes
    # first reaches DIAMOND's (CPU pass-1) count in the same families
    python scripts/l2gpu_sparsify_posthoc.py match --gpu GPU.npz --cpu CPU.npz --families FAM.npy \
        --genes genes.tsv --out params.json
    # apply
    python scripts/l2gpu_sparsify_posthoc.py apply --gpu GPU.npz --families FAM.npy --genes genes.tsv \
        --spec topk:40 --out scores.npz
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.scoring import ScoreTable, ScoringConfig, sparsify_table  # noqa: E402


class GT:
    def __init__(self, path):
        import pandas as pd
        d = pd.read_csv(path, sep="\t")
        _, self.species_of = np.unique(d["species"].to_numpy(), return_inverse=True)
        self.lengths = d["length"].to_numpy().astype(np.int64)
        self.n = len(d)


def per_pair_stats(st, fam, gt, minf):
    sizes = np.bincount(fam)
    q, t = st.q.astype(np.int64), st.t.astype(np.int64)
    idx = np.flatnonzero(sizes[fam[q]] > minf)
    q, t, b = q[idx], t[idx], st.bits[idx].astype(np.float64)
    m = len(q)
    dq, dt, db = np.concatenate([q, t]), np.concatenate([t, q]), np.concatenate([b, b])
    tsp = gt.species_of[dt].astype(np.int64)
    order = np.lexsort((dt, -db, tsp, dq))
    k1, k2 = dq[order], tsp[order]
    start = np.r_[True, (k1[1:] != k1[:-1]) | (k2[1:] != k2[:-1])]
    gid = np.cumsum(start) - 1
    first = np.flatnonzero(start)
    rank = np.empty(2 * m, np.int64)
    ratio = np.empty(2 * m)
    rank[order] = np.arange(len(order)) - first[gid]
    ratio[order] = db[order] / db[order][first][gid]
    selfb = st.self_bits.astype(np.float64).copy()
    miss = ~np.isfinite(selfb)
    selfb[miss] = 2.2 * gt.lengths[miss]
    norm = b / np.maximum(selfb[q], selfb[t])
    return np.minimum(rank[:m], rank[m:]), np.maximum(ratio[:m], ratio[m:]), norm, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["match", "apply"])
    ap.add_argument("--gpu", required=True)
    ap.add_argument("--cpu")
    ap.add_argument("--families", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--spec")
    ap.add_argument("--min-family", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    fam = np.load(a.families)
    gt = GT(a.genes)
    st = ScoreTable.load(a.gpu)
    if a.mode == "match":
        sizes = np.bincount(fam)
        zc = np.load(a.cpu)
        target = int((sizes[fam[zc["q"]]] > a.min_family).sum())
        rank, ratio, norm, m = per_pair_stats(st, fam, gt, a.min_family)
        cnt = np.bincount(np.minimum(rank, 10 ** 6))
        cum = np.cumsum(cnt)                 # cum[k-1] = kept with rank < k
        k = int(np.searchsorted(cum, target) + 1)
        alpha = float(np.sort(ratio)[::-1][target - 1])
        x = float(np.sort(norm)[::-1][target - 1])
        res = dict(target_pairs=target, gpu_pairs_big=m, k=k, kept_k=int(cum[k - 1]),
                   kept_k10=int(cum[9]), alpha=alpha, kept_alpha=int((ratio >= alpha).sum()),
                   x=x, kept_x=int((norm >= x).sum()), rand_fraction=float(cum[k - 1] / m),
                   runtime_s=round(time.time() - t0, 1))
        json.dump(res, open(a.out, "w"), indent=1)
        print(json.dumps(res, indent=1))
        return
    cfg = ScoringConfig(backend="gpusw", gpu_max_family=0, gpu_sparsify=a.spec, rescue_max_family=a.min_family)
    stats: dict = {}
    out = sparsify_table(st, fam, gt, cfg, stats)
    out.save(a.out)
    stats.update(spec=a.spec, n_in=len(st), n_out=len(out), runtime_s=round(time.time() - t0, 1))
    json.dump(stats, open(a.out + ".json", "w"), indent=1)
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
