"""G2: cross-piece rep-vs-rep edges for Leiden-split S0 components (cap unchanged, 7.5k).

No new search: S0's stage-2 search (every gene vs all stage-1 reps, DIAMOND --sensitive, e <= 1e-5, -k 500;
stages/families/s0/reps_hits.tsv) already holds the rep-vs-rep hits. Kept: query and subject are both
stage-1 reps of the same split parent component but in different S2-final families (= different pieces),
e <= 1e-5, aligned span >= 0.3 of the shorter protein (S0's own edge filter); then per query rep the top 2
hits per other species and 20 same-species hits (S2 cap rule). The pairs are appended to the S2 ScoreTable
(bits, coverages, e-value from the S0 search; pident 0), families unchanged. The merge sees them as ordinary
edges, so sibling-piece units can become mutual best at any level (bounded: <= 22 edges per rep).

    python scripts/giant_g2_edges.py {qfo,ob}
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import load_proteomes, read_gene_table  # noqa
from ancestra.scoring import ScoreTable, ScoringConfig, apply_cap  # noqa

CFG = {"qfo": dict(run="work/rescue/qfo/base", prot="data/inputs/qfo2022/proteomes", out="work/giant/qfo"),
       "ob": dict(run="work/core_v0_ob/run_det", prot="data/orthobench/Open_Orthobench/BENCHMARKS/Input",
                  out="work/giant/ob")}


def main():
    bench = sys.argv[1]
    C = CFG[bench]
    t0 = time.time()
    gt = read_gene_table(os.path.join(C["run"], "stages/input/genes.tsv"))
    full = load_proteomes(C["prot"])
    seqs = full.seqs
    order = np.array(sorted(range(gt.n), key=lambda g: (-len(seqs[g]), seqs[g], gt.names[g])), np.int64)
    s0 = os.path.join(C["run"], "stages/families/s0")
    comp = np.load(os.path.join(C["out"], "s0_g7500_component.npy"))
    fam = np.load(os.path.join(C["run"], "stages/scores/families_final.npy")).astype(np.int64)
    rep = np.zeros(gt.n, bool)
    for line in open(os.path.join(s0, "stage1_clusters.tsv")):
        r, m = line.split("\t")[:2]
        if r == m.strip():
            rep[order[int(r[1:])]] = True
    df = pd.read_csv(os.path.join(s0, "reps_hits.tsv"), sep="\t", header=None,
                     names=["q", "s", "ev", "bits", "ql", "sl", "qs", "qe", "ss", "se"],
                     dtype={"q": str, "s": str})
    q = order[df["q"].str[1:].astype(np.int64).to_numpy()]
    t = order[df["s"].str[1:].astype(np.int64).to_numpy()]
    ev, bits = df["ev"].to_numpy(np.float64), df["bits"].to_numpy(np.float64)
    ql, sl = df["ql"].to_numpy(np.float64), df["sl"].to_numpy(np.float64)
    qcov = (np.abs(df["qe"].to_numpy() - df["qs"].to_numpy()) + 1) / ql
    tcov = (np.abs(df["se"].to_numpy() - df["ss"].to_numpy()) + 1) / sl
    del df
    span = np.maximum(qcov * ql, tcov * sl) / np.minimum(ql, sl)
    m = (rep[q] & rep[t] & (comp[q] >= 0) & (comp[q] == comp[t]) & (fam[q] != fam[t]) & (ev <= 1e-5)
         & (span >= 0.3))
    n_cand = int(m.sum())
    q, t, bits, qcov, tcov, ev = q[m], t[m], bits[m], qcov[m], tcov[m], ev[m]
    # orient q<t, dedup (max bits)
    sw = q > t
    q, t = np.where(sw, t, q), np.where(sw, q, t)
    qcov, tcov = np.where(sw, tcov, qcov), np.where(sw, qcov, tcov)
    o = np.lexsort((-bits, t, q))
    q, t, bits, qcov, tcov, ev = (x[o] for x in (q, t, bits, qcov, tcov, ev))
    k = np.r_[True, (q[1:] != q[:-1]) | (t[1:] != t[:-1])]
    q, t, bits, qcov, tcov, ev = (x[k] for x in (q, t, bits, qcov, tcov, ev))
    st = ScoreTable.load(os.path.join(C["run"], "stages/scores/scores.npz"))
    edge = ScoreTable(q=q.astype(np.int32), t=t.astype(np.int32), bits=bits.astype(np.float32),
                      qcov=qcov.astype(np.float32), tcov=tcov.astype(np.float32),
                      pident=np.zeros(len(q), np.float32), evalue=ev.astype(np.float32), self_bits=st.self_bits)
    # cap: all edges treated as one big "family" so the per-query top-2/species + 20 same-species rule applies
    cfg = ScoringConfig(cap_all=0)
    s2 = {}
    edge = apply_cap(edge, np.zeros(gt.n, np.int64), gt.species_of, cfg, s2)
    out = ScoreTable(*(np.concatenate([getattr(st, x), getattr(edge, x)]) for x in
                       ("q", "t", "bits", "qcov", "tcov", "pident", "evalue")), st.self_bits)
    out.save(os.path.join(C["out"], "scores_G2.npz"))
    stats = dict(n_reps_in_split_components=int((rep & (comp >= 0)).sum()), n_cross_piece_candidates=n_cand,
                 n_unique_pairs=int(len(q)), n_edges_after_cap=int(len(edge)), n_pairs_total=int(len(out)),
                 n_pairs_old=int(len(st)), wall_s=round(time.time() - t0, 1))
    json.dump(stats, open(os.path.join(C["out"], "g2_stats.json"), "w"), indent=1)
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
