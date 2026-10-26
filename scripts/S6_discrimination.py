"""S6 diagnostic: does the 3Di+AA score rank orthologs above paralogs better than sequence?
Triplets (a, o, p) from DEV RefOG pair labels: o ortholog and p paralog of a, both in the same
species; correct if S(a,o) > S(a,p) (ties = 0.5). Only triplets where every score compared is
available (both pairs aligned by sequence and by Foldseek). Scores: S_seq (run_v2 norm_max),
S_3di (Foldseek 3Di+AA, bits/max self) for real and shuffled 3Di, and the combined avg/max S
used by the merge arms. Split by the species-pair class of (a, o). Paired bootstrap over
RefOGs (B=2000, seed 11) for accuracy(real 3Di) - accuracy(seq) and - accuracy(shuffled).
    python scripts/S6_discrimination.py -> results/S6_structure/discrimination.json"""
import json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import read_gene_table
from ancestra.scoring import ScoreTable
from orthobench_dev_pairs_eval import load_truth

RUN = os.path.join(ROOT, "work/core_v0_ob/run_v2")
DEEP = {"vert-teleost", "vert-invert", "invert-invert"}
FIT = (0.733, 0.881)   # d_seq ~ a + b d_3di fitted by combine_scores on the real-3Di pairs (merge_runs.jsonl)


def pairdict(q, t, s, n):
    a, b = np.minimum(q, t).astype(np.int64), np.maximum(q, t).astype(np.int64)
    return dict(zip((a * n + b).tolist(), s.tolist()))


def main():
    gt = read_gene_table(os.path.join(RUN, "stages/input/genes.tsv"))
    n = gt.n
    idx = {x: i for i, x in enumerate(gt.names)}
    st = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
    sb = st.self_bits_filled(gt.lengths)
    q, t = st.q.astype(np.int64), st.t.astype(np.int64)
    seq = pairdict(q, t, st.bits / np.maximum(sb[q], sb[t]), n)
    s3 = {}
    for v in ["real", "shuf11", "shuf23", "shuf47"]:
        z = np.load(os.path.join(ROOT, "work/S6/pairs/seq", v, "fs.npz"))
        fq, ft, fb = z["q"], z["t"], z["bits"]
        self_ = np.full(n, np.nan); m = fq == ft; self_[fq[m]] = fb[m]
        fq, ft, fb = fq[~m], ft[~m], fb[~m]
        s3[v] = pairdict(fq, ft, fb / np.maximum(self_[fq], self_[ft]), n)
    truth = load_truth()
    ga = np.array([idx[x] for x in truth.gene_a]); gb = np.array([idx[x] for x in truth.gene_b])
    lab = truth.label.values; cls = truth.cls.values; ref = truth.refog.values
    # per gene a, species -> lists of (partner, label)
    from collections import defaultdict
    part = defaultdict(list)
    for a, b, l, c, r in zip(ga, gb, lab, cls, ref):
        part[(a, gt.species_of[b])].append((b, l, c, r))
        part[(b, gt.species_of[a])].append((a, l, c, r))
    key = lambda a, b: min(a, b) * n + max(a, b)
    trip = []   # (refog, deep, seq_cmp, {v: cmp})
    for (a, s), lst in part.items():
        O = [(b, c, r) for b, l, c, r in lst if l == 1]
        P = [b for b, l, c, r in lst if l == 0]
        for o, c, r in O:
            for p in P:
                ko, kp = key(a, o), key(a, p)
                if ko not in seq or kp not in seq or any(ko not in s3[v] or kp not in s3[v] for v in s3):
                    continue
                cmp = lambda d: 1.0 if d[ko] > d[kp] else (0.5 if d[ko] == d[kp] else 0.0)
                comb = {}
                for k_ in (ko, kp):   # merge-arm combinations with the real-3Di calibration (fit a, b)
                    c3 = np.exp(-(FIT[0] + FIT[1] * -np.log(min(max(s3["real"][k_], 1e-6), 1.0))))
                    comb.setdefault("avg_real", {})[k_] = np.sqrt(max(seq[k_], 1e-6) * c3)
                    comb.setdefault("max_real", {})[k_] = max(seq[k_], c3)
                d_ = {v: cmp(s3[v]) for v in s3}
                d_.update({k_: cmp(dd) for k_, dd in comb.items()})
                trip.append((r, c in DEEP, cmp(seq), d_))
    refs = sorted({x[0] for x in trip})
    ri = {r: i for i, r in enumerate(refs)}
    out = dict(n_triplets=len(trip))
    rng = np.random.default_rng(11)
    W = rng.multinomial(len(refs), np.full(len(refs), 1 / len(refs)), size=2000).astype(float)
    for grp, sel in (("all", lambda x: True), ("deep", lambda x: x[1]), ("shallow", lambda x: not x[1])):
        T = [x for x in trip if sel(x)]
        if not T:
            continue
        cols = {"seq": np.array([x[2] for x in T])}
        for v in list(s3) + ["avg_real", "max_real"]:
            cols[v] = np.array([x[3][v] for x in T])
        r_ = np.array([ri[x[0]] for x in T])
        cnt = np.bincount(r_, minlength=len(refs)).astype(float)
        sums = {k: np.bincount(r_, weights=c, minlength=len(refs)) for k, c in cols.items()}
        acc = {k: float(c.mean()) for k, c in cols.items()}
        boot = {k: (W @ sums[k]) / (W @ cnt) for k in cols}
        sh = (boot["shuf11"] + boot["shuf23"] + boot["shuf47"]) / 3
        ci = lambda d: [round(float(np.percentile(d, 2.5)), 4), round(float(np.percentile(d, 97.5)), 4)]
        out[grp] = dict(n=len(T), n_refogs=int((cnt > 0).sum()), acc={k: round(v, 4) for k, v in acc.items()},
                        real_minus_seq=dict(delta=round(acc["real"] - acc["seq"], 4), ci95=ci(boot["real"] - boot["seq"])),
                        real_minus_shuf=dict(delta=round(acc["real"] - np.mean([acc[s] for s in ("shuf11", "shuf23", "shuf47")]), 4),
                                             ci95=ci(boot["real"] - sh)))
    json.dump(out, open(os.path.join(ROOT, "results/S6_structure/discrimination.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
