#!/usr/bin/env python
"""L2-plm step (a): does a pLM score rank a gene's ortholog above its same-species paralog?

Pre-registered in docs/L2_plm.md §0. Triplets (a, o, p) from the DEV pair labels; accuracy of S(a,o) > S(a,p)
(ties 0.5), groups all / deep / shallow / hard (S5 rule on frozen bits), sets T_seq (both sequence scores exist)
and T_all (missing sequence score = 0). Arms: seq, *_mean (cosine of mean-pooled embeddings), *_aln (embedding
alignment, scripts/l2_plm_alnscore.py), comb_*@0.9 (sequence unless near-tie, then pLM), head_* (idea 2).
Control: embeddings permuted among the RefOG genes of the same species (seeds 11/23/47).
Paired bootstrap over RefOGs, B = 2000, seed 11. Also pooled pair ROC-AUC (RefOG bootstrap CI).

    python scripts/l2_plm_discrim.py --arms e35_mean e650_mean pt5_mean e650_aln pt5_aln --tag main
Runs in conda_envs/orth; the alignment arms call the GPU scorer in conda_envs/gpu.
"""
import argparse, json, os, subprocess, sys, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2_plm_common import DATA, RUN, WORK, RES, SEEDS, DEEP, gene_table, manifest

GPU_PY = os.path.join(DATA, "conda_envs/gpu/bin/python")
HERE_S = os.path.dirname(os.path.abspath(__file__))
B = 2000


def seq_tables(gt):
    from ancestra.scoring import ScoreTable
    st = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
    sb = st.self_bits_filled(gt.lengths)
    q, t = st.q.astype(np.int64), st.t.astype(np.int64)
    n = gt.n
    a, b = np.minimum(q, t), np.maximum(q, t)
    key = a * n + b
    norm = st.bits / np.maximum(sb[q], sb[t])
    sn, sbits = {}, {}
    for k, v, w in zip(key.tolist(), norm.tolist(), st.bits.tolist()):
        if v > sn.get(k, -1):
            sn[k] = v
        if w > sbits.get(k, -1):
            sbits[k] = w
    return sn, sbits


def build(gt):
    from orthobench_dev_pairs_eval import load_truth
    truth = load_truth()
    idx = {x: i for i, x in enumerate(gt.names)}
    ga = np.array([idx[x] for x in truth.gene_a]); gb = np.array([idx[x] for x in truth.gene_b])
    return truth, ga, gb


def permutation(genes, species_of, seed):
    rng = np.random.default_rng(seed)
    pi = {}
    by = defaultdict(list)
    for g in genes:
        by[species_of[g]].append(g)
    for s in sorted(by):
        gs = np.array(sorted(by[s]))
        pg = rng.permutation(gs)
        pi.update(zip(gs.tolist(), pg.tolist()))
    return pi


class MeanArm:
    def __init__(self, model):
        if model == "e35":
            E = np.load(os.path.join(DATA, "work/gpuv2/esm_knn/emb_esm2_t12_35M_UR50D.npy"))
            self.row = lambda g: g
        else:
            E = np.load(os.path.join(WORK, f"emb_{model}_refog.mean.npy"))
            gids = np.load(os.path.join(WORK, f"emb_{model}_refog.gids.npy"))
            m = {int(g): i for i, g in enumerate(gids)}
            self.row = lambda g: m[g]
        self.E = E / np.linalg.norm(E, axis=1, keepdims=True)

    def score(self, pairs):
        r1 = np.array([self.row(a) for a, _ in pairs]); r2 = np.array([self.row(b) for _, b in pairs])
        return np.einsum("ij,ij->i", self.E[r1], self.E[r2])


class AlnArm:
    def __init__(self, model):
        self.model = model
        gids = np.load(os.path.join(WORK, f"emb_{model}_refog.gids.npy"))
        self.m = {int(g): i for i, g in enumerate(gids)}
        self.cache_f = os.path.join(WORK, f"aln_{model}_cache.npz")
        self.cache = {}
        if os.path.exists(self.cache_f):
            z = np.load(self.cache_f)
            self.cache = dict(zip(map(tuple, z["p"].tolist()), z["s"].tolist()))

    def ensure(self, pairs):
        need = sorted({(self.m[a], self.m[b]) for a, b in pairs} | {(self.m[a], self.m[a]) for a, b in pairs}
                      | {(self.m[b], self.m[b]) for a, b in pairs})
        need = [p for p in need if p not in self.cache]
        if not need:
            return
        pf = os.path.join(WORK, f"aln_{self.model}_need.npy"); of = pf.replace("_need", "_out")
        np.save(pf, np.array(need, np.int64))
        r = subprocess.run([GPU_PY, os.path.join(HERE_S, "l2_plm_alnscore.py"), "--emb",
                            os.path.join(WORK, f"emb_{self.model}_refog"), "--pairs", pf, "--out", of],
                           capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(r.stderr[-3000:])
        print(self.model, r.stdout.strip().splitlines()[-1], flush=True)
        self.cache.update(zip(need, np.load(of).tolist()))
        P = np.array(list(self.cache.keys()), np.int64)
        np.savez(self.cache_f, p=P, s=np.array([self.cache[tuple(x)] for x in P.tolist()]))

    def score(self, pairs):
        self.ensure(pairs)
        out = []
        for a, b in pairs:
            x, y = self.m[a], self.m[b]
            out.append(self.cache[(x, y)] / max(min(self.cache[(x, x)], self.cache[(y, y)]), 1e-9))
        return np.array(out)


def roc_auc(s, y):
    from scipy.stats import rankdata
    r = rankdata(s)
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["e35_mean", "e650_mean", "pt5_mean"])
    ap.add_argument("--heads", nargs="*", default=[], help="name=path.npz score files from l2_plm_head.py")
    ap.add_argument("--tag", default="main")
    ap.add_argument("--contrasts", nargs="*", default=[], help="A:B -> paired bootstrap of acc(A real) - acc(B real)")
    a = ap.parse_args()
    t0 = time.time()
    gt = gene_table(); n = gt.n
    sp = np.asarray(gt.species_of)
    sn, sbits = seq_tables(gt)
    truth, ga, gb = build(gt)
    lab = truth.label.values; cls = truth.cls.values; ref = truth.refog.values
    genes = sorted(set(ga.tolist()) | set(gb.tolist()))
    key = lambda x, y: min(x, y) * n + max(x, y)
    part = defaultdict(list)
    for x, y, l, c, r in zip(ga, gb, lab, cls, ref):
        part[(x, sp[y])].append((y, l, c, r)); part[(y, sp[x])].append((x, l, c, r))
    trip = []   # (refog, deep, hard, both_seq, a, o, p)
    for (x, s), lst in part.items():
        O = [(y, c, r) for y, l, c, r in lst if l == 1]
        P = [y for y, l, c, r in lst if l == 0]
        if not O or not P:
            continue
        bo = max(sbits.get(key(x, y), 0.0) for y, _, _ in O)
        bp = max(sbits.get(key(x, y), 0.0) for y in P)
        hard = bp > 0 and bp >= 0.9 * bo
        for o, c, r in O:
            for p in P:
                both = key(x, o) in sn and key(x, p) in sn
                trip.append((r, c in DEEP, hard, both, x, o, p))
    T = np.array([(t[1], t[2], t[3]) for t in trip], bool)
    tri = np.array([(t[4], t[5], t[6]) for t in trip], np.int64)
    refs = sorted(set(ref)); ri = {r: i for i, r in enumerate(refs)}
    tref = np.array([ri[t[0]] for t in trip])
    print(f"triplets {len(trip)}  T_seq {T[:, 2].sum()}  hard {T[:, 1].sum()}  deep {T[:, 0].sum()}", flush=True)

    # scores per arm per variant on (a,o) and (a,p)
    seqv = lambda P: np.array([sn.get(key(x, y), 0.0) for x, y in P])
    pao = list(map(tuple, tri[:, [0, 1]].tolist())); pap = list(map(tuple, tri[:, [0, 2]].tolist()))
    lp = list(zip(ga.tolist(), gb.tolist()))
    S = {"seq": {"real": (seqv(pao), seqv(pap), seqv(lp))}}
    perms = {f"shuf{s}": permutation(genes, sp, s) for s in SEEDS}
    arms = {}
    for arm in a.arms:
        model, kind = arm.split("_")
        arms[arm] = MeanArm(model) if kind == "mean" else AlnArm(model)
    for nm in a.heads:
        name, path = nm.split("=", 1)
        arms[name] = HeadArm(path)
    for arm, obj in arms.items():
        if hasattr(obj, "ensure"):   # one GPU job for all variants
            allp = []
            for v in ["real"] + list(perms):
                mp = (lambda g: g) if v == "real" else (lambda g, pi=perms[v]: pi[g])
                allp += [(mp(x), mp(y)) for x, y in pao + pap + lp]
            obj.ensure(allp)
        S[arm] = {}
        for v in ["real"] + list(perms):
            mp = (lambda g: g) if v == "real" else (lambda g, pi=perms[v]: pi[g])
            m2 = lambda P: [(mp(x), mp(y)) for x, y in P]
            S[arm][v] = (obj.score(m2(pao)), obj.score(m2(pap)), obj.score(m2(lp)))
        print(arm, "scored", round(time.time() - t0), flush=True)
    cmp = lambda so, sp_: np.where(so > sp_, 1.0, np.where(so == sp_, 0.5, 0.0))
    so_seq, sp_seq, _ = S["seq"]["real"]
    near = np.minimum(so_seq, sp_seq) >= 0.9 * np.maximum(so_seq, sp_seq)
    cols = {"seq": cmp(so_seq, sp_seq)}
    for arm in arms:
        for v in S[arm]:
            so, spp, _ = S[arm][v]
            cols[f"{arm}|{v}"] = cmp(so, spp)
            cols[f"comb_{arm}|{v}"] = np.where(near, cmp(so, spp), cols["seq"])
    np.savez_compressed(os.path.join(WORK, f"triplets_{a.tag}.npz"), tri=tri, tref=tref, T=T, near=near,
                        names=np.array(list(cols)), cols=np.array([cols[k] for k in cols], np.float32),
                        so_seq=so_seq, sp_seq=sp_seq,
                        **{f"so|{arm}": S[arm]["real"][0] for arm in arms}, **{f"sp|{arm}": S[arm]["real"][1] for arm in arms})
    rng = np.random.default_rng(11)
    W = rng.multinomial(len(refs), np.full(len(refs), 1 / len(refs)), size=B).astype(float)
    ci = lambda d: [round(float(np.percentile(d, 2.5)), 4), round(float(np.percentile(d, 97.5)), 4)]
    out = dict(n_triplets=len(trip), groups={}, pair_auc={}, near_tie_share={})
    groups = {"all": np.ones(len(trip), bool), "deep": T[:, 0], "shallow": ~T[:, 0], "hard": T[:, 1],
              "hard_deep": T[:, 1] & T[:, 0]}
    for setname, smask in (("T_seq", T[:, 2]), ("T_all", np.ones(len(trip), bool))):
        for g, gm in groups.items():
            m = smask & gm
            if m.sum() == 0:
                continue
            cnt = np.bincount(tref[m], minlength=len(refs)).astype(float)
            acc, boot = {}, {}
            for k, c in cols.items():
                sums = np.bincount(tref[m], weights=c[m], minlength=len(refs))
                acc[k] = float(c[m].mean()); boot[k] = (W @ sums) / np.maximum(W @ cnt, 1e-12)
            res = dict(n=int(m.sum()), n_refogs=int((cnt > 0).sum()), near_tie_share=round(float(near[m].mean()), 4),
                       acc={"seq": round(acc["seq"], 4)}, arms={})
            for arm in arms:
                for pre in ("", "comb_"):
                    k = pre + arm
                    sh = [f"{k}|shuf{s}" for s in SEEDS]
                    bsh = sum(boot[x] for x in sh) / 3
                    res["arms"][k] = dict(
                        acc=round(acc[f"{k}|real"], 4), shuf=[round(acc[x], 4) for x in sh],
                        minus_seq=dict(delta=round(acc[f"{k}|real"] - acc["seq"], 4),
                                       ci95=ci(boot[f"{k}|real"] - boot["seq"])),
                        minus_shuf=dict(delta=round(acc[f"{k}|real"] - float(np.mean([acc[x] for x in sh])), 4),
                                        ci95=ci(boot[f"{k}|real"] - bsh)))
            res["contrasts"] = {}
            for cs in a.contrasts:
                x_, y_ = cs.split(":")
                kx = "seq" if x_ == "seq" else f"{x_}|real"; ky = "seq" if y_ == "seq" else f"{y_}|real"
                res["contrasts"][cs] = dict(delta=round(acc[kx] - acc[ky], 4), ci95=ci(boot[kx] - boot[ky]))
            out["groups"][f"{setname}/{g}"] = res
    # pooled pair ROC-AUC (all labelled within-RefOG pairs), RefOG bootstrap
    pref = np.array([ri[r] for r in ref]); y = lab.astype(int)
    pdeep = np.isin(cls, list(DEEP))
    rng2 = np.random.default_rng(11)
    bidx = [rng2.choice(len(refs), len(refs)) for _ in range(200)]
    rows_by = [np.flatnonzero(pref == i) for i in range(len(refs))]
    for g, gm in (("all", np.ones(len(y), bool)), ("deep", pdeep), ("shallow", ~pdeep)):
        for arm in ["seq"] + list(arms):
            for v in (["real"] if arm == "seq" else ["real", "shuf11"]):
                s = S[arm][v][2]
                auc = roc_auc(s[gm], y[gm])
                bs = []
                for bi in bidx:
                    rr = np.concatenate([rows_by[i] for i in bi]); rr = rr[gm[rr]]
                    bs.append(roc_auc(s[rr], y[rr]))
                out["pair_auc"][f"{g}/{arm}|{v}"] = dict(auc=round(float(auc), 4), ci95=ci(np.array(bs)))
    out["manifest"] = manifest(dict(arms=a.arms, heads=a.heads, seeds=SEEDS, B=B, near_tie=0.9, hard_rule="S5 0.9 bits"),
                               dict(scores=os.path.join(RUN, "stages/scores/scores.npz"), labels="work/S5/refog_pair_labels.parquet"),
                               t0)
    os.makedirs(RES, exist_ok=True)
    f = os.path.join(RES, f"discrimination_{a.tag}.json")
    json.dump(out, open(f, "w"), indent=1)
    for k, r in out["groups"].items():
        line = f"{k:18s} n={r['n']:6d} seq={r['acc']['seq']:.3f} " + " ".join(
            f"{arm}={v['acc']:.3f}[{v['minus_seq']['ci95'][0]:+.3f},{v['minus_seq']['ci95'][1]:+.3f}]/sh{np.mean(v['shuf']):.3f}"
            for arm, v in r["arms"].items())
        print(line)
        for cs, v in r.get("contrasts", {}).items():
            print(f"   {cs}: {v['delta']:+.4f} {v['ci95']}")
    print(f)


class HeadArm:
    """Scores from a trained head: an npz with gids and projected unit vectors (z)."""
    def __init__(self, path):
        z = np.load(path)
        self.m = {int(g): i for i, g in enumerate(z["gids"])}
        self.Z = z["z"]

    def score(self, pairs):
        r1 = np.array([self.m[a] for a, _ in pairs]); r2 = np.array([self.m[b] for _, b in pairs])
        return np.einsum("ij,ij->i", self.Z[r1], self.Z[r2])


if __name__ == "__main__":
    main()
