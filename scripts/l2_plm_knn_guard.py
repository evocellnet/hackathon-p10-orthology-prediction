#!/usr/bin/env python
"""L2-plm idea 4: pLM kNN with a family guard as an S0 completeness rescue (DEV only; standalone S0 metrics).

Pre-registered guard (docs/L2_plm.md §0): a cross-family kNN link (k = 25, ESM-2 35M mean embeddings of all
OrthoBench proteins, work/gpuv2/esm_knn) is accepted only if it is MUTUAL and its cosine is above the 99th
percentile of a null; a family join needs >= 2 accepted links between the two families (union-find).
Null (our operationalisation, fixed before the run): for every mutual cross-family link (i, j), the cosine of
(i, j') with j' a random gene of j's species (seed 11) -- "shuffled within species".
Control: embeddings permuted among the genes of each species (seeds 11/23/47), same procedure.
Metrics on the 34 DEV RefOGs vs the frozen S0 partition (work/ob_final/main_accurate_go): completeness (all RefOG
genes in one family), DEV RefOG pairs newly put in one family (fusions; RefOG010 flagged), families joined,
largest family.     python scripts/l2_plm_knn_guard.py      (conda_envs/gpu)
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2_plm_common import DATA, RES, SEEDS, manifest

K = 25


def uf_join(fam, pairs):
    parent = {}
    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x]); x = parent[x]
        return x
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    return np.array([find(f) for f in fam.tolist()])


def run(E, fam, sp, rng_seed=11):
    import faiss
    res = faiss.StandardGpuResources()
    idx = faiss.GpuIndexFlatIP(res, E.shape[1]); idx.add(E)
    D, I = idx.search(E, K + 1)
    n = len(E)
    nb = np.where(I[:, :1] == np.arange(n)[:, None], I[:, 1:], I[:, :-1])
    cs = np.where(I[:, :1] == np.arange(n)[:, None], D[:, 1:], D[:, :-1])
    rows = np.repeat(np.arange(n), K); cols = nb.reshape(-1); cos = cs.reshape(-1)
    key = set((rows * n + cols).tolist())
    mutual = np.fromiter(((c * n + r) in key for r, c in zip(rows.tolist(), cols.tolist())), bool, len(rows))
    m = mutual & (fam[rows] != fam[cols]) & (rows < cols)
    r, c, co = rows[m], cols[m], cos[m]
    rng = np.random.default_rng(rng_seed)
    by_sp = {s: np.flatnonzero(sp == s) for s in np.unique(sp)}
    cj = np.array([rng.choice(by_sp[sp[j]]) for j in c.tolist()]) if len(c) else np.array([], int)
    null = np.einsum("ij,ij->i", E[r], E[cj]) if len(c) else np.array([0.0])
    tau = float(np.percentile(null, 99))
    acc = co > tau
    fa, fb = np.minimum(fam[r[acc]], fam[c[acc]]), np.maximum(fam[r[acc]], fam[c[acc]])
    u, cnt = np.unique(fa.astype(np.int64) * (fam.max() + 1) + fb, return_counts=True)
    joins = [(int(x // (fam.max() + 1)), int(x % (fam.max() + 1))) for x in u[cnt >= 2]]
    return dict(n_mutual_cross=int(m.sum()), tau=round(tau, 4), n_accepted=int(acc.sum()), n_family_pairs_joined=len(joins)), joins


def evaluate(fam0, fam1, R, gid):
    comp0 = comp1 = 0
    fused = set()
    refs = sorted(R)
    main1 = {}
    for rid in refs:
        ii = [gid[g] for g in R[rid] if g in gid]
        comp0 += len(np.unique(fam0[ii])) == 1; comp1 += len(np.unique(fam1[ii])) == 1
        main1[rid] = set(fam1[ii].tolist())
        main1[rid + "_0"] = set(fam0[ii].tolist())
    for i, x in enumerate(refs):
        for y in refs[i + 1:]:
            if (main1[x] & main1[y]) and not (main1[x + "_0"] & main1[y + "_0"]):
                fused.add((x, y))
    sz = np.bincount(np.unique(fam1, return_inverse=True)[1])
    return dict(completeness_before=round(comp0 / len(refs), 4), completeness_after=round(comp1 / len(refs), 4),
                complete_before=int(comp0), complete_after=int(comp1), new_dev_fusions=sorted(map(list, fused)),
                refog010_fused=any("RefOG010" in p for p in fused), largest_family=int(sz.max()),
                n_families=int(len(sz)))


def main():
    t0 = time.time()
    import pandas as pd
    from ancestra_eval import families_eval as FE
    run_dir = os.path.join(DATA, "work/ob_final/main_accurate_go")
    genes = pd.read_csv(os.path.join(run_dir, "stages/input/genes.tsv"), sep="\t")
    names = genes["protein"].tolist(); gid = {x: i for i, x in enumerate(names)}
    sp = pd.factorize(genes["species"])[0]
    fam = np.load(os.path.join(run_dir, "stages/families/families.npy"))
    split = json.load(open(os.path.join(DATA, "splits/orthobench_split.json")))
    R = FE.read_refogs(os.path.join(DATA, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"), split["dev"])
    E = np.load(os.path.join(DATA, "work/gpuv2/esm_knn/emb_esm2_t12_35M_UR50D.npy")).astype(np.float32)
    E /= np.linalg.norm(E, axis=1, keepdims=True)
    out = {}
    for v in ["real"] + [f"shuf{s}" for s in SEEDS]:
        Ev = E
        if v != "real":
            rng = np.random.default_rng(int(v[4:]))
            perm = np.arange(len(E))
            for s in np.unique(sp):
                ii = np.flatnonzero(sp == s); perm[ii] = rng.permutation(ii)
            Ev = np.ascontiguousarray(E[perm])
        info, joins = run(Ev, fam, sp)
        fam1 = uf_join(fam, joins)
        info.update(evaluate(fam, fam1, R, gid))
        out[v] = info
        print(v, json.dumps(info), flush=True)
    out["manifest"] = manifest(dict(k=K, null="random same-species partner, p99, seed 11", min_links=2, mutual=True),
                               dict(emb="work/gpuv2/esm_knn/emb_esm2_t12_35M_UR50D.npy", s0=run_dir), t0)
    json.dump(out, open(os.path.join(RES, "knn_guard.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
