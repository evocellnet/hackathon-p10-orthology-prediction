#!/usr/bin/env python
"""L2-plm idea 2: self-derived ortholog/paralog ranking labels from the project's gene-order data.

No benchmark labels are used. RefOG membership (dev AND test) is read only to EXCLUDE every S0 family that holds
a RefOG gene. Within the remaining ("training") families of the frozen S0 partition (run_v2):
  anchors   = unambiguous reciprocal best hits between two species (frozen sequence bits; runner-up in the other
              species < 0.9 x best on both sides; S5 rule), computed genome-wide;
  support   = #anchors (x, y) with x within +-10 genes of a and y within +-10 genes of b, anchors whose genes lie in
              the family of a or b ignored (S5 context score, w = 10);
  positive  = anchor (a, b) inside a training family with support >= 1;
  neg (i)   = (a, b') with b' != b in the species of b and in the same family, support(a, b') = 0;
  neg (ii)  = (a, a') with a' != a in the species of a and in the same family.
Families are sampled (seed 11) until ~--max-genes genes are involved, then split 90/10 train/val by family.
Output: work/l2_plm/selflabels.npz (+ .json summary), and the gene list to embed.
    python scripts/l2_plm_selflabels.py [--max-genes 25000]
"""
import argparse, json, os, sys, time
from collections import defaultdict
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2_plm_common import DATA, RUN, WORK, RES, gene_table, manifest

W = 10


def refog_genes():
    from ancestra_eval import families_eval as FE
    split = json.load(open(os.path.join(DATA, "splits/orthobench_split.json")))
    R = FE.read_refogs(os.path.join(DATA, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"),
                       split["dev"] + split["test"], remove_uncertain=False)
    return {g for v in R.values() for g in v}


def gene_order(gt):
    contig = np.full(gt.n, -1, np.int64); order = np.full(gt.n, -1, np.int64)
    idx = {x: i for i, x in enumerate(gt.names)}
    base = 0
    for sp in gt.species:
        f = os.path.join(DATA, "data/genomes/orthobench", sp.split(".")[0] + ".parquet")
        df = pd.read_parquet(f, columns=["protein_id", "contig", "order"])
        gi = df.protein_id.map(idx); ok = gi.notna().values
        codes, uniq = pd.factorize(df.contig[ok])
        contig[gi[ok].astype(int).values] = codes + base; base += len(uniq)
        order[gi[ok].astype(int).values] = df.order[ok].values
    return contig, order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-genes", type=int, default=25000)
    ap.add_argument("--mode", default="prereg", choices=["prereg", "hard"],
                    help="hard = logged deviation D2: synteny decides sequence near-ties (see docs/L2_plm.md)")
    a = ap.parse_args()
    t0 = time.time()
    from ancestra.scoring import ScoreTable
    from ancestra.genomes import neighbour_csr
    gt = gene_table(); n = gt.n; sp = np.asarray(gt.species_of)
    fam = np.load(os.path.join(RUN, "stages/scores/families_final.npy"))
    st = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
    q, t, b = st.q.astype(np.int64), st.t.astype(np.int64), st.bits.astype(np.float64)
    m = (q != t) & (sp[q] != sp[t])
    q, t, b = q[m], t[m], b[m]
    # symmetric max bits per unordered pair
    lo, hi = np.minimum(q, t), np.maximum(q, t)
    d = pd.DataFrame(dict(k=lo * n + hi, b=b)).groupby("k").b.max()
    k = d.index.values; bb = d.values
    x, y = k // n, k % n
    D = pd.DataFrame(dict(q=np.r_[x, y], t=np.r_[y, x], b=np.r_[bb, bb]))
    D["s"] = sp[D.t.values]
    D = D.sort_values(["q", "s", "b"], ascending=[True, True, False])
    g = D.groupby(["q", "s"], sort=False)
    top = g.head(1).set_index(["q", "s"])
    second = g.nth(1).set_index(["q", "s"]).b
    top["b2"] = second.reindex(top.index).fillna(0.0).values
    top = top.reset_index()
    best = dict(zip(zip(top.q.values, top.s.values), zip(top.t.values, top.b.values, top.b2.values)))
    anchors = []
    for (qq, s), (tt, bq, b2) in best.items():
        if qq < tt and b2 < 0.9 * bq:
            r = best.get((tt, sp[qq]))
            if r and r[0] == qq and r[2] < 0.9 * r[1]:
                anchors.append((qq, tt))
    anchors = np.array(anchors, np.int64)
    print("anchors", len(anchors), round(time.time() - t0), flush=True)
    apart = defaultdict(dict)       # x -> {species: y}
    for u, v in anchors.tolist():
        apart[u][sp[v]] = v; apart[v][sp[u]] = u
    contig, order = gene_order(gt)
    ip, ix = neighbour_csr(contig, order, W)
    nb = lambda g_: ix[ip[g_]:ip[g_ + 1]]

    def support(u, v):
        if contig[u] < 0 or contig[v] < 0:
            return 0
        fv = {fam[u], fam[v]}
        Nv = set(nb(v).tolist()); s_v = sp[v]; c = 0
        for xx in nb(u).tolist():
            yy = apart.get(xx, {}).get(s_v)
            if yy is not None and yy in Nv and fam[xx] not in fv and fam[yy] not in fv:
                c += 1
        return c

    rg = refog_genes(); idx = {x_: i for i, x_ in enumerate(gt.names)}
    excl = np.unique(fam[[idx[x_] for x_ in rg if x_ in idx]])
    train_fam = ~np.isin(fam, excl)
    members = defaultdict(list)
    for gg in np.flatnonzero(train_fam):
        members[fam[gg]].append(gg)
    if a.mode == "hard":
        return hard_mode(a, t0, gt, sp, fam, train_fam, members, best, support, excl)
    # candidate positives inside training families
    pos = [(u, v) for u, v in anchors.tolist() if fam[u] == fam[v] and train_fam[u]]
    rng = np.random.default_rng(11)
    fams = sorted({fam[u] for u, _ in pos}); rng.shuffle(fams)
    posf = defaultdict(list)
    for u, v in pos:
        posf[fam[u]].append((u, v))
    rows, used, genes = [], [], set()
    for f in fams:
        if len(members[f]) > 2000:
            continue
        fr = []
        for u0, v0 in posf[f]:
            for u, v in ((u0, v0), (v0, u0)):        # both directions as queries
                if support(u, v) < 1:
                    continue
                n1 = [w for w in members[f] if sp[w] == sp[v] and w != v and support(u, w) == 0]
                n2 = [w for w in members[f] if sp[w] == sp[u] and w != u]
                if n1 or n2:
                    fr.append((u, v, n1, n2))
        if not fr:
            continue
        rows += [(f,) + r for r in fr]; used.append(f)
        for u, v, n1, n2 in fr:
            genes.update([u, v] + n1 + n2)
        if len(genes) >= a.max_genes:
            break
    used = np.array(used)
    val_f = set(rng.choice(used, max(1, len(used) // 10), replace=False).tolist())
    # flatten to (query, pos, neg, negtype, family, is_val)
    Q, P, N, T, F, V = [], [], [], [], [], []
    for f, u, v, n1, n2 in rows:
        for w, ty in [(w, 1) for w in n1] + [(w, 2) for w in n2]:
            Q.append(u); P.append(v); N.append(w); T.append(ty); F.append(f); V.append(f in val_f)
    genes = np.array(sorted(genes), np.int64)
    np.savez(os.path.join(WORK, "selflabels.npz"), q=np.array(Q), p=np.array(P), n=np.array(N), ty=np.array(T),
             fam=np.array(F), val=np.array(V), genes=genes)
    np.save(os.path.join(WORK, "genes_train.npy"), genes)
    T_, V_ = np.array(T), np.array(V)
    info = dict(n_anchors_genome=int(len(anchors)), n_excluded_families=int(len(excl)),
                n_training_family_positives_candidates=len(pos), n_families_used=int(len(used)),
                n_val_families=len(val_f), n_queries=len(rows), n_genes=int(len(genes)),
                n_triplets=len(Q), n_neg_i=int((T_ == 1).sum()), n_neg_ii=int((T_ == 2).sum()),
                n_val_triplets=int(V_.sum()),
                overlap_with_refog_genes=int(len(set(genes.tolist()) & {idx[x_] for x_ in rg if x_ in idx})))
    info["manifest"] = manifest(dict(w=W, amb=0.9, max_genes=a.max_genes, seed=11, max_family=2000),
                                dict(scores=os.path.join(RUN, "stages/scores/scores.npz"),
                                     families=os.path.join(RUN, "stages/scores/families_final.npy"),
                                     gene_order="data/genomes/orthobench/*.parquet"), t0)
    json.dump(info, open(os.path.join(RES, "selflabels.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


def hard_mode(a, t0, gt, sp, fam, train_fam, members, best, support, excl):
    """D2: for query a and species s != sp(a): near-tie set C = family members of species s with bits(a, c) >=
    0.9 x best bits of a in s (and >= 1 competitor). Positive = the unique c in C with synteny support >= 2 while
    every other c' in C has support 0; negatives = the other members of C. Same family sampling / val split."""
    import pandas as pd
    from l2_plm_discrim import seq_tables
    sn, sb = seq_tables(gt); n = gt.n
    key = lambda x, y: min(x, y) * n + max(x, y)
    rng = np.random.default_rng(11)
    fams = sorted(f for f, m in members.items() if 2 < len(m) <= 2000); rng.shuffle(fams)
    rows, used, genes = [], [], set()
    for f in fams:
        mem = members[f]; fr = []
        for u in mem:
            for s in set(sp[mem].tolist()) - {sp[u]}:
                cs = [w for w in mem if sp[w] == s]
                if len(cs) < 2:
                    continue
                bits = np.array([sb.get(key(u, w), 0.0) for w in cs])
                if bits.max() <= 0:
                    continue
                C = [w for w, b_ in zip(cs, bits) if b_ >= 0.9 * bits.max()]
                if len(C) < 2:
                    continue
                sup = [support(u, w) for w in C]
                win = [w for w, x in zip(C, sup) if x >= 2]
                if len(win) == 1 and sum(1 for x in sup if x > 0) == 1:
                    fr.append((u, win[0], [w for w in C if w != win[0]], []))
        if not fr:
            continue
        rows += [(f,) + r for r in fr]; used.append(f)
        for u, v, n1, n2 in fr:
            genes.update([u, v] + n1)
        if len(genes) >= a.max_genes:
            break
    used = np.array(used)
    val_f = set(rng.choice(used, max(1, len(used) // 10), replace=False).tolist())
    Q, P, N, T, F, V = [], [], [], [], [], []
    for f, u, v, n1, n2 in rows:
        for w in n1:
            Q.append(u); P.append(v); N.append(w); T.append(1); F.append(f); V.append(f in val_f)
    genes = np.array(sorted(genes), np.int64)
    np.savez(os.path.join(WORK, "selflabels_hard.npz"), q=np.array(Q), p=np.array(P), n=np.array(N), ty=np.array(T),
             fam=np.array(F), val=np.array(V), genes=genes)
    np.save(os.path.join(WORK, "genes_train_hard.npy"), genes)
    info = dict(mode="hard (deviation D2)", n_families_used=int(len(used)), n_val_families=len(val_f),
                n_queries=len(rows), n_genes=int(len(genes)), n_triplets=len(Q), n_val_triplets=int(np.sum(V)))
    info["manifest"] = manifest(dict(w=W, near_tie=0.9, min_support_winner=2, max_genes=a.max_genes, seed=11, max_family=2000),
                                dict(scores=os.path.join(RUN, "stages/scores/scores.npz")), t0)
    json.dump(info, open(os.path.join(RES, "selflabels_hard.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
