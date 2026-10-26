#!/usr/bin/env python
"""S5 task 3: ancestral adjacency + context_support on toy simulations, at every tree level.

At each internal node X of the true species tree:
* units = true ancestral genes of X's children (truth HOGs at the child levels);
* candidates = cross-child unit pairs in the same root family; label 1 iff the
  two units belong to the same true HOG at X (orthologs at X), else 0 (paralogs);
* unit sequence score = avg-link of leaf DIAMOND bitscores (all-vs-all in-family);
* anchors = pass 1, :func:`ancestra.synteny.select_anchors` (sequence only,
  unambiguous 1:1 mutual best, all families at once) -- or oracle anchors (true
  1:1 orthologous unit pairs) as an upper bound;
* support = :func:`ancestra.synteny.context_support` with the unit adjacency from
  :class:`AncestralAdjacency` (w genes); anchors in the candidate's family ignored.
Reports per node and per depth: AUC(support; ortholog vs paralog) on all and on
hard candidates (query unit with both, best paralog score >= 0.9 x best ortholog),
fraction supported, and the shuffled-order control (3 seeds).
Output: results/S5_synteny/ancestral_toy_<sim>.json
"""
import argparse, glob, json, os, subprocess, sys, time

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.genomes import assign_order, parse_gff3_proteins  # noqa: E402
from ancestra.newick import read_newick  # noqa: E402
from ancestra.orthoxml import hogs_at_level, read_orthoxml  # noqa: E402
from ancestra.simulate import SimParams, simulate, write_simulation  # noqa: E402
from ancestra.synteny import AncestralAdjacency, context_support, select_anchors  # noqa: E402

SEEDS = (11, 23, 47)


def auc(y, s):
    y = np.asarray(y)
    if len(y) == 0 or y.min() == y.max():
        return None
    return round(float(roc_auc_score(y, np.asarray(s, float))), 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="wgd8_s11")
    ap.add_argument("--species", type=int, default=8)
    ap.add_argument("--genes", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--wgd", default="random")
    ap.add_argument("--dup", type=float, default=0.3)
    ap.add_argument("--w", type=int, nargs="+", default=[5, 10])
    ap.add_argument("--amb", type=float, default=0.9)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    d = f"{ROOT}/work/S5/toy/{a.name}"
    if not os.path.exists(f"{d}/truth/hogs.orthoxml"):
        p = SimParams(n_species=a.species, n_genes=a.genes, wgd_branch=a.wgd if a.wgd != "none" else None,
                      dup_rate=a.dup, n_chromosomes=4)
        write_simulation(simulate(p, seed=a.seed), d)
    tree = read_newick(f"{d}/species_tree.nwk")
    sp = sorted(os.path.basename(x)[:-3] for x in glob.glob(f"{d}/proteomes/*.fa"))
    # genes, positions
    names, sp_of, cc, oo = [], [], [], []
    base = 0
    for si, s in enumerate(sp):
        df = assign_order(parse_gff3_proteins(f"{d}/gff/{s}.gff3").assign(primary=True, source="toy"))
        codes, uniq = pd.factorize(df.contig)
        names += df.protein_id.tolist(); sp_of += [si] * len(df)
        cc.append(codes + base); oo.append(df.order.values); base += len(uniq)
    cc, oo = np.concatenate(cc).astype(np.int64), np.concatenate(oo).astype(np.int64)
    sp_of = np.array(sp_of)
    N = len(names)
    idx = {x: i for i, x in enumerate(names)}
    gf = pd.read_csv(f"{d}/truth/gene_families.tsv", sep="\t")
    fam = np.full(N, -1, np.int64)
    fcode = pd.factorize(gf["family"])[0]
    for g, f in zip(gf["protein"], fcode):
        if g in idx:
            fam[idx[g]] = f
    nf = fam < 0
    fam[nf] = fcode.max() + 1 + np.arange(nf.sum())
    # all-vs-all DIAMOND (in-family pairs kept)
    allfa = f"{d}/all.fa"
    hits = f"{d}/all_vs_all.tsv"
    if not os.path.exists(hits):
        subprocess.run(f"cat {d}/proteomes/*.fa > {allfa} && diamond blastp -q {allfa} -d {allfa} -o {hits} "
                       f"--sensitive -k 0 --max-hsps 1 -e 10 -p {a.threads} --quiet --outfmt 6 qseqid sseqid bitscore",
                       shell=True, check=True)
    h = pd.read_csv(hits, sep="\t", header=None, names=["q", "t", "b"])
    q = h.q.map(idx).values.astype(np.int64); t = h.t.map(idx).values.astype(np.int64); b = h.b.values
    k = (q != t) & (fam[q] == fam[t])
    q, t, b = q[k], t[k], b[k]
    lo, hi = np.minimum(q, t), np.maximum(q, t)
    s = pd.Series(b).groupby(lo * N + hi).max()
    pk = s.index.values
    pa_, pb_, ps_ = pk // N, pk % N, s.values
    ox = read_orthoxml(f"{d}/truth/hogs.orthoxml")
    adjs = {(w, None): AncestralAdjacency(cc, oo, w, fam) for w in a.w}
    for w in a.w:
        for sd in SEEDS:
            adjs[(w, sd)] = AncestralAdjacency(cc, oo, w, fam, shuffle_seed=sd)
    depth = {}
    for n in tree.preorder():
        depth[n.name] = 0 if n.parent is None else depth[n.parent.name] + 1
    rows = []
    for X in tree.postorder():
        if not X.children:
            continue
        kids = X.children
        hog_x = hogs_at_level(ox, X.name, species_tree=tree)
        gx = np.full(N, -1, np.int64)
        for i, hg in enumerate(hog_x):
            for g in hg:
                gx[idx[g]] = i
        gene_unit = np.full(N, -1, np.int64)
        unit_child = []
        nu = 0
        for ci, C in enumerate(kids):
            for hg in hogs_at_level(ox, C.name, species_tree=tree):
                ii = [idx[g] for g in hg]
                gene_unit[ii] = nu
                unit_child.append(ci)
                nu += 1
        unit_child = np.array(unit_child)
        # unit -> X-HOG (consistent: all genes of a unit in one X HOG)
        unit_x = np.full(nu, -1, np.int64)
        m = gene_unit >= 0
        unit_x[gene_unit[m]] = gx[m]
        unit_f = np.full(nu, -1, np.int64)
        unit_f[gene_unit[m]] = fam[m]
        # unit-level avg-link scores for cross-child pairs (new pairs at X)
        ua, ub = gene_unit[pa_], gene_unit[pb_]
        ok = (ua >= 0) & (ub >= 0)
        ua, ub, us = ua[ok], ub[ok], ps_[ok]
        ok = unit_child[ua] != unit_child[ub]
        ua, ub, us = ua[ok], ub[ok], us[ok]
        x, y = np.minimum(ua, ub), np.maximum(ua, ub)
        agg = pd.DataFrame({"k": x * nu + y, "s": us}).groupby("k").s.mean()
        ca, cb = agg.index.values // nu, agg.index.values % nu
        cs = agg.values
        lab = (unit_x[ca] == unit_x[cb]).astype(int)
        pa, pb = select_anchors(ca, cb, cs, unit_child, a.amb)
        # oracle anchors: true 1:1 ortholog unit pairs (X-HOG with exactly one unit per child pair)
        o_mask = lab == 1
        cnt = pd.Series(np.r_[ca[o_mask], cb[o_mask]]).value_counts()
        one = cnt[cnt == 1].index.values
        om = o_mask & np.isin(ca, one) & np.isin(cb, one)
        oa, ob = ca[om], cb[om]
        # hard: per query unit with both labels, best paralog >= 0.9 best ortholog
        dfc = pd.DataFrame({"a": ca, "b": cb, "s": cs, "lab": lab})
        hard = np.zeros(len(ca), bool)
        for side in ("a", "b"):
            g = dfc.groupby(side)
            mo = g.apply(lambda z: z.s[z.lab == 1].max() if (z.lab == 1).any() else np.nan)
            mp = g.apply(lambda z: z.s[z.lab == 0].max() if (z.lab == 0).any() else np.nan)
            hq = mo.index[(mp >= 0.9 * mo).fillna(False).values]
            hard |= dfc[side].isin(hq).values
        anchor_precision = float(np.mean(unit_x[pa] == unit_x[pb])) if len(pa) else None
        rec = dict(node=X.name, depth=depth[X.name], n_children=len(kids), n_units=int(nu),
                   n_cand=int(len(ca)), n_orth=int(lab.sum()), n_hard=int(hard.sum()),
                   n_anchors=int(len(pa)), anchor_precision=anchor_precision,
                   n_oracle_anchors=int(len(oa)), seq_auc=auc(lab, cs), seq_auc_hard=auc(lab[hard], cs[hard]))
        for w in a.w:
            nadj = adjs[(w, None)].unit_adjacency(gene_unit, nu)
            sup = context_support((ca, cb), (pa, pb), nadj)
            supo = context_support((ca, cb), (oa, ob), nadj)
            rec[f"w{w}"] = dict(auc=auc(lab, sup), auc_hard=auc(lab[hard], sup[hard]),
                                auc_oracle=auc(lab, supo), auc_hard_oracle=auc(lab[hard], supo[hard]),
                                frac_orth_sup=round(float((sup[lab == 1] > 0).mean()), 4) if lab.any() else None,
                                frac_par_sup=round(float((sup[lab == 0] > 0).mean()), 4) if (lab == 0).any() else None,
                                n_adj_edges=int(nadj.n_edges))
            sh, shh = [], []
            for sd in SEEDS:
                nadj2 = adjs[(w, sd)].unit_adjacency(gene_unit, nu)
                sup2 = context_support((ca, cb), (pa, pb), nadj2)
                sh.append(auc(lab, sup2)); shh.append(auc(lab[hard], sup2[hard]))
            rec[f"w{w}"]["shuf_auc"] = _m(sh)
            rec[f"w{w}"]["shuf_auc_hard"] = _m(shh)
            # pooled arrays for depth summaries
            rec.setdefault("_arr", {})[w] = (lab, sup, hard)
        rows.append(rec)
        print(X.name, {k: v for k, v in rec.items() if k != "_arr"}, flush=True)
    # by depth (pooled)
    by_depth = {}
    for dd in sorted({r["depth"] for r in rows}):
        rr = [r for r in rows if r["depth"] == dd]
        ent = dict(nodes=[r["node"] for r in rr], n_cand=sum(r["n_cand"] for r in rr),
                   n_orth=sum(r["n_orth"] for r in rr), n_hard=sum(r["n_hard"] for r in rr))
        for w in a.w:
            lab = np.concatenate([r["_arr"][w][0] for r in rr]); sup = np.concatenate([r["_arr"][w][1] for r in rr])
            hd = np.concatenate([r["_arr"][w][2] for r in rr])
            ent[f"w{w}_auc"] = auc(lab, sup); ent[f"w{w}_auc_hard"] = auc(lab[hd], sup[hd])
        by_depth[dd] = ent
    for r in rows:
        r.pop("_arr", None)
    out = f"{ROOT}/results/S5_synteny/ancestral_toy_{a.name}.json"
    json.dump(dict(sim=d, params=json.load(open(f"{d}/params.json")), windows=a.w, amb=a.amb,
                   shuffle_seeds=SEEDS, nodes=rows, by_depth=by_depth, wall_s=round(time.time() - t0, 1)),
              open(out, "w"), indent=1, default=str)
    print(json.dumps(by_depth, indent=1))
    print("wrote", out)


def _m(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.mean(xs)), 4) if xs else None


if __name__ == "__main__":
    main()
