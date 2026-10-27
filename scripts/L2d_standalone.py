"""L2-domains step (a): standalone ortholog-vs-paralog signal of Pfam domain architectures (DEV only).

    python scripts/L2d_standalone.py --domtbl work/l2_domains/pfam_stage1/pfam.domtbl --out results/l2_domains

OrthoBench dev (work/S5/refog_pair_labels.parquet: within-RefOG + deep pairs) and TreeFam-A dev reference
pairs. Score J = Jaccard of clan sets, ident = identical collapsed ordered architecture. Hard queries
(docs/S5_synteny.md definition, with the pipeline's norm_max leaf score). Controls: architectures permuted
among annotated genes of the same species (seeds 11/23/47), and within species x S0 family.
Pre-registration: docs/L2_domains.md §0.3. Benchmark-only code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
WT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, WT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.domains import DomainArch, load_clans, parse_domtbl  # noqa: E402
from ancestra.io import read_gene_table  # noqa: E402
from ancestra.merge import leaf_scores  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402

SEEDS = (11, 23, 47)
CLANS = os.path.join(ROOT, "data/extra/l2_domains/Pfam-A.clans.tsv.gz")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def pair_scores(gt, st, ia, ib):
    """norm_max leaf score for gene index pairs (0 if unscored)."""
    s = leaf_scores(gt, st, "norm_max")
    q, t = st.q.astype(np.int64), st.t.astype(np.int64)
    need = np.unique(np.concatenate([ia, ib]))
    m = np.isin(q, need) & np.isin(t, need)
    N = gt.n
    key = np.minimum(q[m], t[m]) * N + np.maximum(q[m], t[m])
    o = np.argsort(key)
    key, val = key[o], s[m][o]
    k2 = np.minimum(ia, ib) * N + np.maximum(ia, ib)
    pos = np.searchsorted(key, k2)
    pos = np.clip(pos, 0, max(len(key) - 1, 0))
    hit = (len(key) > 0) & (key[pos] == k2) if len(key) else np.zeros(len(k2), bool)
    out = np.where(hit, val[pos], 0.0)
    return out


def hard_mask(df):
    """S5 hard queries: (query gene, other species) with >=1 ortholog and >=1 paralog candidate and
    best paralog score >= 0.9 x best ortholog score. Returns (pair mask, list of query groups)."""
    rows = []
    for col_q, col_o, sp_o in (("gene_a", "gene_b", "sp_b"), ("gene_b", "gene_a", "sp_a")):
        rows.append(pd.DataFrame({"q": df[col_q].values, "osp": df[sp_o].values, "i": np.arange(len(df))}))
    r = pd.concat(rows, ignore_index=True)
    r["label"] = df.label.values[r.i]
    r["seq"] = df.seq.values[r.i]
    hard = np.zeros(len(df), bool)
    queries = []
    for (q, osp), g in r.groupby(["q", "osp"], sort=True):
        lab = g.label.values
        if lab.max() < 1 or lab.min() > 0:
            continue
        bo, bp = g.seq.values[lab == 1].max(), g.seq.values[lab == 0].max()
        if bp >= 0.9 * bo:
            hard[g.i.values] = True
            queries.append(g.i.values)
    return hard, queries


def auc(y, s):
    y = np.asarray(y)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, s))


def boot_auc(y, s, grp, B=1000, seed=11):
    """Percentile CI, resampling groups (RefOGs / families)."""
    y, s, grp = np.asarray(y), np.asarray(s), np.asarray(grp)
    ug = np.unique(grp)
    idx = {g: np.flatnonzero(grp == g) for g in ug}
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(B):
        pick = rng.choice(ug, len(ug), replace=True)
        ii = np.concatenate([idx[g] for g in pick])
        v = auc(y[ii], s[ii])
        if np.isfinite(v):
            vals.append(v)
    return [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)]


def rerank(df, queries, J, lam=0.3, window=0.75):
    """Hard-query top-1: sequence vs score*(1+lam*J) among candidates >= window x best."""
    seq_ok, comb_ok, grp = [], [], []
    Jz = np.nan_to_num(J, nan=0.0)
    for ii in queries:
        s = df.seq.values[ii]
        lab = df.label.values[ii]
        o = np.lexsort((ii, -s))           # deterministic tie order
        seq_ok.append(lab[o[0]])
        c = s * (1 + np.where(s >= window * s.max(), lam * Jz[ii], 0.0))
        o2 = np.lexsort((ii, -s, -c))
        comb_ok.append(lab[o2[0]])
        grp.append(df.grp.values[ii[0]])
    return np.array(seq_ok), np.array(comb_ok), np.array(grp)


def boot_diff(a, b, grp, B=2000, seed=11):
    ug = np.unique(grp)
    idx = {g: np.flatnonzero(grp == g) for g in ug}
    rng = np.random.default_rng(seed)
    d = []
    for _ in range(B):
        ii = np.concatenate([idx[g] for g in rng.choice(ug, len(ug), replace=True)])
        d.append(b[ii].mean() - a[ii].mean())
    return [round(float(np.percentile(d, 2.5)), 4), round(float(np.percentile(d, 97.5)), 4)]


def evaluate(df, archs: dict, hard, queries):
    """archs: name -> DomainArch. df has ia, ib, label, grp."""
    out = {}
    for name, da in archs.items():
        J, same = da.gene_sim(df.ia.values, df.ib.values)
        ok = np.isfinite(J)
        r = dict(coverage=round(float(ok.mean()), 4),
                 coverage_orth=round(float(ok[df.label.values == 1].mean()), 4))
        for sub, m in (("all", ok), ("hard", ok & hard)):
            y = df.label.values[m]
            r[f"n_{sub}"] = int(m.sum())
            r[f"n_orth_{sub}"] = int(y.sum())
            r[f"aucJ_{sub}"] = round(auc(y, J[m]), 4)
            r[f"aucIdent_{sub}"] = round(auc(y, same[m].astype(float)), 4)
            r[f"meanJ_orth_{sub}"] = round(float(J[m][y == 1].mean()), 4) if y.sum() else None
            r[f"meanJ_par_{sub}"] = round(float(J[m][y == 0].mean()), 4) if (y == 0).sum() else None
            r[f"ident_orth_{sub}"] = round(float(same[m][y == 1].mean()), 4) if y.sum() else None
            r[f"ident_par_{sub}"] = round(float(same[m][y == 0].mean()), 4) if (y == 0).sum() else None
            if name == "real":
                r[f"aucJ_{sub}_ci"] = boot_auc(y, J[m], df.grp.values[m])
        # sequence AUC on the same pairs (reference)
        if name == "real":
            for sub, m in (("all", ok), ("hard", ok & hard)):
                r[f"aucSeq_{sub}"] = round(auc(df.label.values[m], df.seq.values[m]), 4)
        so, co, grp = rerank(df, queries, J)
        r["hard_queries"] = len(queries)
        r["top1_seq"] = round(float(so.mean()), 4)
        r["top1_comb"] = round(float(co.mean()), 4)
        r["top1_delta"] = round(float(co.mean() - so.mean()), 4)
        r["top1_delta_ci"] = boot_diff(so.astype(float), co.astype(float), grp)
        out[name] = r
    return out


def build_archs(gt, arch_by_name, family):
    archs = {"real": DomainArch(gt.names, gt.species_of, arch_by_name)}
    for s in SEEDS:
        archs[f"shuf{s}"] = DomainArch(gt.names, gt.species_of, arch_by_name, shuffle_seed=s)
    for s in SEEDS:
        archs[f"shufFam{s}"] = DomainArch(gt.names, gt.species_of, arch_by_name, shuffle_seed=s,
                                          scope="family", family=family)
    return archs


def summarize_controls(res):
    for sub in ("all", "hard"):
        sh = [res[f"shuf{s}"][f"aucJ_{sub}"] for s in SEEDS]
        shf = [res[f"shufFam{s}"][f"aucJ_{sub}"] for s in SEEDS]
        res["real"][f"shuf_max_{sub}"] = max(sh)
        res["real"][f"shufFam_max_{sub}"] = max(shf)
        res["real"][f"gate_{sub}"] = bool(res["real"][f"aucJ_{sub}_ci"][0] > max(sh))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    clans = load_clans(CLANS)
    arch_by_name = parse_domtbl(a.domtbl, clans)
    os.makedirs(a.out, exist_ok=True)
    result = {}
    # ---------------- OrthoBench dev
    run = os.path.join(ROOT, "work/core_v0_ob/run_det")
    gt = read_gene_table(os.path.join(run, "stages/input/genes.tsv"))
    st = ScoreTable.load(os.path.join(run, "stages/scores/scores.npz"))
    fam = np.load(os.path.join(run, "stages/scores/families_final.npy"))
    d = pd.read_parquet(os.path.join(ROOT, "work/S5/refog_pair_labels.parquet"))
    idx = gt.index()
    d = d[d.gene_a.isin(idx) & d.gene_b.isin(idx)].reset_index(drop=True)
    d["ia"] = d.gene_a.map(idx).astype(np.int64)
    d["ib"] = d.gene_b.map(idx).astype(np.int64)
    d["seq"] = pair_scores(gt, st, d.ia.values, d.ib.values)
    d["grp"] = d.refog_a
    hard, queries = hard_mask(d)
    archs = build_archs(gt, arch_by_name, fam)
    ob = evaluate(d, archs, hard, queries)
    ob = summarize_controls(ob)
    # within-RefOG only (no deep pairs)
    w = d.kind.values == "within"
    dw = d[w].reset_index(drop=True)
    hw, qw = hard_mask(dw)
    ob_w = summarize_controls(evaluate(dw, archs, hw, qw))
    # family-fusion signal: same RefOG (any label) vs different RefOG (deep), all pairs
    J, same = archs["real"].gene_sim(d.ia.values, d.ib.values)
    ok = np.isfinite(J)
    y_same = (d.kind.values == "within").astype(int)
    fusion = dict(n=int(ok.sum()), aucJ_sameRefOG=round(auc(y_same[ok], J[ok]), 4),
                  aucJ_sameRefOG_ci=boot_auc(y_same[ok], J[ok], d.grp.values[ok]),
                  aucSeq_sameRefOG=round(auc(y_same[ok], d.seq.values[ok]), 4),
                  shuf=[round(auc(y_same[ok], np.nan_to_num(archs[f"shuf{s}"].gene_sim(
                      d.ia.values, d.ib.values)[0][ok], nan=0)), 4) for s in SEEDS])
    result["orthobench_dev"] = dict(n_pairs=len(d), n_hard=int(hard.sum()), n_hard_queries=len(queries),
                                    all_pairs=ob, within_only=ob_w, family_fusion=fusion,
                                    annotated_in_scope=int((archs["real"].code >= 0).sum()))
    print("OB", json.dumps({k: ob["real"][k] for k in ob["real"] if "auc" in k or "top1" in k or "gate" in k
                             or "shuf" in k or "cover" in k}))
    # ---------------- TreeFam-A dev
    from qfo_treefam_dev_fast import dev_refpairs
    refs = dev_refpairs()
    qrun = os.path.join(ROOT, "work/qfo_recall/base")
    gq = read_gene_table(os.path.join(qrun, "stages/input/genes.tsv"))
    sq = ScoreTable.load(os.path.join(qrun, "stages/scores/scores.npz"))
    fq = np.load(os.path.join(qrun, "stages/scores/families_final.npy"))
    iq = gq.index()
    t = pd.DataFrame(refs, columns=["gene_a", "gene_b", "t", "grp"])
    t = t[t.gene_a.isin(iq) & t.gene_b.isin(iq)].reset_index(drop=True)
    t["ia"] = t.gene_a.map(iq).astype(np.int64)
    t["ib"] = t.gene_b.map(iq).astype(np.int64)
    t["label"] = (t.t == "S").astype(int)
    t["sp_a"] = [gq.species[x] for x in gq.species_of[t.ia.values]]
    t["sp_b"] = [gq.species[x] for x in gq.species_of[t.ib.values]]
    t["seq"] = pair_scores(gq, sq, t.ia.values, t.ib.values)
    ht, qt = hard_mask(t)
    archs_q = build_archs(gq, arch_by_name, fq)
    tf = summarize_controls(evaluate(t, archs_q, ht, qt))
    result["treefamA_dev"] = dict(n_pairs=len(t), n_hard=int(ht.sum()), n_hard_queries=len(qt), res=tf,
                                  annotated_in_scope=int((archs_q["real"].code >= 0).sum()))
    print("TF", json.dumps({k: tf["real"][k] for k in tf["real"] if "auc" in k or "top1" in k or "gate" in k
                             or "shuf" in k or "cover" in k}))
    gate = any(result[b]["all_pairs" if b == "orthobench_dev" else "res"]["real"]["gate_hard"]
               for b in ("orthobench_dev", "treefamA_dev"))
    result["gate_step_b"] = bool(gate)
    git = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    result["manifest"] = dict(git_commit=git, script=os.path.relpath(__file__, WT), argv=sys.argv,
                              inputs={"domtbl": a.domtbl, "domtbl_md5": md5(a.domtbl), "clans_md5": md5(CLANS),
                                      "ob_run": run, "qfo_run": qrun},
                              config=dict(seeds=SEEDS, rerank_lam=0.3, window=0.75, hard_ratio=0.9,
                                          boot_B_auc=1000, boot_B_delta=2000, boot_seed=11),
                              runtime_s=round(time.time() - t0, 1))
    json.dump(result, open(os.path.join(a.out, "standalone.json"), "w"), indent=1)
    print("gate_step_b", gate)


if __name__ == "__main__":
    main()
