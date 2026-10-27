#!/usr/bin/env python
"""L2-network step (a): standalone signal of network-neighbourhood conservation (interologs) and
of protein-abundance similarity (PaxDb) for ortholog vs paralog candidates, OrthoBench DEV only.

Truth / candidates / hard cases are reused from S5/S9 (work/S5/refog_pair_labels.parquet; S7
canonical DIAMOND bits as the sequence score; S5 query decision rule). Pre-registration:
docs/L2_network.md §1.

Interolog score of a candidate (a in A, b in B):
  anchors   = 1:1 ortholog pairs of the FROZEN Ancestra accurate run (work/ob_final/
              main_accurate_go/pairs.tsv): x has exactly one predicted ortholog in B and it has
              exactly one in A. Data-derived, never benchmark labels. Anchors whose genes lie in
              the S0 family of a or b are ignored (this also drops the pair being scored).
  count     = #{(x, y) anchors : x ~ a in net_A, y ~ b in net_B}
  cos       = count / sqrt(n_a * n_b), n_a / n_b = anchored (family-filtered) neighbours of a / b
  evaluable = n_a > 0 and n_b > 0
Shuffled control: node labels permuted within species among the channel's network nodes
(degree sequence and coverage kept), seeds 11/23/47.

Abundance score: -|q_a - q_b|, q = within-species percentile rank of PaxDb ppm (integrated
whole-organism set). Shuffled: abundances permuted among a species' covered genes.
"""
import json, os, sys, time

import numpy as np
import pandas as pd

MAIN = os.environ.get("P10_MAIN", "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction")
WT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{MAIN}/scripts")
from S5_pair_benchmark import auc, bootstrap_auc, pooled_queries, query_decisions, query_table  # noqa: E402

W = f"{MAIN}/work/l2_network"
OUT = f"{WT}/results/l2_network"
RUN = f"{MAIN}/work/ob_final/main_accurate_go"
SEEDS = (11, 23, 47)
HARD = 0.9
MIN_NODES = 100          # a species' channel network is used only with >= this many nodes

# channels (pre-registered, docs/L2_network.md §1); 'db' is a flagged circular reference
CHANNELS = {
    "exp150": lambda e: e.exp >= 150,
    "exp400": lambda e: e.exp >= 400,
    "coexp400": lambda e: e.coexp >= 400,
    "own": lambda e: (e.exp >= 150) | (e.coexp >= 400),    # PRIMARY: all own-species evidence
    "db400_CIRCULAR": lambda e: e.db >= 400,
}
PRIMARY = "own"
# secondary run (pre-registration §1.6): human network + hu.MAP 3.0 complex co-membership pairs
# (CF-MS/AP-MS integration); other species keep their own channel. L2_RUN=humap selects it.
RUNMODE = os.environ.get("L2_RUN", "main")
if RUNMODE == "humap":
    CHANNELS = {"humap_only": CHANNELS["own"], "own_humap": CHANNELS["own"]}
    PRIMARY = "own_humap"


def load_species():
    g = pd.read_csv(f"{RUN}/stages/input/genes.tsv", sep="\t")
    g["sp"] = g.species.str.split(".").str[0]
    fam = np.load(f"{RUN}/stages/scores/families_final.npy")
    return dict(zip(g.protein, g.sp)), dict(zip(g.protein, fam))


def anchors(sp_of):
    p = pd.read_csv(f"{RUN}/pairs.tsv", sep="\t", comment="#", header=None, names=["a", "b", "level"],
                    usecols=[0, 1])
    p["sa"] = p.a.map(sp_of); p["sb"] = p.b.map(sp_of)
    sw = (p.sa > p.sb).values                     # orient every pair by species name
    p.loc[sw, ["a", "b", "sa", "sb"]] = p.loc[sw, ["b", "a", "sb", "sa"]].values
    p = p.drop_duplicates(["a", "b"])
    out = {}
    for (sa, sb), d in p.groupby(["sa", "sb"]):
        # 1:1 within this species pair
        ca = d.a.value_counts(); cb = d.b.value_counts()
        d = d[d.a.map(ca).eq(1) & d.b.map(cb).eq(1)]
        out[(sa, sb)] = dict(zip(d.a, d.b)); out[(sb, sa)] = dict(zip(d.b, d.a))
    return out


def networks():
    nets = {}
    for f in sorted(os.listdir(W)):
        if f.startswith("edges_"):
            sp = f[6:-8]
            e = pd.read_parquet(f"{W}/{f}")
            for ch, fn in CHANNELS.items():
                x = e[fn(e)]
                if sp == "Homo_sapiens" and ch in ("own_humap", "humap_only"):
                    hm = pd.read_parquet(f"{W}/humap_Homo_sapiens.parquet")[["a", "b"]]
                    x = hm if ch == "humap_only" else pd.concat([x[["a", "b"]], hm]).drop_duplicates()
                nodes = set(x.a) | set(x.b)
                nets[(sp, ch)] = (x.a.values, x.b.values, nodes)
    return nets


def adj_from(a, b, perm=None):
    adj = {}
    if perm is not None:
        a = [perm[x] for x in a]; b = [perm[x] for x in b]
    for x, y in zip(a, b):
        adj.setdefault(x, set()).add(y); adj.setdefault(y, set()).add(x)
    return adj


def permutation(nodes, seed):
    n = sorted(nodes)
    rng = np.random.default_rng(seed)
    return dict(zip(n, [n[i] for i in rng.permutation(len(n))]))


def canonical_seq_scores():
    h = pd.read_csv(f"{MAIN}/work/S7/hits.tsv", sep="\t", header=None, usecols=[0, 1, 2], names=["q", "t", "bits"])
    q = h.q.str.split("|", n=1, expand=True); t = h.t.str.split("|", n=1, expand=True)
    m = (q[0] == q[1]) & (t[0] == t[1]) & (q[0] != t[0])
    a, b = q[0][m].values, t[0][m].values
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    return pd.Series(h.bits[m].values).groupby([lo, hi]).max()


def candidates(A, B, lab, canon):
    c = lab[((lab.sp_a == A) & (lab.sp_b == B)) | ((lab.sp_a == B) & (lab.sp_b == A))].copy()
    sw = (c.sp_a != A).values
    for x, y in (("gene_a", "gene_b"), ("sp_a", "sp_b"), ("refog_a", "refog_b")):
        c.loc[sw, [x, y]] = c.loc[sw, [y, x]].values
    c = c.drop_duplicates(["gene_a", "gene_b"]).reset_index(drop=True)
    lo = np.where(c.gene_a < c.gene_b, c.gene_a, c.gene_b); hi = np.where(c.gene_a < c.gene_b, c.gene_b, c.gene_a)
    c["seq"] = canon.reindex(pd.MultiIndex.from_arrays([lo, hi])).fillna(0.0).values
    return c[["gene_a", "gene_b", "sp_a", "sp_b", "refog_a", "refog_b", "label", "kind", "seq"]]


def mark_hard(c, mask=None):
    """S5 hard flag (best paralog bits >= 0.9 x best ortholog bits), on the rows in ``mask``."""
    hard = np.zeros(len(c), bool)
    d0 = c if mask is None else c[mask]
    for side in ("gene_a", "gene_b"):
        for g, d in d0.groupby(side):
            o, p = d[d.label == 1], d[d.label == 0]
            if len(o) and len(p) and p.seq.max() >= HARD * o.seq.max():
                hard[d.index.values] = True
    return hard


def interolog(c, adjA, adjB, M, fam):
    cnt = np.zeros(len(c)); cos = np.zeros(len(c)); na_ = np.zeros(len(c), int); nb_ = np.zeros(len(c), int)
    Minv = {v: k for k, v in M.items()}
    for i, (a, b) in enumerate(zip(c.gene_a.values, c.gene_b.values)):
        Na = adjA.get(a); Nb = adjB.get(b)
        if not Na or not Nb:
            continue
        excl = {fam.get(a), fam.get(b)}
        ma = {M[x] for x in Na if x in M and fam.get(x) not in excl and fam.get(M[x]) not in excl}
        nb = {y for y in Nb if y in Minv and fam.get(y) not in excl and fam.get(Minv[y]) not in excl}
        na_[i], nb_[i] = len(ma), len(nb)
        if ma and nb:
            k = len(ma & nb)
            cnt[i] = k; cos[i] = k / np.sqrt(len(ma) * len(nb))
    return cnt, cos, na_, nb_


def abundance_scores(c, ab, A, B):
    qa, qb = ab[A], ab[B]
    xa = c.gene_a.map(qa).values.astype(float); xb = c.gene_b.map(qb).values.astype(float)
    return -np.abs(xa - xb)


def load_abundance():
    ab, raw = {}, {}
    for f in sorted(os.listdir(W)):
        if f.startswith("abund_"):
            sp = f[6:-8]
            d = pd.read_parquet(f"{W}/{f}")
            d = d[d.ppm > 0]
            raw[sp] = d
            ab[sp] = dict(zip(d.ob, d.ppm.rank(pct=True).values))
    return ab


def shuffled_abundance(ab, seed):
    out = {}
    for i, (sp, d) in enumerate(sorted(ab.items())):
        k = sorted(d); v = np.array([d[x] for x in k])
        rng = np.random.default_rng(seed + i)
        out[sp] = dict(zip(k, v[rng.permutation(len(v))]))
    return out


def query_rows(c, col, near):
    c = c.assign(a=c.gene_a.values, b=c.gene_b.values)
    qs = [x + (False,) for x in query_table(c)]
    return query_decisions(c, qs, col, near=near, return_rows=True)


def boot_diff(ra, rb, key, B=1000, seed=11):
    a = np.array([r[key] for r in ra], float); b = np.array([r[key] for r in rb], float)
    if len(a) == 0:
        return None, None
    cl = np.array([str(r["refog"]) for r in ra]); u = np.unique(cl); g = {k: np.flatnonzero(cl == k) for k in u}
    rng = np.random.default_rng(seed); ds = []
    for _ in range(B):
        ix = np.concatenate([g[k] for k in rng.choice(u, len(u), replace=True)])
        ds.append(a[ix].mean() - b[ix].mean())
    return round(float(a.mean() - b.mean()), 4), [round(float(np.percentile(ds, 2.5)), 4), round(float(np.percentile(ds, 97.5)), 4)]


def evaluate(c, score_cols, cov_col, label):
    """AUC (covered pairs) + hard-query decisions (hard flag recomputed on covered candidates)."""
    c = c.copy(); c["cluster"] = c.refog_a.astype(str)
    e = c[c[cov_col]].reset_index(drop=True)
    res = dict(set=label, n_candidates=int(len(c)), n_covered=int(len(e)),
               frac_covered=round(float(len(e) / max(1, len(c))), 4),
               n_orth_covered=int((e.label == 1).sum()), n_orth_total=int((c.label == 1).sum()),
               n_par_covered=int((e.label == 0).sum()))
    if len(e) == 0 or e.label.nunique() < 2:
        return res
    e["hard"] = mark_hard(e)
    h = e.hard.values
    res["n_hard_pairs_covered"] = int(h.sum())
    res["auc"] = {}
    for col in ["seq"] + score_cols:
        res["auc"][col] = dict(all=auc(e.label, e[col]), all_ci95=bootstrap_auc(e, col, np.ones(len(e), bool)),
                               hard=auc(e.label[h], e[col][h]) if h.any() and e.label[h].nunique() == 2 else None,
                               hard_ci95=bootstrap_auc(e, col, h) if h.any() else None)
        if col != "seq":
            res["auc"][col]["shuf_all"] = [auc(e.label, e[f"shuf{s}_{col}"]) for s in SEEDS]
            res["auc"][col]["shuf_hard"] = [auc(e.label[h], e[f"shuf{s}_{col}"][h]) if h.any() and e.label[h].nunique() == 2 else None for s in SEEDS]
    res["queries"] = {}
    for qset, near in (("hard", HARD), ("all", 0.0)):
        allrows = {}
        for col in score_cols + [f"shuf{s}_{m}" for s in SEEDS for m in score_cols]:
            rows = []
            for pr, d in e.groupby("pair", sort=False):
                _, r = query_rows(d.reset_index(drop=True), col, near)
                for x in r:
                    x.update(pair=pr, col=col)
                rows += r
            allrows[col] = rows
        pq = pooled_queries([{"_qrows": sum(allrows.values(), [])}])
        q = {}
        for m in score_cols:
            if m not in pq:
                continue
            d = {k: pq[m][k] for k in ("n", "seq_acc", "syn_correct", "syn_wrong", "syn_abstain",
                                       "syn_precision_when_decides", "seq_acc_when_syn_decides", "comb0.9_acc",
                                       "comb0.9_minus_seq", "comb0.9_minus_seq_ci95", "comb0.75_acc",
                                       "comb0.75_minus_seq", "comb0.75_minus_seq_ci95")}
            d["shuf_precision_when_decides"] = [pq[f"shuf{s}_{m}"]["syn_precision_when_decides"] for s in SEEDS]
            d["shuf_decides"] = [pq[f"shuf{s}_{m}"]["syn_correct"] + pq[f"shuf{s}_{m}"]["syn_wrong"] for s in SEEDS]
            for r in ("0.9", "0.75"):
                d[f"shuf_comb{r}_minus_seq"] = [pq[f"shuf{s}_{m}"][f"comb{r}_minus_seq"] for s in SEEDS]
                d[f"comb{r}_real_minus_shuf"] = [boot_diff(allrows[m], allrows[f"shuf{s}_{m}"], f"comb{r}") for s in SEEDS]
            q[m] = d
        res["queries"][qset] = q
    return res


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    sp_of, fam = load_species()
    anc = anchors(sp_of)
    nets = networks()
    ab = load_abundance()
    ab_shuf = {s: shuffled_abundance(ab, s) for s in SEEDS}
    lab = pd.read_parquet(f"{MAIN}/work/S5/refog_pair_labels.parquet")
    canon = canonical_seq_scores()
    sps = sorted({k[0] for k in nets})
    adj, adj_shuf = {}, {}
    for (sp, ch), (a, b, nodes) in nets.items():
        if len(nodes) < MIN_NODES:
            continue
        adj[(sp, ch)] = adj_from(a, b)
        for i, s in enumerate(SEEDS):
            adj_shuf[(sp, ch, s)] = adj_from(a, b, permutation(nodes, s * 1000 + sps.index(sp)))
    tabs, per_pair = [], {}
    for i, A in enumerate(sps):
        for B in sps[i + 1:]:
            c = candidates(A, B, lab, canon)
            if len(c) == 0:
                continue
            M = anc.get((A, B), {})
            pp = dict(n=int(len(c)), n_orth=int((c.label == 1).sum()), n_anchors=len(M))
            for ch in CHANNELS:
                ok = (A, ch) in adj and (B, ch) in adj
                if ok:
                    cnt, cos, na_, nb_ = interolog(c, adj[(A, ch)], adj[(B, ch)], M, fam)
                else:
                    cnt = cos = np.zeros(len(c)); na_ = nb_ = np.zeros(len(c), int)
                c[f"{ch}_count"] = cnt; c[f"{ch}_cos"] = cos
                c[f"{ch}_eval"] = (na_ > 0) & (nb_ > 0)
                c[f"{ch}_bothnodes"] = (c.gene_a.isin(adj[(A, ch)].keys()) & c.gene_b.isin(adj[(B, ch)].keys())) if ok else False
                for s in SEEDS:
                    if ok:
                        k1, k2, _, _ = interolog(c, adj_shuf[(A, ch, s)], adj_shuf[(B, ch, s)], M, fam)
                    else:
                        k1 = k2 = np.zeros(len(c))
                    c[f"shuf{s}_{ch}_count"] = k1; c[f"shuf{s}_{ch}_cos"] = k2
                pp[ch] = dict(both_nodes=round(float(np.mean(c[f"{ch}_bothnodes"])), 4),
                              evaluable=round(float(c[f"{ch}_eval"].mean()), 4),
                              count_gt0=round(float((c[f"{ch}_count"] > 0).mean()), 4),
                              orth_count_gt0=round(float((c[f"{ch}_count"][c.label == 1] > 0).mean()), 4) if (c.label == 1).any() else None,
                              par_count_gt0=round(float((c[f"{ch}_count"][c.label == 0] > 0).mean()), 4) if (c.label == 0).any() else None)
            # abundance
            if A in ab and B in ab:
                c["abund"] = abundance_scores(c, ab, A, B)
                for s in SEEDS:
                    c[f"shuf{s}_abund"] = abundance_scores(c, ab_shuf[s], A, B)
            else:
                c["abund"] = np.nan
                for s in SEEDS:
                    c[f"shuf{s}_abund"] = np.nan
            c["abund_eval"] = c.abund.notna()
            pp["abund_evaluable"] = round(float(c.abund_eval.mean()), 4)
            c["pair"] = f"{A}:{B}"
            per_pair[f"{A}:{B}"] = pp
            tabs.append(c)
            print(A, B, json.dumps(pp), flush=True)
    allc = pd.concat(tabs, ignore_index=True)
    allc["hard_global"] = mark_hard(allc.assign(gene_a=allc.pair + "|" + allc.gene_a, gene_b=allc.pair + "|" + allc.gene_b))
    allc.to_parquet(f"{W}/standalone_pairs{'' if RUNMODE == 'main' else '_' + RUNMODE}.parquet", index=False)
    res = dict(benchmark="orthobench_dev (34 RefOGs)", seeds=SEEDS, hard_ratio=HARD, min_nodes=MIN_NODES,
               primary_channel=PRIMARY, per_pair=per_pair)
    # coverage over ALL dev candidates and over all dev hard queries (labels only used for counting)
    cov = {}
    for ch in list(CHANNELS) + ["abund"]:
        ev = allc[f"{ch}_eval"].values
        ho = allc.hard_global.values
        cov[ch] = dict(frac_candidates_evaluable=round(float(ev.mean()), 4),
                       frac_orth_pairs_evaluable=round(float(ev[allc.label == 1].mean()), 4),
                       frac_hard_pairs_evaluable=round(float(ev[ho].mean()), 4) if ho.any() else None,
                       n_evaluable=int(ev.sum()))
        if ch != "abund":
            cov[ch]["frac_candidates_count_gt0"] = round(float((allc[f"{ch}_count"] > 0).mean()), 4)
    res["coverage"] = cov
    res["n_candidates"] = int(len(allc)); res["n_hard_pairs_global"] = int(allc.hard_global.sum())
    res["evaluation"] = {}
    for ch in CHANNELS:
        res["evaluation"][ch] = evaluate(allc, [f"{ch}_count", f"{ch}_cos"], f"{ch}_eval", ch)
        print(ch, json.dumps({k: v for k, v in res["evaluation"][ch].items() if k != "queries"})[:1500], flush=True)
    res["evaluation"]["abund"] = evaluate(allc, ["abund"], "abund_eval", "abund")
    res["wall_s"] = round(time.time() - t0, 1)
    res["manifest"] = manifest(t0)
    json.dump(res, open(f"{OUT}/standalone{'' if RUNMODE == 'main' else '_' + RUNMODE}.json", "w"), indent=1, default=str)
    print("wall", res["wall_s"])


def manifest(t0):
    import hashlib, subprocess
    def md5f(p):
        return hashlib.md5(open(p, "rb").read()).hexdigest()
    git = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "-C", WT, "status", "--porcelain", "scripts"], capture_output=True, text=True).stdout.strip())
    return dict(git_commit=git, scripts_dirty=dirty, script=os.path.abspath(__file__), script_md5=md5f(os.path.abspath(__file__)),
                config=dict(channels={k: k for k in CHANNELS}, primary=PRIMARY, seeds=SEEDS, hard=HARD, min_nodes=MIN_NODES),
                inputs={"labels": f"{MAIN}/work/S5/refog_pair_labels.parquet", "seq": f"{MAIN}/work/S7/hits.tsv",
                        "anchors": f"{RUN}/pairs.tsv", "anchors_md5": md5f(f"{RUN}/pairs.tsv"),
                        "families": f"{RUN}/stages/scores/families_final.npy",
                        "edges": {f: md5f(f"{W}/{f}") for f in sorted(os.listdir(W)) if f.startswith(("edges_", "abund_"))},
                        "string_md5s": f"{MAIN}/data/extra/l2_network/string/md5sums.txt",
                        "paxdb_md5s": f"{MAIN}/data/extra/l2_network/paxdb/md5sums.txt"},
                runtime_s=round(time.time() - t0, 1), host=os.uname().nodename)


if __name__ == "__main__":
    main()
