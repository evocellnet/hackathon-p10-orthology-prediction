#!/usr/bin/env python
"""S5 standalone sub-benchmark: does gene-order context separate orthologs from paralogs?

For a species pair (A, B):
* candidates = homologous gene pairs with a truth label
  - OrthoBench: pairs inside DEV RefOGs labelled from the rooted RefOG trees
    (species overlap; scripts/S5_refog_labels.py) + cross-RefOG pairs of the
    same homology component (deep paralogs);
  - toy: same-root-family cross-species pairs, ortholog iff in truth/ortholog_pairs.tsv.
* anchors = genome-wide unambiguous reciprocal best hits (DIAMOND, both proteomes);
  anchors inside the S0 family (toy: true root family) of a or b are ignored.
* context score(a, b) = # anchors (x, y) with x in ±w genes of a, y in ±w genes of b.
* sequence score = DIAMOND bitscore (max of both directions; 0 if not reported).
Controls: shuffled gene order (positions permuted within each genome; seeds 11/23/47);
toy simulation (positive control).
Hard cases: queries a with >=1 ortholog and >=1 paralog candidate in B whose best
paralog bitscore >= 0.9 x best ortholog bitscore (sequence ranks the paralog above or near).

Outputs: results/S5_synteny/pair_benchmark[_toy].json (+ per-pair parquet in work/S5).
"""
import argparse, glob, json, os, subprocess, sys, time

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.genomes import neighbour_csr, shuffle_order  # noqa: E402
from ancestra.synteny import context_support, rbh_anchors, _Node  # noqa: E402

WS = (5, 10, 20)
SEEDS = (11, 23, 47)
HARD_RATIO = 0.9


def read_fasta_ids(path):
    out = []
    op = __import__("gzip").open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                out.append(l[1:].split()[0])
    return out


def load_hits(path, idx):
    h = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["q", "t", "bits"])
    q = h["q"].map(idx); t = h["t"].map(idx)
    ok = q.notna() & t.notna()
    return q[ok].astype(np.int64).values, t[ok].astype(np.int64).values, h["bits"][ok].values


def positions(names, sp_of, species, table_dir):
    n = len(names)
    cc = np.full(n, -1, np.int64); oo = np.zeros(n, np.int64)
    idx = {x: i for i, x in enumerate(names)}
    base = 0
    cov = {}
    for si, sp in enumerate(species):
        df = pd.read_parquet(f"{table_dir}/{sp}.parquet", columns=["protein_id", "contig", "order"])
        gi = df["protein_id"].map(idx)
        ok = gi.notna().values
        codes, uniq = pd.factorize(df["contig"][ok])
        cc[gi[ok].astype(np.int64).values] = codes + base
        oo[gi[ok].astype(np.int64).values] = df["order"].values[ok]
        base += len(uniq)
        cov[sp] = round(float(ok.sum()) / int((sp_of == si).sum()), 5)
    return cc, oo, cov


def auc(y, s):
    y = np.asarray(y); s = np.asarray(s, float)
    if len(y) == 0 or y.min() == y.max():
        return None
    return round(float(roc_auc_score(y, s)), 4)


def ap(y, s):
    y = np.asarray(y)
    if len(y) == 0 or y.min() == y.max():
        return None
    return round(float(average_precision_score(y, np.asarray(s, float))), 4)


def query_table(c):
    """Per (query gene, direction): ortholog/paralog candidate stats."""
    rows = []
    for side, oth in (("a", "b"), ("b", "a")):
        g = c.groupby(side)
        for q, d in g:
            o = d[d.label == 1]; p = d[d.label == 0]
            if len(o) == 0 or len(p) == 0:
                continue
            rows.append((side, q, d.index.values, len(o), len(p), o.seq.max(), p.seq.max()))
    return rows


def evaluate(c, hard_mask, syn_cols):
    """AUC / AP of each score column, on all and hard pairs; query-level decisions."""
    out = {}
    for name in ["seq"] + syn_cols:
        out[name] = dict(auc_all=auc(c.label, c[name]), ap_all=ap(c.label, c[name]),
                         auc_hard=auc(c.label[hard_mask], c[name][hard_mask]),
                         ap_hard=ap(c.label[hard_mask], c[name][hard_mask]))
    return out


NEARS = (0.9, 0.75, 0.5)
MARGINS = (1, 2, 3)


def query_decisions(c, queries, syn_col, near=HARD_RATIO, return_rows=False):
    """Hard queries: top-1 by sequence / by synteny (abstain on ties) / combined.

    combined@r: among candidates with seq >= r * best seq of the query, pick the
    highest synteny score, ties by seq (r = 0.9, 0.75, 0.5).
    """
    res = dict(n=0, seq_correct=0, syn_correct=0, syn_wrong=0, syn_abstain=0,
               n_tandem_competitor=0, syn_abstain_tandem=0, seq_correct_when_syn_decides=0)
    for r in NEARS:
        res[f"comb{r}_correct"] = 0
    rows = []
    for side, q, ix, no, npar, so, sp, tandem in queries:
        d = c.loc[ix]
        if not (sp >= near * so):
            continue
        res["n"] += 1
        top = d.loc[d.seq.idxmax()]
        sc = int(top.label == 1)
        res["seq_correct"] += sc
        mo = d[d.label == 1][syn_col].max(); mp = d[d.label == 0][syn_col].max()
        res["n_tandem_competitor"] += int(tandem)
        dec = 1 if mo > mp else (-1 if mo < mp else 0)
        if dec == 1:
            res["syn_correct"] += 1
        elif dec == -1:
            res["syn_wrong"] += 1
        else:
            res["syn_abstain"] += 1
            res["syn_abstain_tandem"] += int(tandem)
        if dec != 0:
            res["seq_correct_when_syn_decides"] += sc
        row = dict(refog=d.get("refog_a", pd.Series([None])).iloc[0], seq=sc, syn=dec)
        for r in NEARS:
            cand = d[d.seq >= r * d.seq.max()]
            best = cand.sort_values([syn_col, "seq"], ascending=False).iloc[0]
            res[f"comb{r}_correct"] += int(best.label == 1)
            row[f"comb{r}"] = int(best.label == 1)
            # override rule: switch from the sequence top hit only with a support margin >= m
            for m in MARGINS:
                pick = best if best[syn_col] - top[syn_col] >= m else top
                row[f"ovr{m}@{r}"] = int(pick.label == 1)
        rows.append(row)
    if res["n"]:
        for k in list(res):
            if k.endswith("_correct") or k in ("syn_wrong", "syn_abstain"):
                res[k + "_frac"] = round(res[k] / res["n"], 4)
    return (res, rows) if return_rows else res


def run_pair(A, B, names_a, names_b, hits_paths, cand, table_dir, fam_of, amb, log):
    names = list(names_a) + list(names_b)
    sp_of = np.r_[np.zeros(len(names_a), np.int64), np.ones(len(names_b), np.int64)]
    idx = {x: i for i, x in enumerate(names)}
    q, t, bits = [], [], []
    for p in hits_paths:
        a, b, s = load_hits(p, idx)
        q.append(a); t.append(b); bits.append(s)
    q, t, bits = np.concatenate(q), np.concatenate(t), np.concatenate(bits)
    cross = sp_of[q] != sp_of[t]
    q, t, bits = q[cross], t[cross], bits[cross]
    aa, ab = rbh_anchors(q, t, bits, amb=amb)
    aa1, ab1 = rbh_anchors(q, t, bits, amb=1.0)
    # sequence score per unordered pair
    lo, hi = np.minimum(q, t), np.maximum(q, t)
    seq = pd.Series(bits).groupby(lo * len(names) + hi).max()
    cc, oo, cov = positions(names, sp_of, [A, B], table_dir)
    fam = np.array([fam_of.get(x, -1 - i) for i, x in enumerate(names)], dtype=object)
    fam = pd.factorize(pd.Series(fam).astype(str))[0].astype(np.int64)
    # candidates oriented a in A, b in B
    c = cand.copy()
    sw = c.sp_a != A
    c.loc[sw, ["gene_a", "gene_b"]] = c.loc[sw, ["gene_b", "gene_a"]].values
    c["a"] = c.gene_a.map(idx); c["b"] = c.gene_b.map(idx)
    n0 = len(c)
    c = c[c.a.notna() & c.b.notna()].copy()
    c["a"] = c.a.astype(np.int64); c["b"] = c.b.astype(np.int64)
    c = c.drop_duplicates(["a", "b"]).reset_index(drop=True)
    key = np.minimum(c.a, c.b) * len(names) + np.maximum(c.a, c.b)
    c["seq"] = seq.reindex(key.values).fillna(0.0).values
    c["has_coords"] = (cc[c.a] >= 0) & (cc[c.b] >= 0)
    syn_cols = []
    for w in WS:
        ptr, ix = neighbour_csr(cc, oo, w)
        c[f"syn_w{w}"] = context_support((c.a.values, c.b.values), (aa, ab), _Node(ptr, ix), fam)
        syn_cols.append(f"syn_w{w}")
        for sd in SEEDS:
            c2, o2 = shuffle_order(cc, oo, sd)
            ptr, ix = neighbour_csr(c2, o2, w)
            c[f"shuf{sd}_w{w}"] = context_support((c.a.values, c.b.values), (aa, ab), _Node(ptr, ix), fam)
            syn_cols.append(f"shuf{sd}_w{w}")
    # tandem competitor flag: the query's best ortholog and best paralog (by seq) lie within 5 genes
    def tandem_flag(side, other):
        flags = {}
        for g, d in c.groupby(side):
            o = d[d.label == 1]; p = d[d.label == 0]
            if len(o) == 0 or len(p) == 0:
                continue
            xo = o.loc[o.seq.idxmax(), other]; xp = p.loc[p.seq.idxmax(), other]
            flags[g] = bool(cc[xo] >= 0 and cc[xo] == cc[xp] and abs(oo[xo] - oo[xp]) <= 5)
        return flags
    queries = query_table(c)
    hard = np.zeros(len(c), bool)
    hq = [x for x in queries if x[6] >= HARD_RATIO * x[5]]
    for x in hq:
        hard[x[2]] = True
    c["hard"] = hard
    tf_a, tf_b = tandem_flag("a", "b"), tandem_flag("b", "a")
    queries = [x + ((tf_a if x[0] == "a" else tf_b).get(x[1], False),) for x in queries]
    res = dict(pair=f"{A}__{B}", n_candidates_input=int(n0), n_candidates=int(len(c)),
               n_orth=int((c.label == 1).sum()), n_par=int((c.label == 0).sum()),
               n_par_deep=int(((c.label == 0) & (c.get("kind", "within") == "deep")).sum()) if "kind" in c else None,
               n_hard_pairs=int(hard.sum()), n_hard_orth=int((hard & (c.label == 1)).sum()),
               n_queries_both=len(queries), n_queries_hard=len(hq),
               n_anchors=int(len(aa)), n_anchors_plain_rbh=int(len(aa1)),
               coverage=cov, frac_candidates_with_coords=round(float(c.has_coords.mean()), 4),
               frac_cand_seq_missing=round(float((c.seq == 0).mean()), 4))
    res["metrics"] = evaluate(c, hard, syn_cols)
    res["queries"] = {}
    for w in WS:
        for col in [f"syn_w{w}"] + [f"shuf{sd}_w{w}" for sd in SEEDS]:
            res["queries"][col], qrows = query_decisions(c, queries, col, return_rows=True)
            for x in qrows:
                x.update(pair=f"{A}__{B}", col=col)
            res.setdefault("_qrows", []).extend(qrows)
            # all queries with >= 1 ortholog and >= 1 paralog candidate (label-free selection)
            res.setdefault("queries_all", {})[col], qrows = query_decisions(c, queries, col, near=0.0,
                                                                              return_rows=True)
            for x in qrows:
                x.update(pair=f"{A}__{B}", col=col)
            res.setdefault("_qrows_all", []).extend(qrows)
    # summary: real vs shuffled mean
    summ = {}
    for w in WS:
        m = res["metrics"]
        sh = [m[f"shuf{sd}_w{w}"] for sd in SEEDS]
        summ[f"w{w}"] = dict(
            auc_all=m[f"syn_w{w}"]["auc_all"], auc_hard=m[f"syn_w{w}"]["auc_hard"],
            ap_all=m[f"syn_w{w}"]["ap_all"], ap_hard=m[f"syn_w{w}"]["ap_hard"],
            shuf_auc_all=_mean([s["auc_all"] for s in sh]), shuf_auc_hard=_mean([s["auc_hard"] for s in sh]),
            frac_orth_supported=round(float((c[c.label == 1][f"syn_w{w}"] > 0).mean()), 4),
            frac_par_supported=round(float((c[c.label == 0][f"syn_w{w}"] > 0).mean()), 4))
    summ["seq"] = res["metrics"]["seq"]
    res["summary"] = summ
    log(f"{A}-{B}: cand={len(c)} orth={res['n_orth']} par={res['n_par']} hard={hard.sum()} anchors={len(aa)} "
        + " ".join(f"w{w}:AUC={summ[f'w{w}']['auc_all']}/{summ[f'w{w}']['auc_hard']} shuf={summ[f'w{w}']['shuf_auc_all']}" for w in WS)
        + f" seqAUC={summ['seq']['auc_all']}/{summ['seq']['auc_hard']}")
    return res, c


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(float(np.mean(xs)), 4) if xs else None


def bootstrap_auc(c, col, mask, B=1000, seed=11):
    """95% CI of AUC(col) by resampling clusters (RefOG / family) with replacement."""
    cl = c["cluster"].values
    u = np.unique(cl)
    groups = {k: np.flatnonzero(cl == k) for k in u}
    rng = np.random.default_rng(seed)
    vals = []
    y_all, s_all = c.label.values, c[col].values
    for _ in range(B):
        pick = rng.choice(u, len(u), replace=True)
        ix = np.concatenate([groups[k] for k in pick])
        ix = ix[mask[ix]]
        y = y_all[ix]
        if len(y) and y.min() != y.max():
            vals.append(roc_auc_score(y, s_all[ix]))
    return [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)] if vals else None


def pooled_queries(results, B=1000, seed=11, key="_qrows"):
    """Pool hard-query decisions over pairs; bootstrap CI of (combined - seq) accuracy over RefOGs."""
    rows = pd.DataFrame([r for res in results for r in res.get(key, [])])
    out = {}
    if rows.empty:
        return out
    rng = np.random.default_rng(seed)
    for col, d in rows.groupby("col"):
        e = dict(n=len(d), seq_acc=round(d.seq.mean(), 4), syn_correct=int((d.syn == 1).sum()),
                 syn_wrong=int((d.syn == -1).sum()), syn_abstain=int((d.syn == 0).sum()))
        e["syn_precision_when_decides"] = round(e["syn_correct"] / max(1, e["syn_correct"] + e["syn_wrong"]), 4)
        dec = d[d.syn != 0]
        e["seq_acc_when_syn_decides"] = round(dec.seq.mean(), 4) if len(dec) else None
        for r in NEARS:
            for m in MARGINS:
                e[f"ovr{m}@{r}_acc"] = round(d[f"ovr{m}@{r}"].mean(), 4)
                e[f"ovr{m}@{r}_minus_seq"] = round(float(d[f"ovr{m}@{r}"].mean() - d.seq.mean()), 4)
            e[f"comb{r}_acc"] = round(d[f"comb{r}"].mean(), 4)
            cl = d.refog.fillna("na").values
            u = np.unique(cl)
            g = {k: np.flatnonzero(cl == k) for k in u}
            diffs = []
            for _ in range(B):
                ix = np.concatenate([g[k] for k in rng.choice(u, len(u), replace=True)])
                diffs.append(d[f"comb{r}"].values[ix].mean() - d.seq.values[ix].mean())
            e[f"comb{r}_minus_seq"] = round(float(d[f"comb{r}"].mean() - d.seq.mean()), 4)
            e[f"comb{r}_minus_seq_ci95"] = [round(float(np.percentile(diffs, 2.5)), 4),
                                            round(float(np.percentile(diffs, 97.5)), 4)]
        out[col] = e
    return out


def pooled(cs, ws=WS):
    """Pooled metrics over all pairs (concatenated candidates)."""
    c = pd.concat(cs, ignore_index=True)
    c["cluster"] = c["refog_a"].astype(str)
    out = dict(n=len(c), n_orth=int((c.label == 1).sum()), n_hard=int(c.hard.sum()))
    h = c.hard.values
    out["seq"] = dict(auc_all=auc(c.label, c.seq), auc_hard=auc(c.label[h], c.seq[h]))
    for w in ws:
        col = f"syn_w{w}"
        out[col] = dict(auc_all=auc(c.label, c[col]), auc_hard=auc(c.label[h], c[col][h]),
                        ap_all=ap(c.label, c[col]), ap_hard=ap(c.label[h], c[col][h]),
                        shuf_auc_all=_mean([auc(c.label, c[f"shuf{s}_w{w}"]) for s in SEEDS]),
                        shuf_auc_hard=_mean([auc(c.label[h], c[f"shuf{s}_w{w}"][h]) for s in SEEDS]),
                        shuf_ap_hard=_mean([ap(c.label[h], c[f"shuf{s}_w{w}"][h]) for s in SEEDS]),
                        auc_all_ci95=bootstrap_auc(c, col, np.ones(len(c), bool)),
                        auc_hard_ci95=bootstrap_auc(c, col, h))
    out["seq"]["auc_all_ci95"] = bootstrap_auc(c, "seq", np.ones(len(c), bool))
    out["seq"]["auc_hard_ci95"] = bootstrap_auc(c, "seq", h)
    return out


# ----------------------------------------------------------------------------
def main_orthobench(a, log):
    IN = f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/Input"
    lab = pd.read_parquet(f"{ROOT}/work/S5/refog_pair_labels.parquet")
    fam = pd.read_csv(a.families, sep="\t", comment="#", header=None, names=["g", "f"])
    fam_of = dict(zip(fam.g, fam.f))
    results, cs = [], []
    for pr in a.pairs.split():
        A, B = pr.split(":")
        na = read_fasta_ids(glob.glob(f"{IN}/{A}.*.pep.all.fa")[0])
        nb = read_fasta_ids(glob.glob(f"{IN}/{B}.*.pep.all.fa")[0])
        cand = lab[((lab.sp_a == A) & (lab.sp_b == B)) | ((lab.sp_a == B) & (lab.sp_b == A))]
        hp = [f"{ROOT}/work/S5/diamond/{A}__{B}.tsv", f"{ROOT}/work/S5/diamond/{B}__{A}.tsv"]
        r, c = run_pair(A, B, na, nb, hp, cand, a.tables, fam_of, a.amb, log)
        c["pair"] = pr
        results.append(r); cs.append(c)
        c.to_parquet(f"{ROOT}/work/S5/pair_scores_{A}__{B}.parquet", index=False)
    return results, cs


def main_toy(a, log):
    from ancestra.genomes import parse_gff3_proteins, assign_order, write_gene_order
    from ancestra.simulate import SimParams, simulate, write_simulation
    d = f"{ROOT}/work/S5/toy/{a.toy_name}"
    if not os.path.exists(f"{d}/truth/ortholog_pairs.tsv"):
        p = SimParams(n_species=a.toy_species, n_genes=a.toy_genes, wgd_branch="random",
                      dup_rate=0.3, n_chromosomes=4)
        write_simulation(simulate(p, seed=a.toy_seed), d)
    sp = sorted(os.path.basename(x)[:-3] for x in glob.glob(f"{d}/proteomes/*.fa"))
    tdir = f"{d}/genomes"
    for s in sp:
        if not os.path.exists(f"{tdir}/{s}.parquet"):
            df = parse_gff3_proteins(f"{d}/gff/{s}.gff3")
            df["primary"] = True; df["source"] = "toy"
            write_gene_order(assign_order(df), f"{tdir}/{s}.parquet")
    gf = pd.read_csv(f"{d}/truth/gene_families.tsv", sep="\t")
    fam_of = dict(zip(gf["protein"], gf["family"]))
    tp = pd.read_csv(f"{d}/truth/ortholog_pairs.tsv", sep="\t")
    truth = set(zip(tp.iloc[:, 0], tp.iloc[:, 1])) | set(zip(tp.iloc[:, 1], tp.iloc[:, 0]))
    # species pairs: the most distant pair, a close pair, and pairs across the WGD
    tree = open(f"{d}/species_tree.nwk").read()
    params = json.load(open(f"{d}/params.json"))
    pairs = [(sp[i], sp[j]) for i in range(len(sp)) for j in range(i + 1, len(sp))]
    rng = np.random.default_rng(a.toy_seed)
    pairs = [pairs[i] for i in rng.choice(len(pairs), min(a.toy_pairs, len(pairs)), replace=False)]
    results, cs = [], []
    os.makedirs(f"{d}/diamond", exist_ok=True)
    for A, B in pairs:
        hp = []
        for x, y in ((A, B), (B, A)):
            out = f"{d}/diamond/{x}__{y}.tsv"
            if not os.path.exists(out):
                subprocess.run(f"diamond blastp -q {d}/proteomes/{x}.fa -d {d}/proteomes/{y}.fa -o {out} --sensitive "
                               f"-k 100 --max-hsps 1 -e 1e-3 -p {a.threads} --quiet --outfmt 6 qseqid sseqid bitscore",
                               shell=True, check=True)
            hp.append(out)
        na = read_fasta_ids(f"{d}/proteomes/{A}.fa"); nb = read_fasta_ids(f"{d}/proteomes/{B}.fa")
        fa = pd.Series(fam_of)
        ga = [g for g in na if g in fam_of]; gb = [g for g in nb if g in fam_of]
        da = pd.DataFrame({"gene_a": ga, "f": [fam_of[g] for g in ga]})
        db = pd.DataFrame({"gene_b": gb, "f": [fam_of[g] for g in gb]})
        cand = da.merge(db, on="f")
        cand["sp_a"], cand["sp_b"] = A, B
        cand["label"] = [int((x, y) in truth) for x, y in zip(cand.gene_a, cand.gene_b)]
        cand["kind"] = "within"
        cand["refog_a"] = cand["f"]
        r, c = run_pair(A, B, na, nb, hp, cand, tdir, fam_of, a.amb, log)
        c["pair"] = f"{A}:{B}"
        results.append(r); cs.append(c)
    return results, cs, dict(sim=d, params=params.get("params", params), tree=tree.strip())


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--toy", action="store_true")
    ap_.add_argument("--pairs", default="Homo_sapiens:Mus_musculus Homo_sapiens:Gallus_gallus Homo_sapiens:Danio_rerio "
                                       "Danio_rerio:Tetraodon_nigroviridis Homo_sapiens:Drosophila_melanogaster "
                                       "Homo_sapiens:Tetraodon_nigroviridis")
    ap_.add_argument("--tables", default=f"{ROOT}/data/genomes/orthobench")
    ap_.add_argument("--families", default=f"{ROOT}/work/S0/2sA_dc30_cc_sh30_RA/families.tsv")
    ap_.add_argument("--amb", type=float, default=0.9)
    ap_.add_argument("--out", default=None)
    ap_.add_argument("--threads", type=int, default=4)
    ap_.add_argument("--toy-name", default="wgd8_s11")
    ap_.add_argument("--toy-species", type=int, default=8)
    ap_.add_argument("--toy-genes", type=int, default=1500)
    ap_.add_argument("--toy-seed", type=int, default=11)
    ap_.add_argument("--toy-pairs", type=int, default=6)
    a = ap_.parse_args()
    t0 = time.time()
    outdir = f"{ROOT}/results/S5_synteny"
    os.makedirs(outdir, exist_ok=True)
    msgs = []

    def log(m):
        print(m, flush=True); msgs.append(m)
    extra = {}
    if a.toy:
        results, cs, extra = main_toy(a, log)
        out = a.out or f"{outdir}/pair_benchmark_toy_{a.toy_name}.json"
    else:
        results, cs = main_orthobench(a, log)
        out = a.out or f"{outdir}/pair_benchmark_orthobench_dev.json"
    pool = pooled(cs)
    pq = pooled_queries(results)
    pqa = pooled_queries(results, key="_qrows_all")
    for r in results:
        r.pop("_qrows", None)
        r.pop("_qrows_all", None)
    log("POOLED_ALLQ " + json.dumps({k: v for k, v in pqa.items() if not k.startswith("shuf") or k.startswith("shuf11")}))
    log("POOLED " + json.dumps(pool))
    log("POOLED_QUERIES " + json.dumps({k: v for k, v in pq.items() if not k.startswith("shuf") or k.startswith("shuf11")}))
    json.dump(dict(benchmark="toy" if a.toy else "orthobench_dev", amb=a.amb, windows=WS, shuffle_seeds=SEEDS,
                   hard_ratio=HARD_RATIO, pooled=pool, pooled_hard_queries=pq, pooled_all_queries=pqa, pairs=results, log=msgs, extra=extra,
                   wall_s=round(time.time() - t0, 1)), open(out, "w"), indent=1, default=str)
    print("wrote", out)


if __name__ == "__main__":
    main()
