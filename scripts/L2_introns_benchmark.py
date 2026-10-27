#!/usr/bin/env python
"""L2-genome idea 3 step (a): do shared intron positions separate orthologs from paralogs?

Pre-registered in docs/L2_genome.md §1.4. Benchmark-only (reads dev RefOG labels).

* Pairs: all labelled OrthoBench DEV pairs (work/S5/refog_pair_labels.parquet; within-RefOG
  + deep paralogs). Introns from scripts/L2_introns_table.py.
* Alignment: Biopython PairwiseAligner, global, BLOSUM62, open -10, extend -0.5.
* Intron of a shared with b: an intron of b with the same phase at the aligned position
  (phase 1/2: the codon residue; phase 0: the residue after the junction, else the one before + 1).
  Dice = 2 x shared / (n_a + n_b) over introns inside the aligned span of each protein;
  undefined when both counts are 0.
* Shuffled control: per gene, the same number of junctions drawn without replacement over the
  3L-1 interior CDS nt positions (uniform position, random phase); seeds 11/23/47, per-gene rng
  seeded by (seed, crc32(gene)).
* Metrics: ROC-AUC (all / within-RefOG / hard pairs / per class), bootstrap CI over RefOGs
  (B=1000, seed 11); hard queries (best paralog seq >= 0.9 x best ortholog seq, seq =
  canonical DIAMOND bits, work/S7/hits.tsv) with the S5 rule (top-1 among candidates with
  seq >= 0.9 x best by intron score, ties by seq) vs sequence top-1, paired RefOG bootstrap.

Output: results/l2_genome/introns_standalone.json, work/l2_genome/intron_pairs.parquet.
"""
import hashlib, json, os, subprocess, sys, time, zlib
from multiprocessing import Pool

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W = f"{ROOT}/work/l2_genome"
OUT = f"{ROOT}/results/l2_genome"
INP = f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/Input"
SEEDS = (11, 23, 47)
HARD = 0.9
NEAR = 0.9
MAMMALS = {"Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Rattus_norvegicus", "Canis_familiaris",
           "Monodelphis_domestica"}
TETRA = MAMMALS | {"Gallus_gallus"}
TELEO = {"Danio_rerio", "Tetraodon_nigroviridis"}
INV = {"Drosophila_melanogaster", "Caenorhabditis_elegans", "Ciona_intestinalis"}

_G = {}


def pair_class(a, b):
    if a in MAMMALS and b in MAMMALS:
        return "mammal-mammal"
    if a in TETRA and b in TETRA:
        return "tetrapod-tetrapod"
    if a in TELEO and b in TELEO:
        return "teleost-teleost"
    if a in INV and b in INV:
        return "invert-invert"
    if a in INV or b in INV:
        return "vert-invert"
    return "vert-teleost"


def shuffled_introns(gene, n, L, seed):
    rng = np.random.default_rng([seed, zlib.crc32(gene.encode())])
    J = rng.choice(np.arange(1, max(2, 3 * L)), size=min(n, max(1, 3 * L - 1)), replace=False)
    J = np.sort(J)
    return (J // 3).tolist(), (J % 3).tolist()


def _init(seqs, intr):
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10
    al.extend_gap_score = -0.5
    _G.update(al=al, seqs=seqs, intr=intr)


def _key_map(i_pos, i_ph, m, L):
    """aligned key (position in the other protein) of each intron; -1 if unaligned."""
    out = []
    for r, p in zip(i_pos, i_ph):
        k = -1
        if p == 0:
            if r < L and m[r] >= 0:
                k = m[r]
            elif 0 <= r - 1 < L and m[r - 1] >= 0:
                k = m[r - 1] + 1
        elif r < L and m[r] >= 0:
            k = m[r]
        out.append(k)
    return out


def _dice(ia, ib, ma, sa, sb):
    """ia/ib = (pos, phase) lists; ma = map a->b residue; sa/sb = aligned spans (lo, hi)."""
    La = len(ma)
    ka = _key_map(ia[0], ia[1], ma, La)
    in_a = [(k, p) for r, p, k in zip(ia[0], ia[1], ka) if sa[0] <= r <= sa[1] + 1]
    in_b = {(r, p) for r, p in zip(ib[0], ib[1]) if sb[0] <= r <= sb[1] + 1}
    if not in_a and not in_b:
        return np.nan, 0, 0, 0
    shared = len({x for x in in_a if x[0] >= 0} & in_b)
    return 2 * shared / (len(in_a) + len(in_b)), shared, len(in_a), len(in_b)


def _work(chunk):
    al, seqs, intr = _G["al"], _G["seqs"], _G["intr"]
    out = []
    for a, b in chunk:
        sa_, sb_ = seqs[a], seqs[b]
        aln = al.align(sa_, sb_)[0]
        blocks_a, blocks_b = aln.aligned
        ma = np.full(len(sa_), -1, np.int64)
        for (x0, x1), (y0, y1) in zip(blocks_a, blocks_b):
            ma[x0:x1] = np.arange(y0, y1)
        span_a = (int(blocks_a[0][0]), int(blocks_a[-1][1]) - 1) if len(blocks_a) else (0, -2)
        span_b = (int(blocks_b[0][0]), int(blocks_b[-1][1]) - 1) if len(blocks_b) else (0, -2)
        ident = sum(sum(1 for i in range(x1 - x0) if sa_[x0 + i] == sb_[y0 + i])
                    for (x0, x1), (y0, y1) in zip(blocks_a, blocks_b))
        row = dict(gene_a=a, gene_b=b, aln_score=float(aln.score), aln_ident=ident / max(1, min(len(sa_), len(sb_))))
        for tag in ("real",) + tuple(f"shuf{s}" for s in SEEDS):
            ia, ib = intr[tag][a], intr[tag][b]
            d, sh, na, nb = _dice(ia, ib, ma, span_a, span_b)
            row[f"dice_{tag}"] = d
            if tag == "real":
                row.update(shared=sh, n_in_a=na, n_in_b=nb)
        out.append(row)
    return out


def canonical_seq_scores():
    h = pd.read_csv(f"{ROOT}/work/S7/hits.tsv", sep="\t", header=None, usecols=[0, 1, 2], names=["q", "t", "bits"])
    q = h.q.str.split("|", n=1, expand=True); t = h.t.str.split("|", n=1, expand=True)
    m = (q[0] == q[1]) & (t[0] == t[1]) & (q[0] != t[0])
    a, b = q[0][m].values, t[0][m].values
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    return pd.Series(h.bits[m].values).groupby([lo, hi]).max()


def read_fasta(path, want):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if name in want:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name in want:
        seqs[name] = "".join(buf)
    return seqs


def auc(y, s):
    y = np.asarray(y); s = np.asarray(s)
    return round(float(roc_auc_score(y, s)), 4) if len(y) and y.min() != y.max() else None


def boot_auc(d, col, B=1000, seed=11):
    cl = d.refog_a.values; u = np.unique(cl); g = {k: np.flatnonzero(cl == k) for k in u}
    rng = np.random.default_rng(seed); y, s = d.label.values, d[col].values; v = []
    for _ in range(B):
        ix = np.concatenate([g[k] for k in rng.choice(u, len(u), replace=True)])
        if y[ix].min() != y[ix].max():
            v.append(roc_auc_score(y[ix], s[ix]))
    return [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]


def queries(c):
    """(query gene, other species) groups with >= 1 ortholog and >= 1 paralog candidate."""
    rows = []
    for side, oth in (("gene_a", "sp_b"), ("gene_b", "sp_a")):
        for (q, s), d in c.groupby([side, oth]):
            o, p = d[d.label == 1], d[d.label == 0]
            if len(o) and len(p):
                rows.append((q, d.index.values, o.seq.max(), p.seq.max()))
    return rows


def decisions(c, qs, col):
    rows = []
    for q, ix, so, sp in qs:
        if not sp >= HARD * so:
            continue
        d = c.loc[ix]
        top = d.loc[d.seq.idxmax()]
        cand = d[d.seq >= NEAR * d.seq.max()]
        best = cand.sort_values([col, "seq"], ascending=False).iloc[0]
        mo, mp = d[d.label == 1][col].max(), d[d.label == 0][col].max()
        rows.append(dict(refog=d.refog_a.iloc[0], seq=int(top.label == 1), comb=int(best.label == 1),
                         dec=1 if mo > mp else (-1 if mo < mp else 0)))
    return pd.DataFrame(rows)


def boot_gain(r, B=1000, seed=11):
    cl = r.refog.values; u = np.unique(cl); g = {k: np.flatnonzero(cl == k) for k in u}
    rng = np.random.default_rng(seed); diff = (r.comb - r.seq).values; v = []
    for _ in range(B):
        ix = np.concatenate([g[k] for k in rng.choice(u, len(u), replace=True)])
        v.append(diff[ix].mean())
    return [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def main():
    t0, c0 = time.time(), os.times()
    os.makedirs(OUT, exist_ok=True)
    lab = pd.read_parquet(f"{ROOT}/work/S5/refog_pair_labels.parquet")
    lo = np.where(lab.gene_a < lab.gene_b, lab.gene_a, lab.gene_b)
    hi = np.where(lab.gene_a < lab.gene_b, lab.gene_b, lab.gene_a)
    lab = lab.assign(lo=lo, hi=hi).drop_duplicates(["lo", "hi"]).reset_index(drop=True)
    it = pd.read_parquet(f"{W}/introns.parquet").set_index("protein")
    genes = set(lab.gene_a) | set(lab.gene_b)
    seqs = {}
    for f in os.listdir(INP):
        seqs.update(read_fasta(f"{INP}/{f}", genes))
    seqs = {k: v.rstrip("*") for k, v in seqs.items()}
    intr = {"real": {g: (list(r.pos), list(r.phase)) for g, r in it.iterrows()}}
    for s in SEEDS:
        intr[f"shuf{s}"] = {g: shuffled_introns(g, int(r.n_introns), int(r.prot_len), s) if r.n_introns else ([], [])
                            for g, r in it.iterrows()}
    pairs = [(a, b) for a, b in zip(lab.gene_a, lab.gene_b) if a in intr["real"] and b in intr["real"]]
    cache = f"{W}/intron_pairs.parquet"
    if os.path.exists(cache):
        sc = pd.read_parquet(cache)
    else:
        chunks = [pairs[i:i + 200] for i in range(0, len(pairs), 200)]
        ta = time.time()
        with Pool(3, initializer=_init, initargs=(seqs, intr)) as pool:
            res = pool.map(_work, chunks, chunksize=1)
        sc = pd.DataFrame([r for x in res for r in x])
        sc.attrs["align_wall_s"] = time.time() - ta
        sc.to_parquet(cache, index=False)
        print(f"aligned {len(sc)} pairs in {time.time() - ta:.0f}s", flush=True)
    c = lab.merge(sc, on=["gene_a", "gene_b"], how="left")
    canon = canonical_seq_scores()
    c["seq"] = canon.reindex(pd.MultiIndex.from_arrays([c.lo.values, c.hi.values])).fillna(0.0).values
    c["cls"] = [pair_class(a, b) for a, b in zip(c.sp_a, c.sp_b)]
    tags = ["real"] + [f"shuf{s}" for s in SEEDS]
    cov = c.dice_real.notna()
    out = dict(n_pairs=int(len(c)), n_orth=int((c.label == 1).sum()), n_par=int((c.label == 0).sum()),
               coverage=dict(pairs_with_dice=round(float(cov.mean()), 4),
                             orth=round(float(cov[c.label == 1].mean()), 4), par=round(float(cov[c.label == 0].mean()), 4),
                             seq_missing=round(float((c.seq == 0).mean()), 4)),
               mean_dice=dict(orth=round(float(c[c.label == 1].dice_real.mean()), 4),
                              par_within=round(float(c[(c.label == 0) & (c.kind == "within")].dice_real.mean()), 4),
                              par_deep=round(float(c[(c.label == 0) & (c.kind == "deep")].dice_real.mean()), 4)))
    e = c[cov].reset_index(drop=True)
    # shuffled introns can all fall outside the aligned span where the real ones do not:
    # nothing shared -> Dice 0 (coverage is defined by the real introns only)
    out["shuf_nan_filled"] = {s: int(e[f"dice_shuf{s}"].isna().sum()) for s in SEEDS}
    for s in SEEDS:
        e[f"dice_shuf{s}"] = e[f"dice_shuf{s}"].fillna(0.0)
    # hard queries / pairs recomputed on the intron-covered candidates
    qs = queries(e)
    hard = np.zeros(len(e), bool)
    for q, ix, so, sp in qs:
        if sp >= HARD * so:
            hard[ix] = True
    e["hard"] = hard
    subsets = {"all": np.ones(len(e), bool), "within": (e.kind == "within").values, "hard": hard}
    out["auc"] = {}
    for name, m in subsets.items():
        d = e[m]
        r = dict(n=int(len(d)), n_orth=int((d.label == 1).sum()), seq_auc=auc(d.label, d.seq))
        for tg in tags:
            r[tg] = auc(d.label, d[f"dice_{tg}"])
        r["real_ci"] = boot_auc(d, "dice_real")
        r["shuf_max"] = max(r[f"shuf{s}"] for s in SEEDS)
        r["shuf_ci"] = {f"shuf{s}": boot_auc(d, f"dice_shuf{s}") for s in SEEDS}
        out["auc"][name] = r
    out["auc_by_class"] = {}
    for cl, d in e.groupby("cls"):
        out["auc_by_class"][cl] = dict(n=int(len(d)), n_orth=int((d.label == 1).sum()), real=auc(d.label, d.dice_real),
                                       shuf_mean=round(float(np.mean([auc(d.label, d[f"dice_shuf{s}"]) or np.nan for s in SEEDS])), 4),
                                       real_hard=auc(d.label[d.hard], d.dice_real[d.hard]),
                                       n_hard=int(d.hard.sum()))
    out["hard_queries"] = {}
    for tg in tags:
        r = decisions(e, qs, f"dice_{tg}")
        out["hard_queries"][tg] = dict(n=int(len(r)), seq_acc=round(float(r.seq.mean()), 4),
                                       comb_acc=round(float(r.comb.mean()), 4),
                                       gain=round(float((r.comb - r.seq).mean()), 4), gain_ci=boot_gain(r),
                                       intron_correct=int((r.dec == 1).sum()), intron_wrong=int((r.dec == -1).sum()),
                                       intron_abstain=int((r.dec == 0).sum()))
    a_h = out["auc"]["hard"]; g = out["hard_queries"]
    shuf_gain_max = max(g[f"shuf{s}"]["gain"] for s in SEEDS)
    pass_auc = a_h["real_ci"][0] > a_h["shuf_max"]
    pass_gain = g["real"]["gain_ci"][0] > shuf_gain_max
    out["stop_rule"] = dict(hard_auc_ci_lo=a_h["real_ci"][0], hard_auc_shuf_max=a_h["shuf_max"], pass_auc=bool(pass_auc),
                            gain_ci_lo=g["real"]["gain_ci"][0], gain_shuf_max=shuf_gain_max, pass_gain=bool(pass_gain),
                            verdict="PASS" if (pass_auc or pass_gain) else "STOP")
    c1 = os.times()
    try:
        commit = subprocess.check_output(["git", "-C", ROOT, "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = None
    out["manifest"] = dict(
        git_commit=commit, script="scripts/L2_introns_benchmark.py",
        config=dict(aligner="PairwiseAligner global BLOSUM62 open -10 extend -0.5", score="Dice shared introns (same aligned pos + phase) within aligned spans",
                    shuffle_seeds=list(SEEDS), hard_ratio=HARD, near=NEAR, bootstrap="RefOG clusters B=1000 seed 11", workers=3),
        inputs={p: md5(f"{ROOT}/{p}") for p in ("work/S5/refog_pair_labels.parquet", "work/l2_genome/introns.parquet", "work/S7/hits.tsv")},
        introns_qc=json.load(open(f"{W}/introns_qc.json")).get("_total"),
        wall_s=round(time.time() - t0, 1),
        cpu_s=round((c1.user - c0.user) + (c1.children_user - c0.children_user) + (c1.system - c0.system), 1),
        align_ms_per_pair_cpu=None)
    json.dump(out, open(f"{OUT}/introns_standalone.json", "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("coverage", "mean_dice", "auc", "hard_queries", "stop_rule")}, indent=1))


if __name__ == "__main__":
    main()
