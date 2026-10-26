"""S6 -- structure (3Di) evidence (DESIGN.md §1.3, §3 row S6).

3Di strings are read from a FASTA (``>protein_id|...`` header, one 3Di string of
the same length as the protein; e.g. exact-sequence AFDB 3Di). Proteins without a
3Di string fall back to sequence everywhere (they are simply absent from the
Foldseek DB).

Building blocks
---------------
``load_3di``             protein id -> 3Di string (length-checked against the AA sequence).
``shuffled_3di``         pre-registered negative control: 3Di strings permuted among proteins of
                         similar length (length-sorted bins), resampled to the recipient length.
``build_db``             Foldseek DB from AA + 3Di strings without structures
                         (``base:createdb`` of the AA FASTA and of the 3Di FASTA with
                         ``--shuffle 0`` so keys agree; the second becomes ``<db>_ss``).
                         Headers are integer gene ids.
``search``               ``foldseek search`` (3Di+AA, ``--alignment-type 2``) -> hit arrays.
``align_pairs``          3Di+AA alignment of an explicit pair list (``tsv2db`` prefilter DB +
                         ``structurealign``), self pairs give self-scores.
``link_families``        S0 rescue: union of families linked by 3Di+AA hits (e-value, coverage
                         of the shorter protein, size cap, optional small-family-only rule).
``combine_scores``       structure-aware ScoreTable for the merge: per pair, the normalised
                         3Di+AA score is calibrated to the sequence scale (linear fit in
                         -ln S on pairs that have both) and averaged with the sequence score in
                         distance space; pairs without 3Di keep the sequence score; an optional
                         node filter applies structure only to pairs whose species LCA is a deep
                         node. The merge reads only ``bits``/``self_bits`` via ``norm_max``, so the
                         combined S is written as bits' = S * max(self_a, self_b).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import time

import numpy as np

from .io import GeneTable, read_fasta
from .scoring import ScoreTable

__all__ = ["load_3di", "shuffled_3di", "build_db", "search", "align_pairs", "link_families",
           "combine_scores", "pair_lca_nodes", "FOLDSEEK"]

FOLDSEEK = os.environ.get("ANCESTRA_FOLDSEEK", "foldseek")
THREEDI = set("ACDEFGHIKLMNPQRSTVWY")


def _run(cmd, log=None):
    cmd = [str(c) for c in cmd]
    if log is not None:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
    r = subprocess.run(cmd, stdout=log if log is not None else subprocess.DEVNULL,
                       stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(cmd)}")


# ----------------------------------------------------------------------------
# 3Di input
# ----------------------------------------------------------------------------

def load_3di(path: str, gt: GeneTable | None = None) -> dict:
    """protein id -> 3Di string. The id is the first ``|`` field of the header. With ``gt``,
    only proteins of the gene table are kept and the 3Di length must equal the protein
    length (entries with a different length are dropped and counted)."""
    out = {}
    lens = None
    if gt is not None:
        lens = {n: int(l) for n, l in zip(gt.names, gt.lengths)}
    bad = 0
    for hid, s in read_fasta(path):
        pid = hid.split("|")[0]
        s = s.strip().upper()
        if lens is not None:
            if pid not in lens:
                continue
            if len(s) != lens[pid]:
                bad += 1
                continue
        out[pid] = s
    load_3di.n_length_mismatch = bad
    return out


def shuffled_3di(ss: dict, seed: int, bin_size: int = 20) -> dict:
    """Negative control: each protein gets the 3Di string of another protein of similar
    length. Proteins are sorted by 3Di length, cut into consecutive bins of ``bin_size``
    and permuted within each bin (``np.random.default_rng(seed)``; a derangement, so no
    protein keeps its own string; a last bin of one protein keeps it). The donor string is resampled to the
    recipient's length by nearest index (lengths differ by a few residues within a bin)."""
    rng = np.random.default_rng(seed)
    ids = sorted(ss, key=lambda k: (len(ss[k]), k))
    out = {}
    for s in range(0, len(ids), bin_size):
        b = ids[s:s + bin_size]
        perm = rng.permutation(len(b))
        while len(b) > 1 and (perm == np.arange(len(b))).any():
            perm = rng.permutation(len(b))
        for i, j in zip(range(len(b)), perm):
            rec, don = b[i], ss[b[j]]
            L, Ld = len(ss[rec]), len(don)
            if Ld == L:
                out[rec] = don
            else:
                idx = np.minimum((np.arange(L) * Ld) // L, Ld - 1)
                out[rec] = "".join(don[k] for k in idx)
    return out


# ----------------------------------------------------------------------------
# Foldseek DB / search / pair alignment
# ----------------------------------------------------------------------------

def build_db(gt: GeneTable, seqs: list, ss: dict, out_prefix: str, log=None,
             genes: np.ndarray | None = None) -> np.ndarray:
    """Foldseek DB (AA + 3Di, no C-alpha) for the genes that have a 3Di string.
    ``seqs``: AA sequence per gene id. Returns the gene ids in DB-key order (key k = genes[k]);
    headers are the integer gene ids."""
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    cand = np.arange(gt.n) if genes is None else np.asarray(genes)
    keep = [int(g) for g in cand if gt.names[int(g)] in ss]
    aa_fa, ss_fa = out_prefix + ".aa.fasta", out_prefix + ".3di.fasta"
    with open(aa_fa, "w") as fa, open(ss_fa, "w") as fs:
        for g in keep:
            a, s = seqs[g], ss[gt.names[g]]
            if len(a) != len(s):
                raise ValueError(f"3Di length {len(s)} != AA length {len(a)} for {gt.names[g]}")
            fa.write(f">{g}\n{a}\n")
            fs.write(f">{g}\n{s}\n")
    tmp_ss = out_prefix + "_sstmp"
    _run([FOLDSEEK, "base:createdb", aa_fa, out_prefix, "--shuffle", "0", "-v", "1"], log)
    _run([FOLDSEEK, "base:createdb", ss_fa, tmp_ss, "--shuffle", "0", "-v", "1"], log)
    for suf in ("", ".index", ".dbtype"):
        os.replace(tmp_ss + suf, out_prefix + "_ss" + suf)
    for f in os.listdir(os.path.dirname(os.path.abspath(out_prefix))):
        p = os.path.join(os.path.dirname(os.path.abspath(out_prefix)), f)
        if f.startswith(os.path.basename(tmp_ss)):
            os.remove(p)
    os.remove(aa_fa)
    os.remove(ss_fa)
    return np.asarray(keep, dtype=np.int64)


def _read_hits(path):
    if os.path.getsize(path) == 0:
        z = np.zeros(0)
        return dict(q=z.astype(np.int64), t=z.astype(np.int64), bits=z, evalue=z, qcov=z, tcov=z,
                    pident=z)
    a = np.loadtxt(path, dtype=np.float64, ndmin=2)
    return dict(q=a[:, 0].astype(np.int64), t=a[:, 1].astype(np.int64), bits=a[:, 2],
                evalue=a[:, 3], qcov=a[:, 4], tcov=a[:, 5], pident=a[:, 6])


FMT = "query,target,bits,evalue,qcov,tcov,fident"


def search(qdb: str, tdb: str, wd: str, threads: int = 10, evalue: float = 1e-2,
           max_seqs: int = 1000, sens: float = 9.5, extra=(), log=None) -> dict:
    """``foldseek search`` (3Di+AA) of qdb vs tdb -> dict of arrays (gene ids)."""
    os.makedirs(wd, exist_ok=True)
    res, tmp = os.path.join(wd, "res"), os.path.join(wd, "tmp")
    _run([FOLDSEEK, "search", qdb, tdb, res, tmp, "--alignment-type", "2", "-e", evalue,
          "--max-seqs", max_seqs, "-s", sens, "--threads", threads, "-v", "1", *extra], log)
    out = os.path.join(wd, "hits.m8")
    _run([FOLDSEEK, "convertalis", qdb, tdb, res, out, "--format-output", FMT,
          "--threads", threads, "-v", "1"], log)
    h = _read_hits(out)
    shutil.rmtree(tmp, ignore_errors=True)
    return h


def align_pairs(db: str, keys_of_gene: dict, q: np.ndarray, t: np.ndarray, wd: str,
                threads: int = 10, log=None) -> dict:
    """3Di+AA local alignment (``structurealign --alignment-type 2``) of an explicit gene-pair
    list, both genes in ``db`` (``keys_of_gene``: gene id -> DB key). E-value threshold off,
    so every pair is reported unless the aligner returns nothing."""
    os.makedirs(wd, exist_ok=True)
    qk = np.fromiter((keys_of_gene[int(x)] for x in q), np.int64, len(q))
    tk = np.fromiter((keys_of_gene[int(x)] for x in t), np.int64, len(t))
    o = np.lexsort((tk, qk))
    qk, tk = qk[o], tk[o]
    pref_tsv = os.path.join(wd, "pairs.tsv")
    with open(pref_tsv, "w") as fh:
        blk = 2_000_000
        for s in range(0, len(qk), blk):
            fh.write("".join(f"{a}\t{b}\t0\t0\n" for a, b in zip(qk[s:s + blk].tolist(),
                                                               tk[s:s + blk].tolist())))
    pref = os.path.join(wd, "pref")
    _run([FOLDSEEK, "base:tsv2db", pref_tsv, pref, "--output-dbtype", "7", "-v", "1"], log)
    os.remove(pref_tsv)
    aln = os.path.join(wd, "aln")
    _run([FOLDSEEK, "structurealign", db, db, pref, aln, "--alignment-type", "2", "-e", "inf",
          "--tmscore-threshold", "0", "--lddt-threshold", "0", "--threads", threads, "-a",
          "-v", "1"], log)
    out = os.path.join(wd, "aln.m8")
    _run([FOLDSEEK, "convertalis", db, db, aln, out, "--format-output", FMT,
          "--threads", threads, "-v", "1"], log)
    return _read_hits(out)


# ----------------------------------------------------------------------------
# S0 rescue
# ----------------------------------------------------------------------------

def link_families(fam: np.ndarray, hits: dict, lengths: np.ndarray, evalue: float = 1e-3,
                  cov_shorter: float = 0.3, max_size: int = 7500, small_max: int | None = None,
                  cov_min: float = 0.0, mode: str = "union") -> tuple[np.ndarray, dict]:
    """Join S0 families linked by structure hits.

    A hit (q, t) links fam[q] and fam[t] if e <= ``evalue`` and the alignment covers >=
    ``cov_shorter`` of the shorter protein (qcov/tcov are fractions of q/t). With
    ``small_max``, only links where at least one side has <= small_max genes (singletons
    and small fragments of a split family) are used. Links are applied in increasing
    e-value; a join that would exceed ``max_size`` genes is skipped (no giant families).
    ``cov_min``: additionally both qcov and tcov >= cov_min (bidirectional guard).
    ``mode="attach"`` (needs ``small_max``): each small family (<= small_max genes) is joined
    only to the family of its single best link (lowest e-value); families above small_max
    have no outgoing link. This is a functional graph whose sinks are the large families,
    so two large families are never joined (no chaining through small units).
    Returns new labels (old family id of the representative) and stats."""
    fam = np.asarray(fam).copy()
    q, t, e = hits["q"], hits["t"], hits["evalue"]
    L = lengths.astype(np.float64)
    alen_q = hits["qcov"] * L[q]
    alen_t = hits["tcov"] * L[t]
    cov_sh = np.maximum(alen_q, alen_t) / np.minimum(L[q], L[t])   # aligned span / shorter length
    fq, ft = fam[q], fam[t]
    ok = (e <= evalue) & (cov_sh >= cov_shorter) & (fq != ft)
    if cov_min > 0:
        ok &= (hits["qcov"] >= cov_min) & (hits["tcov"] >= cov_min)
    uniq, inv, size0 = np.unique(fam, return_inverse=True, return_counts=True)
    a, b, ee = inv[q[ok]], inv[t[ok]], e[ok]
    if mode == "attach":
        if small_max is None:
            raise ValueError("attach mode needs small_max")
        src = np.concatenate([a, b]); dst = np.concatenate([b, a]); e2 = np.concatenate([ee, ee])
        m = size0[src] <= small_max
        src, dst, e2 = src[m], dst[m], e2[m]
        o = np.lexsort((-size0[dst], e2, src))          # per source: lowest e, then larger partner
        src, dst, e2 = src[o], dst[o], e2[o]
        first = np.ones(len(src), bool)
        first[1:] = src[1:] != src[:-1]
        a, b, ee = src[first], dst[first], e2[first]
    elif small_max is not None:
        sm = (size0[a] <= small_max) | (size0[b] <= small_max)
        a, b, ee = a[sm], b[sm], ee[sm]
    parent = np.arange(len(uniq))
    size = size0.astype(np.int64).copy()

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    n_join = n_capped = 0
    for i in np.argsort(ee, kind="stable"):
        ra, rb = find(a[i]), find(b[i])
        if ra == rb:
            continue
        if size[ra] + size[rb] > max_size:
            n_capped += 1
            continue
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]
        n_join += 1
    root = np.array([find(x) for x in range(len(uniq))])
    new = uniq[root][inv]
    return new, dict(n_links=int(ok.sum()), n_links_used=int(len(a)), n_joins=n_join,
                     n_capped=n_capped, n_families_before=int(len(uniq)),
                     n_families_after=int(len(np.unique(new))))


# ----------------------------------------------------------------------------
# structure-aware score table for the merge
# ----------------------------------------------------------------------------

def pair_lca_nodes(gt: GeneTable, tree, q: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, list]:
    """Name index of the species-tree LCA of each gene pair, and the node-name list."""
    names, sidx = [], {s: i for i, s in enumerate(gt.species)}
    nodes = list(tree.postorder())
    nid = {nd.name: i for i, nd in enumerate(nodes)}
    S = len(gt.species)
    lca = np.full((S, S), -1, dtype=np.int32)
    clade = {}
    for nd in nodes:
        names.append(nd.name)
        if not nd.children:
            clade[nd.name] = [sidx[nd.name]]
            lca[sidx[nd.name], sidx[nd.name]] = nid[nd.name]
            continue
        cl = [clade[c.name] for c in nd.children]
        for i in range(len(cl)):
            for j in range(i + 1, len(cl)):
                lca[np.ix_(cl[i], cl[j])] = nid[nd.name]
                lca[np.ix_(cl[j], cl[i])] = nid[nd.name]
        clade[nd.name] = sum(cl, [])
    return lca[gt.species_of[q], gt.species_of[t]], names


def _pair_key(q, t, n):
    a, b = np.minimum(q, t).astype(np.int64), np.maximum(q, t).astype(np.int64)
    return a * n + b


def combine_scores(gt: GeneTable, st: ScoreTable, fs: dict, fs_self: np.ndarray,
                   mode: str = "avg", node_mask: np.ndarray | None = None,
                   fs_genes: np.ndarray | None = None, add_missing: bool = True,
                   fit: tuple | None = None, miss_evalue: float = 1e-3) -> tuple[ScoreTable, dict]:
    """Structure-aware ScoreTable.

    ``fs``: 3Di+AA pair hits (gene ids, bits); ``fs_self``: 3Di+AA self bits per gene (NaN =
    no 3Di). S_seq = bits/max(self) (the merge's norm_max); S_3di likewise from Foldseek.
    Calibration (``fit`` or fitted here): d_3di' = a + b * d_3di in d = -ln S, least squares on
    pairs with both scores. ``mode``: ``avg`` -> d = (d_seq + d_3di')/2; ``struct`` -> d =
    d_3di'; ``max`` -> S = max(S_seq, S_3di'). Pairs without a 3Di score keep S_seq; pairs with
    only a 3Di score are added when ``add_missing`` (they were unscored by sequence).
    ``node_mask``: per pair of the *output* table, True where structure may be used
    (e.g. deep LCA); computed by the caller via :func:`pair_lca_nodes` -- pass a callable
    ``f(q, t) -> bool array`` instead of an array. Self-bits are unchanged, so the merge's
    norm_max of the output equals the combined S."""
    n = gt.n
    sb = st.self_bits_filled(gt.lengths)
    q0, t0 = st.q.astype(np.int64), st.t.astype(np.int64)
    s_seq = st.bits.astype(np.float64) / np.maximum(sb[q0], sb[t0])
    fq, ft = fs["q"].astype(np.int64), fs["t"].astype(np.int64)
    m = fq != ft
    fq, ft, fb = fq[m], ft[m], fs["bits"][m].astype(np.float64)
    fe = fs["evalue"][m].astype(np.float64) if "evalue" in fs else np.zeros(len(fq))
    den = np.maximum(fs_self[fq], fs_self[ft])
    s3 = fb / den
    ok = np.isfinite(s3) & (s3 > 0)
    fq, ft, s3, fe = fq[ok], ft[ok], np.minimum(s3[ok], 1.0), fe[ok]
    k3 = _pair_key(fq, ft, n)
    o = np.argsort(k3)
    k3, s3, fq, ft, fe = k3[o], s3[o], fq[o], ft[o], fe[o]
    # duplicates (both directions): keep max
    u, first = np.unique(k3, return_index=True)
    s3m = np.maximum.reduceat(s3, first) if len(s3) else s3
    fem = np.minimum.reduceat(fe, first) if len(fe) else fe
    k3, s3, fq, ft, fe = u, s3m, fq[first], ft[first], fem
    k0 = _pair_key(q0, t0, n)
    pos = np.searchsorted(k3, k0)
    pos = np.minimum(pos, max(len(k3) - 1, 0))
    has = (len(k3) > 0) & (k3[pos] == k0) if len(k3) else np.zeros(len(k0), bool)
    s3_on_seq = np.where(has, s3[pos] if len(k3) else 0.0, np.nan)
    d_seq = -np.log(np.clip(s_seq, 1e-6, 1.0))
    both = has & (s_seq > 0)
    if fit is None:
        x = -np.log(np.clip(s3_on_seq[both], 1e-6, 1.0))
        y = d_seq[both]
        if len(x) >= 10:
            b, a = np.polyfit(x, y, 1)
        else:
            a, b = 0.0, 1.0
        fit = (float(a), float(b))
    a, b = fit

    def cal(s):
        return np.exp(-(a + b * (-np.log(np.clip(s, 1e-6, 1.0)))))

    # merged pair list: sequence pairs + structure-only pairs
    only3 = ~np.isin(k3, k0) & (fe <= miss_evalue)   # structure-only pairs need e <= miss_evalue
    if add_missing:
        q = np.concatenate([q0, fq[only3]])
        t = np.concatenate([t0, ft[only3]])
        sseq_all = np.concatenate([s_seq, np.full(only3.sum(), np.nan)])
        s3_all = np.concatenate([s3_on_seq, s3[only3]])
    else:
        q, t, sseq_all, s3_all = q0, t0, s_seq, s3_on_seq
    use = np.isfinite(s3_all)
    if node_mask is not None:
        nm = node_mask(q, t) if callable(node_mask) else np.asarray(node_mask)
        use &= nm
    s3c = cal(np.where(use, s3_all, 1.0))
    out = np.where(np.isfinite(sseq_all), sseq_all, 0.0)
    seq_ok = np.isfinite(sseq_all)
    if mode == "avg":
        comb = np.where(seq_ok, np.sqrt(np.clip(out, 1e-6, 1) * s3c), s3c)
    elif mode == "struct":
        comb = s3c
    elif mode == "max":
        comb = np.where(seq_ok, np.maximum(out, s3c), s3c)
    else:
        raise ValueError(mode)
    S = np.where(use, comb, out)
    keep = seq_ok | use
    q, t, S = q[keep], t[keep], S[keep]
    nseq = len(q0)
    pad = len(q) - int(keep[:nseq].sum())

    def ext(arr, fill):
        return np.concatenate([arr[keep[:nseq]], np.full(pad, fill, dtype=arr.dtype)])

    qq, tt = np.minimum(q, t), np.maximum(q, t)
    bits = (S * np.maximum(sb[qq], sb[tt])).astype(np.float32)
    out_st = ScoreTable(q=qq.astype(st.q.dtype), t=tt.astype(st.t.dtype), bits=bits,
                        qcov=ext(st.qcov, np.float32(1.0)), tcov=ext(st.tcov, np.float32(1.0)),
                        pident=ext(st.pident, np.float32(0.0)), evalue=ext(st.evalue, np.float32(1.0)),
                        self_bits=st.self_bits)
    stats = dict(mode=mode, fit_a=a, fit_b=b, n_seq_pairs=int(nseq), n_3di_pairs=int(len(k3)),
                 n_both=int(both.sum()), n_struct_only_added=int(pad),
                 n_pairs_changed=int(use[:nseq][keep[:nseq]].sum()) if nseq else 0,
                 n_out=int(len(qq)))
    return out_st, stats
