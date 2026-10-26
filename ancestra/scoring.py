"""S2 -- within-family pairwise similarity (DESIGN.md §1.2, §3 row S2).

Only pairs of genes in the same S0 family are scored. Ancestral genes are
never searched again: the merge (S3) aggregates these leaf-pair scores.

Backends (``ScoringConfig.backend``)
-----------------------------------
``mmseqs`` (default)
    Exact targeted alignment: one ``mmseqs createdb`` of all proteins, a
    *custom prefilter DB* (``tsv2db --output-dbtype 7``) that lists exactly the
    unordered within-family pairs ``q <= t`` (self pairs included, they give the
    self-score), then ``mmseqs align -a --alignment-mode 3`` (full SW with
    backtrace -> bitscore, coverages, identity) and ``convertalis``. Work is
    exactly sum_f g_f(g_f+1)/2 alignments and nothing else; no seeding, so no
    sensitivity loss inside families.
``diamond``
    Families are packed (largest first) into batches of <= ``batch_genes``
    proteins; per batch one ``diamond blastp`` batch-vs-batch with ``-k`` =
    batch size, then hits between different families are discarded. Wasted
    work = cross-family pairs inside a batch; DIAMOND's seeding makes this
    cheap but also loses weak within-family pairs (sensitivity flag).
``foldseek3di``
    Hook only (DESIGN.md S6): same interface, 3Di+AA scores. Not implemented.

Large families (> ``cap`` members): all-pairs scoring would be quadratic, so
each member is scored only against its ``knn`` nearest family members, found
by a DIAMOND search restricted to that family (``-k knn``, sensitive). Those
pairs then go through the same backend output format. The merge sees a
sparser graph for such families (documented limitation; sub-clustering is the
planned alternative).

Output :class:`ScoreTable` (``scores.npz``): unordered pairs ``q < t`` with
``bits``, ``qcov`` (fraction of q covered), ``tcov``, ``pident``, ``evalue``;
plus ``self_bits`` per gene (self-alignment bitscore; NaN if unavailable).
Directed duplicates (DIAMOND reports both directions) are reduced to the
higher-bitscore record.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass

import numpy as np

from .io import GeneTable, write_fasta

__all__ = ["ScoringConfig", "ScoreTable", "score_families", "BACKENDS", "within_family_pairs"]


@dataclass
class ScoringConfig:
    backend: str = "mmseqs"         # mmseqs | diamond | foldseek3di (hook)
    evalue: float = 1e-3
    cap: int = 2000                 # families larger than this: kNN pairs only
    knn: int = 200
    batch_genes: int = 6000         # diamond backend: proteins per batch
    diamond_sens: str = "--sensitive"
    threads: int = 24

    def key(self) -> dict:
        d = asdict(self)
        d.pop("threads")
        return d


@dataclass
class ScoreTable:
    q: np.ndarray        # int32, q < t
    t: np.ndarray
    bits: np.ndarray     # float32
    qcov: np.ndarray     # float32 fraction of q covered
    tcov: np.ndarray
    pident: np.ndarray
    evalue: np.ndarray
    self_bits: np.ndarray  # float32 per gene

    def __len__(self):
        return len(self.q)

    def save(self, path: str) -> None:
        tmp = path + ".tmp.npz"
        np.savez(tmp, **{k: getattr(self, k) for k in self.__dataclass_fields__})
        os.replace(tmp, path)

    @classmethod
    def load(cls, path: str) -> "ScoreTable":
        z = np.load(path)
        return cls(**{k: z[k] for k in cls.__dataclass_fields__})

    def self_bits_filled(self, lengths: np.ndarray) -> np.ndarray:
        """Self-scores with gaps filled by max(best hit, ~2.2 bits/residue)."""
        sb = self.self_bits.astype(np.float64).copy()
        miss = ~np.isfinite(sb)
        if miss.any():
            best = np.zeros(len(sb))
            np.maximum.at(best, self.q, self.bits)
            np.maximum.at(best, self.t, self.bits)
            sb[miss] = np.maximum(best[miss], 2.2 * lengths[miss])
        return sb


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

def _run(cmd, log):
    log.write("$ " + " ".join(map(str, cmd)) + "\n")
    log.flush()
    r = subprocess.run([str(c) for c in cmd], stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}")


def family_members(fam: np.ndarray) -> list[np.ndarray]:
    order = np.argsort(fam, kind="stable")
    sf = fam[order]
    cuts = np.flatnonzero(np.diff(sf)) + 1
    return np.split(order, cuts)


def within_family_pairs(members: np.ndarray, self_pairs: bool = True):
    """All unordered pairs (q <= t) of one family."""
    m = np.sort(members)
    i, j = np.triu_indices(len(m), k=0 if self_pairs else 1)
    return m[i], m[j]


def _reduce_records(q, t, bits, qcov, tcov, pid, ev, n_genes):
    """Self records -> self_bits; others oriented q<t and deduplicated (max bits)."""
    self_bits = np.full(n_genes, np.nan, dtype=np.float32)
    s = q == t
    if s.any():
        np.fmax.at(self_bits, q[s], bits[s])
    q, t, bits, qcov, tcov, pid, ev = (x[~s] for x in (q, t, bits, qcov, tcov, pid, ev))
    sw = q > t
    q, t = np.where(sw, t, q), np.where(sw, q, t)
    qcov, tcov = np.where(sw, tcov, qcov), np.where(sw, qcov, tcov)
    order = np.lexsort((-bits, t, q))
    q, t, bits, qcov, tcov, pid, ev = (x[order] for x in (q, t, bits, qcov, tcov, pid, ev))
    keep = np.ones(len(q), bool)
    if len(q):
        keep[1:] = (q[1:] != q[:-1]) | (t[1:] != t[:-1])
    return ScoreTable(q[keep].astype(np.int32), t[keep].astype(np.int32),
                      bits[keep].astype(np.float32), qcov[keep].astype(np.float32),
                      tcov[keep].astype(np.float32), pid[keep].astype(np.float32),
                      ev[keep].astype(np.float32), self_bits)


def _read_m8(path, ncol):
    if os.path.getsize(path) == 0:
        return [np.zeros(0)] * ncol
    a = np.loadtxt(path, dtype=np.float64, ndmin=2)
    return [a[:, i] for i in range(ncol)]


def _knn_pairs(gt: GeneTable, members: np.ndarray, cfg: ScoringConfig, wd: str, log):
    """Pairs (q<=t) of each member with its knn nearest family members (DIAMOND)."""
    fa = os.path.join(wd, "big.fa")
    write_fasta(gt, fa, members)
    _run(["diamond", "makedb", "--in", fa, "-d", fa + ".dmnd", "--threads", cfg.threads, "--quiet"], log)
    out = os.path.join(wd, "big.m8")
    _run(["diamond", "blastp", "-q", fa, "-d", fa + ".dmnd", "-o", out, "-k", cfg.knn,
          "-e", cfg.evalue, cfg.diamond_sens, "--threads", cfg.threads, "--quiet",
          "--outfmt", "6", "qseqid", "sseqid"], log)
    q, t = _read_m8(out, 2)
    q, t = q.astype(np.int64), t.astype(np.int64)
    lo, hi = np.minimum(q, t), np.maximum(q, t)
    lo = np.concatenate([lo, members])        # self pairs
    hi = np.concatenate([hi, members])
    key = np.unique(lo * (gt.n + 1) + hi)
    return key // (gt.n + 1), key % (gt.n + 1)


# ----------------------------------------------------------------------------
# backends
# ----------------------------------------------------------------------------

def _pair_lists(gt, fam, cfg, wd, log, stats):
    qs, ts = [], []
    n_big = 0
    for m in family_members(fam):
        if len(m) > cfg.cap:
            n_big += 1
            q, t = _knn_pairs(gt, m, cfg, wd, log)
        else:
            q, t = within_family_pairs(m)
        qs.append(q)
        ts.append(t)
    stats["n_big_families"] = n_big
    q = np.concatenate(qs).astype(np.int64)
    t = np.concatenate(ts).astype(np.int64)
    return q, t


def _backend_mmseqs(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, wd: str, log, stats):
    fa = os.path.join(wd, "genes.fa")
    write_fasta(gt, fa)
    db = os.path.join(wd, "db")
    _run(["mmseqs", "createdb", fa, db, "-v", "1"], log)
    # db key -> gene id (keys are assigned in input order; verify through .lookup)
    key2gid = np.empty(gt.n, dtype=np.int64)
    with open(db + ".lookup") as fh:
        for line in fh:
            k, nm, _ = line.split("\t")
            key2gid[int(k)] = int(nm)
    gid2key = np.empty(gt.n, dtype=np.int64)
    gid2key[key2gid] = np.arange(gt.n)
    t0 = time.time()
    q, t = _pair_lists(gt, fam, cfg, wd, log, stats)
    qk, tk = gid2key[q], gid2key[t]
    order = np.lexsort((tk, qk))
    qk, tk = qk[order], tk[order]
    stats["n_pairs_aligned"] = int(len(qk))
    pref_tsv = os.path.join(wd, "pairs.tsv")
    with open(pref_tsv, "w") as fh:
        # prefilter entry line: targetKey  prefScore  diagonal
        blk = 2_000_000
        for s in range(0, len(qk), blk):
            fh.write("".join(f"{a}\t{b}\t0\t0\n" for a, b in zip(qk[s:s + blk].tolist(),
                                                               tk[s:s + blk].tolist())))
    stats["pairs_write_s"] = round(time.time() - t0, 2)
    pref = os.path.join(wd, "pref")
    _run(["mmseqs", "tsv2db", pref_tsv, pref, "--output-dbtype", "7", "-v", "1"], log)
    os.remove(pref_tsv)
    aln = os.path.join(wd, "aln")
    t0 = time.time()
    _run(["mmseqs", "align", db, db, pref, aln, "-a", "-e", cfg.evalue, "--alignment-mode", "3",
          "--threads", cfg.threads, "-v", "1"], log)
    stats["align_s"] = round(time.time() - t0, 2)
    out = os.path.join(wd, "aln.m8")
    _run(["mmseqs", "convertalis", db, db, aln, out, "--format-output",
          "query,target,bits,qcov,tcov,pident,evalue", "--threads", cfg.threads, "-v", "1"], log)
    qq, tt, bits, qcov, tcov, pid, ev = _read_m8(out, 7)
    return _reduce_records(qq.astype(np.int64), tt.astype(np.int64), bits, qcov, tcov,
                           pid / 100.0, ev, gt.n)


def _backend_diamond(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, wd: str, log, stats):
    mems = family_members(fam)
    small = [m for m in mems if len(m) <= cfg.cap]
    big = [m for m in mems if len(m) > cfg.cap]
    small.sort(key=len, reverse=True)
    batches, cur, size = [], [], 0
    for m in small:
        if len(m) == 1:
            continue                        # singletons: self-score only (filled later)
        if size + len(m) > cfg.batch_genes and cur:
            batches.append(np.concatenate(cur))
            cur, size = [], 0
        cur.append(m)
        size += len(m)
    if cur:
        batches.append(np.concatenate(cur))
    cols = []
    t0 = time.time()
    for bi, genes in enumerate(batches):
        fa = os.path.join(wd, f"b{bi}.fa")
        write_fasta(gt, fa, genes)
        _run(["diamond", "makedb", "--in", fa, "-d", fa + ".dmnd", "--threads", cfg.threads,
              "--quiet"], log)
        out = os.path.join(wd, f"b{bi}.m8")
        _run(["diamond", "blastp", "-q", fa, "-d", fa + ".dmnd", "-o", out, "-k", len(genes),
              "-e", cfg.evalue, cfg.diamond_sens, "--threads", cfg.threads, "--quiet",
              "--outfmt", "6", "qseqid", "sseqid", "bitscore", "qcovhsp", "scovhsp", "pident",
              "evalue"], log)
        c = _read_m8(out, 7)
        q, t = c[0].astype(np.int64), c[1].astype(np.int64)
        keep = fam[q] == fam[t]
        cols.append([x[keep] for x in c])
        os.remove(fa)
        os.remove(fa + ".dmnd")
        os.remove(out)
    stats["n_batches"] = len(batches)
    stats["search_s"] = round(time.time() - t0, 2)
    stats["n_big_families"] = len(big)
    for m in big:  # kNN pairs found by DIAMOND directly carry the scores we need
        fa = os.path.join(wd, "big.fa")
        write_fasta(gt, fa, m)
        _run(["diamond", "makedb", "--in", fa, "-d", fa + ".dmnd", "--threads", cfg.threads,
              "--quiet"], log)
        out = os.path.join(wd, "big.m8")
        _run(["diamond", "blastp", "-q", fa, "-d", fa + ".dmnd", "-o", out, "-k", cfg.knn,
              "-e", cfg.evalue, cfg.diamond_sens, "--threads", cfg.threads, "--quiet",
              "--outfmt", "6", "qseqid", "sseqid", "bitscore", "qcovhsp", "scovhsp", "pident",
              "evalue"], log)
        cols.append(_read_m8(out, 7))
    if not cols:
        cols = [[np.zeros(0)] * 7]
    c = [np.concatenate([x[i] for x in cols]) for i in range(7)]
    return _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), c[2], c[3] / 100.0,
                           c[4] / 100.0, c[5] / 100.0, c[6], gt.n)


def _backend_foldseek3di(*a, **k):
    raise NotImplementedError("3Di backend is a DESIGN.md S6 hook; not part of core v0")


BACKENDS = {"mmseqs": _backend_mmseqs, "diamond": _backend_diamond,
            "foldseek3di": _backend_foldseek3di}


def score_families(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, workdir: str,
                   keep_workdir: bool = False) -> tuple[ScoreTable, dict]:
    """Run the configured backend; returns (ScoreTable, stats)."""
    os.makedirs(workdir, exist_ok=True)
    stats: dict = {"backend": cfg.backend}
    sizes = np.bincount(fam)
    stats["n_families"] = int((sizes > 0).sum())
    stats["sum_g2"] = int((sizes.astype(np.int64) ** 2).sum())
    stats["max_family"] = int(sizes.max()) if len(sizes) else 0
    t0 = time.time()
    with open(os.path.join(workdir, "commands.log"), "a") as log:
        st = BACKENDS[cfg.backend](gt, fam, cfg, workdir, log, stats)
    stats["wall_s"] = round(time.time() - t0, 2)
    stats["n_scored_pairs"] = int(len(st))
    stats["n_self_scores"] = int(np.isfinite(st.self_bits).sum())
    if not keep_workdir:
        for f in os.listdir(workdir):
            if f != "commands.log":
                p = os.path.join(workdir, f)
                shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    return st, stats
