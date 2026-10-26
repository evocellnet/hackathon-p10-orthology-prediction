"""S2 -- within-family pairwise similarity (DESIGN.md §1.2, §3 row S2).

Only pairs of genes in the same S0 family are scored. Ancestral genes are
never searched again: the merge (S3) aggregates these leaf-pair scores.

Backends (``ScoringConfig.backend``)
-----------------------------------
``diamond`` (default; PI decision after r2_algorithms_compute)
    Whole families are packed (largest first) into blocks of <= ``block_genes``
    (100k) proteins; a family larger than a block is its own block. Per block
    one ``diamond blastp`` block-vs-block with ``-k 0 --max-hsps 1 -b 4 -c 1
    --sensitive``; hits between different families are discarded (their count
    is reported: a free S0-fission signal). Cost ~ N*B, not sum g_f^2.
``mmseqs`` (exhaustive reference)
    One ``mmseqs createdb``, a *custom prefilter DB* (``tsv2db --output-dbtype
    7``) listing exactly the unordered within-family pairs ``q <= t`` (self pairs
    included -> self-scores), then ``mmseqs align -a --alignment-mode 3`` (full
    SW) and ``convertalis``. Work = sum_f g_f(g_f+1)/2 alignments, no seeding
    loss. Families > ``cap_all`` go through the DIAMOND block path.
``foldseek3di``
    Hook only (DESIGN.md S6): same interface, 3Di+AA scores. Not implemented.

Large families (r2 M3): g <= ``cap_all`` (5k): all in-family hits kept.
5k < g: per query keep the top ``top_per_species`` (2) hits per other species
plus <= ``top_same`` (20) same-species hits (a pair survives if either direction
keeps it); dropped edges are logged per family. g > ``g_max`` (20k) should have
been split by S0; such families are counted (``n_families_over_gmax``) and get
the same reduction.

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
    backend: str = "diamond"        # diamond | mmseqs | foldseek3di (hook)
    evalue: float = 10.0            # within-family: homology already established by S0
    block_genes: int = 100_000      # diamond: proteins per block (whole families)
    diamond_sens: str = "--sensitive"
    diamond_b: float = 4.0
    cap_all: int = 5000             # families <= cap_all keep all in-family hits
    g_max: int = 20_000             # families > g_max should be split by S0 (flagged here)
    top_per_species: int = 2        # cap rule for families > cap_all
    top_same: int = 20
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


# ----------------------------------------------------------------------------
# large-family cap (r2_algorithms_compute M3)
# ----------------------------------------------------------------------------

def apply_cap(st: ScoreTable, fam: np.ndarray, species_of: np.ndarray, cfg: "ScoringConfig",
              stats: dict) -> ScoreTable:
    """For families with > ``cap_all`` members keep, per query, the top ``top_per_species``
    hits per other species plus the top ``top_same`` same-species hits (a pair is kept if
    either direction keeps it). Families <= cap_all keep all hits."""
    sizes = np.bincount(fam)
    big = sizes > cfg.cap_all
    stats["n_families_over_cap"] = int(big.sum())
    stats["n_families_over_gmax"] = int((sizes > cfg.g_max).sum())
    sel = big[fam[st.q]]
    if not sel.any():
        stats["n_pairs_dropped_by_cap"] = 0
        return st
    idx = np.flatnonzero(sel)
    q = np.concatenate([st.q[idx], st.t[idx]]).astype(np.int64)
    t = np.concatenate([st.t[idx], st.q[idx]]).astype(np.int64)
    b = np.concatenate([st.bits[idx], st.bits[idx]])
    e = np.concatenate([idx, idx])
    tsp = species_of[t].astype(np.int64)
    samesp = tsp == species_of[q]
    order = np.lexsort((-b, tsp, q))
    q, tsp, samesp, e = q[order], tsp[order], samesp[order], e[order]
    grp_start = np.r_[True, (q[1:] != q[:-1]) | (tsp[1:] != tsp[:-1])]
    gid = np.cumsum(grp_start) - 1
    first = np.flatnonzero(grp_start)
    rank = np.arange(len(q)) - first[gid]
    keep_dir = np.where(samesp, rank < cfg.top_same, rank < cfg.top_per_species)
    keep_e = np.zeros(len(st.q), bool)
    keep_e[~sel] = True
    keep_e[e[keep_dir]] = True
    per_fam = {}
    dropped = np.flatnonzero(~keep_e)
    if len(dropped):
        f, c = np.unique(fam[st.q[dropped]], return_counts=True)
        per_fam = {int(a): int(n) for a, n in zip(f, c)}
    stats["n_pairs_dropped_by_cap"] = int(len(dropped))
    stats["pairs_dropped_per_family"] = per_fam
    return ScoreTable(*(getattr(st, k)[keep_e] for k in ("q", "t", "bits", "qcov", "tcov",
                                                         "pident", "evalue")), st.self_bits)


# ----------------------------------------------------------------------------
# backends
# ----------------------------------------------------------------------------

def pack_blocks(fam: np.ndarray, block_genes: int, families=None) -> list[np.ndarray]:
    """Pack whole families (largest first) into blocks of <= block_genes proteins;
    a family larger than the block is its own block. Singleton families are skipped."""
    mems = [m for m in family_members(fam) if len(m) > 1]
    if families is not None:
        mems = [m for m in mems if fam[m[0]] in families]
    mems.sort(key=len, reverse=True)
    blocks, cur, size = [], [], 0
    for m in mems:
        if size + len(m) > block_genes and cur:
            blocks.append(np.concatenate(cur))
            cur, size = [], 0
        cur.append(m)
        size += len(m)
    if cur:
        blocks.append(np.concatenate(cur))
    return blocks


def _diamond_blocks(gt, fam, blocks, cfg, wd, log, stats):
    cols = []
    t0 = time.time()
    for bi, genes in enumerate(blocks):
        fa = os.path.join(wd, f"b{bi}.fa")
        write_fasta(gt, fa, genes)
        _run(["diamond", "makedb", "--in", fa, "-d", fa + ".dmnd", "--threads", cfg.threads,
              "--quiet"], log)
        out = os.path.join(wd, f"b{bi}.m8")
        cmd = ["diamond", "blastp", "-q", fa, "-d", fa + ".dmnd", "-o", out, "-k", "0",
               "--max-hsps", "1", "-b", cfg.diamond_b, "-c", "1", "-e", cfg.evalue,
               "--threads", cfg.threads, "--quiet",
               "--outfmt", "6", "qseqid", "sseqid", "bitscore", "qcovhsp", "scovhsp", "pident",
               "evalue"]
        if cfg.diamond_sens:
            cmd.append(cfg.diamond_sens)
        _run(cmd, log)
        c = _read_m8(out, 7)
        q, t = c[0].astype(np.int64), c[1].astype(np.int64)
        keep = fam[q] == fam[t]
        stats["n_cross_family_hits"] = stats.get("n_cross_family_hits", 0) + int((~keep).sum())
        cols.append([x[keep] for x in c])
        for f in (fa, fa + ".dmnd", out):
            os.remove(f)
    stats["n_blocks"] = stats.get("n_blocks", 0) + len(blocks)
    stats["diamond_s"] = round(stats.get("diamond_s", 0) + time.time() - t0, 2)
    if not cols:
        return [np.zeros(0)] * 7
    c = [np.concatenate([x[i] for x in cols]) for i in range(7)]
    c[3], c[4], c[5] = c[3] / 100.0, c[4] / 100.0, c[5] / 100.0
    return c


def _backend_mmseqs(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, wd: str, log, stats):
    """Exhaustive: mmseqs align on the explicit in-family pair list (families <= cap_all);
    larger families go through the DIAMOND block path (then the cap rule)."""
    fa = os.path.join(wd, "genes.fa")
    write_fasta(gt, fa)
    db = os.path.join(wd, "db")
    _run(["mmseqs", "createdb", fa, db, "-v", "1"], log)
    key2gid = np.empty(gt.n, dtype=np.int64)
    with open(db + ".lookup") as fh:
        for line in fh:
            k, nm, _ = line.split("\t")
            key2gid[int(k)] = int(nm)
    gid2key = np.empty(gt.n, dtype=np.int64)
    gid2key[key2gid] = np.arange(gt.n)
    t0 = time.time()
    qs, ts, big = [], [], set()
    for m in family_members(fam):
        if len(m) > cfg.cap_all:
            big.add(int(fam[m[0]]))
            continue
        q, t = within_family_pairs(m)
        qs.append(q)
        ts.append(t)
    q = np.concatenate(qs).astype(np.int64)
    t = np.concatenate(ts).astype(np.int64)
    qk, tk = gid2key[q], gid2key[t]
    order = np.lexsort((tk, qk))
    qk, tk = qk[order], tk[order]
    stats["n_pairs_aligned"] = int(len(qk))
    pref_tsv = os.path.join(wd, "pairs.tsv")
    with open(pref_tsv, "w") as fh:
        # prefilter entry line: queryKey  targetKey  prefScore  diagonal
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
    c = _read_m8(out, 7)
    c[5] = c[5] / 100.0
    if big:
        c2 = _diamond_blocks(gt, fam, pack_blocks(fam, cfg.block_genes, big), cfg, wd, log, stats)
        c = [np.concatenate([a, b]) for a, b in zip(c, c2)]
    return _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), *c[2:], gt.n)


def _backend_diamond(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, wd: str, log, stats):
    """Whole families packed into blocks of <= block_genes; DIAMOND block-vs-block
    (-k 0 --max-hsps 1 -b 4 -c 1 --sensitive), cross-family hits discarded."""
    c = _diamond_blocks(gt, fam, pack_blocks(fam, cfg.block_genes), cfg, wd, log, stats)
    return _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), *c[2:], gt.n)


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
    st = apply_cap(st, fam, gt.species_of, cfg, stats)
    stats["wall_s"] = round(time.time() - t0, 2)
    stats["n_scored_pairs"] = int(len(st))
    stats["n_self_scores"] = int(np.isfinite(st.self_bits).sum())
    if not keep_workdir:
        for f in os.listdir(workdir):
            if f != "commands.log":
                p = os.path.join(workdir, f)
                shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    return st, stats
