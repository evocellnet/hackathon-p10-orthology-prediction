"""S2 -- within-family pairwise similarity (DESIGN.md §1.2, §3 row S2).

Only pairs of genes in the same S0 family are scored. Ancestral genes are
never searched again: the merge (S3) aggregates these leaf-pair scores.

Backends (``ScoringConfig.backend``)
-----------------------------------
``diamond`` (PI choice after r2_algorithms_compute; the ``hybrid`` default adds a rescue)
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
``hybrid`` (default since the core-v0 toy grid, docs/core_v0.md)
    DIAMOND blocks, then ``mmseqs align`` only on the in-family pairs (families
    <= cap_all) DIAMOND did not report: exhaustive recall at ~DIAMOND cost plus
    the missed fraction (~7-10% of pairs on the toy).
``foldseek3di``
    Hook only (DESIGN.md S6): same interface, 3Di+AA scores. Not implemented.

Large families (r2 M3): g <= ``cap_all`` (5k): all in-family hits kept.
5k < g: per query keep the top ``top_per_species`` (2) hits per other species
plus <= ``top_same`` (20) same-species hits (a pair survives if either direction
keeps it); dropped edges are logged per family. g > ``g_max`` (20k) should have
been split by S0; such families are counted (``n_families_over_gmax``) and get
the same reduction.

Cross-family rescue (``family_rescue``, default on; r2 external-PI M1): the
DIAMOND pass also sees hits *between* S0 families (same block). Families linked
by a hit with min(qcov, tcov) >= 0.5 and e <= 1e-5 are joined and only the
joined families are re-scored (``score_with_rescue``). On the toy this repairs
about half of the S0 fissions (docs/core_v0.md).

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

__all__ = ["ScoringConfig", "ScoreTable", "score_families", "score_with_rescue", "join_families",
           "BACKENDS", "within_family_pairs", "apply_cap"]


@dataclass
class ScoringConfig:
    backend: str = "hybrid"         # hybrid | diamond | mmseqs | foldseek3di (hook)
    evalue: float = 10.0            # within-family: homology already established by S0
    block_genes: int = 100_000      # diamond: proteins per block (whole families)
    diamond_sens: str = "--sensitive"
    diamond_b: float = 4.0
    cap_all: int = 5000             # families <= cap_all keep all in-family hits
    g_max: int = 20_000             # families > g_max should be split by S0 (flagged here)
    top_per_species: int = 2        # cap rule for families > cap_all
    top_same: int = 20
    family_rescue: bool = True      # join S0 families linked by a cross-family hit (S0 fission repair)
    rescue_cov: float = 0.5         # ... with min(qcov, tcov) >= rescue_cov
    rescue_max_join: int = 7500     # never join families into one larger than this (S0 s2_gmax)
    rescue_max_family: int = 1000   # hybrid: MMseqs2 re-alignment of DIAMOND-missed pairs only in
                                    # families <= this size (bounded cost; big families get many hits)
    rescue_evalue: float = 1e-5     # ... and evalue <= rescue_evalue
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
    cross = None           # optional dict of cross-family hits (diamond/hybrid), saved separately

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
    """Read a tab-separated hit table (ids first two columns) with the pandas C parser."""
    if os.path.getsize(path) == 0:
        return [np.zeros(0)] * ncol
    import pandas as pd
    df = pd.read_csv(path, sep="\t", header=None, engine="c",
                     dtype={i: (np.int64 if i < 2 else np.float32) for i in range(ncol)})
    return [df[i].to_numpy() for i in range(ncol)]


def _iter_m8(path, ncol, chunksize=20_000_000):
    """Chunked version of _read_m8 for very large DIAMOND outputs."""
    if os.path.getsize(path) == 0:
        return
    import pandas as pd
    for df in pd.read_csv(path, sep="\t", header=None, engine="c", chunksize=chunksize,
                          dtype={i: (np.int64 if i < 2 else np.float32) for i in range(ncol)}):
        yield [df[i].to_numpy() for i in range(ncol)]


def _orient_dedup(c):
    """Orient non-self records to q<t (swapping coverages) and keep the best record per pair.
    Self records (q == t) are kept unchanged."""
    q, t = c[0], c[1]
    selfm = q == t
    out_self = [x[selfm] for x in c]
    q, t, bits, qcov, tcov, pid, ev = (x[~selfm] for x in c)
    sw = q > t
    q, t = np.where(sw, t, q), np.where(sw, q, t)
    qcov, tcov = np.where(sw, tcov, qcov), np.where(sw, qcov, tcov)
    order = np.lexsort((-bits, t, q))
    q, t, bits, qcov, tcov, pid, ev = (x[order] for x in (q, t, bits, qcov, tcov, pid, ev))
    keep = np.ones(len(q), bool)
    if len(q):
        keep[1:] = (q[1:] != q[:-1]) | (t[1:] != t[:-1])
    red = [x[keep] for x in (q, t, bits, qcov, tcov, pid, ev)]
    return [np.concatenate([a, b]) for a, b in zip(out_self, red)]


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
    out = ScoreTable(*(getattr(st, k)[keep_e] for k in ("q", "t", "bits", "qcov", "tcov",
                                                        "pident", "evalue")), st.self_bits)
    return out


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
        bcols = []
        for c in _iter_m8(out, 7):
            q, t = c[0], c[1]
            keep = fam[q] == fam[t]
            stats["n_cross_family_hits"] = stats.get("n_cross_family_hits", 0) + int((~keep).sum())
            stats["n_raw_hits"] = stats.get("n_raw_hits", 0) + int(len(q))
            cr = [x[~keep] for x in c]
            # cross-family hits: only those that could join families are kept (S0-fission signal)
            good = np.minimum(cr[3], cr[4]) >= 100 * cfg.rescue_cov * 0.6
            stats.setdefault("_cross", []).append([x[good] for x in cr])
            bcols.append([x[keep] for x in c])
        if bcols:
            cols.append(_orient_dedup([np.concatenate([x[i] for x in bcols]) for i in range(7)]))
        del bcols
        for f in (fa, fa + ".dmnd", out):
            os.remove(f)
    stats["n_blocks"] = stats.get("n_blocks", 0) + len(blocks)
    stats["diamond_s"] = round(stats.get("diamond_s", 0) + time.time() - t0, 2)
    if not cols:
        return [np.zeros(0)] * 7
    c = [np.concatenate([x[i] for x in cols]) for i in range(7)]
    c[3], c[4], c[5] = c[3] / 100.0, c[4] / 100.0, c[5] / 100.0
    return c


def _mmseqs_align_pairs(gt: GeneTable, q: np.ndarray, t: np.ndarray, cfg, wd, log, stats):
    """Align exactly the listed pairs (q, t) with mmseqs align on a custom prefilter DB."""
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
    qk, tk = gid2key[q], gid2key[t]
    order = np.lexsort((tk, qk))
    qk, tk = qk[order], tk[order]
    stats["n_pairs_aligned_mmseqs"] = int(len(qk))
    if len(qk) == 0:
        return [np.zeros(0)] * 7
    pref_tsv = os.path.join(wd, "pairs.tsv")
    with open(pref_tsv, "w") as fh:
        # tsv2db line: queryKey  targetKey  prefScore  diagonal
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
    return c


def _all_small_family_pairs(fam, cap_all):
    qs, ts, big = [], [], set()
    for m in family_members(fam):
        if len(m) > cap_all:
            big.add(int(fam[m[0]]))
            continue
        q, t = within_family_pairs(m)
        qs.append(q)
        ts.append(t)
    q = np.concatenate(qs).astype(np.int64) if qs else np.zeros(0, np.int64)
    t = np.concatenate(ts).astype(np.int64) if ts else np.zeros(0, np.int64)
    return q, t, big


def _backend_mmseqs(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, wd: str, log, stats):
    """Exhaustive: mmseqs align on the explicit in-family pair list (families <= cap_all);
    larger families go through the DIAMOND block path (then the cap rule)."""
    q, t, big = _all_small_family_pairs(fam, cfg.cap_all)
    c = _mmseqs_align_pairs(gt, q, t, cfg, wd, log, stats)
    if big:
        c2 = _diamond_blocks(gt, fam, pack_blocks(fam, cfg.block_genes, big), cfg, wd, log, stats)
        c = [np.concatenate([a, b]) for a, b in zip(c, c2)]
    return _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), *c[2:], gt.n)


def _backend_hybrid(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, wd: str, log, stats):
    """DIAMOND blocks first; then mmseqs align only the in-family pairs (families <= cap_all)
    that DIAMOND did not report (targeted rescue of divergent pairs)."""
    c = _diamond_blocks(gt, fam, pack_blocks(fam, cfg.block_genes), cfg, wd, log, stats)
    q, t, _ = _all_small_family_pairs(fam, min(cfg.cap_all, cfg.rescue_max_family))
    M = gt.n + 1
    dq, dt = c[0].astype(np.int64), c[1].astype(np.int64)
    found = np.unique(np.minimum(dq, dt) * M + np.maximum(dq, dt))
    want = q * M + t
    miss = ~np.isin(want, found)
    stats["n_pairs_diamond"] = int(len(found))
    stats["n_pairs_rescue_candidates"] = int(miss.sum())
    c2 = _mmseqs_align_pairs(gt, q[miss], t[miss], cfg, wd, log, stats)
    c = [np.concatenate([a, b]) for a, b in zip(c, c2)]
    return _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), *c[2:], gt.n)


def _backend_diamond(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, wd: str, log, stats):
    """Whole families packed into blocks of <= block_genes; DIAMOND block-vs-block
    (-k 0 --max-hsps 1 -b 4 -c 1 --sensitive), cross-family hits discarded."""
    c = _diamond_blocks(gt, fam, pack_blocks(fam, cfg.block_genes), cfg, wd, log, stats)
    return _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), *c[2:], gt.n)


def _backend_foldseek3di(*a, **k):
    raise NotImplementedError("3Di backend is a DESIGN.md S6 hook; not part of core v0")


BACKENDS = {"mmseqs": _backend_mmseqs, "diamond": _backend_diamond, "hybrid": _backend_hybrid,
            "foldseek3di": _backend_foldseek3di}


def join_families(fam: np.ndarray, cross: dict, cov: float, evalue: float, max_size: int | None = None):
    """Union families linked by a cross-family hit passing the coverage/e-value guard
    (strongest hits first; a join that would exceed ``max_size`` genes is skipped, so
    S0's deliberate splits of giant components are not undone).
    Returns (new family array relabelled 0.., number of joins)."""
    parent = np.arange(int(fam.max()) + 1)
    fsize = np.bincount(fam, minlength=len(parent)).astype(np.int64)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    m = np.flatnonzero((np.minimum(cross["qcov"], cross["tcov"]) >= cov) & (cross["evalue"] <= evalue))
    m = m[np.argsort(-np.asarray(cross["bits"])[m], kind="stable")] if "bits" in cross else m
    n = 0
    for a, b in zip(fam[cross["q"][m]].tolist(), fam[cross["t"][m]].tolist()):
        ra, rb = find(a), find(b)
        if ra != rb:
            if max_size is not None and fsize[ra] + fsize[rb] > max_size:
                continue
            lo_, hi_ = min(ra, rb), max(ra, rb)
            parent[hi_] = lo_
            fsize[lo_] += fsize[hi_]
            n += 1
    roots = np.array([find(f) for f in range(len(parent))])
    _, new = np.unique(roots[fam], return_inverse=True)
    return new.astype(np.int32), n


def score_with_rescue(gt: GeneTable, fam: np.ndarray, cfg: ScoringConfig, workdir: str,
                      pass1_path: str | None = None):
    """score_families + cross-family rescue: families joined by ``join_families`` are
    re-scored (only those); returns (ScoreTable, final families, stats)."""
    st, stats = score_families(gt, fam, cfg, workdir)
    if not cfg.family_rescue or st.cross is None or len(st.cross["q"]) == 0:
        stats["rescue_joins"] = 0
        return st, fam, stats
    fam2, nj = join_families(fam, st.cross, cfg.rescue_cov, cfg.rescue_evalue, cfg.rescue_max_join)
    if pass1_path:
        st.save(pass1_path)
    stats["rescue_joins"] = nj
    if nj == 0:
        return st, fam, stats
    # families whose membership changed: re-score them only (others -> singletons, skipped)
    sizes_old = np.bincount(fam)
    changed_new = np.zeros(fam2.max() + 1, bool)
    first_old = {}
    for g, (fo, fn) in enumerate(zip(fam.tolist(), fam2.tolist())):
        if first_old.setdefault(fn, fo) != fo:
            changed_new[fn] = True
    ch = changed_new[fam2]
    sub = np.where(ch, fam2, fam2.max() + 1 + np.arange(gt.n)).astype(np.int64)
    _, sub = np.unique(sub, return_inverse=True)
    st2, stats2 = score_families(gt, sub, cfg, workdir)
    keep_old = ~(ch[st.q] | ch[st.t])
    cols = {k: np.concatenate([getattr(st, k)[keep_old], getattr(st2, k)[~(~ch[st2.q] & ~ch[st2.t])]])
            for k in ("q", "t", "bits", "qcov", "tcov", "pident", "evalue")}
    out = ScoreTable(**cols, self_bits=st.self_bits)
    out.cross = st.cross
    stats.update(rescue_rescored_genes=int(ch.sum()), rescue_wall_s=stats2["wall_s"],
                 n_scored_pairs=int(len(out)), n_families_after_rescue=int(fam2.max() + 1),
                 n_families_before_rescue=int((sizes_old > 0).sum()))
    return out, fam2, stats


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
    cr = stats.pop("_cross", [])
    if cr:
        c = [np.concatenate([x[i] for x in cr]) for i in range(7)]
        cross = dict(q=c[0].astype(np.int32), t=c[1].astype(np.int32), bits=c[2].astype(np.float32),
                     qcov=(c[3] / 100).astype(np.float32), tcov=(c[4] / 100).astype(np.float32),
                     evalue=c[6].astype(np.float32))
    else:
        cross = {k: np.zeros(0, np.float32) for k in ("q", "t", "bits", "qcov", "tcov", "evalue")}
    st.cross = cross     # cross-family DIAMOND hits (not part of the table; see pipeline)
    stats["wall_s"] = round(time.time() - t0, 2)
    stats["n_scored_pairs"] = int(len(st))
    stats["n_self_scores"] = int(np.isfinite(st.self_bits).sum())
    if not keep_workdir:
        for f in os.listdir(workdir):
            if f != "commands.log":
                p = os.path.join(workdir, f)
                shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    return st, stats
