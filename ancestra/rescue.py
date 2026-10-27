"""S0 family rescue after S2 scoring (`--family-rescue {none,3di,profile,both}`).

Problem: S0 (sequence, DIAMOND e <= 1e-5) leaves remote homologs in singletons and
small families; their orthologs in other families can then never be paired. This
module attaches small families to a larger family with remote-homology evidence and
re-scores only the new pairs, so it runs on cached S0 + S2 outputs.

Units (deterministic; no random step)
-------------------------------------
* *small* family: <= ``small_max`` genes (10). Queries = one gene per S0 stage-1
  cluster inside small families (the DIAMOND cluster representative when it lies in
  the family, else the longest member; without stage-1 info every gene is a unit).
* *target* family: > ``small_max`` genes. Targets = representatives of the
  ``reps_per_family`` largest stage-1 clusters of the family (ties: rep gene id).

Evidence (``mode``)
-------------------
* ``3di``     Foldseek 3Di+AA search (``--alignment-type 2``) of small-family queries
              against target representatives (exact-sequence AFDB 3Di; genes without
              3Di are skipped). Reuses :mod:`ancestra.structure`.
* ``profile`` MMseqs2 profiles: one profile per target representative (centre), built
              with ``result2profile`` from the centre's alignment to <= ``profile_members``
              genes of its family (the family's other representatives first, then its
              own cluster), searched against the small-family queries (profile query vs
              sequence target, ``-s`` ``mm_sens``).
* ``both``    union of the two hit lists.

Attach rule (S6 ``attach``, restricted to small -> target)
---------------------------------------------------------
A hit links a small family S to target family T if e <= ``evalue`` and the aligned
span covers >= ``cov_shorter`` of the shorter protein. Each S joins only the family of
its single best link (lowest e, then higher bits, then larger T, then lower T label).
Joins are applied in increasing e-value; a join that would make T exceed
``max_size`` genes is skipped. Two target families are never joined, and small
families never chain.

Re-scoring (only new pairs)
---------------------------
New pairs = attached gene x (old members of T + other genes attached to T). Same
backend as the hybrid S2: DIAMOND (``--sensitive -k 0 --max-hsps 1 -b 4 -c 1``,
e <= 10) of the attached genes against the receiving families, then MMseqs2 SW of
the in-family pairs DIAMOND missed in families <= ``rescue_max_family`` (1000) genes,
plus self-scores of attached genes. For families > ``cap_all`` (5000) only the
attached gene's top 2 hits per other species and 20 same-species hits are kept.

Merge-side guard (``guard``; R4)
--------------------------------
After a first merge, an attached small family is kept only if at least one of its
genes is predicted orthologous to an *original* member t of its target family and
the pair is a reciprocal best hit (norm_max score) between their two species within
the new family. Otherwise it is detached (its genes get their old family back, all
pairs to other families are dropped) and the merge is run again
(:func:`guard_detach`).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass

import numpy as np

from .io import GeneTable, read_fasta, write_fasta
from .scoring import ScoreTable, ScoringConfig, _iter_m8, _mmseqs_align_pairs, _reduce_records

__all__ = ["RescueConfig", "RESCUE_MODES", "stage1_rep_of", "select_units", "load_3di_any",
           "search_3di", "search_profiles", "attach_small", "rescore_new_pairs", "run_rescue",
           "guard_detach", "revert_families", "write_families"]

RESCUE_MODES = ("none", "3di", "profile", "both")


@dataclass
class RescueConfig:
    mode: str = "none"
    small_max: int = 10
    evalue: float = 1e-3
    cov_shorter: float = 0.3
    max_size: int = 7500            # = S0 s2_gmax / S2 rescue_max_join
    reps_per_family: int = 20
    fs_sens: float = 7.5
    fs_max_seqs: int = 100
    mm_sens: float = 7.5
    mm_max_seqs: int = 1000
    profile_members: int = 100
    guard: bool = False
    threads: int = 10

    def key(self) -> dict:
        d = asdict(self)
        d.pop("threads")
        if self.mode == "profile":
            for k in ("fs_sens", "fs_max_seqs"):
                d.pop(k)
        if self.mode == "3di":
            for k in ("mm_sens", "mm_max_seqs", "profile_members"):
                d.pop(k)
        return d


def _run(cmd, log):
    cmd = [str(c) for c in cmd]
    log.write("$ " + " ".join(cmd) + "\n")
    log.flush()
    r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(cmd)}")


# ----------------------------------------------------------------------------
# units
# ----------------------------------------------------------------------------

def stage1_rep_of(gt: GeneTable, s0_dir: str | None) -> np.ndarray:
    """gene -> gene id of its S0 stage-1 cluster representative. The S0 run clusters a
    canonical copy of the input (sequences sorted by (length desc, sequence, id), renamed
    ``c<rank>``, ``families._canonical_fasta``); the same order is rebuilt here from the
    gene table's sequences. Without ``s0_dir`` (or stage-1 file) every gene is its own rep."""
    rep = np.arange(gt.n, dtype=np.int64)
    p = os.path.join(s0_dir, "stage1_clusters.tsv") if s0_dir else None
    if not p or not os.path.exists(p):
        return rep
    if gt.seqs is None:
        raise ValueError("stage1_rep_of needs sequences (gt.seqs)")
    order = sorted(range(gt.n), key=lambda g: (-len(gt.seqs[g]), gt.seqs[g], gt.names[g]))
    canon = np.asarray(order, dtype=np.int64)          # canonical rank -> gene id
    with open(p) as fh:
        for line in fh:
            r, m = line.rstrip("\n").split("\t")[:2]
            rep[canon[int(m[1:])]] = canon[int(r[1:])]
    return rep


def select_units(fam: np.ndarray, rep_of: np.ndarray, lengths: np.ndarray, small_max: int,
                 reps_per_family: int):
    """(query genes of small families, target genes of large families, family sizes)."""
    sizes = np.bincount(fam)
    n = len(fam)
    # cluster key inside a family: (family, rep) -- members whose rep lies in another family
    # (S0 reassignment) form their own unit keyed by themselves
    rep_ok = fam[rep_of] == fam
    key_rep = np.where(rep_ok, rep_of, np.arange(n))
    order = np.lexsort((np.arange(n), -lengths.astype(np.int64), key_rep, fam))
    f_s, k_s = fam[order], key_rep[order]
    start = np.r_[True, (f_s[1:] != f_s[:-1]) | (k_s[1:] != k_s[:-1])]
    gidx = np.cumsum(start) - 1
    first = order[start]                      # longest member of each unit
    ufam, ukey = f_s[start], k_s[start]
    usize = np.bincount(gidx)
    unit_gene = np.where(rep_ok[ukey] & (fam[ukey] == ufam) & (rep_of[ukey] == ukey), ukey, first)
    small = sizes[ufam] <= small_max
    queries = np.sort(unit_gene[small])
    big = np.flatnonzero(~small)
    o = np.lexsort((unit_gene[big], -usize[big], ufam[big]))
    big = big[o]
    bf = ufam[big]
    st = np.r_[True, bf[1:] != bf[:-1]]
    rank = np.arange(len(big)) - np.flatnonzero(st)[np.cumsum(st) - 1]
    targets = np.sort(unit_gene[big[rank < reps_per_family]])
    return queries.astype(np.int64), targets.astype(np.int64), sizes


# ----------------------------------------------------------------------------
# evidence
# ----------------------------------------------------------------------------

def load_3di_any(path: str, gt: GeneTable) -> dict:
    """protein id -> 3Di string; the id is whichever ``|`` field of the header is a gene of
    ``gt`` (OrthoBench headers start with the protein id, QfO ones with the proteome id).
    Strings whose length differs from the protein are dropped."""
    lens = {n: int(l) for n, l in zip(gt.names, gt.lengths)}
    out = {}
    for hid, s in read_fasta(path):
        for pid in hid.split()[0].split("|"):
            if pid in lens:
                s = s.strip().upper()
                if len(s) == lens[pid]:
                    out[pid] = s
                break
    return out


def search_3di(gt: GeneTable, queries, targets, ss: dict, wd: str, cfg: RescueConfig, log) -> dict:
    from . import structure as S
    os.makedirs(wd, exist_ok=True)
    qk = S.build_db(gt, gt.seqs, ss, os.path.join(wd, "qdb"), log, genes=queries)
    tk = S.build_db(gt, gt.seqs, ss, os.path.join(wd, "tdb"), log, genes=targets)
    h = S.search(os.path.join(wd, "qdb"), os.path.join(wd, "tdb"), os.path.join(wd, "search"),
                 threads=cfg.threads, evalue=cfg.evalue, max_seqs=cfg.fs_max_seqs, sens=cfg.fs_sens,
                 log=log)
    h["n_queries_3di"], h["n_targets_3di"] = len(qk), len(tk)
    return h


def _mm_lookup(db):
    keys, gids = [], []
    with open(db + ".lookup") as fh:
        for line in fh:
            k, nm = line.split("\t")[:2]
            keys.append(int(k))
            gids.append(int(nm))
    return np.asarray(keys, np.int64), np.asarray(gids, np.int64)


def search_profiles(gt: GeneTable, fam: np.ndarray, rep_of: np.ndarray, queries, targets, wd: str,
                    cfg: RescueConfig, log) -> dict:
    """One MMseqs2 profile per target representative, searched against the small-family
    queries. Returns hits with q = profile centre gene, t = query gene."""
    os.makedirs(wd, exist_ok=True)
    tmpd = os.path.join(wd, "tmp")
    genes_needed = set(targets.tolist()) | set(queries.tolist())
    # members for each centre: other target reps of its family, then its own cluster members
    tset_by_fam: dict = {}
    for g in targets.tolist():
        tset_by_fam.setdefault(int(fam[g]), []).append(g)
    tcent = set(targets.tolist())
    clus: dict = {}
    for g in np.flatnonzero(np.isin(rep_of, targets) & (fam == fam[rep_of])).tolist():
        clus.setdefault(int(rep_of[g]), []).append(g)
    pairs_q, pairs_t = [], []
    P = cfg.profile_members
    for c in targets.tolist():
        mem = [x for x in tset_by_fam[int(fam[c])] if x != c][:P]
        if len(mem) < P:
            own = [x for x in sorted(clus.get(c, [])) if x != c and x not in tcent]
            mem += own[:P - len(mem)]
        pairs_q += [c] * (len(mem) + 1)
        pairs_t += [c] + mem
        genes_needed.update(mem)
    genes = np.asarray(sorted(genes_needed), np.int64)
    fa = os.path.join(wd, "genes.fa")
    write_fasta(gt, fa, genes)
    db = os.path.join(wd, "db")
    _run(["mmseqs", "createdb", fa, db, "--shuffle", "0", "-v", "1"], log)
    os.remove(fa)
    keys, gids = _mm_lookup(db)
    g2k = dict(zip(gids.tolist(), keys.tolist()))
    pq = np.fromiter((g2k[x] for x in pairs_q), np.int64, len(pairs_q))
    pt = np.fromiter((g2k[x] for x in pairs_t), np.int64, len(pairs_t))
    o = np.lexsort((pt, pq))
    pq, pt = pq[o], pt[o]
    with open(os.path.join(wd, "pairs.tsv"), "w") as fh:
        fh.write("".join(f"{a}\t{b}\t0\t0\n" for a, b in zip(pq.tolist(), pt.tolist())))
    _run(["mmseqs", "tsv2db", os.path.join(wd, "pairs.tsv"), os.path.join(wd, "pref"),
          "--output-dbtype", "7", "-v", "1"], log)
    # centre vs members; members must align with e <= 1e-3 to enter the profile
    _run(["mmseqs", "align", db, db, os.path.join(wd, "pref"), os.path.join(wd, "aln"), "-a",
          "-e", 1e-3, "--threads", cfg.threads, "-v", "1"], log)
    _run(["mmseqs", "result2profile", db, db, os.path.join(wd, "aln"), os.path.join(wd, "prof"),
          "--threads", cfg.threads, "-v", "1"], log)
    qlist = os.path.join(wd, "qkeys.txt")
    with open(qlist, "w") as fh:
        fh.write("".join(f"{g2k[g]}\n" for g in queries.tolist()))
    _run(["mmseqs", "createsubdb", qlist, db, os.path.join(wd, "qdb"), "-v", "1"], log)
    res = os.path.join(wd, "res")
    _run(["mmseqs", "search", os.path.join(wd, "prof"), os.path.join(wd, "qdb"), res, tmpd,
          "-s", cfg.mm_sens, "-e", cfg.evalue, "--max-seqs", cfg.mm_max_seqs, "-a",
          "--threads", cfg.threads, "-v", "1"], log)
    out = os.path.join(wd, "hits.m8")
    _run(["mmseqs", "convertalis", os.path.join(wd, "prof"), os.path.join(wd, "qdb"), res, out,
          "--format-output", "query,target,bits,evalue,qcov,tcov,fident", "--threads", cfg.threads,
          "-v", "1"], log)
    shutil.rmtree(tmpd, ignore_errors=True)
    h = _read_hits_by_name(out)
    h["n_profiles"] = len(targets)
    h["n_profile_pairs"] = int(len(pq))
    return h


def _read_hits_by_name(path):
    z = np.zeros(0)
    if os.path.getsize(path) == 0:
        return dict(q=z.astype(np.int64), t=z.astype(np.int64), bits=z, evalue=z, qcov=z, tcov=z)
    a = np.loadtxt(path, dtype=np.float64, ndmin=2)
    return dict(q=a[:, 0].astype(np.int64), t=a[:, 1].astype(np.int64), bits=a[:, 2], evalue=a[:, 3],
                qcov=a[:, 4], tcov=a[:, 5])


# ----------------------------------------------------------------------------
# attach rule
# ----------------------------------------------------------------------------

def attach_small(fam: np.ndarray, hits: list, lengths: np.ndarray, cfg: RescueConfig):
    """Attach each small family to the target family of its best link (see module doc).
    ``hits``: list of (source name, hit dict) with gene ids q/t, evalue, bits, qcov, tcov
    (fractions or percent-free fractions). Returns (new families, attach table dict, stats)."""
    sizes = np.bincount(fam)
    srcs, a_s, a_t, e_s, b_s = [], [], [], [], []
    L = lengths.astype(np.float64)
    for si, (name, h) in enumerate(hits):
        q, t = h["q"].astype(np.int64), h["t"].astype(np.int64)
        if len(q) == 0:
            continue
        cov_sh = np.maximum(h["qcov"] * L[q], h["tcov"] * L[t]) / np.minimum(L[q], L[t])
        ok = (h["evalue"] <= cfg.evalue) & (cov_sh >= cfg.cov_shorter) & (fam[q] != fam[t])
        q, t, e, b = q[ok], t[ok], h["evalue"][ok], h["bits"][ok]
        # orient: source = gene of the small family, dest = gene of a target family
        sq, st_ = sizes[fam[q]] <= cfg.small_max, sizes[fam[t]] <= cfg.small_max
        src = np.where(sq & ~st_, q, np.where(st_ & ~sq, t, -1))
        dst = np.where(sq & ~st_, t, np.where(st_ & ~sq, q, -1))
        m = src >= 0
        a_s.append(src[m]); a_t.append(dst[m]); e_s.append(e[m]); b_s.append(b[m])
        srcs.append(np.full(m.sum(), si, np.int64))
    stats = {"n_links": int(sum(len(x) for x in a_s))}
    if not a_s or stats["n_links"] == 0:
        return fam.copy(), {k: np.zeros(0, np.int64) for k in ("small", "target", "source")} | \
            {"evalue": np.zeros(0)}, dict(stats, n_attached=0, n_capped=0, n_genes_attached=0)
    s, t, e, b, so = (np.concatenate(x) for x in (a_s, a_t, e_s, b_s, srcs))
    fs, ft = fam[s], fam[t]
    o = np.lexsort((ft, -sizes[ft], -b, e, fs))        # per small family: best link first
    fs, ft, e, b, so = fs[o], ft[o], e[o], b[o], so[o]
    first = np.r_[True, fs[1:] != fs[:-1]]
    fs, ft, e, so = fs[first], ft[first], e[first], so[first]
    size = sizes.astype(np.int64).copy()
    new_label = np.arange(len(sizes))
    keep = np.zeros(len(fs), bool)
    for i in np.lexsort((fs, e)):                      # increasing e (ties: small label)
        if size[ft[i]] + sizes[fs[i]] > cfg.max_size:
            continue
        size[ft[i]] += sizes[fs[i]]
        new_label[fs[i]] = ft[i]
        keep[i] = True
    newfam = new_label[fam]
    att = dict(small=fs[keep], target=ft[keep], evalue=e[keep], source=so[keep])
    stats.update(n_attached=int(keep.sum()), n_capped=int((~keep).sum()),
                 n_genes_attached=int(sizes[fs[keep]].sum()),
                 n_target_families=int(len(np.unique(ft[keep]))),
                 by_source={hits[k][0]: int((so[keep] == k).sum()) for k in range(len(hits))})
    return newfam, att, stats


# ----------------------------------------------------------------------------
# re-scoring of the new pairs only
# ----------------------------------------------------------------------------

def rescore_new_pairs(gt: GeneTable, fam_old: np.ndarray, fam_new: np.ndarray, st: ScoreTable,
                      scfg: ScoringConfig, wd: str, log) -> tuple[ScoreTable, dict]:
    os.makedirs(wd, exist_ok=True)
    moved = np.flatnonzero(fam_old != fam_new)         # attached genes (labels of their family change)
    stats = {"n_attached_genes": int(len(moved))}
    if len(moved) == 0:
        return st, stats
    recv = np.unique(fam_new[moved])
    dbg = np.flatnonzero(np.isin(fam_new, recv))
    t0 = time.time()
    qfa, dfa = os.path.join(wd, "q.fa"), os.path.join(wd, "d.fa")
    write_fasta(gt, qfa, moved)
    write_fasta(gt, dfa, dbg)
    _run(["diamond", "makedb", "--in", dfa, "-d", dfa + ".dmnd", "--threads", scfg.threads, "--quiet"], log)
    out = os.path.join(wd, "new.m8")
    cmd = ["diamond", "blastp", "-q", qfa, "-d", dfa + ".dmnd", "-o", out, "-k", "0", "--max-hsps", "1",
           "-b", scfg.diamond_b, "-c", "1", "-e", scfg.evalue, "--threads", scfg.threads, "--quiet",
           "--outfmt", "6", "qseqid", "sseqid", "bitscore", "qcovhsp", "scovhsp", "pident", "evalue"]
    if scfg.diamond_sens:
        cmd.append(scfg.diamond_sens)
    _run(cmd, log)
    cols = []
    for c in _iter_m8(out, 7):
        q, t = c[0].astype(np.int64), c[1].astype(np.int64)
        keep = (fam_new[q] == fam_new[t]) & ((fam_old[q] != fam_old[t]) | (q == t))
        cols.append([x[keep] for x in c])
    c = [np.concatenate([x[i] for x in cols]) for i in range(7)] if cols else [np.zeros(0)] * 7
    c[0], c[1] = c[0].astype(np.int64), c[1].astype(np.int64)
    c[3], c[4], c[5] = c[3] / 100.0, c[4] / 100.0, c[5] / 100.0
    stats["diamond_s"] = round(time.time() - t0, 1)
    stats["n_diamond_new"] = int((c[0] != c[1]).sum())
    for f in (qfa, dfa, dfa + ".dmnd", out):
        os.remove(f)
    # MMseqs2 SW of missed new pairs in families <= rescue_max_family (hybrid S2 rule) + self pairs
    sizes_new = np.bincount(fam_new)
    M = gt.n + 1
    found = np.unique(np.minimum(c[0], c[1]) * M + np.maximum(c[0], c[1]))
    qs, ts = [moved], [moved]                           # self pairs (self-scores)
    small_recv = recv[sizes_new[recv] <= min(scfg.cap_all, scfg.rescue_max_family)]
    mem = {int(f): np.flatnonzero(fam_new == f) for f in small_recv}
    for g in moved.tolist():
        f = int(fam_new[g])
        if f in mem:
            o = mem[f]
            o = o[fam_old[o] != fam_old[g]]
            qs.append(np.full(len(o), g, np.int64))
            ts.append(o)
    q, t = np.concatenate(qs), np.concatenate(ts)
    q, t = np.minimum(q, t), np.maximum(q, t)
    k = np.unique(q * M + t)
    k = k[~np.isin(k, found) | ((k // M) == (k % M))]
    q, t = k // M, k % M
    t1 = time.time()
    os.makedirs(os.path.join(wd, "mm"), exist_ok=True)
    c2 = _mmseqs_align_pairs(gt, q, t, scfg, os.path.join(wd, "mm"), log, stats)
    shutil.rmtree(os.path.join(wd, "mm"), ignore_errors=True)
    stats["mmseqs_s"] = round(time.time() - t1, 1)
    c = [np.concatenate([a, np.asarray(b)]) for a, b in zip(c, c2)]
    new = _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), *c[2:], gt.n)
    # cap rule for families > cap_all: keep the attached gene's top hits per species
    big = sizes_new[fam_new[new.q]] > scfg.cap_all
    if big.any():
        isq = np.isin(new.q, moved)
        anc = np.where(isq, new.q, new.t).astype(np.int64)
        oth = np.where(isq, new.t, new.q).astype(np.int64)
        sp_o = gt.species_of[oth].astype(np.int64)
        same = sp_o == gt.species_of[anc]
        idx = np.flatnonzero(big)
        o = idx[np.lexsort((-new.bits[idx], sp_o[idx], anc[idx]))]
        st0 = np.r_[True, (anc[o][1:] != anc[o][:-1]) | (sp_o[o][1:] != sp_o[o][:-1])]
        rank = np.arange(len(o)) - np.flatnonzero(st0)[np.cumsum(st0) - 1]
        keep = np.ones(len(new.q), bool)
        keep[o] = np.where(same[o], rank < scfg.top_same, rank < scfg.top_per_species)
        stats["n_new_pairs_dropped_by_cap"] = int((~keep).sum())
        new = ScoreTable(*(getattr(new, x)[keep] for x in ("q", "t", "bits", "qcov", "tcov", "pident",
                                                           "evalue")), new.self_bits)
    sb = st.self_bits.copy()
    fill = ~np.isfinite(sb) & np.isfinite(new.self_bits)
    sb[fill] = new.self_bits[fill]
    out_st = ScoreTable(*(np.concatenate([getattr(st, x), getattr(new, x).astype(getattr(st, x).dtype)])
                          for x in ("q", "t", "bits", "qcov", "tcov", "pident", "evalue")), sb)
    stats.update(n_new_pairs=int(len(new.q)), n_self_filled=int(fill.sum()),
                 rescore_s=round(time.time() - t0, 1))
    return out_st, stats


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------

def run_rescue(gt: GeneTable, fam: np.ndarray, st: ScoreTable, cfg: RescueConfig, wd: str,
               s0_dir: str | None = None, threedi: str | None = None,
               scoring: ScoringConfig | None = None, log_path: str | None = None,
               hits_cache: str | None = None, input_key: str = ""):
    """fam (S2 final families) + scores -> (new families, new ScoreTable, attach table, stats).
    ``gt.seqs`` must be loaded. ``threedi``: 3Di FASTA (modes 3di/both). ``hits_cache``: optional
    directory; search hits are stored as ``<dir>/<source>_<hash>.npz`` (hash over ``input_key``
    -- the caller's key of fam/scores --, the source's search fields and, for 3di, the 3Di
    file md5) and reused when present, so arms that share a search do not repeat it."""
    if cfg.mode not in RESCUE_MODES:
        raise ValueError(f"unknown rescue mode {cfg.mode!r}")
    if cfg.mode == "none":
        return fam, st, None, {"mode": "none"}
    if cfg.mode in ("3di", "both") and not threedi:
        raise ValueError("family rescue mode %r needs a 3Di FASTA (--structure-3di)" % cfg.mode)
    os.makedirs(wd, exist_ok=True)
    scoring = scoring or ScoringConfig()
    scoring.threads = cfg.threads
    stats = {"mode": cfg.mode, "config": asdict(cfg)}
    t0 = time.time()
    with open(log_path or os.path.join(wd, "commands.log"), "a") as log:
        rep_of = stage1_rep_of(gt, s0_dir)
        queries, targets, sizes = select_units(fam, rep_of, gt.lengths, cfg.small_max, cfg.reps_per_family)
        stats.update(n_queries=int(len(queries)), n_targets=int(len(targets)),
                     n_small_families=int(((sizes > 0) & (sizes <= cfg.small_max)).sum()),
                     n_target_families=int((sizes > cfg.small_max).sum()),
                     units_s=round(time.time() - t0, 1))
        hits = []
        def cached(src, fn):
            p = None
            if hits_cache:
                import hashlib
                import json
                fields = ["small_max", "reps_per_family", "evalue"] + (
                    ["fs_sens", "fs_max_seqs"] if src == "3di" else ["mm_sens", "mm_max_seqs", "profile_members"])
                kk = {"input": input_key, "src": src, **{f: getattr(cfg, f) for f in fields}}
                if src == "3di":
                    h5 = hashlib.md5()
                    with open(threedi, "rb") as fh:
                        for chunk in iter(lambda: fh.read(1 << 22), b""):
                            h5.update(chunk)
                    kk["3di"] = h5.hexdigest()
                hh = hashlib.sha256(json.dumps(kk, sort_keys=True).encode()).hexdigest()[:16]
                p = os.path.join(hits_cache, f"{src}_{hh}.npz")
            if p and os.path.exists(p):
                z = np.load(p)
                stats[f"{src}_from_cache"] = p
                return {k: z[k] for k in z.files}
            t1 = time.time()
            h = fn()
            h = {k: np.asarray(v) for k, v in h.items()}
            h["search_s"] = np.asarray(round(time.time() - t1, 1))
            if p:
                os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
                np.savez(p + ".tmp.npz", **h)
                os.replace(p + ".tmp.npz", p)
            return h
        if cfg.mode in ("3di", "both"):
            def f3():
                ss = load_3di_any(threedi, gt)
                h = search_3di(gt, queries, targets, ss, os.path.join(wd, "fs"), cfg, log)
                shutil.rmtree(os.path.join(wd, "fs"), ignore_errors=True)
                h["n_3di_genes"] = len(ss)
                h.pop("pident", None)
                return h
            h = cached("3di", f3)
            stats.update({k: (h.pop(k).item()) for k in ("n_3di_genes", "n_queries_3di", "n_targets_3di")})
            stats.update(n_3di_hits=int(len(h["q"])), search_3di_s=float(h.pop("search_s")))
            hits.append(("3di", h))
        if cfg.mode in ("profile", "both"):
            def fp():
                h = search_profiles(gt, fam, rep_of, queries, targets, os.path.join(wd, "prof"), cfg, log)
                shutil.rmtree(os.path.join(wd, "prof"), ignore_errors=True)
                return h
            h = cached("profile", fp)
            stats.update({k: (h.pop(k).item()) for k in ("n_profiles", "n_profile_pairs")})
            stats.update(n_profile_hits=int(len(h["q"])), search_profile_s=float(h.pop("search_s")))
            hits.append(("profile", h))
        fam_new, att, ast = attach_small(fam, hits, gt.lengths, cfg)
        stats["attach"] = ast
        st_new, rst = rescore_new_pairs(gt, fam, fam_new, st, scoring, os.path.join(wd, "rescore"), log)
        stats["rescore"] = rst
    shutil.rmtree(os.path.join(wd, "rescore"), ignore_errors=True)
    stats["wall_s"] = round(time.time() - t0, 1)
    return fam_new, st_new, att, stats


def write_families(gt: GeneTable, fam: np.ndarray, path: str) -> None:
    """families.tsv (gene id, family label F<id>) in the S0 format."""
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("#gene_id\tfamily_id\n")
        for nm, f in zip(gt.names, fam.tolist()):
            fh.write(f"{nm}\tF{f:07d}\n")
    os.replace(tmp, path)


# ----------------------------------------------------------------------------
# merge-side guard (R4)
# ----------------------------------------------------------------------------

def _norm_scores(gt: GeneTable, st: ScoreTable) -> np.ndarray:
    sb = st.self_bits_filled(gt.lengths)
    return st.bits.astype(np.float64) / np.maximum(sb[st.q], sb[st.t])


def guard_detach(res, gt: GeneTable, tree, fam_old: np.ndarray, fam_new: np.ndarray,
                 st: ScoreTable) -> tuple[np.ndarray, dict]:
    """Small families (old labels) whose attachment is NOT confirmed by the merge: no gene
    of the small family is predicted orthologous to an original member of the target family
    with a reciprocal best hit (norm_max) between the two species within the new family.
    Returns (array of detached old family labels, stats)."""
    from .merge import to_orthoxml
    from .orthoxml import Group, OrthoXML, hogs_to_pairs
    moved = fam_old != fam_new
    if not moved.any():
        return np.zeros(0, np.int64), {"n_attached": 0, "n_detached": 0}
    att_fams = np.unique(fam_old[moved])
    ox = to_orthoxml(res, gt, tree)
    mnames = {gt.names[g] for g in np.flatnonzero(moved)}
    keep = [g for g in ox.groups if isinstance(g, Group) and mnames & set(g.genes())]
    idx = gt.index()
    pred = hogs_to_pairs(OrthoXML(genes=ox.genes, groups=keep, species_tree=ox.species_tree))
    cand = set()
    for a, b in pred:
        ia, ib = idx[a], idx[b]
        if moved[ia] and not moved[ib]:
            cand.add((ia, ib))
        elif moved[ib] and not moved[ia]:
            cand.add((ib, ia))
    # reciprocal best hit within the new family (norm_max), per species pair
    S = _norm_scores(gt, st)
    q, t = st.q.astype(np.int64), st.t.astype(np.int64)
    sp = gt.species_of.astype(np.int64)
    nsp = int(sp.max()) + 1
    inv = (fam_new[q] == fam_new[t])
    rel = inv & (moved[q] | moved[t] | np.isin(q, [c[1] for c in cand]) | np.isin(t, [c[1] for c in cand]))
    rq, rt, rs = q[rel], t[rel], S[rel]
    A = np.concatenate([rq, rt]); B = np.concatenate([rt, rq]); W = np.concatenate([rs, rs])
    key = A * nsp + sp[B]
    o = np.lexsort((B, -W, key))
    first = np.r_[True, key[o][1:] != key[o][:-1]]
    best = dict(zip(key[o][first].tolist(), B[o][first].tolist()))
    ok_fams = set()
    for s_, t_ in cand:
        if best.get(s_ * nsp + sp[t_]) == t_ and best.get(t_ * nsp + sp[s_]) == s_:
            ok_fams.add(int(fam_old[s_]))
    det = np.array(sorted(set(att_fams.tolist()) - ok_fams), dtype=np.int64)
    return det, {"n_attached": int(len(att_fams)), "n_confirmed": len(ok_fams),
                 "n_detached": int(len(det)), "n_candidate_pairs": len(cand),
                 "n_genes_detached": int(np.isin(fam_old, det).sum())}


def revert_families(fam_old: np.ndarray, fam_new: np.ndarray, st: ScoreTable, detached: np.ndarray):
    """Give detached small families their old labels back and drop their cross-family pairs."""
    fam = fam_new.copy()
    m = np.isin(fam_old, detached)
    # old labels of small families are never reused as targets, so they are free
    fam[m] = fam_old[m]
    keep = fam[st.q] == fam[st.t]
    st2 = ScoreTable(*(getattr(st, x)[keep] for x in ("q", "t", "bits", "qcov", "tcov", "pident",
                                                      "evalue")), st.self_bits)
    return fam, st2
