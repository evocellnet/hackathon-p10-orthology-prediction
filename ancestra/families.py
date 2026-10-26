"""S0 -- de novo family partition of all input proteins (DESIGN.md §1, §3 row S0).

Near-linear clustering of every protein into homologous families, without a
reference database and without an all-vs-all search over all N proteins.

Backends
--------
``diamond``  DIAMOND ``cluster`` / ``deepclust`` / ``linclust`` (cascaded, greedy
             vertex cover).  Knobs: ``approx_id``, ``member_cover``, optional
             sensitivity flag and round schedules.
``mmseqs``   MMseqs2 ``cluster`` (cascaded) or ``linclust``. Knobs: ``min_seq_id``,
             ``cov`` + ``cov_mode``, ``sens`` (-s), ``cluster_mode``.
two-stage    any of the above as stage 1, then a representative-only
             all-vs-all (DIAMOND blastp on the R << N stage-1 representatives),
             filtered by e-value and coverage, merged by connected components
             (``cc``) or MCL (``mcl``) on the reps graph.  Families = union of the
             stage-1 clusters whose reps fall in the same reps-graph component.

Output: TSV ``gene_id<TAB>family_id`` (every input gene, singletons included;
family ids ``F0000001``... ordered by decreasing size, ties by smallest gene id)
plus a stats JSON (size distribution, Σ g_f², top-10 sizes).

CLI
---
    python -m ancestra.families cluster --fasta in.fa --outdir D --config cfg.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from typing import Iterable

__all__ = [
    "UnionFind", "read_fasta_ids", "read_pairs_tsv", "partition_from_pairs",
    "label_families", "family_stats", "read_families_tsv", "write_families_tsv",
    "FamilyConfig", "run_families",
]

# --------------------------------------------------------------------------
# generic helpers (pure python, unit-tested in tests/test_families.py)
# --------------------------------------------------------------------------


class UnionFind:
    """Union-find over arbitrary hashable keys (path halving, union by size)."""

    def __init__(self, items: Iterable = ()):
        self.parent: dict = {}
        self.size: dict = {}
        for x in items:
            self.add(x)

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.size[x] = 1

    def find(self, x):
        p = self.parent
        if x not in p:
            self.add(x)
            return x
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        return ra

    def groups(self) -> dict:
        out = defaultdict(list)
        for x in self.parent:
            out[self.find(x)].append(x)
        return out


def read_fasta_ids(path: str) -> list[str]:
    """Return sequence ids (first whitespace token of each header) in file order."""
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                tok = line[1:].split(None, 1)
                ids.append(tok[0] if tok else "")
    return ids


def _canonical_fasta(fasta: str, outdir: str):
    """Write ``outdir/canonical_input.fa``: sequences sorted by (length desc, sequence,
    original id), renamed ``c0000000``.. by rank. Returns (path, {canonical id: original id})."""
    recs, name, buf = [], None, []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                tok = line[1:].split(None, 1)
                name, buf = (tok[0] if tok else ""), []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    # same normalisation as ancestra.io.load_proteomes (terminal '*' stripped, internal '*' -> X)
    recs = [(n, sq.rstrip("*").replace("*", "X").upper()) for n, sq in recs]
    recs.sort(key=lambda r: (-len(r[1]), r[1], r[0]))
    width = max(7, len(str(len(recs))))
    out = os.path.join(outdir, "canonical_input.fa")
    c2o = {}
    with open(out, "w") as fh:
        for i, (nm, sq) in enumerate(recs):
            cid = f"c{i:0{width}d}"
            c2o[cid] = nm
            fh.write(f">{cid}\n{sq}\n")
    return out, c2o


def read_pairs_tsv(path: str, cols=(0, 1)) -> list[tuple[str, str]]:
    """Read 2-column (representative, member) cluster TSVs (DIAMOND / MMseqs createtsv).

    Lines starting with '#' and blank lines are skipped. Extra columns ignored.
    """
    a, b = cols
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            out.append((t[a], t[b]))
    return out


def partition_from_pairs(gene_ids: Iterable[str], pairs: Iterable[tuple[str, str]],
                         strict: bool = True) -> dict[str, str]:
    """Connected components of ``pairs`` over ``gene_ids`` -> {gene: component root}.

    Genes absent from all pairs become singletons. With ``strict``, a pair that
    names an unknown gene raises ``KeyError`` (catches id-parsing bugs).
    """
    uf = UnionFind(gene_ids)
    known = uf.parent
    for x, y in pairs:
        if strict and (x not in known or y not in known):
            raise KeyError(f"pair ({x!r}, {y!r}) names a gene absent from the input")
        uf.union(x, y)
    return {g: uf.find(g) for g in known}


def label_families(gene2root: dict[str, str], prefix: str = "F") -> dict[str, str]:
    """Relabel components deterministically: largest first, ties by smallest gene id."""
    members = defaultdict(list)
    for g, r in gene2root.items():
        members[r].append(g)
    order = sorted(members.values(), key=lambda m: (-len(m), min(m)))
    width = max(7, len(str(len(order))))
    out = {}
    for i, m in enumerate(order, 1):
        fid = f"{prefix}{i:0{width}d}"
        for g in m:
            out[g] = fid
    return out


def family_stats(gene2fam: dict[str, str], top: int = 10) -> dict:
    sizes = sorted(Counter(gene2fam.values()).values(), reverse=True)
    n = sum(sizes)
    sum_sq = sum(s * s for s in sizes)
    multi = [s for s in sizes if s > 1]
    return {
        "n_genes": n,
        "n_families": len(sizes),
        "n_singletons": sum(1 for s in sizes if s == 1),
        "n_families_ge2": len(multi),
        "frac_genes_in_singletons": (sum(1 for s in sizes if s == 1) / n) if n else 0.0,
        "sum_g2": sum_sq,
        "sum_g2_over_N2": sum_sq / (n * n) if n else 0.0,
        "max_family": sizes[0] if sizes else 0,
        "max_family_over_N": sizes[0] / n if n else 0.0,
        "top_sizes": sizes[:top],
        "size_quantiles_ge2": _quantiles(multi, (0.5, 0.9, 0.99, 0.999)),
    }


def _quantiles(xs, qs):
    if not xs:
        return {}
    s = sorted(xs)
    return {str(q): s[min(len(s) - 1, int(q * len(s)))] for q in qs}


def write_families_tsv(gene2fam: dict[str, str], path: str):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("#gene_id\tfamily_id\n")
        for g, f in sorted(gene2fam.items(), key=lambda kv: (kv[1], kv[0])):
            fh.write(f"{g}\t{f}\n")
    os.replace(tmp, path)


def read_families_tsv(path: str) -> dict[str, str]:
    return dict(read_pairs_tsv(path))


# --------------------------------------------------------------------------
# backends
# --------------------------------------------------------------------------


@dataclass
class FamilyConfig:
    name: str = "default"
    backend: str = "diamond"          # diamond | mmseqs
    # diamond
    dmnd_cmd: str = "cluster"         # cluster | deepclust | linclust
    approx_id: float = 30.0
    member_cover: float = 50.0
    mutual_cover: float | None = None
    dmnd_sens: str | None = None      # e.g. "--sensitive" (passed verbatim if set)
    dmnd_extra: list = field(default_factory=list)
    # mmseqs
    mm_cmd: str = "cluster"           # cluster | linclust
    min_seq_id: float = 0.0
    cov: float = 0.5
    cov_mode: int = 1
    sens: float = 4.0
    cluster_mode: int = 0
    evalue: float = 1e-3
    mm_extra: list = field(default_factory=list)
    # two-stage
    stage2: str | None = None         # None | cc | mcl
    s2_sens: str = "--sensitive"
    s2_evalue: float = 1e-5
    s2_cov: float = 0.0               # coverage threshold (fraction 0..1) for reps edges
    s2_cov_mode: str = "min"          # min: both qcov,scov >= c ; max: either ; shorter: cov of shorter seq
    s2_max_target: int = 500
    s2_block: float = 4.0             # DIAMOND -b for the stage-2 search (pinned: thread-independent)
    s2_chunks: int = 1                # DIAMOND -c
    s2_bitscore: float = 0.0
    s2_query: str = "reps"            # reps: reps vs reps ; all: every gene vs reps (N x R)
    s2_gmax: float | None = None      # split connected components larger than this many genes with Leiden;
                                      # a value in (0, 1) means that fraction of N (e.g. 0.03 -> 3% of N)
    s2_reassign: bool = False         # after a Leiden split: move each gene of a split part to the part
                                      # holding its best-scoring (filter-passing) rep hit
    mcl_I: float = 1.4
    s1_reuse: str | None = None       # outdir of an earlier run with identical stage-1 params: reuse its clusters
    seed: int = 11
    # resources
    threads: int = 24
    memory_gb: int = 64

    def hash(self) -> str:
        d = asdict(self)
        d.pop("name"); d.pop("threads"); d.pop("memory_gb"); d.pop("s1_reuse")
        return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()[:12]


def _run(cmd: list[str], log, cwd=None):
    log.write("$ " + " ".join(map(str, cmd)) + "\n"); log.flush()
    r = subprocess.run([str(c) for c in cmd], stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}")


def _diamond_cluster(fasta, outdir, cfg: FamilyConfig, log) -> str:
    out = os.path.join(outdir, "stage1_clusters.tsv")
    tmp = os.path.join(outdir, "tmp"); os.makedirs(tmp, exist_ok=True)
    cmd = ["diamond", cfg.dmnd_cmd, "-d", fasta, "-o", out, "--threads", cfg.threads,
           "-M", f"{cfg.memory_gb}G", "--tmpdir", tmp, "--approx-id", cfg.approx_id,
           "--member-cover", cfg.member_cover]
    if cfg.mutual_cover is not None:
        cmd += ["--mutual-cover", cfg.mutual_cover]
    if cfg.dmnd_sens:
        cmd += [cfg.dmnd_sens]
    cmd += list(cfg.dmnd_extra)
    _run(cmd, log)
    shutil.rmtree(tmp, ignore_errors=True)
    return out


def _mmseqs_cluster(fasta, outdir, cfg: FamilyConfig, log) -> str:
    db = os.path.join(outdir, "mm", "db"); os.makedirs(os.path.dirname(db), exist_ok=True)
    clu = os.path.join(outdir, "mm", "clu")
    tmp = os.path.join(outdir, "mm", "tmp")
    out = os.path.join(outdir, "stage1_clusters.tsv")
    _run(["mmseqs", "createdb", fasta, db, "-v", "1"], log)
    cmd = ["mmseqs", cfg.mm_cmd, db, clu, tmp, "--threads", cfg.threads, "--min-seq-id", cfg.min_seq_id,
           "-c", cfg.cov, "--cov-mode", cfg.cov_mode, "--cluster-mode", cfg.cluster_mode, "-e", cfg.evalue,
           "-v", "2"]
    if cfg.mm_cmd == "cluster":
        cmd += ["-s", cfg.sens]
    cmd += list(cfg.mm_extra)
    _run(cmd, log)
    _run(["mmseqs", "createtsv", db, db, clu, out, "--threads", cfg.threads, "-v", "1"], log)
    shutil.rmtree(os.path.join(outdir, "mm"), ignore_errors=True)
    return out


def _write_subset_fasta(fasta, keep: set[str], out):
    w = False
    with open(fasta) as fi, open(out, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                w = line[1:].split(None, 1)[0] in keep
            if w:
                fo.write(line)


def _passes_cov(qcov, scov, qlen, slen, c, mode):
    if c <= 0:
        return True
    if mode == "min":
        return qcov >= c and scov >= c
    if mode == "max":
        return qcov >= c or scov >= c
    if mode == "shorter":
        return (qcov if qlen <= slen else scov) >= c
    raise ValueError(mode)


def filter_hits(rows: Iterable[list[str]], cfg: FamilyConfig) -> list[tuple[str, str, float]]:
    """rows: qseqid sseqid evalue bitscore qlen slen qstart qend sstart send.

    Returns undirected edges (a<b, bitscore) passing e-value, bitscore and coverage filters.
    """
    best = {}
    for t in rows:
        q, s = t[0], t[1]
        if q == s:
            continue
        ev, bs = float(t[2]), float(t[3])
        ql, sl = int(t[4]), int(t[5])
        qcov = (abs(int(t[7]) - int(t[6])) + 1) / ql
        scov = (abs(int(t[9]) - int(t[8])) + 1) / sl
        if ev > cfg.s2_evalue or bs < cfg.s2_bitscore:
            continue
        if not _passes_cov(qcov, scov, ql, sl, cfg.s2_cov, cfg.s2_cov_mode):
            continue
        k = (q, s) if q < s else (s, q)
        if bs > best.get(k, -1.0):
            best[k] = bs
    return [(a, b, w) for (a, b), w in best.items()]


def best_hits(rows: Iterable[list[str]], cfg: FamilyConfig) -> dict[str, tuple[str, float]]:
    """Directed best hit per query (highest bitscore; ties -> smaller subject id) among rows that pass
    the same e-value / bitscore / coverage filters as :func:`filter_hits`. Self hits are ignored."""
    best = {}
    for t in rows:
        q, s = t[0], t[1]
        if q == s:
            continue
        ev, bs = float(t[2]), float(t[3])
        ql, sl = int(t[4]), int(t[5])
        if ev > cfg.s2_evalue or bs < cfg.s2_bitscore:
            continue
        qcov = (abs(int(t[7]) - int(t[6])) + 1) / ql
        scov = (abs(int(t[9]) - int(t[8])) + 1) / sl
        if not _passes_cov(qcov, scov, ql, sl, cfg.s2_cov, cfg.s2_cov_mode):
            continue
        cur = best.get(q)
        if cur is None or bs > cur[1] or (bs == cur[1] and s < cur[0]):
            best[q] = (s, bs)
    return best


def _mcl(edges, outdir, cfg: FamilyConfig, log) -> list[tuple[str, str]]:
    abc = os.path.join(outdir, "reps_graph.abc")
    with open(abc, "w") as fh:
        for a, b, w in edges:
            fh.write(f"{a}\t{b}\t{w}\n")
    out = os.path.join(outdir, "reps_mcl.txt")
    _run(["mcl", abc, "--abc", "-I", cfg.mcl_I, "-te", cfg.threads, "-o", out], log)
    pairs = []
    for line in open(out):
        t = line.split()
        pairs += [(t[0], x) for x in t[1:]]
    return pairs


STAGE1_FIELDS = ("backend", "dmnd_cmd", "approx_id", "member_cover", "mutual_cover", "dmnd_sens", "dmnd_extra",
                 "mm_cmd", "min_seq_id", "cov", "cov_mode", "sens", "cluster_mode", "evalue", "mm_extra")


def leiden_split(nodes: list, edges: list, node_w: dict, gmax: int, seed: int = 11):
    """Recursively split a connected node set whose total weight exceeds ``gmax``.

    Uses Leiden (RBConfiguration) with the smallest resolution in a fixed ladder that
    yields more than one community, so splits are as coarse as possible; parts still
    above ``gmax`` are split again. Returns a list of node lists.
    """
    import igraph as ig
    import leidenalg
    ladder = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0)
    out, queue = [], [list(nodes)]
    while queue:
        ns = queue.pop()
        if len(ns) == 1 or sum(node_w[n] for n in ns) <= gmax:
            out.append(ns)
            continue
        idx = {n: i for i, n in enumerate(ns)}
        es = [(idx[a], idx[b], w) for a, b, w in edges if a in idx and b in idx]
        g = ig.Graph(n=len(ns), edges=[(a, b) for a, b, _ in es])
        g.es["weight"] = [w for _, _, w in es]
        parts = None
        for gamma in ladder:
            p = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights="weight",
                                         resolution_parameter=gamma, seed=seed)
            if len(p) > 1:
                parts = [[ns[i] for i in c] for c in p]
                break
        if parts is None:  # cannot split further (e.g. a clique): keep as is
            out.append(ns)
            continue
        queue.extend(parts)
    return out


def run_families(fasta: str, outdir: str, cfg: FamilyConfig) -> dict:
    """Run S0 with ``cfg``; write families.tsv + stats.json in ``outdir``; return stats."""
    os.makedirs(outdir, exist_ok=True)
    fasta = os.path.abspath(fasta)
    t0 = time.time()
    timings = {}
    orig_ids = read_fasta_ids(fasta)
    if len(set(orig_ids)) != len(orig_ids):
        raise ValueError("duplicate sequence ids in input FASTA")
    # Canonical input (determinism): DIAMOND's greedy clustering breaks ties by input order
    # and IDs, so the same sequences in another order / with other names gave a different
    # partition (ARI 0.78 on OrthoBench). We cluster a canonical copy -- sorted by
    # (length desc, sequence, original id), IDs = canonical rank -- and map back at the end.
    fasta, canon2orig = _canonical_fasta(fasta, outdir)
    gene_ids = sorted(canon2orig)
    flagged = {}
    with open(os.path.join(outdir, "commands.log"), "a") as log:
        t = time.time()
        if cfg.s1_reuse:
            prev = json.load(open(os.path.join(cfg.s1_reuse, "stats.json")))
            for k in STAGE1_FIELDS:
                if prev["config"].get(k) != getattr(cfg, k):
                    raise ValueError(f"s1_reuse: stage-1 field {k} differs from {cfg.s1_reuse}")
            if not prev.get("canonical_input"):
                raise ValueError(f"s1_reuse: {cfg.s1_reuse} was not run on canonical input (pre-fix S0)")
            s1 = os.path.join(outdir, "stage1_clusters.tsv")
            shutil.copyfile(os.path.join(cfg.s1_reuse, "stage1_clusters.tsv"), s1)
            timings["stage1_reused_s"] = prev["timings"].get("stage1_s", prev["timings"].get("stage1_reused_s"))
            log.write(f"# stage 1 reused from {cfg.s1_reuse}\n")
        elif cfg.backend == "diamond":
            s1 = _diamond_cluster(fasta, outdir, cfg, log)
        elif cfg.backend == "mmseqs":
            s1 = _mmseqs_cluster(fasta, outdir, cfg, log)
        else:
            raise ValueError(cfg.backend)
        timings["stage1_s"] = time.time() - t
        pairs = read_pairs_tsv(s1)
        g2r = partition_from_pairs(gene_ids, pairs)
        n_stage1 = len(set(g2r.values()))
        extra = {"n_stage1_clusters": n_stage1}
        if cfg.stage2:
            t = time.time()
            # use the tool's own representative (first column) where available
            rep_of = {m: r for r, m in pairs}
            root2rep = {}
            for g in gene_ids:
                if rep_of.get(g) == g:
                    root2rep[g2r[g]] = g
            for g in gene_ids:  # clusters without explicit rep (singletons)
                root2rep.setdefault(g2r[g], g2r[g])
            rep_set = set(root2rep.values())
            gene2rep = {g: root2rep[g2r[g]] for g in gene_ids}
            rep_size = Counter(gene2rep.values())
            rfa = os.path.join(outdir, "reps.fa")
            _write_subset_fasta(fasta, rep_set, rfa)
            dmnd = os.path.join(outdir, "reps.dmnd")
            hits = os.path.join(outdir, "reps_hits.tsv")
            query = rfa if cfg.s2_query == "reps" else fasta
            if cfg.s2_query not in ("reps", "all"):
                raise ValueError(cfg.s2_query)
            _run(["diamond", "makedb", "--in", rfa, "-d", dmnd, "--threads", cfg.threads], log)
            # -b/-c pinned: DIAMOND's automatic block/chunk sizing depends on --threads and
            # changed ~0.02% of hits between 12 and 24 threads (-> different families)
            cmd = ["diamond", "blastp", "-q", query, "-d", dmnd, "-o", hits, "--threads", cfg.threads,
                   "-b", cfg.s2_block, "-c", cfg.s2_chunks,
                   "-k", cfg.s2_max_target, "-e", cfg.s2_evalue,
                   "--outfmt", "6", "qseqid", "sseqid", "evalue", "bitscore", "qlen", "slen",
                   "qstart", "qend", "sstart", "send"]
            if cfg.s2_sens:
                cmd += [cfg.s2_sens]
            _run(cmd, log)
            timings["stage2_search_s"] = time.time() - t
            with open(hits) as fh:
                gene_edges = filter_hits((l.rstrip("\n").split("\t") for l in fh), cfg)
            best = {}
            for a, b, w in gene_edges:   # lift gene-level edges to stage-1 cluster (rep) edges
                ra, rb = gene2rep[a], gene2rep[b]
                if ra == rb:
                    continue
                k = (ra, rb) if ra < rb else (rb, ra)
                if w > best.get(k, -1.0):
                    best[k] = w
            edges = sorted((a, b, w) for (a, b), w in best.items())   # order-independent
            if cfg.stage2 == "cc":
                rep_pairs = [(a, b) for a, b, _ in edges]
            elif cfg.stage2 == "mcl":
                rep_pairs = _mcl(edges, outdir, cfg, log)
            else:
                raise ValueError(cfg.stage2)
            rep2root = partition_from_pairs(sorted(rep_set), rep_pairs)
            n_split = 0
            gmax = cfg.s2_gmax * len(gene_ids) if cfg.s2_gmax and cfg.s2_gmax < 1 else cfg.s2_gmax
            extra["s2_gmax_genes"] = gmax
            if gmax:
                comp = defaultdict(list)
                for r, root in rep2root.items():
                    comp[root].append(r)
                big = {root for root, rs in comp.items() if sum(rep_size[r] for r in rs) > gmax}
                if big:
                    by_comp = defaultdict(list)
                    for a, b, w in edges:   # sorted -> Leiden sees a fixed edge order
                        if rep2root[a] in big and rep2root[a] == rep2root[b]:
                            by_comp[rep2root[a]].append((a, b, w))
                    for root in sorted(big):
                        parts = leiden_split(sorted(comp[root]), by_comp[root], rep_size, gmax, cfg.seed)
                        n_split += 1
                        for part in parts:
                            proot = min(part)
                            for r in part:
                                rep2root[r] = proot
                                flagged[r] = root
            g2r = {g: rep2root[gene2rep[g]] for g in gene_ids}
            n_moved = 0
            if cfg.s2_reassign and flagged:
                with open(hits) as fh:
                    bh = best_hits((l.rstrip("\n").split("\t") for l in fh), cfg)
                new = dict(g2r)
                for g, (s_, _) in bh.items():
                    if gene2rep[g] in flagged and s_ in flagged and g2r[s_] != g2r[g]:
                        new[g] = g2r[s_]
                        n_moved += 1
                g2r = new
            extra["n_genes_reassigned"] = n_moved
            timings["stage2_s"] = time.time() - t
            extra.update(n_reps=len(rep_set), n_gene_edges=len(gene_edges), n_rep_edges=len(edges),
                         n_components_split=n_split)
            for f in (dmnd, rfa, os.path.join(outdir, "reps_graph.abc")):
                if os.path.exists(f):
                    os.remove(f)
        gene2fam = label_families(g2r)
    canon_fam = gene2fam
    gene2fam = {canon2orig[g]: f for g, f in canon_fam.items()}
    write_families_tsv(gene2fam, os.path.join(outdir, "families.tsv"))
    os.remove(fasta)
    if flagged:
        fam_from = {}
        for g, f in canon_fam.items():
            r = gene2rep[g]
            if r in flagged:
                fam_from.setdefault(f, flagged[r])
            elif g in flagged:
                fam_from.setdefault(f, flagged[g])
        with open(os.path.join(outdir, "families_split_flags.tsv"), "w") as fh:
            fh.write("#family_id\tsplit_from_component\n")
            for f in sorted(fam_from):
                fh.write(f"{f}\t{fam_from[f]}\n")
        extra["n_split_families"] = len(fam_from)
        extra["genes_in_split_families"] = sum(1 for g, f in gene2fam.items() if f in fam_from)
    timings["total_s"] = time.time() - t0
    stats = family_stats(gene2fam)
    stats.update(extra)
    stats["timings"] = timings
    stats["config"] = asdict(cfg)
    stats["config_hash"] = cfg.hash()
    stats["canonical_input"] = True
    with open(os.path.join(outdir, "stats.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    return stats


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ancestra.families", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cluster", help="partition a FASTA into families")
    c.add_argument("--fasta", required=True)
    c.add_argument("--outdir", required=True)
    c.add_argument("--config", help="JSON file with FamilyConfig fields")
    c.add_argument("--set", nargs="*", default=[], help="override field=value (JSON-parsed)")
    a = ap.parse_args(argv)
    if a.cmd == "cluster":
        d = json.load(open(a.config)) if a.config else {}
        for kv in a.set:
            k, v = kv.split("=", 1)
            try:
                d[k] = json.loads(v)
            except json.JSONDecodeError:
                d[k] = v
        cfg = FamilyConfig(**d)
        st = run_families(a.fasta, a.outdir, cfg)
        print(json.dumps({k: st[k] for k in ("n_families", "max_family", "sum_g2_over_N2", "timings")}))


if __name__ == "__main__":
    sys.exit(main())
