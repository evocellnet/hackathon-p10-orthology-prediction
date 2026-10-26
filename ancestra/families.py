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

Evaluation (``evaluate`` sub-command / :func:`evaluate_refogs`) scores a
partition against OrthoBench RefOGs restricted to one split subset (dev or
test) -- see docs/S0_families.md for the metric definitions.

CLI
---
    python -m ancestra.families cluster --fasta in.fa --outdir D --config cfg.json
    python -m ancestra.families evaluate --families D/families.tsv --subset dev
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
    "read_refogs", "evaluate_refogs", "FamilyConfig", "run_families",
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
# RefOG evaluation
# --------------------------------------------------------------------------


def read_refogs(refog_dir: str, ids: Iterable[str] | None = None,
                remove_uncertain: bool = True) -> dict[str, set[str]]:
    """Read OrthoBench RefOG gene sets (optionally minus low_certainty_assignments)."""
    if ids is None:
        ids = sorted(f[:-4] for f in os.listdir(refog_dir) if f.startswith("RefOG") and f.endswith(".txt"))
    out = {}
    for rid in ids:
        genes = {l.strip() for l in open(os.path.join(refog_dir, rid + ".txt")) if l.strip()}
        unc = os.path.join(refog_dir, "low_certainty_assignments", rid + ".txt")
        if remove_uncertain and os.path.exists(unc):
            genes -= {l.strip() for l in open(unc) if l.strip()}
        out[rid] = genes
    return out


def evaluate_refogs(gene2fam: dict[str, str], refogs: dict[str, set[str]]) -> dict:
    """Score a family partition against reference groups.

    Per RefOG R (genes missing from ``gene2fam`` count as singleton families):
      complete     -- all of R lies in one family
      n_fam        -- number of distinct families R is spread over (fission; 1 = none)
      best_frac    -- largest fraction of R in a single family
      best_fam_size-- size of that family; size_ratio = best_fam_size / |R|
      foreign_refog_frac -- fraction of the best family's genes that belong to *other*
                       RefOGs in ``refogs`` (only the RefOGs passed in -- dev only during tuning)
      n_refogs_in_best -- number of RefOGs (of those passed) touching the best family
    Aggregate: completeness (= mean complete), mean / median fission, mean best_frac,
    median / max size_ratio, fusion rate (fraction of RefOGs whose best family also
    contains genes of another RefOG), and the OrthoBench-style pairwise
    recall / precision (per-RefOG normalised by |R|-1, as benchmark.py q_even).
    """
    fam_members = defaultdict(set)
    for g, f in gene2fam.items():
        fam_members[f].add(g)
    gene2ref = {}
    for rid, genes in refogs.items():
        for g in genes:
            gene2ref[g] = rid
    per = {}
    tp = fp = fn = 0.0
    for rid, R in refogs.items():
        if not R:
            continue
        cnt = Counter(gene2fam.get(g, "__single__" + g) for g in R)
        best_f, best_n = max(cnt.items(), key=lambda kv: (kv[1], kv[0]))
        fam = fam_members.get(best_f, {best_f[len("__single__"):]} if best_f.startswith("__single__") else set())
        others = Counter(gene2ref[g] for g in fam if g in gene2ref and gene2ref[g] != rid)
        n = len(R)
        # OrthoBench-style pair counts
        t_tp = t_fp = t_fn = 0.0
        for f, ov in cnt.items():
            fs = len(fam_members[f]) if not f.startswith("__single__") else 1
            t_tp += ov * (ov - 1) / 2
            t_fp += ov * (fs - ov)
            t_fn += (n - ov) * ov / 2
        if n > 1:
            tp += t_tp / (n - 1); fp += t_fp / (n - 1); fn += t_fn / (n - 1)
        per[rid] = {
            "size": n,
            "complete": len(cnt) == 1,
            "n_fam": len(cnt),
            "best_frac": best_n / n,
            "best_fam_size": len(fam),
            "size_ratio": len(fam) / n,
            "foreign_refog_frac": sum(others.values()) / len(fam) if fam else 0.0,
            "n_refogs_in_best": 1 + len(others),
        }
    vals = list(per.values())
    k = len(vals)

    def med(xs):
        s = sorted(xs)
        return s[len(s) // 2] if len(s) % 2 else 0.5 * (s[len(s) // 2 - 1] + s[len(s) // 2])

    agg = {
        "n_refogs": k,
        "completeness": sum(v["complete"] for v in vals) / k,
        "n_complete": sum(v["complete"] for v in vals),
        "mean_n_fam": sum(v["n_fam"] for v in vals) / k,
        "median_n_fam": med([v["n_fam"] for v in vals]),
        "max_n_fam": max(v["n_fam"] for v in vals),
        "mean_best_frac": sum(v["best_frac"] for v in vals) / k,
        "median_size_ratio": med([v["size_ratio"] for v in vals]),
        "mean_size_ratio": sum(v["size_ratio"] for v in vals) / k,
        "max_size_ratio": max(v["size_ratio"] for v in vals),
        "fusion_rate": sum(v["n_refogs_in_best"] > 1 for v in vals) / k,
        "mean_foreign_refog_frac": sum(v["foreign_refog_frac"] for v in vals) / k,
        "pair_recall": tp / (tp + fn) if tp + fn else 0.0,
        "pair_precision": tp / (tp + fp) if tp + fp else 0.0,
    }
    return {"aggregate": agg, "per_refog": per}


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
    s2_bitscore: float = 0.0
    mcl_I: float = 1.4
    # resources
    threads: int = 24
    memory_gb: int = 64

    def hash(self) -> str:
        d = asdict(self)
        d.pop("name"); d.pop("threads"); d.pop("memory_gb")
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


def run_families(fasta: str, outdir: str, cfg: FamilyConfig) -> dict:
    """Run S0 with ``cfg``; write families.tsv + stats.json in ``outdir``; return stats."""
    os.makedirs(outdir, exist_ok=True)
    fasta = os.path.abspath(fasta)
    t0 = time.time()
    timings = {}
    gene_ids = read_fasta_ids(fasta)
    if len(set(gene_ids)) != len(gene_ids):
        raise ValueError("duplicate sequence ids in input FASTA")
    with open(os.path.join(outdir, "commands.log"), "a") as log:
        t = time.time()
        if cfg.backend == "diamond":
            s1 = _diamond_cluster(fasta, outdir, cfg, log)
        elif cfg.backend == "mmseqs":
            s1 = _mmseqs_cluster(fasta, outdir, cfg, log)
        else:
            raise ValueError(cfg.backend)
        timings["stage1_s"] = time.time() - t
        pairs = read_pairs_tsv(s1)
        g2r = partition_from_pairs(gene_ids, pairs)
        reps = sorted(set(g2r.values()))
        n_stage1 = len(reps)
        extra = {"n_stage1_clusters": n_stage1}
        if cfg.stage2:
            t = time.time()
            # use the tool's own representative (first column) where available
            rep_of = {}
            for r, m in pairs:
                rep_of[m] = r
            root2rep = {}
            for g in gene_ids:
                root = g2r[g]
                if g in rep_of and rep_of[g] == g:
                    root2rep[root] = g
            for g in gene_ids:  # clusters without explicit rep (singletons)
                root2rep.setdefault(g2r[g], g2r[g])
            rep_set = set(root2rep.values())
            rfa = os.path.join(outdir, "reps.fa")
            _write_subset_fasta(fasta, rep_set, rfa)
            dmnd = os.path.join(outdir, "reps.dmnd")
            hits = os.path.join(outdir, "reps_hits.tsv")
            _run(["diamond", "makedb", "--in", rfa, "-d", dmnd, "--threads", cfg.threads], log)
            cmd = ["diamond", "blastp", "-q", rfa, "-d", dmnd, "-o", hits, "--threads", cfg.threads,
                   "-k", cfg.s2_max_target, "-e", cfg.s2_evalue, "-M", f"{cfg.memory_gb}G",
                   "--outfmt", "6", "qseqid", "sseqid", "evalue", "bitscore", "qlen", "slen",
                   "qstart", "qend", "sstart", "send"]
            if cfg.s2_sens:
                cmd += [cfg.s2_sens]
            _run(cmd, log)
            timings["stage2_search_s"] = time.time() - t
            with open(hits) as fh:
                edges = filter_hits((l.rstrip("\n").split("\t") for l in fh), cfg)
            if cfg.stage2 == "cc":
                rep_pairs = [(a, b) for a, b, _ in edges]
            elif cfg.stage2 == "mcl":
                rep_pairs = _mcl(edges, outdir, cfg, log)
            else:
                raise ValueError(cfg.stage2)
            rep2root = partition_from_pairs(sorted(rep_set), rep_pairs)
            g2r = {g: rep2root[root2rep[g2r[g]]] for g in gene_ids}
            timings["stage2_s"] = time.time() - t
            extra.update(n_reps=len(rep_set), n_rep_edges=len(edges))
            for f in (dmnd, rfa, os.path.join(outdir, "reps_graph.abc")):
                if os.path.exists(f):
                    os.remove(f)
        gene2fam = label_families(g2r)
    write_families_tsv(gene2fam, os.path.join(outdir, "families.tsv"))
    timings["total_s"] = time.time() - t0
    stats = family_stats(gene2fam)
    stats.update(extra)
    stats["timings"] = timings
    stats["config"] = asdict(cfg)
    stats["config_hash"] = cfg.hash()
    with open(os.path.join(outdir, "stats.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    return stats


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_REFOGS = os.path.join(_ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs")
_SPLIT = os.path.join(_ROOT, "splits/orthobench_split.json")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ancestra.families", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cluster", help="partition a FASTA into families")
    c.add_argument("--fasta", required=True)
    c.add_argument("--outdir", required=True)
    c.add_argument("--config", help="JSON file with FamilyConfig fields")
    c.add_argument("--set", nargs="*", default=[], help="override field=value (JSON-parsed)")
    e = sub.add_parser("evaluate", help="score families.tsv vs OrthoBench RefOGs of one split subset")
    e.add_argument("--families", required=True)
    e.add_argument("--subset", choices=["dev", "test"], required=True)
    e.add_argument("--split", default=_SPLIT)
    e.add_argument("--refogs", default=_REFOGS)
    e.add_argument("--keep-uncertain", action="store_true")
    e.add_argument("--out")
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
    else:
        split = json.load(open(a.split))
        refogs = read_refogs(a.refogs, split[a.subset], remove_uncertain=not a.keep_uncertain)
        res = evaluate_refogs(read_families_tsv(a.families), refogs)
        res["subset"] = a.subset
        res["remove_uncertain"] = not a.keep_uncertain
        s = json.dumps(res, indent=1)
        if a.out:
            open(a.out, "w").write(s)
        print(json.dumps(res["aggregate"], indent=1))


if __name__ == "__main__":
    sys.exit(main())
