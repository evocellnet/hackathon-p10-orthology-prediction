"""S8: validate apparent gene losses with genome DNA (demo layer).

An *apparent loss* is a HOG at a clade level (default ``Theria``) that has
members in every species of the clade but one (the *target*). For each one:

1. **Syntenic locus** (:func:`syntenic_locus`). In a *source* species that has
   the HOG (closest to the target first), take the flanking genes (up to ``K``
   on each side, primary-assembly genes only), map each one through its own
   clade-level HOG to a *single-copy* member of the target, and keep the
   nearest left/right pair of target anchors that lie on one contig within
   ``max_span``. The locus is that span padded by ``pad`` (200 kb). If only
   one side maps, the locus is the anchor +- ``pad`` ("one_sided"); with no
   anchor there is no locus.
2. **Protein-to-genome alignment** (miniprot) of the HOG's proteins from the
   present species against the locus DNA (:func:`run_miniprot`).
3. **Classification** (:func:`classify_locus`):

   * ``A_annotated``: intact ORF that overlaps an annotated protein of the
     target that sits in a target-only HOG -> the gene is present, the HOG
     inference split it off (orthology miss, not a loss);
   * ``A_gap``: intact ORF with no annotated protein-coding gene of the input
     proteome -> unannotated gene / annotation gap (not a loss);
   * ``B_pseudo``: >= ``min_cov_hit`` of the query aligned but with frameshifts
     or in-frame stops -> pseudogene (true loss, pseudogenized);
   * ``B_partial``: >= ``min_cov_hit`` aligned without disruptions but not
     intact (truncated / degraded remnant; loss or partial annotation gap);
   * ``C_deleted``: no hit in a bracketed locus with < ``max_nfrac`` Ns
     (flanks conserved, gene sequence gone -> true loss by deletion);
   * ``D_undetermined``: no bracketed locus, or locus mostly assembly gap.

   Hits that overlap an annotated target gene belonging to a *different* HOG
   that has members in other species (i.e. a paralog with its own orthologs)
   are excluded before classifying ("paralog hits").
4. Genome-wide miniprot for B/C/D detects relocated intact copies.

The module only reads the core outputs (OrthoXML, gene-order tables); it never
modifies them. Driver: ``scripts/run_S8_geneloss.py``; summary: ``scripts/summarize_S8.py``.
"""
from __future__ import annotations

import gzip
import json
import os
import random
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd

__all__ = ["short", "HogIndex", "GeneOrderIndex", "FastaIndex", "Locus", "MiniprotHit",
           "select_apparent_losses", "select_complete_single_copy", "syntenic_locus",
           "parse_miniprot_paf", "CdsIndex", "query_votes", "classify_hits", "run_miniprot", "INTACT_COV", "MAMMALS",
           "SOURCE_ORDER", "LocusSearcher", "random_window", "read_fasta", "write_fasta", "hit_record"]

INTACT_COV = 0.8        # query fraction aligned for an "intact ORF"
MIN_COV_HIT = 0.3       # query fraction aligned for "a hit"


def short(species: str) -> str:
    """'Homo_sapiens.GRCh38.pep.all' -> 'Homo_sapiens'."""
    return species.split(".")[0]


# ----------------------------------------------------------------------------
# HOGs
# ----------------------------------------------------------------------------

class HogIndex:
    """Clade-level HOGs from an OrthoXML file: gene -> hog id, hog -> genes."""

    def __init__(self, hogs: list, species_of: dict, outgroup_hogs: list | None = None):
        self.members = [sorted(h) for h in hogs]
        self.species_of = species_of            # protein -> short species name
        self.hog_of = {}
        for i, h in enumerate(self.members):
            for g in h:
                self.hog_of[g] = i
        self.out_hog_of = {}
        self.out_members = []
        if outgroup_hogs is not None:
            self.out_members = [sorted(h) for h in outgroup_hogs]
            for i, h in enumerate(self.out_members):
                for g in h:
                    self.out_hog_of[g] = i

    @classmethod
    def from_orthoxml(cls, path: str, tree_path: str, level: str = "Theria",
                      outgroup_level: str | None = "Amniota"):
        from .orthoxml import read_orthoxml, hogs_at_level
        from .newick import read_newick
        tree = read_newick(tree_path)
        ox = read_orthoxml(path)
        sp = {g: short(e.species) for g, e in ox.genes.items()}
        hogs = hogs_at_level(ox, level, species_tree=tree, include_singletons=True)
        out = (hogs_at_level(ox, outgroup_level, species_tree=tree, include_singletons=True)
               if outgroup_level else None)
        return cls(hogs, sp, out)

    def species_counts(self, i: int) -> Counter:
        return Counter(self.species_of[g] for g in self.members[i])

    def in_species(self, i: int, sp: str) -> list:
        return [g for g in self.members[i] if self.species_of[g] == sp]

    def outgroup_species(self, i: int) -> set:
        """Species of the enclosing outgroup-level HOG (e.g. Amniota)."""
        if not self.out_members:
            return set()
        g = self.members[i][0]
        j = self.out_hog_of.get(g)
        return set() if j is None else {self.species_of[x] for x in self.out_members[j]}


def select_apparent_losses(hi: HogIndex, clade: list, max_copies: int = 2,
                           outgroup: list | None = None, root_species: set | None = None) -> pd.DataFrame:
    """HOGs present in all species of ``clade`` but exactly one (the target).

    ``max_copies``: at most this many genes per present species (keeps the
    synteny mapping simple). If the target is in ``root_species`` (the species
    whose absence alone could mean a gain on the sister branch, e.g. opossum
    for Theria) the HOG must also have an ``outgroup`` member at the
    outgroup level, otherwise absence there is ambiguous (gain vs loss).
    """
    rows = []
    root_species = root_species or set()
    for i, h in enumerate(hi.members):
        c = hi.species_counts(i)
        present = [s for s in clade if c.get(s, 0) > 0]
        absent = [s for s in clade if c.get(s, 0) == 0]
        if len(absent) != 1 or len(present) != len(clade) - 1:
            continue
        if max(c[s] for s in present) > max_copies:
            continue
        t = absent[0]
        og = sorted(hi.outgroup_species(i) & set(outgroup or []))
        if t in root_species and not og:
            continue
        rows.append({"hog": i, "target": t, "n_genes": len(h),
                     "max_copies": max(c[s] for s in present), "outgroup": ",".join(og)})
    return pd.DataFrame(rows)


def select_complete_single_copy(hi: HogIndex, clade: list) -> list:
    """HOG ids with exactly one gene in every clade species (positive controls)."""
    out = []
    for i in range(len(hi.members)):
        c = hi.species_counts(i)
        if len(h := hi.members[i]) == len(clade) and all(c.get(s, 0) == 1 for s in clade):
            out.append(i)
    return out


# ----------------------------------------------------------------------------
# gene order
# ----------------------------------------------------------------------------

class GeneOrderIndex:
    """Primary-assembly gene coordinates per species (from S5 gene-order tables).

    Alt-haplotype genes projected onto a chromosome (``projected``) and
    non-primary contigs are dropped, and ranks are recomputed on the rest.
    """

    def __init__(self, tables: dict):
        self.tab = {}
        self.by_contig = {}
        for sp, df in tables.items():
            d = df
            if "primary" in d:
                d = d[d["primary"].astype(bool).to_numpy()]
            if "projected" in d:
                d = d[~d["projected"].astype(bool).to_numpy()]
            d = d.sort_values(["contig", "start", "end", "protein_id"]).reset_index(drop=True)
            d["rank"] = d.groupby("contig").cumcount()
            self.tab[sp] = d.set_index("protein_id", drop=False)
            self.by_contig[sp] = {c: g.reset_index(drop=True) for c, g in d.groupby("contig")}

    @classmethod
    def from_dir(cls, table_dir: str, species: list):
        return cls({sp: pd.read_parquet(os.path.join(table_dir, f"{sp}.parquet")) for sp in species})

    def row(self, sp, prot):
        t = self.tab[sp]
        return t.loc[prot] if prot in t.index else None

    def neighbours(self, sp, prot, k: int):
        """(left, right) protein lists, nearest first, up to k each."""
        r = self.row(sp, prot)
        if r is None:
            return [], []
        g = self.by_contig[sp][r["contig"]]
        i = int(r["rank"])
        left = list(g["protein_id"].iloc[max(0, i - k):i])[::-1]
        right = list(g["protein_id"].iloc[i + 1:i + 1 + k])
        return left, right

    def overlapping(self, sp, contig, start, end) -> pd.DataFrame:
        g = self.by_contig.get(sp, {}).get(contig)
        if g is None:
            return g
        return g[(g["start"] <= end) & (g["end"] >= start)]


@dataclass
class Locus:
    target: str
    contig: str
    start: int            # 1-based inclusive
    end: int
    kind: str             # bracketed | one_sided
    source: str           # source species
    source_gene: str
    left_anchor: str | None = None
    right_anchor: str | None = None
    k_left: int | None = None
    k_right: int | None = None
    n_sources: int = 1      # source species giving a bracketed locus
    n_agree: int = 1        # of which agree with this one


def syntenic_locus(hi: HogIndex, go: GeneOrderIndex, hog: int, target: str, sources: list,
                   k: int = 5, pad: int = 200_000, max_span: int = 2_000_000,
                   contig_len: dict | None = None, vote: bool = True) -> Locus | None:
    """Expected locus, with synteny voting across source species.

    Every source species (and copy) proposes a locus (:func:`_locus_one_source`);
    bracketed proposals are clustered (same contig, overlapping) and the
    cluster backed by the most source species wins (ties: closest source).
    This protects against one source whose gene is itself an odd model (e.g.
    an annotated retrocopy) pointing to the wrong place. ``vote=False``
    returns the first bracketed proposal (closest source).
    """
    props = []
    for src in sources:
        loc = _locus_one_source(hi, go, hog, target, [src], k, pad, max_span, contig_len)
        if loc is not None:
            props.append(loc)
            if not vote and loc.kind == "bracketed":
                return loc
    br = [p for p in props if p.kind == "bracketed"]
    if not br:
        return props[0] if props else None
    clusters = []
    for p in br:
        for c in clusters:
            if any(q.contig == p.contig and q.start <= p.end and p.start <= q.end for q in c):
                c.append(p)
                break
        else:
            clusters.append([p])
    best = max(clusters, key=lambda c: (len({q.source for q in c}), -sources.index(c[0].source)))
    loc = best[0]
    loc.n_sources = len(br)
    loc.n_agree = len({q.source for q in best})
    return loc


def _locus_one_source(hi: HogIndex, go: GeneOrderIndex, hog: int, target: str, sources: list,
                      k: int = 5, pad: int = 200_000, max_span: int = 2_000_000,
                      contig_len: dict | None = None, one_sided_pad: int = 1_000_000) -> Locus | None:
    """Expected locus of HOG ``hog`` in ``target`` from flanking orthologs.

    ``sources`` are tried in order (closest species first); a bracketed locus
    from any source beats a one-sided one. Members of the HOG itself in the
    target are never used as anchors (positive controls hide them).
    """
    own = set(hi.members[hog])
    best_one = None
    for src in sources:
        for sg in hi.in_species(hog, src):
            left, right = go.neighbours(src, sg, k)
            if not left and not right:
                continue

            def anchors(side):
                out = []
                for kk, n in enumerate(side, 1):
                    j = hi.hog_of.get(n)
                    if j is None or j == hog:
                        continue
                    tm = [g for g in hi.in_species(j, target) if g not in own]
                    if len(tm) != 1:
                        continue
                    r = go.row(target, tm[0])
                    if r is None:
                        continue
                    out.append((kk, tm[0], r))
                return out
            la, ra = anchors(left), anchors(right)
            pairs = []
            for (kl, gl, rl) in la:
                for (kr, gr, rr) in ra:
                    if rl["contig"] != rr["contig"]:
                        continue
                    lo = min(rl["start"], rr["start"])
                    hi_ = max(rl["end"], rr["end"])
                    if hi_ - lo > max_span:
                        continue
                    pairs.append((kl + kr, kl, kr, gl, gr, rl, rr, lo, hi_))
            if pairs:
                pairs.sort(key=lambda x: (x[0], x[7]))
                _, kl, kr, gl, gr, rl, rr, lo, hi_ = pairs[0]
                s, e = max(1, lo - pad), hi_ + pad
                if contig_len and rl["contig"] in contig_len:
                    e = min(e, contig_len[rl["contig"]])
                return Locus(target, rl["contig"], int(s), int(e), "bracketed", src, sg, gl, gr, kl, kr)
            if best_one is None and (la or ra):
                cand = sorted(la + ra, key=lambda x: x[0])[0]
                kk, g, r = cand
                s, e = max(1, int(r["start"]) - one_sided_pad), int(r["end"]) + one_sided_pad
                if contig_len and r["contig"] in contig_len:
                    e = min(e, contig_len[r["contig"]])
                side_left = any(g == x[1] for x in la)
                best_one = Locus(target, r["contig"], s, int(e), "one_sided", src, sg,
                                 g if side_left else None, None if side_left else g,
                                 kk if side_left else None, None if side_left else kk)
    return best_one


# ----------------------------------------------------------------------------
# genome FASTA access (.fai)
# ----------------------------------------------------------------------------

class FastaIndex:
    """Random access to an uncompressed FASTA with a samtools ``.fai`` index."""

    def __init__(self, path: str):
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai):
            subprocess.run(["samtools", "faidx", path], check=True)
        self.idx = {}
        with open(fai) as fh:
            for line in fh:
                name, ln, off, lb, lbytes = line.split("\t")[:5]
                self.idx[name] = (int(ln), int(off), int(lb), int(lbytes))

    def length(self, name):
        return self.idx[name][0]

    def lengths(self):
        return {k: v[0] for k, v in self.idx.items()}

    def fetch(self, name: str, start: int, end: int) -> str:
        """1-based inclusive."""
        ln, off, lb, lbytes = self.idx[name]
        start = max(1, start)
        end = min(ln, end)
        s0 = start - 1
        b0 = off + (s0 // lb) * lbytes + s0 % lb
        e0 = end
        b1 = off + (e0 // lb) * lbytes + e0 % lb
        with open(self.path, "rb") as fh:
            fh.seek(b0)
            raw = fh.read(b1 - b0)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode()


def read_fasta(path: str) -> dict:
    seqs, name, buf = {}, None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


# ----------------------------------------------------------------------------
# miniprot
# ----------------------------------------------------------------------------

@dataclass
class MiniprotHit:
    query: str
    qlen: int
    qs: int
    qe: int
    strand: str
    contig: str
    ts: int               # 0-based start in the searched sequence
    te: int               # 0-based exclusive end
    n_match: int          # matching nt (PAF col 10)
    aln_len: int          # alignment length in nt (PAF col 11)
    score: int
    n_pos: int
    fs: int               # frameshifts
    st: int               # in-frame stop codons
    da: int = 0           # non-canonical acceptor sites
    do: int = 0           # non-canonical donor sites
    n_introns: int = 0
    blocks: list = field(default_factory=list)   # CDS blocks, 1-based inclusive, searched-seq coords

    @property
    def qcov(self) -> float:
        return (self.qe - self.qs) / max(self.qlen, 1)

    @property
    def identity(self) -> float:
        return self.n_match / max(self.aln_len, 1)

    @property
    def splice_nc(self) -> bool:
        return self.da > 0 or self.do > 0

    @property
    def disrupted(self) -> bool:
        return self.fs > 0 or self.st > 0

    @property
    def intact(self) -> bool:
        return self.qcov >= INTACT_COV and not self.disrupted


def _paf_hit(f: list) -> MiniprotHit:
    tags = {}
    for t in f[12:]:
        p = t.split(":", 2)
        if len(p) == 3 and p[1] == "i":
            tags[p[0]] = int(p[2])
        elif len(p) == 3 and p[0] == "cg":
            tags["n_introns"] = sum(p[2].count(c) for c in "NUV")
    return MiniprotHit(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]),
                       int(f[9]), int(f[10]), tags.get("AS", 0), tags.get("np", 0),
                       tags.get("fs", 0), tags.get("st", 0), max(tags.get("da", 0), 0),
                       max(tags.get("do", 0), 0), tags.get("n_introns", 0))


def parse_miniprot_paf(text: str) -> list:
    """Parse miniprot output: plain PAF, or ``--gff`` output (``##PAF`` lines
    followed by the hit's mRNA/CDS lines, whose CDS blocks are attached).
    Tags: ``fs`` frameshifts, ``st`` in-frame stops, ``da``/``do``
    non-canonical splice sites, ``cg`` CIGAR (introns = N/U/V ops)."""
    hits = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith("##PAF\t"):
            f = line.rstrip("\n").split("\t")[1:]
        elif line.startswith("#"):
            continue
        else:
            f = line.rstrip("\n").split("\t")
            if len(f) == 9:          # GFF line
                if f[2] == "CDS" and hits:
                    hits[-1].blocks.append((int(f[3]), int(f[4])))
                continue
        if len(f) < 12 or f[5] == "*":
            continue
        hits.append(_paf_hit(f))
    for h in hits:
        if not h.blocks:
            h.blocks = [(h.ts + 1, h.te)]
    return hits


def run_miniprot(target_fa: str, query_fa: str, threads: int = 1, extra=None) -> list:
    cmd = ["miniprot", "--gff", "-t", str(threads), "--outs=0.5", "-p", "0.5", "-N", "20", "--outc=0.1"]
    cmd += list(extra or [])
    cmd += [target_fa, query_fa]
    r = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return parse_miniprot_paf(r.stdout)


class CdsIndex:
    """CDS segments of the annotated (input-proteome) proteins of one species."""

    def __init__(self, df: pd.DataFrame):
        self.by_contig = {}
        for c, g in df.groupby("contig"):
            g = g.sort_values("start")
            self.by_contig[c] = (g["start"].to_numpy(), g["end"].to_numpy(), g["protein_id"].to_numpy())

    @classmethod
    def from_parquet(cls, path: str):
        return cls(pd.read_parquet(path))

    def overlap(self, contig: str, blocks: list) -> dict:
        """protein -> bp of ``blocks`` (genomic, 1-based) inside its CDS."""
        t = self.by_contig.get(contig)
        out = defaultdict(int)
        if t is None:
            return out
        st, en, pid = t
        for (s, e) in blocks:
            m = (st <= e) & (en >= s)
            for i in np.nonzero(m)[0]:
                out[pid[i]] += min(e, en[i]) - max(s, st[i]) + 1
        return out


def write_fasta(path: str, recs: dict, width: int = 80):
    with open(path, "w") as fh:
        for k, s in recs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


# ----------------------------------------------------------------------------
# classification
# ----------------------------------------------------------------------------

def classify_hits(hits: list, bracketed: bool, nfrac: float, max_nfrac: float = 0.2,
                  annotated_overlap=None, majority: bool = True) -> dict:
    """Classify one locus from its (paralog/retrocopy-filtered) miniprot hits.

    ``annotated_overlap(hit) -> str | None`` returns the overlapped annotated
    target protein (this HOG or a target-only HOG) or None.

    ``majority``: an *unannotated* intact ORF is accepted (``A_gap``) only if
    more query proteins give an intact best hit than a disrupted one;
    otherwise the locus is ``B_pseudo`` with ``conflict=True``. This guards
    against a pseudogene modelled as protein-coding in one source species
    (its own model aligns intact; the other species' proteins expose the
    disruptions).
    """
    hits = [h for h in hits if h.qcov >= 0.1]
    intact = [h for h in hits if h.intact]
    votes = query_votes(hits, {h.query for h in hits})
    if intact:
        best = max(intact, key=lambda h: (h.qcov, h.score))
        ov = None
        if annotated_overlap:
            for h in sorted(intact, key=lambda h: -h.score):
                ov = annotated_overlap(h)
                if ov:
                    best = h
                    break
        if ov:
            return {"cls": "A_annotated", "best": best, "overlap_gene": ov, "conflict": False}
        if not majority or votes["nq_intact"] > votes["nq_disrupted"]:
            return {"cls": "A_gap", "best": best, "overlap_gene": None, "conflict": False}
        dis = [h for h in hits if h.disrupted and h.qcov >= MIN_COV_HIT]
        return {"cls": "B_pseudo", "best": max(dis, key=lambda h: (h.qcov, h.score)), "overlap_gene": None,
                "conflict": True}
    big = [h for h in hits if h.qcov >= MIN_COV_HIT]
    if big:
        dis = [h for h in big if h.disrupted]
        if dis:
            best = max(dis, key=lambda h: (h.qcov, h.score))
            return {"cls": "B_pseudo", "best": best, "overlap_gene": None, "conflict": False}
        best = max(big, key=lambda h: (h.qcov, h.score))
        return {"cls": "B_partial", "best": best, "overlap_gene": None, "conflict": False}
    best = max(hits, key=lambda h: h.score) if hits else None
    if bracketed and nfrac < max_nfrac:
        return {"cls": "C_deleted", "best": best, "overlap_gene": None, "conflict": False}
    return {"cls": "D_undetermined", "best": best, "overlap_gene": None, "conflict": False}


def query_votes(hits: list, queries) -> dict:
    """Per-query status of the best hit (by score): intact / disrupted
    (qcov >= MIN_COV_HIT with fs/st) / partial / none, plus the species-level
    detail. Used for the majority rule (a pseudogene modelled as protein-coding
    in one source species gives one intact query against several disrupted)."""
    best = {}
    for h in hits:
        if h.query not in best or h.score > best[h.query].score:
            best[h.query] = h
    n = Counter()
    for qn in queries:
        h = best.get(qn)
        if h is None or h.qcov < MIN_COV_HIT:
            n["none"] += 1
        elif h.intact:
            n["intact"] += 1
        elif h.disrupted:
            n["disrupted"] += 1
        else:
            n["partial"] += 1
    return {"nq_intact": n["intact"], "nq_disrupted": n["disrupted"], "nq_partial": n["partial"], "nq_none": n["none"]}


def hit_record(h: MiniprotHit | None, prefix: str = "") -> dict:
    keys = ("query", "qcov", "identity", "fs", "st", "splice_nc", "n_introns", "score")
    if h is None:
        return {f"{prefix}{k}": None for k in keys}
    return {f"{prefix}query": h.query, f"{prefix}qcov": round(h.qcov, 3), f"{prefix}identity": round(h.identity, 3),
            f"{prefix}fs": h.fs, f"{prefix}st": h.st, f"{prefix}splice_nc": int(h.da + h.do),
            f"{prefix}n_introns": h.n_introns, f"{prefix}score": h.score}


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------

MAMMALS = ["Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Rattus_norvegicus",
           "Canis_familiaris", "Monodelphis_domestica"]
# closest-first source order per target (from the OrthoBench species tree)
SOURCE_ORDER = {
    "Homo_sapiens": ["Pan_troglodytes", "Mus_musculus", "Rattus_norvegicus", "Canis_familiaris", "Monodelphis_domestica"],
    "Pan_troglodytes": ["Homo_sapiens", "Mus_musculus", "Rattus_norvegicus", "Canis_familiaris", "Monodelphis_domestica"],
    "Mus_musculus": ["Rattus_norvegicus", "Homo_sapiens", "Pan_troglodytes", "Canis_familiaris", "Monodelphis_domestica"],
    "Rattus_norvegicus": ["Mus_musculus", "Homo_sapiens", "Pan_troglodytes", "Canis_familiaris", "Monodelphis_domestica"],
    "Canis_familiaris": ["Homo_sapiens", "Mus_musculus", "Pan_troglodytes", "Rattus_norvegicus", "Monodelphis_domestica"],
    "Monodelphis_domestica": ["Homo_sapiens", "Canis_familiaris", "Mus_musculus", "Pan_troglodytes", "Rattus_norvegicus"],
}


class LocusSearcher:
    """Runs the locus method for (hog, target) jobs of one target species."""

    def __init__(self, hi: HogIndex, go: GeneOrderIndex, target: str, genome: FastaIndex,
                 proteins: dict, workdir: str, pad: int = 200_000, k: int = 5, max_span: int = 2_000_000,
                 cds: CdsIndex | None = None, n_exons: dict | None = None):
        self.hi, self.go, self.target, self.genome = hi, go, target, genome
        self.cds = cds
        self.n_exons = n_exons or {}       # protein -> number of CDS segments (all species)
        self.proteins = proteins           # protein id -> sequence (all species)
        self.workdir = workdir
        self.pad, self.k, self.max_span = pad, k, max_span
        self.clen = genome.lengths()
        os.makedirs(workdir, exist_ok=True)

    def queries(self, hog: int, hide: set = frozenset()) -> dict:
        return {g: self.proteins[g].replace("*", "") for g in self.hi.members[hog]
                if self.hi.species_of[g] != self.target and g not in hide and g in self.proteins}

    def paralog_or_fission(self, contig, blocks, hog, min_frac: float = 0.2):
        """Annotated target proteins whose CDS covers >= ``min_frac`` of the hit's
        CDS blocks (genomic, 1-based). Returns (fission_gene, paralog_gene):
        a protein in this HOG or in a target-only HOG is a *fission* (the gene
        is annotated but the HOG inference split it off); a protein in another
        HOG with members in other species is a *paralog* hit."""
        ov = self.cds.overlap(contig, blocks) if self.cds is not None else {}
        tot = sum(e - s + 1 for s, e in blocks) or 1
        fission, paralog = None, None
        for g, bp in sorted(ov.items(), key=lambda x: -x[1]):
            if bp < min_frac * tot:
                continue
            j = self.hi.hog_of.get(g)
            if j == hog:
                fission = fission or g
                continue
            spp = {self.hi.species_of[x] for x in self.hi.members[j]} if j is not None else {self.target}
            if spp <= {self.target}:
                fission = fission or g
            else:
                paralog = paralog or g
        return fission, paralog

    def search(self, job_id: str, hog: int, locus: Locus | None, hide: set = frozenset(),
               mask: tuple | None = None, keep_hidden_as_fission: str | None = None) -> dict:
        """Search one locus. ``mask=(start,end)`` (genomic) is replaced by N.

        ``keep_hidden_as_fission``: for positive controls, the hidden target
        gene counts as the expected hit (reported as overlap).
        """
        rec = {"job": job_id, "hog": hog, "target": self.target}
        if locus is None:
            rec.update({"cls": "D_undetermined", "locus_kind": "none"})
            rec.update(hit_record(None))
            return rec
        rec.update({"locus_kind": locus.kind, "contig": locus.contig, "start": locus.start, "end": locus.end,
                    "source": locus.source, "source_gene": locus.source_gene,
                    "left_anchor": locus.left_anchor, "right_anchor": locus.right_anchor,
                    "k_left": locus.k_left, "k_right": locus.k_right,
                    "n_sources": locus.n_sources, "n_agree": locus.n_agree})
        seq = self.genome.fetch(locus.contig, locus.start, locus.end)
        if mask is not None:
            ms, me = max(mask[0], locus.start), min(mask[1], locus.end)
            if me >= ms:
                a, b = ms - locus.start, me - locus.start + 1
                seq = seq[:a] + "N" * (b - a) + seq[b:]
        nfrac = (seq.count("N") + seq.count("n")) / max(len(seq), 1)
        if mask is not None:
            nfrac = 0.0  # the mask is the experiment, not an assembly gap
        rec["nfrac"] = round(nfrac, 4)
        q = self.queries(hog, hide)
        rec["n_queries"] = len(q)
        if not q:
            rec.update({"cls": "D_undetermined"})
            rec.update(hit_record(None))
            return rec
        d = tempfile.mkdtemp(dir=self.workdir)
        tf, qf = os.path.join(d, "t.fa"), os.path.join(d, "q.faa")
        write_fasta(tf, {"locus": seq})
        write_fasta(qf, q)
        hits = run_miniprot(tf, qf, threads=1)
        for f in (tf, qf):
            os.remove(f)
        os.rmdir(d)
        kept, par_hits, retro = [], [], []
        ex = [self.n_exons.get(g, 1) for g in q]
        multi_exon_hog = bool(ex) and float(np.median(ex)) >= 2
        ann = {}
        for h in hits:
            gb = [(locus.start + s0 - 1, locus.start + e0 - 1) for s0, e0 in h.blocks]
            fis, par = self.paralog_or_fission(locus.contig, gb, hog)
            if keep_hidden_as_fission is not None and self.cds is not None:
                if self.cds.overlap(locus.contig, gb).get(keep_hidden_as_fission, 0) > 0:
                    fis, par = keep_hidden_as_fission, None
            if par is not None and fis is None:
                par_hits.append((h, par))
                continue
            if fis is None and h.n_introns == 0 and multi_exon_hog:
                retro.append(h)      # intronless copy of a multi-exon gene: retrocopy
                continue
            ann[id(h)] = fis
            h.genomic_blocks = gb
            kept.append(h)
        rec["n_retro_hits"] = len(retro)
        rec["retro_intact"] = any(h.intact for h in retro)
        n_par = len(par_hits)
        if par_hits:
            hp, gp = max(par_hits, key=lambda x: x[0].score)
            rec.update({"paralog_gene": gp, "paralog_identity": round(hp.identity, 3), "paralog_query": hp.query,
                        "paralog_qcov": round(hp.qcov, 3), "paralog_intact": hp.intact})
        rec["n_hits"] = len(hits)
        rec["any_intact_raw"] = any(h.intact for h in hits)
        rec["n_paralog_hits"] = n_par
        c = classify_hits(kept, locus.kind == "bracketed", nfrac, annotated_overlap=lambda h: ann.get(id(h)))
        rec.update(query_votes(kept, q))
        rec["cls"] = c["cls"]
        rec["conflict"] = c["conflict"]
        rec["overlap_gene"] = c["overlap_gene"]
        rec.update(hit_record(c["best"]))
        if c["best"] is not None:
            rec["hit_start"] = locus.start + c["best"].ts
            rec["hit_end"] = locus.start + c["best"].te - 1
            rec["hit_strand"] = c["best"].strand
            rec["hit_blocks"] = ",".join(f"{x}-{y}" for x, y in c["best"].genomic_blocks)
        # best intact hit regardless of class (for controls)
        intact = [h for h in kept if h.intact]
        rec["any_intact"] = bool(intact)
        return rec


def random_window(genome: FastaIndex, size: int, rng: random.Random, min_contig: int = 5_000_000,
                  allowed: set | None = None) -> tuple:
    lens = {k: v for k, v in genome.lengths().items() if v >= max(min_contig, size + 1)
            and (allowed is None or k in allowed)}
    names = sorted(lens)
    w = np.array([lens[n] for n in names], float)
    w /= w.sum()
    c = names[int(rng.choices(range(len(names)), weights=w)[0])]
    s = rng.randint(1, lens[c] - size)
    return c, s, s + size - 1
