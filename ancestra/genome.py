"""Gene order and ancestral adjacency data for the synteny layer (DESIGN.md S5).

This module only *prepares* the data; synteny scoring lives in the context-evidence
function passed to :func:`ancestra.merge.run_merge` (see docs/core_v0.md §3.4).

* :class:`GeneOrder` -- per gene (integer gene ID of :class:`ancestra.io.GeneTable`):
  contig index (global), rank along the contig among *proteome genes* (genes without a
  proteome entry are skipped, so the order is gap-tolerant w.r.t. unannotated loci),
  start coordinate and strand. Genes without coordinates have contig -1.
* :func:`load_gene_order_gff3` -- Ensembl-style GFF3 (gene / mRNA / CDS with
  ``protein_id``), as in the Ensembl 99 files for OrthoBench and the simulator's
  ``gff/<species>.gff3``. Protein -> transcript -> gene; WormBase-style IDs
  (``<gene_id>.1``) fall back to the gene ID. For human/mouse/zebrafish the
  ``chr_patch_hapl_scaff`` file is preferred (OrthoBench contains alt-haplotype genes;
  those genes sit on their own contigs, so they only get neighbours there).
* :func:`neighbor_pairs` -- all same-contig gene pairs within ``w`` ranks (w=5 per DESIGN).
* :func:`unit_adjacency` -- maps neighbour pairs through a gene -> ancestral-gene map:
  two ancestral genes are adjacent if some of their member genes are neighbours in at
  least one species; weight = number of such species, plus the minimal rank distance.
* :meth:`GeneOrder.shuffled` -- the pre-registered shuffled-order control (ranks permuted
  within each species, contigs kept).
"""
from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass

import numpy as np

from .io import GeneTable

__all__ = ["GeneOrder", "load_gene_order_gff3", "neighbor_pairs", "unit_adjacency"]


@dataclass
class GeneOrder:
    contig: np.ndarray     # int64 per gene, -1 = unknown
    rank: np.ndarray       # int64 per gene: order along its contig (0..), -1 = unknown
    start: np.ndarray      # int64 per gene
    strand: np.ndarray     # int8 per gene: +1 / -1 / 0
    contig_names: list     # global contig index -> "species:seqid"
    coverage: dict         # species -> fraction of proteome genes with coordinates

    def shuffled(self, seed: int = 11) -> "GeneOrder":
        """Shuffled-order control: permute the gene ranks within each species' set of
        located genes (contig sizes kept; adjacency becomes random)."""
        rng = np.random.default_rng(seed)
        contig, rank = self.contig.copy(), self.rank.copy()
        sp_of_contig = np.array([n.split(":", 1)[0] for n in self.contig_names])
        located = np.flatnonzero(contig >= 0)
        sps = sp_of_contig[contig[located]]
        for sp in np.unique(sps):
            g = located[sps == sp]
            perm = rng.permutation(len(g))
            contig[g] = self.contig[g[perm]]
            rank[g] = self.rank[g[perm]]
        return GeneOrder(contig, rank, self.start.copy(), self.strand.copy(), self.contig_names,
                         dict(self.coverage))


def _find_gff(gff_dir: str, species: str) -> str | None:
    files = sorted(os.listdir(gff_dir))
    for f in files:                      # exact: simulator naming <species>.gff3
        if f in (species + ".gff3", species + ".gff3.gz", species + ".gff"):
            return os.path.join(gff_dir, f)
    stem = species.split(".")[0]         # Ensembl: Genus_species.Assembly...
    cand = [f for f in files if f.startswith(stem + ".") and ".gff" in f]
    cand.sort(key=lambda f: ("chr_patch_hapl_scaff" not in f, f))
    return os.path.join(gff_dir, cand[0]) if cand else None


_attr = re.compile(r"([^=;]+)=([^;]*)")


def _parse_gff3(path: str):
    op = gzip.open if path.endswith(".gz") else open
    genes, tx2gene, prot2tx, geneid2gene = {}, {}, {}, {}
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            typ = p[2]
            if typ in ("exon", "five_prime_UTR", "three_prime_UTR", "chromosome", "region",
                       "scaffold", "biological_region"):
                continue
            a = dict(_attr.findall(p[8]))
            fid = a.get("ID", "")
            if typ == "CDS":
                pid = a.get("protein_id")
                if pid and "Parent" in a:
                    prot2tx.setdefault(pid, a["Parent"])
            elif fid.startswith("gene:") or typ in ("gene", "ncRNA_gene", "pseudogene"):
                if fid:
                    genes[fid] = (p[0], int(p[3]), p[6])
                    if "gene_id" in a:
                        geneid2gene[a["gene_id"]] = fid
                    if "Name" in a:
                        geneid2gene.setdefault(a["Name"], fid)
            elif fid and "Parent" in a:          # any transcript type (mRNA, NMD_..., etc.)
                tx2gene[fid] = a["Parent"]
    return genes, tx2gene, prot2tx, geneid2gene


def load_gene_order_gff3(gff_dir: str, gt: GeneTable) -> GeneOrder:
    N = gt.n
    contig = np.full(N, -1, dtype=np.int64)
    start = np.zeros(N, dtype=np.int64)
    strand = np.zeros(N, dtype=np.int8)
    names: list = []
    cidx: dict = {}
    coverage = {}
    for si, (sp, genes_sp) in enumerate(zip(gt.species, gt.genes_of_species())):
        path = _find_gff(gff_dir, sp)
        if path is None:
            coverage[sp] = 0.0
            continue
        genes, tx2gene, prot2tx, geneid2gene = _parse_gff3(path)
        n_ok = 0
        for g in genes_sp.tolist():
            pid = gt.names[g]
            gene = None
            tx = prot2tx.get(pid) or prot2tx.get(pid.split(".")[0])
            if tx is not None:
                gene = tx2gene.get(tx)
            if gene is None:
                base = pid[:-2] if pid.endswith(".1") else pid
                gene = geneid2gene.get(base) or geneid2gene.get(pid)
            if gene is None or gene not in genes:
                continue
            seqid, st, sd = genes[gene]
            key = f"{sp}:{seqid}"
            if key not in cidx:
                cidx[key] = len(names)
                names.append(key)
            contig[g] = cidx[key]
            start[g] = st
            strand[g] = 1 if sd == "+" else (-1 if sd == "-" else 0)
            n_ok += 1
        coverage[sp] = n_ok / max(1, len(genes_sp))
    rank = np.full(N, -1, dtype=np.int64)
    loc = np.flatnonzero(contig >= 0)
    order = loc[np.lexsort((start[loc], contig[loc]))]
    c = contig[order]
    first = np.r_[True, c[1:] != c[:-1]]
    grp_start = np.maximum.accumulate(np.where(first, np.arange(len(order)), 0))
    rank[order] = np.arange(len(order)) - grp_start
    return GeneOrder(contig, rank, start, strand, names, coverage)


def neighbor_pairs(go: GeneOrder, w: int = 5):
    """Same-contig gene pairs (a, b) with 1 <= rank_b - rank_a <= w.
    Returns (a, b, dist) int64 arrays (a, b = gene IDs)."""
    loc = np.flatnonzero(go.contig >= 0)
    order = loc[np.lexsort((go.rank[loc], go.contig[loc]))]
    c = go.contig[order]
    A, B, D = [], [], []
    for k in range(1, w + 1):
        same = c[k:] == c[:-k]
        A.append(order[:-k][same])
        B.append(order[k:][same])
        D.append(np.full(int(same.sum()), k, dtype=np.int64))
    return np.concatenate(A), np.concatenate(B), np.concatenate(D)


def unit_adjacency(na: np.ndarray, nb: np.ndarray, nd: np.ndarray, unit_of_gene: np.ndarray,
                   species_of: np.ndarray, n_units: int):
    """Ancestral adjacency between units.

    ``unit_of_gene[g]`` = local unit index of gene g, or -1 (gene not in scope).
    Returns (ua, ub, n_species, min_dist) with ua < ub: units adjacent through member
    genes that are neighbours in >= 1 species (n_species = number of such species).
    """
    ua, ub = unit_of_gene[na], unit_of_gene[nb]
    m = (ua >= 0) & (ub >= 0) & (ua != ub)
    ua, ub, d, sp = ua[m], ub[m], nd[m], species_of[na[m]].astype(np.int64)
    lo, hi = np.minimum(ua, ub), np.maximum(ua, ub)
    if len(lo) == 0:
        z = np.zeros(0, dtype=np.int64)
        return z, z, z, z
    S = int(species_of.max()) + 1
    key = (lo * n_units + hi) * S + sp
    order = np.argsort(key, kind="stable")
    key, d = key[order], d[order]
    st = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    kd = key[st]
    dmin = np.minimum.reduceat(d, st)
    pair = kd // S
    st2 = np.flatnonzero(np.r_[True, pair[1:] != pair[:-1]])
    nsp = np.diff(np.r_[st2, len(pair)])
    pmin = np.minimum.reduceat(dmin, st2)
    p = pair[st2]
    return p // n_units, p % n_units, nsp.astype(np.int64), pmin
