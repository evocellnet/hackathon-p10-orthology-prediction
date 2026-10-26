"""Gene-order tables (S5 synteny input): protein -> gene -> genomic position.

One table per species, one row per *input protein* (the method's gene unit),
with the columns in :data:`COLUMNS`:

=============  ======  =======================================================
protein_id     str     ID as in the proteome FASTA (the method's gene name)
gene_id        str     annotation gene ID (Ensembl gene, UniProt gene ...)
contig         str     sequence region (chromosome / scaffold / alt haplotype)
start, end     int64   1-based inclusive gene coordinates
strand         int8    +1 / -1 (0 unknown)
order          int32   0-based rank of the gene on its contig among the genes
                       *of this table* (sorted by start, end, gene_id). Gene
                       neighbourhoods (``±w genes``) are defined on this rank,
                       so unannotated/non-input genes are skipped (gap-tolerant)
contig_n       int32   number of table genes on the contig
primary        bool    contig belongs to the primary assembly (False: alt
                       haplotype / patch -- kept, but neighbours only on it)
source         str     provenance tag (``ensembl99_gff3``, ``qfo2022_gcoord``,
                       ``qfo2026_gff_md5`` ...)
=============  ======  =======================================================

Proteins without coordinates are *absent* from the table (coverage is reported
by the builders). Two proteins of the same gene (should not occur with
one-protein-per-gene proteomes) get consecutive ranks.

Builders: :func:`build_ensembl_species` (OrthoBench, Ensembl 99 GFF3) here;
QfO builders in :mod:`ancestra.genomes_qfo`. CLI: ``scripts/build_gene_order.sh``.
"""
from __future__ import annotations

import gzip
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

__all__ = ["COLUMNS", "assign_order", "parse_gff3", "parse_gff3_proteins", "gff3_seqids", "project_non_primary",
           "build_ensembl_species", "write_gene_order", "read_gene_order",
           "neighbour_csr", "shuffle_order", "to_gene_order"]

COLUMNS = ["protein_id", "gene_id", "contig", "start", "end", "strand", "order",
           "contig_n", "primary", "source"]
_DTYPES = {"start": "int64", "end": "int64", "strand": "int8", "order": "int32",
           "contig_n": "int32", "primary": "bool"}


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def assign_order(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``order`` and ``contig_n`` (rank on contig by start, end, gene_id)."""
    df = df.sort_values(["contig", "start", "end", "gene_id", "protein_id"],
                        kind="mergesort").reset_index(drop=True)
    df["order"] = df.groupby("contig", sort=False).cumcount().astype("int32")
    df["contig_n"] = df.groupby("contig", sort=False)["order"].transform("size").astype("int32")
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = {"primary": True, "source": ""}.get(c, 0)
    return df[COLUMNS + [c for c in df.columns if c not in COLUMNS]].astype(_DTYPES)


_ATTR = re.compile(r"([^=;]+)=([^;]*)")


def _attrs(s: str) -> dict:
    return dict(_ATTR.findall(s.strip()))


def gff3_seqids(path: str) -> set:
    """Sequence regions declared (``##sequence-region``) or used in a GFF3."""
    out = set()
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                out.add(line.split()[1])
            elif line and not line.startswith("#"):
                out.add(line.split("\t", 1)[0])
    return out


def parse_gff3(path: str, want: set | None = None):
    """Parse an Ensembl-style GFF3 (also the in-house simulator GFF).

    Returns ``(genes, prot)``: ``genes`` DataFrame (gene_id, name, contig, start,
    end, strand) of every feature whose ID has a ``gene:`` prefix or whose type
    ends in ``gene``; ``prot`` DataFrame (protein_id, gene_id) from CDS
    ``protein_id=`` (fallback ``ID=CDS:X``) -> transcript -> gene.
    """
    feat = {}            # ID -> (seqid, start, end, strand, parent, name, is_gene)
    cds_parent = {}      # protein -> transcript ID
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ = f[2]
            if typ in ("exon", "five_prime_UTR", "three_prime_UTR", "chromosome", "scaffold",
                       "region", "biological_region", "supercontig"):
                continue
            a = _attrs(f[8])
            if typ == "CDS":
                pid = a.get("protein_id") or a.get("ID", "").split(":", 1)[-1]
                if want is not None and pid not in want:
                    continue
                if pid and pid not in cds_parent:
                    cds_parent[pid] = a.get("Parent", "").split(",")[0]
                continue
            fid = a.get("ID")
            if fid:
                is_gene = fid.startswith("gene:") or typ.endswith("gene") or typ == "pseudogene"
                feat[fid] = (f[0], int(f[3]), int(f[4]), f[6], a.get("Parent", "").split(",")[0],
                             a.get("Name", ""), is_gene)
    genes = [(fid.split(":", 1)[-1], r[5], r[0], r[1], r[2], 1 if r[3] == "+" else (-1 if r[3] == "-" else 0))
             for fid, r in feat.items() if r[6]]
    genes = pd.DataFrame(genes, columns=["gene_id", "name", "contig", "start", "end", "strand"])
    rows = []
    for pid, tr in cds_parent.items():
        node = tr
        for _ in range(4):      # climb transcript -> gene
            rec = feat.get(node)
            if rec is None or rec[6] or not rec[4]:
                break
            node = rec[4]
        rec = feat.get(node)
        if rec is not None:
            rows.append((pid, node.split(":", 1)[-1]))
    return genes, pd.DataFrame(rows, columns=["protein_id", "gene_id"])


def parse_gff3_proteins(path: str, want: set | None = None) -> pd.DataFrame:
    """protein_id, gene_id, contig, start, end, strand (gene coordinates)."""
    genes, prot = parse_gff3(path, want)
    return prot.merge(genes.drop(columns="name"), on="gene_id", how="inner")


def _fasta_ids(path: str) -> list:
    ids = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].split()[0]
                if h.startswith(("sp|", "tr|")):
                    h = h.split("|")[1]
                ids.append(h)
    return ids


def project_non_primary(df: pd.DataFrame, genes: pd.DataFrame, primary: set) -> pd.DataFrame:
    """Place alt-haplotype / patch genes on the primary chromosome they replace.

    Ensembl alt loci (e.g. ``CHR_HSCHR6_MHC_COX_CTG1``, zebrafish ``CHR_ALT_CTG*``)
    carry their own gene models, usually with the same gene symbol as the
    primary-assembly copy. Per alt contig: match gene symbols to *unique*
    primary genes, take the majority chromosome and the median offset
    (primary start - alt start); each table gene on the contig is projected to
    its own symbol's primary start if matched, else to alt start + offset.
    Adds ``proj_contig``/``proj_start`` (NaN-free only where projected).
    Contigs with < 2 matches are left as they are (neighbours on the alt contig).
    """
    gp = genes[genes["contig"].isin(primary) & (genes["name"] != "")]
    uniq = gp.drop_duplicates("name", keep=False).set_index("name")
    alt_genes = genes[~genes["contig"].isin(primary) & (genes["name"] != "")]
    m = alt_genes.join(uniq[["contig", "start"]], on="name", rsuffix="_p", how="inner")
    df = df.copy()
    df["proj_contig"] = df["contig"]
    df["proj_start"] = df["start"]
    df["projected"] = False
    if len(m) == 0:
        return df
    name_of = genes.set_index("gene_id")["name"]
    for c, grp in m.groupby("contig"):
        if len(grp) < 2:
            continue
        chrom = grp["contig_p"].mode().iloc[0]
        g = grp[grp["contig_p"] == chrom]
        off = float(np.median(g["start_p"] - g["start"]))
        sel = (df["contig"] == c) & (~df["primary"])
        if not sel.any():
            continue
        direct = g.drop_duplicates("name").set_index("name")["start_p"]
        nm = name_of.reindex(df.loc[sel, "gene_id"]).values
        ps = np.where(pd.Series(nm).isin(direct.index).values,
                      direct.reindex(nm).values, df.loc[sel, "start"].values + off)
        df.loc[sel, "proj_contig"] = chrom
        df.loc[sel, "proj_start"] = ps.astype(np.int64)
        df.loc[sel, "projected"] = True
    return df


def build_ensembl_species(fasta: str, gff_full: str, gff_primary: str | None = None,
                          source: str = "ensembl99_gff3") -> tuple[pd.DataFrame, dict]:
    """Gene-order table for one Ensembl species.

    ``gff_full``: the GFF3 with all sequence regions (``chr_patch_hapl_scaff``
    where it exists, needed because OrthoBench includes alt-haplotype genes);
    ``gff_primary``: the primary-assembly GFF3 whose sequence regions define
    ``primary`` (None: all contigs primary).

    Ranks (``order``) are computed on primary contigs among the table's primary
    genes only. Alt-haplotype/patch genes are projected onto the primary
    chromosome (:func:`project_non_primary`): ``contig`` = that chromosome and
    ``order`` = the rank position of the projected start among primary genes
    (ties with primary genes allowed), ``raw_contig`` keeps the alt contig.
    Unprojectable alt genes stay on their own contig (ranked there).
    Protein IDs not found as CDS ``protein_id`` are matched to gene IDs after
    stripping a ``.N`` suffix (C. elegans OrthoBench uses ``WBGene...1``).
    """
    ids = _fasta_ids(fasta)
    want = set(ids)
    genes, prot = parse_gff3(gff_full, None)
    prot = prot[prot["protein_id"].isin(want)]
    miss = [i for i in ids if i not in set(prot["protein_id"])]
    gset = set(genes["gene_id"])
    extra = [(i, i.rsplit(".", 1)[0]) for i in miss if i.rsplit(".", 1)[0] in gset]
    n_by_gene = len(extra)
    if extra:
        prot = pd.concat([prot, pd.DataFrame(extra, columns=["protein_id", "gene_id"])])
    df = prot.merge(genes.drop(columns="name"), on="gene_id", how="inner")
    prim = gff3_seqids(gff_primary) if gff_primary else set(df["contig"])
    df["primary"] = df["contig"].isin(prim)
    df["source"] = source
    df = project_non_primary(df, genes, prim)
    df["raw_contig"] = df["contig"]
    # ranks: primary genes on their contig; projected genes by searchsorted of proj_start
    pr = assign_order(df[df["primary"]].copy())
    parts = [pr]
    alt = df[~df["primary"]].copy()
    projd = alt[alt["projected"]].copy()
    if len(projd):
        pos = {c: g["start"].values for c, g in pr.groupby("contig")}   # sorted by start
        projd["order"] = [int(np.searchsorted(pos.get(c, np.zeros(0)), s))
                          for c, s in zip(projd["proj_contig"], projd["proj_start"])]
        projd["contig"] = projd["proj_contig"]
        projd["start"] = projd["proj_start"]
        projd["end"] = projd["proj_start"]
        n_on = pr.groupby("contig").size()
        projd["contig_n"] = n_on.reindex(projd["contig"]).fillna(0).astype(int).values
        parts.append(projd)
    rest = alt[~alt["projected"]].copy()
    if len(rest):
        parts.append(assign_order(rest))
    df = pd.concat(parts, ignore_index=True)
    df = df[COLUMNS + ["raw_contig", "projected"]].astype(_DTYPES)
    df = df.sort_values(["contig", "order", "protein_id"], kind="mergesort").reset_index(drop=True)
    n_alt = int((~df["primary"]).sum())
    cov = dict(n_proteins=len(ids), n_mapped=int(len(df)),
               coverage=round(len(df) / max(len(ids), 1), 6),
               n_mapped_via_gene_id=n_by_gene,
               n_non_primary=n_alt, n_projected=int(df["projected"].sum()),
               n_non_primary_unprojected=int((~df["primary"] & ~df["projected"]).sum()),
               n_contigs=int(df["contig"].nunique()),
               n_genes_on_contigs_ge10=int((df["contig_n"] >= 10).sum()),
               missing_examples=sorted(want - set(df["protein_id"]))[:10])
    return df, cov


def write_gene_order(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)


def read_gene_order(path: str) -> pd.DataFrame:
    return pd.read_parquet(path)


# ----------------------------------------------------------------------------
# neighbourhoods on integer gene IDs
# ----------------------------------------------------------------------------

def neighbour_csr(contig_code: np.ndarray, order: np.ndarray, w: int):
    """Neighbour lists within ``±w`` ranks on the same contig.

    ``contig_code[i]``, ``order[i]`` for genes i = 0..n-1 (contig_code < 0 =
    no coordinates -> no neighbours). Neighbours of i: genes j != i with the
    same contig and |order_j - order_i| <= w (ties allowed, e.g. projected
    alt-haplotype genes). Returns CSR (indptr, indices), rows sorted. O(n·w).
    """
    n = len(order)
    ok = contig_code >= 0
    idx = np.flatnonzero(ok)
    key = contig_code[idx].astype(np.int64) * (1 << 32) + order[idx].astype(np.int64)
    srt = np.argsort(key, kind="stable")
    ks, gs = key[srt], idx[srt]
    lo = np.searchsorted(ks, ks - w, "left")
    hi = np.searchsorted(ks, ks + w, "right")
    cnt = hi - lo
    rows = np.repeat(gs, cnt)
    offs = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    cols = gs[np.repeat(lo, cnt) + offs]
    keep = rows != cols
    rows, cols = rows[keep], cols[keep]
    o = np.lexsort((cols, rows))
    rows, cols = rows[o], cols[o]
    indptr = np.zeros(n + 1, np.int64)
    np.add.at(indptr, rows + 1, 1)
    return np.cumsum(indptr), cols.astype(np.int64)


def shuffle_order(contig_code: np.ndarray, order: np.ndarray, seed: int):
    """Negative control: permute gene positions within the genome.

    The multiset of (contig, rank) slots is kept, genes are assigned to slots
    uniformly at random (genes without coordinates stay without).
    """
    rng = np.random.default_rng(seed)
    has = np.flatnonzero(contig_code >= 0)
    perm = rng.permutation(has)
    cc, oo = contig_code.copy(), order.copy()
    cc[perm], oo[perm] = contig_code[has], order[has]
    return cc, oo


def to_gene_order(gt, table_dir: str, contig_col: str = "contig"):
    """Gene-order tables -> the core's :class:`ancestra.genome.GeneOrder` for ``gt``.

    Uses ``table_dir/<species>.parquet`` (this module's schema; alt-haplotype genes
    already projected onto their primary chromosome). ``rank`` = ``order``;
    ``start`` = (projected) start. Genes without a row get contig -1.
    """
    from .genome import GeneOrder
    N = gt.n
    contig = np.full(N, -1, np.int64)
    rank = np.full(N, -1, np.int64)
    start = np.zeros(N, np.int64)
    strand = np.zeros(N, np.int8)
    names, coverage = [], {}
    idx = gt.index()
    for si, sp in enumerate(gt.species):
        p = os.path.join(table_dir, f"{sp}.parquet")
        nsp = int((gt.species_of == si).sum())
        if not os.path.exists(p):
            coverage[sp] = 0.0
            continue
        df = pd.read_parquet(p, columns=["protein_id", contig_col, "order", "start", "strand"])
        gi = df["protein_id"].map(idx)
        ok = gi.notna().values
        df, gi = df[ok], gi[ok].astype(np.int64).values
        codes, uniq = pd.factorize(df[contig_col])
        contig[gi] = codes + len(names)
        names += [f"{sp}:{c}" for c in uniq]
        rank[gi] = df["order"].values
        start[gi] = df["start"].values
        strand[gi] = df["strand"].values
        coverage[sp] = len(gi) / max(nsp, 1)
    return GeneOrder(contig, rank, start, strand, names, coverage)
