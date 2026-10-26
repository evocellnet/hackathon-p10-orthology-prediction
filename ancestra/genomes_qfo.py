"""Gene-order tables for the QfO 2022_02 (updated zebrafish) proteomes.

Schema: :data:`ancestra.genomes.COLUMNS` (one row per input protein with
coordinates; see :mod:`ancestra.genomes`). Two coordinate sources, in order:

1. ``qfo2022_gcoord`` -- the 2022 ``UPxxx_tax_gcoord.xml`` (19 eukaryotes).
   Format (UniProt "proteins API" coordinates): ``<gnEntries>`` of
   ``<gnEntry>`` (``accession``, ``sequence``, ...), each with >= 0
   ``<gnCoordinate ensembl_gene_id= ensembl_transcript_id= ensembl_translation_id=>``
   holding one ``<genomicLocation chromosome= start= end= reverse_strand=>``
   (start > end on the reverse strand) with per-exon ``genomeLocation``s. The
   files cover *all* UniProt entries of the organism (e.g. 86k for human), so
   only accessions of the canonical 2022 FASTA are kept. A protein's gene is
   the Ensembl gene of its gnCoordinates (several genes, e.g. reference +
   patch copies: primary contig first, then most transcripts, then gene ID);
   its span = min/max over that gene's genomicLocations. Entries with an empty
   ``chromosome`` (Helobdella, Phytophthora) are treated as uncoordinated.
   RefSeq-based 2026 files (zebrafish, Xenopus, amphioxus) have no gene IDs
   (``refseq_protein_id``) but a ``nucleotideId`` (NC_/NW_), used as contig so
   that they agree with the RefSeq GFF of the same assembly.
2. ``qfo2026_gff_md5`` -- the 2026_02 per-species GFF (ENA/EMBL, RefSeq or
   Ensembl), *only* for proteins whose 2022 sequence MD5 equals the 2026
   canonical sequence MD5 of the same accession (2026 FASTA from
   ``QfO_release_2026_02.tar.gz``). GFF CDS features are keyed by
   ``protein_id`` / ``Name`` / ``ID`` (version-stripped) and by
   ``Dbxref=UniProtKB/*:ACC``; an accession is matched directly via the
   UniProtKB Dbxref or via the 2026 ``.idmapping`` (EMBL-CDS, RefSeq,
   Ensembl_PRO, EnsemblGenome_PRO ...). The gene = top of the CDS's ``Parent``
   chain. For gcoord species this only fills gaps, and only if the 2026 GFF
   agrees with gcoord on the proteins both place (same assembly; see
   ``MIN_CONCORDANCE``); otherwise mixing assemblies would corrupt gene order.
   Where the 2026 release has no GFF but a GenBank ``.gbff`` (Nitrosopumilus,
   Phytophthora, Helobdella, Anopheles), the same rule applies with
   source ``qfo2026_gbff_md5``. The 2026 proteome is the same UP id, else the
   one with the same taxid (dog, Ixodes, Fusobacterium, Phaeosphaeria).
3. ``qfo2026_gcoord_md5`` -- the 2026 ``_gcoord.xml`` (21 species), for
   accessions whose gcoord ``<sequence>`` MD5 equals the 2022 MD5 (covers
   accessions that are no longer canonical in 2026, e.g. zebrafish).
4. ``qfo2026_gff_vxref`` -- 2026 GFF/GBFF CDS whose *versioned* protein ID
   (RefSeq/EMBL-CDS/Ensembl ``X.N``) equals a versioned xref of the accession
   in the 2022 ``.idmapping`` (a versioned ID pins the sequence, so this is an
   MD5-equivalent identity; covers accessions dropped in 2026).
5. ``eg52_gff3`` -- optional 2022-era GFF3 (Ensembl Genomes 52, same release
   basis as QfO 2022) in ``data/qfo/annotation_2022/UPxxx.gff3.gz``, mapped via
   the 2022 ``.idmapping``; ranks right after the 2022 gcoord (S. pombe).

Layers 2/3 never mix assemblies: see :func:`build_species`.

``primary`` is False for contigs that look like alt haplotypes / patches
(:func:`is_primary_contig`).

CLI::

    python -m ancestra.genomes_qfo --out data/genomes/qfo2022 [--species UP000005640 ...]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from .genomes import COLUMNS, assign_order, write_gene_order

__all__ = ["read_fasta", "seq_md5", "is_primary_contig", "parse_gcoord",
           "parse_gff_cds_genes", "parse_gbff_cds_genes", "index_annotation", "read_idmapping", "map_gff", "build_species",
           "discover_species", "COVERAGE_COLUMNS"]

QFO2022 = "data/qfo/proteomes_2022_02_updZEBRAFISH"
QFO2026 = "data/qfo/proteomes_current"
ANNOT2022 = "data/qfo/annotation_2022"     # optional 2022-era GFF3 (Ensembl Genomes 52)
DOMAINS = ("Archaea", "Bacteria", "Eukaryota")
MIN_CONCORDANCE = 0.95
COVERAGE_COLUMNS = ["species", "domain", "n_proteins", "n_gcoord2022", "n_gff2026_md5",
                    "n_md5_mismatch", "n_no_coord", "coverage", "n_contigs", "n_non_primary",
                    "notes", "n_not_in_2026", "n_md5_ok_no_coord", "n_gcoord2026_md5", "n_gff2026_vxref", "n_eg52_gff3",
                    "gff2026"]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# ----------------------------------------------------------------------------
# sequences
# ----------------------------------------------------------------------------

def _acc(header: str) -> str:
    h = header[1:].split()[0] if header.startswith(">") else header.split()[0]
    if h.startswith(("sp|", "tr|")):
        h = h.split("|")[1]
    return h


def read_fasta(path: str) -> dict:
    """accession -> sequence (uppercase, whitespace removed); file order kept."""
    out, acc, buf = {}, None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    out[acc] = "".join(buf).upper()
                acc, buf = _acc(line), []
            else:
                buf.append(line.strip())
    if acc is not None:
        out[acc] = "".join(buf).upper()
    return out


def seq_md5(seq: str) -> str:
    return hashlib.md5(seq.strip().upper().encode()).hexdigest()


_NONPRIMARY = re.compile(r"(^CHR_|^HSCHR|^HG\d+_|PATCH|_alt$|_ALT$|_hap\d|_HAP\d|NOVEL|"
                         r"^alt[-_]|_fix$|_FIX$)", re.I)


def is_primary_contig(name: str) -> bool:
    """False for contig names that denote alt haplotypes / patches.

    ``CHR_HG1_PATCH``, ``CHR_HSCHR6_MHC_COX_CTG1``, ``HSCHR1_1_CTG3``,
    ``chr6_GL000250v2_alt``, ``*_fix``, ``*_hap1`` -> False; chromosomes,
    unplaced scaffolds (``JH584299.1``, ``KZ626836.1``) -> True.
    """
    return not bool(_NONPRIMARY.search(str(name)))


# ----------------------------------------------------------------------------
# 2022 gcoord XML
# ----------------------------------------------------------------------------

def _local(tag) -> str:
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else ""


def parse_gcoord(path: str, want: set | None = None) -> pd.DataFrame:
    """Gene coordinates per accession from a QfO ``_gcoord.xml`` (streaming).

    Returns protein_id, gene_id, contig, start, end, strand, n_genes (number of
    distinct genes the entry maps to), seq_md5 (of the gcoord ``<sequence>``).
    Entries without a gnCoordinate or with an empty chromosome are dropped.
    """
    from lxml import etree

    rows = []
    src = gzip.open(path, "rb") if str(path).endswith(".gz") else path
    ctx = etree.iterparse(src, events=("end",), tag="gnEntry", huge_tree=True)
    for _, el in ctx:
        acc = el.findtext("accession")
        if acc and (want is None or acc in want):
            # (gene key, contig) -> [starts, ends, strand Counter, n_locations, label]
            genes = {}
            for gc in el.iterfind("gnCoordinate"):
                # Ensembl-based files carry gene IDs; RefSeq-based 2026 files only
                # refseq_protein_id -> all their locations form one per-accession gene
                gkey = gc.get("ensembl_gene_id") or ""
                label = gkey or gc.get("refseq_protein_id") or gc.get("ensembl_transcript_id") or acc
                for gl in gc.iterfind("genomicLocation"):
                    chrom = (gl.get("nucleotideId") or gl.get("chromosome") or "").strip()
                    if not chrom or gl.get("start") is None or gl.get("end") is None:
                        continue
                    s, e = int(gl.get("start")), int(gl.get("end"))
                    rev = (gl.get("reverse_strand") or "").lower() == "true"
                    g = genes.setdefault((gkey, chrom), [[], [], Counter(), 0, label])
                    g[0].append(min(s, e)); g[1].append(max(s, e))
                    g[2][-1 if rev else 1] += 1
                    g[3] += 1
            if genes:
                cand = sorted((not is_primary_contig(chrom), -g[3], g[4], chrom, min(g[0]),
                               max(g[1]), g[2].most_common(1)[0][0])
                              for (gkey, chrom), g in genes.items())
                _, _, gid, chrom, s, e, strand = cand[0]
                seq = el.findtext("sequence") or ""
                rows.append((acc, gid, chrom, s, e, strand, len(cand),
                             seq_md5(seq) if seq else ""))
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
    del ctx
    if src is not path:
        src.close()
    return pd.DataFrame(rows, columns=["protein_id", "gene_id", "contig", "start", "end",
                                       "strand", "n_genes", "seq_md5"])


# ----------------------------------------------------------------------------
# 2026 GFF
# ----------------------------------------------------------------------------

_ATTR = re.compile(r"([^=;]+)=([^;]*)")
_RANK = {"pid": 0, "pidv": 1, "uniprot": 2}
_SKIP = {"exon", "five_prime_UTR", "three_prime_UTR", "chromosome", "scaffold", "region",
         "biological_region", "cDNA_match", "match", "sequence_feature", "start_codon",
         "stop_codon", "intron", "repeat_region", "source", "contig"}
_GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene", "protein_coding_gene", "transposable_element_gene"}


_VER = re.compile(r"\.\d+$")


def _nover(x: str) -> str:
    """Strip a trailing ``.N`` version (ENSP0001.2 -> ENSP0001, AAC73112.1 -> AAC73112)."""
    return re.sub(r"\.\d+$", "", x)


def _strip_prefix(x: str) -> str:
    for p in ("cds-", "CDS:", "cds:", "cds."):
        if x.startswith(p):
            return x[len(p):]
    return x


def parse_gff_cds_genes(path: str) -> tuple[dict, dict]:
    """Index a GFF(3): CDS protein keys -> gene record.

    Returns (key2genes, genes): ``genes[gene_id] = (contig, start, end, strand)``;
    ``key2genes[key] = {gene_id: kind}`` where key is a version-stripped
    protein_id / Name / ID of a CDS (kind 'pid') or a UniProtKB accession from
    its Dbxref (kind 'uniprot'). Gene = top of the CDS ``Parent`` chain.
    """
    feat = {}               # ID -> (seqid, start, end, strand, parent, type)
    cds = {}                # cds key tuple -> parent
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ = f[2]
            if typ in _SKIP:
                continue
            a = dict(_ATTR.findall(f[8].strip()))
            parent = a.get("Parent", "").split(",")[0]
            if typ == "CDS":
                keys = set()
                for k in ("protein_id", "Name", "ID"):
                    v = a.get(k)
                    if v:
                        keys.add(("pid", _nover(_strip_prefix(v))))
                pv = a.get("protein_id", "")
                if pv and not _VER.search(pv) and a.get("version", "").isdigit():
                    pv = f"{pv}.{a['version']}"          # Ensembl: protein_id + version=
                if _VER.search(pv):
                    keys.add(("pidv", pv))
                for x in a.get("Dbxref", "").split(","):
                    if x.startswith(("UniProtKB/Swiss-Prot:", "UniProtKB/TrEMBL:", "UniProtKB:")):
                        keys.add(("uniprot", x.split(":", 1)[1]))
                    elif x.startswith(("NCBI_GP:", "GenBank:", "EMBL:", "RefSeq:")):
                        keys.add(("pid", _nover(x.split(":", 1)[1])))
                        if _VER.search(x):
                            keys.add(("pidv", x.split(":", 1)[1]))
                if not parent and a.get("ID"):
                    # CDS without parent (some ENA files): the CDS is its own gene
                    fid = "cds_gene:" + a["ID"]
                    old = feat.get(fid)
                    s, e = int(f[3]), int(f[4])
                    if old:
                        s, e = min(s, old[1]), max(e, old[2])
                    feat[fid] = (f[0], s, e, f[6], "", "gene")
                    parent = fid
                cds.setdefault(frozenset(keys), parent)
                continue
            fid = a.get("ID")
            if fid:
                feat[fid] = (f[0], int(f[3]), int(f[4]), f[6], parent, typ)
    genes, key2genes = {}, defaultdict(dict)
    for keys, par in cds.items():
        node = par
        for _ in range(6):
            rec = feat.get(node)
            if rec is None or not rec[4] or rec[5] in _GENE_TYPES:
                break
            node = rec[4]
        rec = feat.get(node)
        if rec is None:
            continue
        seqid, s, e, st, _, _ = rec
        gid = re.sub(r"^(gene[-:]|cds_gene:)", "", node)
        genes[gid] = (seqid, s, e, 1 if st == "+" else (-1 if st == "-" else 0))
        for kind, k in keys:
            if _RANK[kind] >= _RANK.get(key2genes[k].get(gid), -1):
                key2genes[k][gid] = kind
    return dict(key2genes), genes


_QUAL = re.compile(r'^\s{21}/(locus_tag|gene|protein_id|db_xref|old_locus_tag)="([^"]*)"')


def parse_gbff_cds_genes(path: str) -> tuple[dict, dict]:
    """Same output as :func:`parse_gff_cds_genes` from a GenBank flat file
    (used where the 2026 release ships a ``.gbff`` but no GFF).

    Contig = ``VERSION`` accession (fallback LOCUS name). Gene = ``gene``
    feature with the CDS's ``locus_tag`` (else ``gene``) on the same record;
    without one, the union of that gene key's CDS spans. CDS keys:
    ``protein_id`` (version-stripped) and ``db_xref`` UniProtKB accessions.
    """
    gene_feat = {}                  # (contig, gkey) -> [s, e, strand]
    cds_feat = {}                   # (contig, gkey) -> [s, e, strand]
    cds_keys = defaultdict(set)     # (contig, gkey) -> {(kind, key)}

    def flush(rec, contig):
        if rec is None or rec[0] not in ("gene", "CDS"):
            return
        typ, loc, q = rec
        nums = [int(x) for x in re.findall(r"\d+", loc)]
        if not nums:
            return
        gkey = q.get("locus_tag") or q.get("gene") or q.get("protein_id")
        if not gkey:
            return
        st = -1 if loc.lstrip("<>").startswith("complement") else 1
        k = (contig, gkey)
        tgt = gene_feat if typ == "gene" else cds_feat
        old = tgt.get(k)
        s, e = min(nums), max(nums)
        tgt[k] = [min(s, old[0]), max(e, old[1]), st] if old else [s, e, st]
        if typ == "CDS":
            if q.get("protein_id"):
                cds_keys[k].add(("pid", _nover(q["protein_id"])))
                if _VER.search(q["protein_id"]):
                    cds_keys[k].add(("pidv", q["protein_id"]))
            for x in q.get("db_xref", []):
                if x.startswith(("UniProtKB/Swiss-Prot:", "UniProtKB/TrEMBL:")):
                    cds_keys[k].add(("uniprot", x.split(":", 1)[1]))

    contig, rec, in_feat, in_loc = None, None, False, False
    with _open(path) as fh:
        for line in fh:
            if not in_feat:
                if line.startswith("LOCUS"):
                    contig = line.split()[1]
                elif line.startswith("VERSION") and len(line.split()) > 1:
                    contig = line.split()[1]
                elif line.startswith("FEATURES"):
                    in_feat = True
                continue
            if line.startswith(("ORIGIN", "CONTIG", "//")):
                flush(rec, contig)
                rec, in_feat, in_loc = None, False, False
                continue
            if len(line) > 5 and line[5] != " " and line.startswith("     "):
                flush(rec, contig)
                p = line.split()
                rec = (p[0], p[1] if len(p) > 1 else "", {"db_xref": []})
                in_loc = True
                continue
            if rec is None:
                continue
            body = line[21:].rstrip("\n")
            if body.startswith("/"):
                in_loc = False
                m = _QUAL.match(line)
                if m:
                    if m.group(1) == "db_xref":
                        rec[2]["db_xref"].append(m.group(2))
                    else:
                        rec[2].setdefault(m.group(1), m.group(2))
            elif in_loc:
                rec = (rec[0], rec[1] + body.strip(), rec[2])
    genes, key2genes = {}, defaultdict(dict)
    for k, keys in cds_keys.items():
        contig, gkey = k
        s, e, st = gene_feat.get(k) or cds_feat[k]
        gid = gkey if gkey not in genes or genes[gkey][0] == contig else f"{contig}:{gkey}"
        genes[gid] = (contig, s, e, st)
        for kind, key in keys:
            if _RANK[kind] >= _RANK.get(key2genes[key].get(gid), -1):
                key2genes[key][gid] = kind
    return dict(key2genes), genes


def index_annotation(path: str) -> tuple[dict, dict]:
    """GFF(3) or GenBank flat file -> (key2genes, genes)."""
    if re.search(r"\.(gbff|gbk|gb)(\.gz)?$", path):
        return parse_gbff_cds_genes(path)
    return parse_gff_cds_genes(path)


_IDMAP_TYPES = {"EMBL-CDS", "RefSeq", "Ensembl_PRO", "EnsemblGenome_PRO", "EnsemblGenome_TRS",
                "Ensembl_TRS", "GeneID", "EnsemblGenome", "Ensembl", "Gene_OrderedLocusName",
                "Gene_ORFName", "Gene_Name"}


def read_idmapping(path: str, accs: set | None = None, keep_version: bool = False) -> dict:
    """UniProt ``.idmapping`` (acc, type, id) -> acc -> {(type, id)}; ids
    version-stripped unless ``keep_version`` (then only versioned protein IDs)."""
    out = defaultdict(set)
    if not path or not os.path.exists(path):
        return out
    with _open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[1] not in _IDMAP_TYPES or p[2] in ("", "-"):
                continue
            acc = p[0].split("-")[0] if p[0].count("-") == 1 and p[0].split("-")[1].isdigit() else p[0]
            if accs is not None and acc not in accs:
                continue
            if keep_version:
                if p[1] in _PROTEIN_TYPES and _VER.search(p[2]):
                    out[acc].add((p[1], p[2]))
            else:
                out[acc].add((p[1], _nover(p[2])))
    return out


def map_versioned_xref(accs, key2genes: dict, genes: dict, idmap_v: dict) -> pd.DataFrame:
    """Accession -> GFF gene via a *versioned* protein xref of the 2022
    release (RefSeq NP_/XP_x.N, EMBL-CDS ABC123.N, Ensembl_PRO ENSP.N) equal
    to a versioned CDS protein_id of the 2026 annotation. Versioned INSDC /
    RefSeq / Ensembl protein IDs change whenever the sequence changes, so an
    exact versioned match implies the 2022 sequence is the annotated one."""
    rows = []
    for acc in accs:
        votes = Counter()
        for _, x in idmap_v.get(acc, ()):
            for gid, kind in key2genes.get(x, {}).items():
                if kind == "pidv":
                    votes[gid] += 1
        if votes:
            best = sorted(votes, key=lambda g: (not is_primary_contig(genes[g][0]), -votes[g], g))[0]
            c, s0, e, st = genes[best]
            rows.append((acc, best, c, s0, e, st))
    return pd.DataFrame(rows, columns=["protein_id", "gene_id", "contig", "start", "end", "strand"])


_PROTEIN_TYPES = ("EMBL-CDS", "RefSeq", "Ensembl_PRO", "EnsemblGenome_PRO")
_GENE_ID_TYPES = ("Gene_OrderedLocusName", "EnsemblGenome", "Ensembl", "Gene_ORFName")


def map_gff(accs, key2genes: dict, genes: dict, idmap: dict) -> tuple[pd.DataFrame, Counter]:
    """Accession -> GFF gene. Priority: UniProtKB Dbxref on the CDS, then
    protein xrefs (EMBL-CDS/RefSeq/Ensembl_PRO/EnsemblGenome_PRO), then
    transcript xrefs, then gene-level xrefs (locus tag / Ensembl gene) matched
    against GFF gene IDs. Ties: primary contig, most votes, gene ID."""
    stats = Counter()
    gene_by_name = {}
    for gid in genes:
        gene_by_name.setdefault(_nover(gid), gid)
    rows = []
    for acc in accs:
        votes = Counter()
        hit = key2genes.get(acc)
        if hit:
            for gid, kind in hit.items():
                if kind == "uniprot":
                    votes[gid] += 100
        how = "uniprot" if votes else None
        if not votes:
            for typ, x in idmap.get(acc, ()):
                if typ in _PROTEIN_TYPES or typ.endswith("_TRS"):
                    for gid in key2genes.get(x, {}):
                        votes[gid] += 10 if typ in _PROTEIN_TYPES else 1
            how = "xref" if votes else None
        if not votes:
            for typ, x in idmap.get(acc, ()):
                if typ in _GENE_ID_TYPES and x in gene_by_name:
                    votes[gene_by_name[x]] += 1
            how = "gene_xref" if votes else None
        if not votes:
            stats["unmapped"] += 1
            continue
        stats[how] += 1
        if len(votes) > 1:
            stats["multi_gene"] += 1
        best = sorted(votes, key=lambda g: (not is_primary_contig(genes[g][0]), -votes[g], g))[0]
        c, s, e, st = genes[best]
        rows.append((acc, best, c, s, e, st))
    return pd.DataFrame(rows, columns=["protein_id", "gene_id", "contig", "start", "end",
                                       "strand"]), stats


# ----------------------------------------------------------------------------
# species discovery and build
# ----------------------------------------------------------------------------

def _gff_for(root2026: str, upid: str):
    """2026 annotation for a UP id: ``gff/UPxxx.gff(3)(.gz)``, else ``gbff/UPxxx.gbff.gz``."""
    for sub, exts in (("gff", (".gff3.gz", ".gff.gz", ".gff3", ".gff")),
                      ("gbff", (".gbff.gz", ".gbff"))):
        for ext in exts:
            p = os.path.join(root2026, sub, upid + ext)
            if os.path.exists(p):
                return p
    return None


def discover_species(root2022: str = QFO2022, root2026: str = QFO2026,
                     annot2022: str = ANNOT2022) -> list[dict]:
    """One dict per 2022 species: paths of the 2022 FASTA/gcoord and the
    2026 FASTA/idmapping/GFF (2026 proteome = same UP id, else same taxid)."""
    fasta_re = re.compile(r"^(UP\d+)_(\d+)\.fasta$")
    fasta26_re = re.compile(r"^(UP\d+)_(\d+)\.fasta(\.gz)?$")
    idx26 = {}
    for dom in DOMAINS:
        d = os.path.join(root2026, dom)
        if os.path.isdir(d):
            for fn in os.listdir(d):
                m = fasta26_re.match(fn)
                if m:
                    idx26[m.group(1)] = (dom, m.group(1), m.group(2), os.path.join(d, fn))
    by_tax = {v[2]: v for v in idx26.values()}
    out = []
    for dom in DOMAINS:
        d = os.path.join(root2022, dom)
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            m = fasta_re.match(fn)
            if not m:
                continue
            up, tax = m.groups()
            sp = dict(species=f"{up}_{tax}", upid=up, taxid=tax, domain=dom,
                      fasta2022=os.path.join(d, fn), gcoord=None, fasta2026=None,
                      idmap2026=None, gff2026=None, upid2026=None, gcoord2026=None,
                      gff2022=None)
            gc = os.path.join(d, f"{up}_{tax}_gcoord.xml")
            if os.path.exists(gc):
                sp["gcoord"] = gc
            hit = idx26.get(up) or by_tax.get(tax)
            if hit:
                dom26, up26, tax26, fa26 = hit
                sp["upid2026"] = up26
                sp["fasta2026"] = fa26
                im = os.path.join(root2026, dom26, f"{up26}_{tax26}.idmapping")
                sp["idmap2026"] = next((x for x in (im + ".gz", im) if os.path.exists(x)), None)
                g26 = os.path.join(root2026, dom26, f"{up26}_{tax26}_gcoord.xml")
                sp["gcoord2026"] = next((x for x in (g26 + ".gz", g26) if os.path.exists(x)), None)
            g22 = os.path.join(annot2022, up + ".gff3.gz")
            sp["gff2022"] = g22 if os.path.exists(g22) else None
            up_gff = sp["upid2026"] or up
            sp["gff2026"] = _gff_for(root2026, up_gff)
            out.append(sp)
    return out


def _concordance(a: pd.DataFrame, b: pd.DataFrame) -> tuple[float, int]:
    """Fraction of shared proteins placed on the same contig with overlapping spans."""
    m = a.merge(b, on="protein_id", suffixes=("_a", "_b"))
    if len(m) == 0:
        return float("nan"), 0
    ok = ((m["contig_a"].astype(str) == m["contig_b"].astype(str))
          & (m["start_a"] <= m["end_b"]) & (m["start_b"] <= m["end_a"]))
    return float(ok.mean()), int(len(m))


MIN_SHARED = 20     # shared proteins needed to verify that two layers use one assembly


def build_species(sp: dict) -> tuple[pd.DataFrame, dict]:
    """Gene-order table + coverage row for one species dict (see discover_species).

    Layers: 2022 gcoord (all its proteins), then the 2026 layers (GFF/GBFF via
    the 2026 canonical FASTA MD5; 2026 gcoord via its own ``<sequence>`` MD5),
    largest first. A later layer only fills proteins not yet placed, and only if
    it agrees with what is already placed (>= MIN_CONCORDANCE same contig and
    overlapping span on >= MIN_SHARED shared proteins) -- or nothing is placed
    yet. This keeps every table on a single assembly.
    """
    seq22 = read_fasta(sp["fasta2022"])
    accs = list(seq22)
    want = set(accs)
    md5_22 = {a: seq_md5(x) for a, x in seq22.items()}
    notes = []
    empty = pd.DataFrame(columns=["protein_id", "gene_id", "contig", "start", "end", "strand"])
    gc = empty.copy()
    if sp.get("gcoord"):
        gc = parse_gcoord(sp["gcoord"], want)
        if len(gc):
            mism = int((gc["seq_md5"] != gc["protein_id"].map(md5_22)).sum())
            if mism:
                notes.append(f"gcoord_seq_differs={mism}")
            amb = int((gc["n_genes"] > 1).sum())
            if amb:
                notes.append(f"gcoord_multi_gene={amb}")
        else:
            notes.append("gcoord2022_unusable(empty chromosome)")
        gc = gc.drop(columns=["n_genes", "seq_md5"])
    if sp.get("upid2026") and sp["upid2026"] != sp["upid"]:
        notes.append(f"2026_proteome={sp['upid2026']}")

    # --- 2026 candidate layers (each over all md5-matching proteins) ---
    seq26 = read_fasta(sp["fasta2026"]) if sp.get("fasta2026") else {}
    md5_26 = {a: seq_md5(x) for a, x in seq26.items() if a in want}
    same_fa = [a for a in accs if md5_26.get(a) == md5_22[a]]
    gc26_md5 = {}
    layers = []
    if sp.get("gff2026") and same_fa:
        key2genes, genes = index_annotation(sp["gff2026"])
        idmap = read_idmapping(sp.get("idmap2026"), set(same_fa))
        g26, st = map_gff(same_fa, key2genes, genes, idmap)
        tag = "qfo2026_gbff_md5" if ".gbff" in sp["gff2026"] else "qfo2026_gff_md5"
        notes.append(f"{tag.split('_')[1]}_match:" + ",".join(
            f"{k}={v}" for k, v in sorted(st.items())))
        layers.append((tag, g26))
    if sp.get("gff2026"):
        # versioned-xref layer for accessions without an MD5-verified 2026 sequence
        if not same_fa:
            key2genes, genes = index_annotation(sp["gff2026"])
        im22 = sp["fasta2022"][:-len(".fasta")] + ".idmapping"
        idv = read_idmapping(im22, want, keep_version=True)
        gv = map_versioned_xref(accs, key2genes, genes, idv)
        if len(gv):
            layers.append(("qfo2026_gff_vxref", gv))
    else:
        notes.append("no_2026_gff")
    if sp.get("gcoord2026"):
        g = parse_gcoord(sp["gcoord2026"], want)
        gc26_md5 = dict(zip(g["protein_id"], g["seq_md5"]))
        g = g[g["seq_md5"] == g["protein_id"].map(md5_22)]
        if len(g):
            layers.append(("qfo2026_gcoord_md5", g.drop(columns=["n_genes", "seq_md5"])))
        else:
            notes.append("gcoord2026_unusable")
    layers.sort(key=lambda t: -len(t[1]))

    tier22 = []
    if sp.get("gff2022"):
        k22, g22 = index_annotation(sp["gff2022"])
        im22 = read_idmapping(sp["fasta2022"][:-len(".fasta")] + ".idmapping", want)
        e52, st = map_gff(accs, k22, g22, im22)
        notes.append("eg52_match:" + ",".join(f"{k}={v}" for k, v in sorted(st.items())))
        tier22.append(("eg52_gff3", e52))

    placed = gc.assign(source="qfo2022_gcoord")
    counts = Counter()
    for tag, lay in tier22 + layers:
        if len(placed):
            conc, n_sh = _concordance(placed, lay)
            notes.append(f"{tag}_vs_placed_concordance={conc:.3f}(n={n_sh})")
            if n_sh < MIN_SHARED or not (conc >= MIN_CONCORDANCE):
                notes.append(f"{tag}_skipped(assembly_not_verified)")
                continue
        add = lay[~lay["protein_id"].isin(set(placed["protein_id"]))].assign(source=tag)
        counts[tag] = len(add)
        placed = pd.concat([placed, add], ignore_index=True) if len(placed) else add

    # --- why the unplaced proteins have no coordinates (partition of n_no_coord) ---
    done = set(placed["protein_id"])
    rest = [a for a in accs if a not in done]
    in26 = [a for a in rest if a in md5_26 or a in gc26_md5]
    ok26 = {a for a in in26 if md5_26.get(a) == md5_22[a] or gc26_md5.get(a) == md5_22[a]}
    n_absent = len(rest) - len(in26)      # accession not in the 2026 FASTA / gcoord
    n_mism = len(in26) - len(ok26)        # 2026 sequence differs -> excluded
    n_ok_nocoord = len(ok26)              # same sequence, but no usable 2026 coordinate

    df = placed if len(placed) else empty.assign(source=[])
    df = df.copy()
    df["strand"] = df["strand"].astype("int8")
    df["primary"] = df["contig"].map(is_primary_contig).astype(bool) if len(df) else \
        pd.Series([], dtype=bool)
    df = assign_order(df)
    dup = int(df["gene_id"].duplicated().sum())
    if dup:
        notes.append(f"proteins_sharing_gene={dup}")
    nonprim = sorted(df.loc[~df["primary"], "contig"].unique())
    if nonprim:
        notes.append("non_primary_contigs=" + ",".join(nonprim[:5])
                     + (f",...({len(nonprim)})" if len(nonprim) > 5 else ""))
    n = len(accs)
    cov = dict(species=sp["species"], domain=sp["domain"], n_proteins=n,
               n_gcoord2022=int(len(gc)),
               n_gff2026_md5=counts["qfo2026_gff_md5"] + counts["qfo2026_gbff_md5"],
               n_md5_mismatch=n_mism, n_no_coord=n - len(df),
               coverage=round(len(df) / max(n, 1), 4), n_contigs=int(df["contig"].nunique()),
               n_non_primary=int((~df["primary"]).sum()), notes="; ".join(notes),
               n_not_in_2026=n_absent, n_md5_ok_no_coord=n_ok_nocoord,
               n_gcoord2026_md5=counts["qfo2026_gcoord_md5"],
               n_gff2026_vxref=counts["qfo2026_gff_vxref"],
               n_eg52_gff3=counts["eg52_gff3"],
               gff2026=os.path.basename(sp["gff2026"]) if sp.get("gff2026") else "")
    return df, cov


def _run_one(args):
    sp, out = args
    df, cov = build_species(sp)
    write_gene_order(df, os.path.join(out, sp["species"] + ".parquet"))
    print(f"{sp['species']}\t{cov['coverage']:.3f}\t{cov['notes']}", file=sys.stderr, flush=True)
    return cov


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="data/genomes/qfo2022")
    ap.add_argument("--qfo2022", default=QFO2022)
    ap.add_argument("--qfo2026", default=QFO2026)
    ap.add_argument("--species", nargs="*", help="UP ids or UPid_taxid (default: all)")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    sps = discover_species(a.qfo2022, a.qfo2026)
    if a.species:
        sel = set(a.species)
        sps = [s for s in sps if s["upid"] in sel or s["species"] in sel]
    # big gcoord files first for load balance
    sps.sort(key=lambda s: -(os.path.getsize(s["gcoord"]) if s["gcoord"] else 0))
    os.makedirs(a.out, exist_ok=True)
    with ProcessPoolExecutor(max_workers=max(1, a.threads)) as ex:
        covs = list(ex.map(_run_one, [(s, a.out) for s in sps]))
    cov = pd.DataFrame(covs, columns=COVERAGE_COLUMNS).sort_values(["domain", "species"])
    path = os.path.join(a.out, "coverage.tsv")
    if a.species and os.path.exists(path):
        old = pd.read_csv(path, sep="\t")
        cov = pd.concat([old[~old["species"].isin(cov["species"])], cov]).sort_values(
            ["domain", "species"])
    cov.to_csv(path, sep="\t", index=False)
    summ = cov.groupby("domain")["coverage"].agg(
        n="size", mean="mean", median="median",
        n_ge90=lambda c: int((c >= 0.9).sum()))
    print(summ.round(3).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
