import gzip
import os

import pandas as pd
import pytest

from ancestra import genomes, genomes_qfo as gq

GCOORD = """<?xml version="1.0" encoding="UTF-8"?>
<gnEntries>
<gnEntry xmlns:ns2="http://www.ebi.ac.uk/proteins/api/doc/xsd/feature">
    <accession>P00001</accession>
    <name>A_HUMAN</name>
    <sequence>MKV</sequence>
    <gnCoordinate ensembl_gene_id="G1" ensembl_transcript_id="T1" ensembl_translation_id="P1">
        <genomicLocation chromosome="1" start="500" end="100" reverse_strand="true">
            <exon id="E1"><genomeLocation><ns2:begin position="500"/><ns2:end position="100"/></genomeLocation></exon>
        </genomicLocation>
        <feature type="domain"><ns2:description>x</ns2:description></feature>
    </gnCoordinate>
    <gnCoordinate ensembl_gene_id="G1" ensembl_transcript_id="T2" ensembl_translation_id="P2">
        <genomicLocation chromosome="1" start="600" end="200" reverse_strand="true"/>
    </gnCoordinate>
    <gnCoordinate ensembl_gene_id="G1P" ensembl_transcript_id="T3" ensembl_translation_id="P3">
        <genomicLocation chromosome="CHR_HG1_PATCH" start="10" end="20" reverse_strand="false"/>
    </gnCoordinate>
</gnEntry>
<gnEntry xmlns:ns2="http://www.ebi.ac.uk/proteins/api/doc/xsd/feature">
    <accession>P00002</accession>
    <sequence>MAA</sequence>
    <gnCoordinate ensembl_gene_id="G2" ensembl_transcript_id="T4" ensembl_translation_id="P4">
        <genomicLocation chromosome="CHR_HG1_PATCH" start="30" end="90" reverse_strand="false"/>
    </gnCoordinate>
</gnEntry>
<gnEntry xmlns:ns2="http://www.ebi.ac.uk/proteins/api/doc/xsd/feature">
    <accession>P00003</accession>
    <sequence>MCC</sequence>
    <gnCoordinate ensembl_gene_id="G3" ensembl_transcript_id="T5" ensembl_translation_id="P5">
        <genomicLocation chromosome="" start="30" end="90" reverse_strand="false"/>
    </gnCoordinate>
</gnEntry>
<gnEntry>
    <accession>P00005</accession>
    <sequence>MFF</sequence>
    <gnCoordinate refseq_nucleotide_id="NC_1.1" refseq_protein_id="NP_1">
        <genomicLocation chromosome="20" start="900" end="700" reverse_strand="true" nucleotideId="NC_1.1"/>
    </gnCoordinate>
    <gnCoordinate refseq_nucleotide_id="NC_1.1" refseq_protein_id="XP_2">
        <genomicLocation chromosome="20" start="950" end="800" reverse_strand="true" nucleotideId="NC_1.1"/>
    </gnCoordinate>
</gnEntry>
<gnEntry xmlns:ns2="http://www.ebi.ac.uk/proteins/api/doc/xsd/feature">
    <accession>Q99999</accession>
    <sequence>MDD</sequence>
    <gnCoordinate ensembl_gene_id="G9" ensembl_transcript_id="T9" ensembl_translation_id="P9">
        <genomicLocation chromosome="2" start="1" end="9" reverse_strand="false"/>
    </gnCoordinate>
</gnEntry>
</gnEntries>
"""

# ENA-style (CDS -> gene, UniProt Dbxref) + Ensembl-style (CDS -> transcript -> gene)
GFF = """##gff-version 3
U00096.3\tGenbank\tregion\t1\t5000\t.\t+\t.\tID=U00096.3:1..5000
U00096.3\tGenbank\tgene\t190\t255\t.\t+\t.\tID=gene-b0001;locus_tag=b0001
U00096.3\tGenbank\tCDS\t190\t255\t.\t+\t0\tID=cds-AAC73112.1;Parent=gene-b0001;Dbxref=UniProtKB/Swiss-Prot:P0AD86,NCBI_GP:AAC73112.1;protein_id=AAC73112.1
U00096.3\tGenbank\tgene\t337\t2799\t.\t-\t.\tID=gene-b0002;locus_tag=b0002
U00096.3\tGenbank\tCDS\t337\t2799\t.\t-\t0\tID=cds-AAC73113.1;Parent=gene-b0002;protein_id=AAC73113.1
3\tensembl\tgene\t1000\t3000\t.\t+\t.\tID=gene:ENSG7;biotype=protein_coding
3\tensembl\tmRNA\t1000\t3000\t.\t+\t.\tID=transcript:ENST7;Parent=gene:ENSG7
3\tensembl\texon\t1000\t1100\t.\t+\t.\tParent=transcript:ENST7
3\tensembl\tCDS\t1050\t1100\t.\t+\t0\tID=CDS:ENSP7;Parent=transcript:ENST7;protein_id=ENSP7
3\tensembl\tCDS\t2000\t2100\t.\t+\t0\tID=CDS:ENSP7;Parent=transcript:ENST7;protein_id=ENSP7
"""

IDMAP = """P0AD86\tEMBL-CDS\tAAC73112.1
P00561\tEMBL-CDS\tAAC73113.1
P00561\tGene_OrderedLocusName\tb0002
Q00007\tEnsembl_PRO\tENSP7.3
Q00008\tGene_OrderedLocusName\tb0001
"""


def _w(path, text, gz=False):
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return str(path)


def test_is_primary_contig():
    for c in ["1", "X", "JH584299.1", "KZ626836.1", "U00096.3", "NC_133176.1", "2R", "LG3"]:
        assert gq.is_primary_contig(c), c
    for c in ["CHR_HG1_PATCH", "CHR_HSCHR6_MHC_COX_CTG1", "HSCHR1_1_CTG3", "chr6_GL000250v2_alt",
              "CHR_HG2_NOVEL_TEST"]:
        assert not gq.is_primary_contig(c), c


def test_read_fasta_and_md5(tmp_path):
    p = _w(tmp_path / "a.fasta", ">sp|P00001|A_HUMAN x\nMK\nv\n>tr|Q2|B_HUMAN\nMAA\n>plain desc\nMC\n")
    s = gq.read_fasta(p)
    assert list(s) == ["P00001", "Q2", "plain"] and s["P00001"] == "MKV"
    assert gq.seq_md5("mkv\n") == gq.seq_md5("MKV")


def test_parse_gcoord(tmp_path):
    p = _w(tmp_path / "x_gcoord.xml", GCOORD)
    df = gq.parse_gcoord(p, {"P00001", "P00002", "P00003", "P00005"}).set_index("protein_id")
    assert set(df.index) == {"P00001", "P00002", "P00005"}
    r5 = df.loc["P00005"]   # RefSeq style: no gene id, nucleotideId as contig, one gene
    assert (r5.gene_id, r5.contig, r5.start, r5.end, r5.strand, r5.n_genes) == \
        ("NP_1", "NC_1.1", 700, 950, -1, 1)      # empty chromosome dropped, Q99999 not wanted
    r = df.loc["P00001"]
    # primary gene preferred over the patch copy; span over both transcripts; reverse strand
    assert (r.gene_id, r.contig, r.start, r.end, r.strand, r.n_genes) == ("G1", "1", 100, 600, -1, 2)
    assert r.seq_md5 == gq.seq_md5("MKV")
    r2 = df.loc["P00002"]
    assert (r2.contig, r2.start, r2.end, r2.strand) == ("CHR_HG1_PATCH", 30, 90, 1)


def test_gff_index_and_map(tmp_path):
    g = _w(tmp_path / "x.gff.gz", GFF, gz=True)
    key2genes, genes = gq.parse_gff_cds_genes(g)
    assert genes["b0001"] == ("U00096.3", 190, 255, 1)
    assert genes["b0002"] == ("U00096.3", 337, 2799, -1)
    assert genes["ENSG7"] == ("3", 1000, 3000, 1)
    assert key2genes["P0AD86"] == {"b0001": "uniprot"}
    assert "AAC73113" in key2genes and "ENSP7" in key2genes
    im = gq.read_idmapping(_w(tmp_path / "x.idmapping", IDMAP))
    assert ("Ensembl_PRO", "ENSP7") in im["Q00007"]
    df, st = gq.map_gff(["P0AD86", "P00561", "Q00007", "Q00008", "NOPE"], key2genes, genes, im)
    d = df.set_index("protein_id")["gene_id"].to_dict()
    assert d == {"P0AD86": "b0001", "P00561": "b0002", "Q00007": "ENSG7", "Q00008": "b0001"}
    assert st["uniprot"] == 1 and st["xref"] == 2 and st["gene_xref"] == 1 and st["unmapped"] == 1


def test_build_species_md5_rule(tmp_path):
    r22, r26 = tmp_path / "q22", tmp_path / "q26"
    for d in (r22 / "Bacteria", r26 / "Bacteria", r26 / "gff"):
        d.mkdir(parents=True)
    # P0AD86 identical; P00561 changed in 2026 -> excluded; Q00008 absent in 2026
    _w(r22 / "Bacteria" / "UP000000625_83333.fasta",
       ">sp|P0AD86|THRL\nMKR\n>sp|P00561|AK1H\nMRV\n>tr|Q00008|X\nMQQ\n")
    _w(r26 / "Bacteria" / "UP000000625_83333.fasta", ">sp|P0AD86|THRL\nMKR\n>sp|P00561|AK1H\nMRVL\n")
    _w(r26 / "Bacteria" / "UP000000625_83333.idmapping", IDMAP)
    _w(r26 / "gff" / "UP000000625.gff.gz", GFF, gz=True)
    sps = gq.discover_species(str(r22), str(r26))
    assert len(sps) == 1 and sps[0]["gff2026"].endswith("UP000000625.gff.gz")
    df, cov = gq.build_species(sps[0])
    assert list(df.columns[:len(genomes.COLUMNS)]) == genomes.COLUMNS
    assert df["protein_id"].tolist() == ["P0AD86"]
    assert df["source"].tolist() == ["qfo2026_gff_md5"]
    assert df["strand"].dtype == "int8" and df["primary"].dtype == bool
    assert (cov["n_proteins"], cov["n_gff2026_md5"], cov["n_md5_mismatch"], cov["n_not_in_2026"],
            cov["n_no_coord"]) == (3, 1, 1, 1, 2)
    out = tmp_path / "o" / "x.parquet"
    genomes.write_gene_order(df, str(out))
    assert genomes.read_gene_order(str(out)).shape == df.shape


def test_build_species_gcoord_first(tmp_path):
    r22, r26 = tmp_path / "q22", tmp_path / "q26"
    (r22 / "Eukaryota").mkdir(parents=True)
    (r26 / "gff").mkdir(parents=True)
    _w(r22 / "Eukaryota" / "UP000005640_9606.fasta",
       ">sp|P00001|A\nMKV\n>sp|P00002|B\nMAA\n>sp|P00003|C\nMCC\n")
    _w(r22 / "Eukaryota" / "UP000005640_9606_gcoord.xml", GCOORD)
    df, cov = gq.build_species(gq.discover_species(str(r22), str(r26))[0])
    d = df.set_index("protein_id")
    assert set(d.index) == {"P00001", "P00002"}
    assert (d["source"] == "qfo2022_gcoord").all()
    assert bool(d.loc["P00001", "primary"]) and not bool(d.loc["P00002", "primary"])
    assert cov["n_gcoord2022"] == 2 and cov["n_non_primary"] == 1 and cov["n_no_coord"] == 1
    assert "no_2026_gff" in cov["notes"]


GBFF = """LOCUS       CP000866             2000 bp    DNA     circular BCT 01-JAN-2020
DEFINITION  test.
VERSION     CP000866.1
FEATURES             Location/Qualifiers
     source          1..2000
                     /organism="x"
     gene            complement(100..400)
                     /locus_tag="Nmar_0001"
     CDS             complement(join(100..200,
                     300..390))
                     /locus_tag="Nmar_0001"
                     /protein_id="ABX00001.1"
                     /db_xref="UniProtKB/TrEMBL:A9A001"
                     /translation="MKVLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
                     LLLL"
     CDS             <500..>900
                     /locus_tag="Nmar_0002"
                     /protein_id="ABX00002.1"
ORIGIN
        1 acgt
//
"""


def test_parse_gbff(tmp_path):
    p = _w(tmp_path / "UP1.gbff.gz", GBFF, gz=True)
    key2genes, genes = gq.index_annotation(p)
    assert genes["Nmar_0001"] == ("CP000866.1", 100, 400, -1)   # gene feature span
    assert genes["Nmar_0002"] == ("CP000866.1", 500, 900, 1)    # CDS span (no gene feature)
    assert key2genes["A9A001"] == {"Nmar_0001": "uniprot"}
    assert key2genes["ABX00002"] == {"Nmar_0002": "pid"}


def test_layers_md5_and_assembly_check(tmp_path):
    """2026 gcoord fills a gap only when its sequence MD5 matches and it agrees
    with the 2022 placements; a disagreeing 2026 layer is skipped."""
    r22, r26 = tmp_path / "q22", tmp_path / "q26"
    (r22 / "Eukaryota").mkdir(parents=True)
    (r26 / "Eukaryota").mkdir(parents=True)
    _w(r22 / "Eukaryota" / "UP000005640_9606.fasta",
       ">sp|P00001|A\nMKV\n>sp|P00002|B\nMAA\n>sp|P00003|C\nMCC\n>sp|P00004|D\nMEE\n")
    _w(r22 / "Eukaryota" / "UP000005640_9606_gcoord.xml", GCOORD)
    _w(r26 / "Eukaryota" / "UP000005640_9606.fasta.gz", ">sp|P00001|A\nMKV\n", gz=True)

    def entry(acc, seq, chrom, s, e):
        return (f'<gnEntry><accession>{acc}</accession><sequence>{seq}</sequence>'
                f'<gnCoordinate ensembl_gene_id="N{acc}"><genomicLocation chromosome="{chrom}" '
                f'start="{s}" end="{e}" reverse_strand="false"/></gnCoordinate></gnEntry>\n')

    def g26(chrom1):
        body = entry("P00001", "MKV", chrom1, 150, 550)          # shared with 2022 gcoord
        body += "".join(entry(f"X{i}", "MXX", "5", i, i + 1) for i in range(30))
        body += entry("P00004", "MEE", "7", 10, 20)               # md5 ok -> fill
        body += entry("P00003", "MCCX", "7", 30, 40)              # md5 differs -> excluded
        return f"<gnEntries>\n{body}</gnEntries>\n"

    gq.MIN_SHARED, old = 1, gq.MIN_SHARED
    try:
        _w(r26 / "Eukaryota" / "UP000005640_9606_gcoord.xml.gz", g26("1"), gz=True)
        sp = gq.discover_species(str(r22), str(r26))[0]
        assert sp["gcoord2026"].endswith("_gcoord.xml.gz")
        df, cov = gq.build_species(sp)
        d = df.set_index("protein_id")
        assert d.loc["P00004", "source"] == "qfo2026_gcoord_md5"
        assert d.loc["P00001", "source"] == "qfo2022_gcoord"
        assert "P00003" not in d.index
        assert (cov["n_gcoord2026_md5"], cov["n_md5_mismatch"], cov["n_not_in_2026"]) == (1, 1, 0)
        # different assembly (P00001 elsewhere) -> 2026 layer skipped entirely
        _w(r26 / "Eukaryota" / "UP000005640_9606_gcoord.xml.gz", g26("9"), gz=True)
        df, cov = gq.build_species(sp)
        assert set(df["source"]) == {"qfo2022_gcoord"} and "skipped" in cov["notes"]
    finally:
        gq.MIN_SHARED = old


def test_versioned_xref_layer(tmp_path):
    """Accession absent in 2026 but its versioned EMBL-CDS id is annotated in the
    2026 GFF -> placed (qfo2026_gff_vxref); a different version is not."""
    r22, r26 = tmp_path / "q22", tmp_path / "q26"
    for d in (r22 / "Bacteria", r26 / "Bacteria", r26 / "gff"):
        d.mkdir(parents=True)
    _w(r22 / "Bacteria" / "UP000000625_83333.fasta", ">tr|OLD1|X\nMKR\n>tr|OLD2|Y\nMRV\n")
    _w(r22 / "Bacteria" / "UP000000625_83333.idmapping",
       "OLD1\tEMBL-CDS\tAAC73112.1\nOLD2\tEMBL-CDS\tAAC73113.2\n")
    _w(r26 / "Bacteria" / "UP000000625_83333.fasta", ">sp|P0AD86|THRL\nMKR\n")
    _w(r26 / "gff" / "UP000000625.gff.gz", GFF, gz=True)
    df, cov = gq.build_species(gq.discover_species(str(r22), str(r26))[0])
    assert df["protein_id"].tolist() == ["OLD1"] and df["source"].tolist() == ["qfo2026_gff_vxref"]
    assert (cov["n_gff2026_vxref"], cov["n_not_in_2026"], cov["n_no_coord"]) == (1, 1, 1)
    k2g, _ = gq.parse_gff_cds_genes(str(r26 / "gff" / "UP000000625.gff.gz"))
    assert k2g["ENSP7"] == {"ENSG7": "pid"}            # Ensembl CDS without version attr
