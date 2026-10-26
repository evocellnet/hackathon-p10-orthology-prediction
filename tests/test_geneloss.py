"""Tests for ancestra.geneloss (S8 gene-loss validation)."""
import os
import random
import shutil

import pandas as pd
import pytest

from ancestra.geneloss import (CdsIndex, FastaIndex, GeneOrderIndex, HogIndex, LocusSearcher, MiniprotHit,
                               classify_hits, parse_miniprot_paf, query_votes, select_apparent_losses,
                               select_complete_single_copy, syntenic_locus)

CLADE = ["A", "B", "C"]


def toy_hogs():
    # genes are named <species>_<locus index>; HOG i = genes at index i in each species
    # HOG 2 missing in C (apparent loss); HOG 5 has 2 copies in A.
    hogs = []
    for i in range(6):
        h = {f"A_{i}", f"B_{i}", f"C_{i}"}
        if i == 2:
            h.discard("C_2")
        hogs.append(h)
    hogs[5].add("A_5b")
    sp = {g: g.split("_")[0] for h in hogs for g in h}
    return HogIndex(hogs, sp, outgroup_hogs=[set(h) for h in hogs])


def toy_order():
    rows = []
    for s in CLADE:
        for i in range(6):
            if s == "C" and i == 2:
                continue
            rows.append({"protein_id": f"{s}_{i}", "gene_id": f"g{s}{i}", "contig": "1",
                         "start": 10_000 * (i + 1), "end": 10_000 * (i + 1) + 2_000, "strand": 1,
                         "primary": True, "projected": False})
    rows.append({"protein_id": "A_5b", "gene_id": "gA5b", "contig": "2", "start": 1, "end": 100, "strand": 1,
                 "primary": True, "projected": False})
    df = pd.DataFrame(rows)
    return GeneOrderIndex({s: df[df.protein_id.str.startswith(s + "_")].copy() for s in CLADE})


def test_select_apparent_losses():
    hi = toy_hogs()
    c = select_apparent_losses(hi, CLADE, max_copies=2)
    assert list(c.hog) == [2] and list(c.target) == ["C"]
    # root-species rule: C is ambiguous (gain vs loss) without an outgroup member
    assert len(select_apparent_losses(hi, CLADE, root_species={"C"}, outgroup=["X"])) == 0
    assert set(select_complete_single_copy(hi, CLADE)) == {0, 1, 3, 4}


def test_syntenic_locus_bracketed_and_voting():
    hi, go = toy_hogs(), toy_order()
    loc = syntenic_locus(hi, go, 2, "C", ["B", "A"], k=2, pad=1_000)
    assert loc.kind == "bracketed" and loc.contig == "1"
    assert loc.left_anchor == "C_1" and loc.right_anchor == "C_3"
    assert loc.start == 20_000 - 1_000 and loc.end == 42_000 + 1_000
    assert loc.n_sources == 2 and loc.n_agree == 2


def test_syntenic_locus_hides_own_target_gene():
    hi, go = toy_hogs(), toy_order()
    # positive control: HOG 3 in C -- its own C member is never an anchor
    # k=1: left neighbour B_2 has no C ortholog -> one-sided
    loc1 = syntenic_locus(hi, go, 3, "C", ["B"], k=1, pad=0)
    assert loc1.kind == "one_sided" and loc1.right_anchor == "C_4"
    loc = syntenic_locus(hi, go, 3, "C", ["B"], k=2, pad=0)
    assert loc.kind == "bracketed" and loc.left_anchor == "C_1" and loc.k_left == 2 and loc.right_anchor == "C_4"
    assert loc.start <= 40_000 <= loc.end    # the hidden gene C_3 lies inside


PAF_GFF = """##gff-version 3
##PAF\tq1\t100\t0\t100\t+\tlocus\t5000\t100\t400\t290\t300\t0\tAS:i:400\tms:i:450\tnp:i:98\tfs:i:0\tst:i:0\tda:i:0\tdo:i:0\tcg:Z:50M200N50M
locus\tminiprot\tmRNA\t101\t400\t450\t+\t.\tID=MP000001
locus\tminiprot\tCDS\t101\t250\t200\t+\t0\tParent=MP000001
locus\tminiprot\tCDS\t251\t400\t200\t+\t0\tParent=MP000001
##PAF\tq2\t100\t10\t60\t-\tlocus\t5000\t1000\t1200\t140\t150\t0\tAS:i:100\tms:i:120\tnp:i:40\tfs:i:2\tst:i:1\tda:i:-1\tdo:i:0\tcg:Z:50M
"""


def test_parse_miniprot():
    h = parse_miniprot_paf(PAF_GFF)
    assert len(h) == 2
    assert h[0].intact and h[0].n_introns == 1 and h[0].blocks == [(101, 250), (251, 400)]
    assert abs(h[0].identity - 290 / 300) < 1e-9
    assert h[1].disrupted and not h[1].intact and h[1].qcov == 0.5 and h[1].da == 0
    assert h[1].blocks == [(1001, 1200)]


def _hit(q, qcov=1.0, fs=0, st=0, score=100):
    return MiniprotHit(q, 100, 0, int(100 * qcov), "+", "locus", 0, 300, 290, 300, score, 90, fs, st)


def test_classify_rules():
    # intact, unannotated, majority intact -> annotation gap
    c = classify_hits([_hit("a"), _hit("b"), _hit("c", fs=1)], True, 0.0)
    assert c["cls"] == "A_gap"
    # one intact (e.g. a pseudogene modelled as coding in one species) vs two disrupted -> pseudogene
    c = classify_hits([_hit("a"), _hit("b", st=1), _hit("c", fs=1)], True, 0.0)
    assert c["cls"] == "B_pseudo" and c["conflict"]
    # the same with an annotated overlap -> annotated (fission)
    c = classify_hits([_hit("a"), _hit("b", st=1)], True, 0.0, annotated_overlap=lambda h: "G" if h.query == "a" else None)
    assert c["cls"] == "A_annotated" and c["overlap_gene"] == "G"
    assert classify_hits([_hit("a", qcov=0.5)], True, 0.0)["cls"] == "B_partial"
    assert classify_hits([], True, 0.0)["cls"] == "C_deleted"
    assert classify_hits([], False, 0.0)["cls"] == "D_undetermined"
    assert classify_hits([], True, 0.5)["cls"] == "D_undetermined"
    v = query_votes([_hit("a"), _hit("a", fs=1, score=50), _hit("b", st=1)], ["a", "b", "c"])
    assert v == {"nq_intact": 1, "nq_disrupted": 1, "nq_partial": 0, "nq_none": 1}


def _write_fa(path, recs, width=60):
    with open(path, "w") as fh:
        for k, s in recs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")
    off = 0
    with open(path + ".fai", "w") as fai:
        for k, s in recs.items():
            off += len(f">{k}\n")
            fai.write(f"{k}\t{len(s)}\t{off}\t{width}\t{width + 1}\n")
            off += len(s) + (len(s) + width - 1) // width


def test_fasta_index(tmp_path):
    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(1000))
    p = str(tmp_path / "g.fa")
    _write_fa(p, {"c1": seq, "c2": seq[::-1]})
    fx = FastaIndex(p)
    assert fx.fetch("c1", 1, 10) == seq[:10]
    assert fx.fetch("c1", 55, 130) == seq[54:130]
    assert fx.fetch("c2", 991, 2000) == seq[::-1][990:]


CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
         "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
         "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}


@pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot not installed")
def test_end_to_end_locus_search(tmp_path):
    """Synthetic genome: intact gene -> A, frameshifted -> B_pseudo, deleted -> C."""
    rng = random.Random(7)
    prot = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(299))
    cds = "".join(CODON[a] for a in prot) + "TAA"
    bg = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    flank1, flank2 = bg(20_000), bg(20_000)
    broken = cds[:300] + "A" + cds[300:600] + "TAG" + cds[603:]
    genome = {"1": flank1 + cds + flank2, "2": flank1 + broken + flank2, "3": flank1 + flank2}
    gpath = str(tmp_path / "g.fa")
    _write_fa(gpath, genome)
    hi = HogIndex([{"S_q"}], {"S_q": "S"})
    go = GeneOrderIndex({"T": pd.DataFrame(columns=["protein_id", "contig", "start", "end", "primary", "projected"])})
    S = LocusSearcher(hi, go, "T", FastaIndex(gpath), {"S_q": prot}, str(tmp_path / "w"),
                      cds=CdsIndex(pd.DataFrame(columns=["protein_id", "contig", "start", "end"])))
    from ancestra.geneloss import Locus
    res = {}
    for c in "123":
        loc = Locus("T", c, 1, len(genome[c]), "bracketed", "S", "S_q")
        res[c] = S.search(c, 0, loc)
    assert res["1"]["cls"] == "A_gap" and res["1"]["fs"] == 0 and res["1"]["st"] == 0
    assert res["2"]["cls"] == "B_pseudo" and (res["2"]["fs"] + res["2"]["st"]) >= 1
    assert res["3"]["cls"] == "C_deleted"
    # the same intact gene, annotated in the target in a target-only HOG -> A_annotated
    s0 = len(flank1) + 1
    cds_df = pd.DataFrame([{"protein_id": "T_x", "contig": "1", "start": s0, "end": s0 + len(cds) - 1}])
    hi2 = HogIndex([{"S_q"}, {"T_x"}], {"S_q": "S", "T_x": "T"})
    S2 = LocusSearcher(hi2, go, "T", FastaIndex(gpath), {"S_q": prot}, str(tmp_path / "w2"), cds=CdsIndex(cds_df))
    r = S2.search("1b", 0, Locus("T", "1", 1, len(genome["1"]), "bracketed", "S", "S_q"))
    assert r["cls"] == "A_annotated" and r["overlap_gene"] == "T_x"
