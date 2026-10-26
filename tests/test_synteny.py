"""Tests for ancestra.genomes (gene-order tables) and ancestra.synteny (S5)."""
import gzip
import os

import numpy as np
import pytest

from ancestra.genomes import (assign_order, build_ensembl_species, neighbour_csr,
                              parse_gff3_proteins, shuffle_order)
from ancestra.synteny import (AncestralAdjacency, context_support, pair_context_scores,
                              rbh_anchors, select_anchors, SyntenyContext)

GFF = """##gff-version 3
##sequence-region   1 1 100000
##sequence-region   2 1 100000
1\tens\tchromosome\t1\t100000\t.\t.\t.\tID=chromosome:1
1\tens\tgene\t5000\t6000\t.\t-\t.\tID=gene:G2;Name=B
1\tens\tmRNA\t5000\t6000\t.\t-\t.\tID=transcript:T2;Parent=gene:G2
1\tens\tCDS\t5100\t5900\t.\t-\t0\tID=CDS:P2;Parent=transcript:T2;protein_id=P2
1\tens\tgene\t100\t900\t.\t+\t.\tID=gene:G1;Name=A
1\tens\tmRNA\t100\t900\t.\t+\t.\tID=transcript:T1;Parent=gene:G1
1\tens\texon\t100\t900\t.\t+\t.\tParent=transcript:T1
1\tens\tCDS\t150\t300\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1;protein_id=P1
1\tens\tCDS\t400\t800\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1;protein_id=P1
1\tens\tgene\t9000\t9500\t.\t+\t.\tID=gene:G3;Name=C
1\tens\tmRNA\t9000\t9500\t.\t+\t.\tID=transcript:T3;Parent=gene:G3
1\tens\tCDS\t9000\t9500\t.\t+\t0\tID=CDS:P3;Parent=transcript:T3;protein_id=P3
2\tens\tgene\t10\t50\t.\t+\t.\tID=gene:G4;Name=D
2\tens\tmRNA\t10\t50\t.\t+\t.\tID=transcript:T4;Parent=gene:G4
2\tens\tCDS\t10\t50\t.\t+\t0\tID=CDS:P4;Parent=transcript:T4;protein_id=P4
ALT1\tens\tgene\t10\t60\t.\t+\t.\tID=gene:G5;Name=B
ALT1\tens\tmRNA\t10\t60\t.\t+\t.\tID=transcript:T5;Parent=gene:G5
ALT1\tens\tCDS\t10\t60\t.\t+\t0\tID=CDS:P5;Parent=transcript:T5;protein_id=P5
ALT1\tens\tgene\t4000\t4100\t.\t+\t.\tID=gene:G6;Name=C
ALT1\tens\tmRNA\t4000\t4100\t.\t+\t.\tID=transcript:T6;Parent=gene:G6
ALT1\tens\tCDS\t4000\t4100\t.\t+\t0\tID=CDS:P6;Parent=transcript:T6;protein_id=P6
ALT1\tens\tgene\t2000\t2100\t.\t+\t.\tID=gene:G7;Name=Z
ALT1\tens\tmRNA\t2000\t2100\t.\t+\t.\tID=transcript:T7;Parent=gene:G7
ALT1\tens\tCDS\t2000\t2100\t.\t+\t0\tID=CDS:P7;Parent=transcript:T7;protein_id=P7
"""


@pytest.fixture
def tiny(tmp_path):
    full = tmp_path / "sp.full.gff3.gz"
    with gzip.open(full, "wt") as fh:
        fh.write(GFF)
    prim = tmp_path / "sp.gff3"
    prim.write_text("\n".join(l for l in GFF.splitlines() if not l.startswith("ALT1")) + "\n")
    fa = tmp_path / "sp.fa"
    fa.write_text("".join(f">{p}\nMA\n" for p in ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "PX"]))
    return str(fa), str(full), str(prim)


def test_parse_gff3_proteins(tiny):
    _, full, _ = tiny
    df = parse_gff3_proteins(full).set_index("protein_id")
    assert df.loc["P1", "gene_id"] == "G1"
    assert (df.loc["P1", "start"], df.loc["P1", "end"]) == (100, 900)  # gene span, not CDS
    assert df.loc["P2", "strand"] == -1
    assert len(df) == 7


def test_assign_order_and_build(tiny):
    fa, full, prim = tiny
    df, cov = build_ensembl_species(fa, full, prim)
    d = df.set_index("protein_id")
    assert cov["n_proteins"] == 8 and cov["n_mapped"] == 7 and cov["missing_examples"] == ["PX"]
    # primary ranks on chr 1 by start: P1 < P2 < P3
    assert [d.loc[p, "order"] for p in ["P1", "P2", "P3"]] == [0, 1, 2]
    assert d.loc["P4", "order"] == 0 and d.loc["P4", "contig_n"] == 1
    # alt contig ALT1: B and C match unique primary genes on chr 1 -> projected
    assert not d.loc["P5", "primary"] and d.loc["P5", "projected"]
    assert d.loc["P5", "contig"] == "1" and d.loc["P5", "raw_contig"] == "ALT1"
    assert d.loc["P5", "order"] == 1          # at G2's position (direct symbol match)
    assert d.loc["P6", "order"] == 2          # at G3's position
    # Z has no primary match: projected by median offset (start 2000 + offset ~ between)
    assert d.loc["P7", "contig"] == "1" and 0 <= d.loc["P7", "order"] <= 3


def test_neighbour_csr_window_and_ties():
    cc = np.array([0, 0, 0, 0, 1, 1, -1, 0])
    oo = np.array([0, 1, 2, 3, 0, 5, 0, 1])     # gene 7 ties with gene 1
    ptr, idx = neighbour_csr(cc, oo, 1)
    nb = {i: set(idx[ptr[i]:ptr[i + 1]].tolist()) for i in range(len(cc))}
    assert nb[0] == {1, 7}
    assert nb[1] == {0, 2, 7}
    assert nb[3] == {2}
    assert nb[4] == set() and nb[5] == set() and nb[6] == set()
    ptr2, idx2 = neighbour_csr(cc, oo, 5)
    assert set(idx2[ptr2[4]:ptr2[5]].tolist()) == {5}


def test_shuffle_keeps_slots():
    cc = np.array([0, 0, 0, 1, 1, -1])
    oo = np.array([0, 1, 2, 0, 1, 0])
    c2, o2 = shuffle_order(cc, oo, 11)
    assert sorted(zip(c2[:5], o2[:5])) == sorted(zip(cc[:5], oo[:5]))
    assert c2[5] == -1


def test_rbh_anchors_unambiguous():
    # species A genes 0,1,2 ; species B genes 10,11,12
    q = [0, 10, 1, 11, 1, 12, 2, 12]
    t = [10, 0, 11, 1, 12, 1, 12, 2]
    b = [100, 100, 80, 80, 78, 78, 50, 50]
    a1, b1 = rbh_anchors(q, t, b, amb=1.0)
    assert set(zip(a1.tolist(), b1.tolist())) == {(0, 10), (1, 11)}
    a2, b2 = rbh_anchors(q, t, b, amb=0.9)   # 1 has 11 (80) and 12 (78): ambiguous
    assert set(zip(a2.tolist(), b2.tolist())) == {(0, 10)}


def _two_genomes():
    """Species A genes 0..9 on one contig, species B genes 10..19 same order;
    a tandem/dispersed paralog: A gene 5 has orthologue 15 and paralogue 19
    whose neighbourhood is unrelated (B order shuffled at the end)."""
    cc = np.r_[np.zeros(10, int), np.ones(10, int)]
    oo = np.r_[np.arange(10), np.arange(10)]
    return cc, oo


def test_pair_context_and_family_exclusion():
    cc, oo = _two_genomes()
    oo = oo.copy()
    oo[19] = 50                          # paralog 19 far away on B
    anchors = (np.array([3, 4, 6, 7]), np.array([13, 14, 16, 17]))
    sc = pair_context_scores([5, 5], [15, 19], *anchors, cc, oo, w=2)
    assert sc.tolist() == [4, 0]
    fam = np.arange(20)
    fam[[4, 14, 5, 15, 19]] = 99        # anchor (4,14) in the candidate's family
    sc = pair_context_scores([5], [15], *anchors, cc, oo, w=2, fam=fam)
    assert sc.tolist() == [3]


def test_ancestral_adjacency_units_and_shuffle_control():
    cc, oo = _two_genomes()
    adj = AncestralAdjacency(cc, oo, w=1)
    # units: A genes (5,6) form unit 0 (in-paralogs), others singletons
    gu = np.full(20, -1)
    gu[[5, 6]] = 0
    gu[4], gu[7], gu[15], gu[14], gu[16] = 1, 2, 3, 4, 5
    nadj = adj.unit_adjacency(gu, 6)
    nb0 = set(nadj.adj_idx[nadj.adj_ptr[0]:nadj.adj_ptr[1]].tolist())
    assert nb0 == {1, 2}
    sup = context_support((np.array([0]), np.array([3])), (np.array([1, 2]), np.array([4, 5])), nadj)
    assert sup.tolist() == [2]
    # shuffled positions destroy support in expectation
    s = []
    for seed in range(30):
        a2 = AncestralAdjacency(cc, oo, w=1, shuffle_seed=seed)
        s.append(context_support((np.array([0]), np.array([3])),
                                 (np.array([1, 2]), np.array([4, 5])), a2.unit_adjacency(gu, 6))[0])
    assert np.mean(s) < 1.0


def test_select_anchors():
    # units 0,1 child 0 ; units 2,3 child 1; 0-2 clear 1:1, 1-3 ambiguous with 1-2
    ca, cb, cs = np.array([0, 1, 1]), np.array([2, 3, 2]), np.array([0.9, 0.5, 0.48])
    child = np.array([0, 0, 1, 1])
    pa, pb = select_anchors(ca, cb, cs, child, amb=0.9)
    assert list(zip(pa.tolist(), pb.tolist())) == [(0, 2)]
    pa, pb = select_anchors(ca, cb, cs, child, amb=1.0)
    assert set(zip(pa.tolist(), pb.tolist())) == {(0, 2), (1, 3)}


def _node_ctx(units_genes, child, ca, cb, cs, anchors):
    from ancestra.merge import NodeContext
    g_sorted = np.concatenate(units_genes)
    b = np.r_[0, np.cumsum([len(u) for u in units_genes])]
    return NodeContext(node="X", children=["c0", "c1"], units=np.arange(len(units_genes)) + 100,
                       child=np.asarray(child), cand_a=np.asarray(ca), cand_b=np.asarray(cb),
                       cand_score=np.asarray(cs, float), anchors=np.asarray(anchors, bool),
                       adjacency=None, members=lambda i, _g=g_sorted, _b=b: _g[_b[i]:_b[i + 1]])


def test_synteny_context_node_context():
    """Units: child 0 = genes 3,4,5 (A); child 1 = genes 13,14,15 (B); collinear.
    Candidates: (u0,u3) (u1,u4) (u2,u5) clear, plus (u1,u5) paralog-like competitor."""
    cc, oo = _two_genomes()
    ug = [np.array([3]), np.array([4]), np.array([5]), np.array([13]), np.array([14]), np.array([15])]
    child = [0, 0, 0, 1, 1, 1]
    ca, cb, cs = [0, 1, 2, 1], [3, 4, 5, 5], [0.9, 0.5, 0.9, 0.49]
    ctx = _node_ctx(ug, child, ca, cb, cs, [True, False, True, False])
    fam = np.arange(20)
    fam[[4, 14, 15, 5]] = 77          # u1, u4, u5, u2 one family
    sc = SyntenyContext(cc, oo, family=fam, w=1, mode="tiebreak", anchor="core")
    sup = sc.support(ctx)
    # (u1,u4): anchors (u0,u3) adjacent on both sides -> 1; (u2,u5) anchor (u2,u5) itself excluded
    assert sup[1] == 1
    assert sup[3] == 0                # (u1,u5): u5's neighbours (14) are not partners of u1's
    b = sc(ctx)
    assert b[1] > 0 and b[3] == 0
    # same through the core's adjacency (ctx.adjacency) instead of own positions
    from ancestra.genome import unit_adjacency
    ug_arr = np.full(20, -1)
    for i, g in enumerate(ug):
        ug_arr[g] = i
    sp_of = (np.arange(20) >= 10).astype(int)
    na = np.array([3, 4, 13, 14]); nb = np.array([4, 5, 14, 15]); nd = np.ones(4, int)
    ctx.adjacency = unit_adjacency(na, nb, nd, ug_arr, sp_of, 6)
    sc2 = SyntenyContext(None, None, family=fam, mode="tiebreak", anchor="core")
    assert sc2.support(ctx).tolist() == sup.tolist()


def test_to_gene_order(tmp_path):
    import pandas as pd
    from ancestra.genomes import to_gene_order, write_gene_order
    from ancestra.io import GeneTable
    gt = GeneTable(names=["a1", "a2", "b1"], species=["A", "B"], species_of=np.array([0, 0, 1]),
                   lengths=np.array([10, 10, 10]))
    df = assign_order(pd.DataFrame(dict(protein_id=["a2", "a1"], gene_id=["g2", "g1"], contig=["1", "1"],
                                        start=[500, 100], end=[600, 200], strand=[1, -1],
                                        primary=[True, True], source="t")))
    write_gene_order(df, str(tmp_path / "A.parquet"))
    go = to_gene_order(gt, str(tmp_path))
    assert go.rank.tolist() == [0, 1, -1] and go.contig.tolist() == [0, 0, -1]
    assert go.coverage == {"A": 1.0, "B": 0.0}


def test_override_bonus():
    from ancestra.synteny import override_bonus
    # unit 0 (child 0) has candidates to 1 (seq 1.0, sup 0) and 2 (seq 0.95, sup 3): override -> 2
    ca, cb, cs = np.array([0, 0]), np.array([1, 2]), np.array([1.0, 0.95])
    child = np.array([0, 1, 1])
    b = override_bonus(ca, cb, cs, child, np.array([0, 3]), near=0.9, margin=1)
    assert b[0] == 0 and cs[1] + b[1] > cs[0]
    # too far in sequence (0.8 < 0.9 x 1.0): no override
    b = override_bonus(ca, cb, np.array([1.0, 0.8]), child, np.array([0, 3]), near=0.9, margin=1)
    assert (b == 0).all()
    # equal support: no override
    b = override_bonus(ca, cb, cs, child, np.array([2, 2]))
    assert (b == 0).all()
