"""Tests for ancestra.expression (S9) and ancestra.isoforms (S7)."""
import gzip

import numpy as np
import pandas as pd
import pytest

from ancestra.expression import (build_profiles, pair_similarity, read_bgee, shared_columns,
                                 shuffled_profiles, spearman)
from ancestra.isoforms import decoy_isoform_sets, gene_pair_scores, parse_ensembl_pep, set_pair_scores

BGEE_HDR = ('Gene ID\t"Gene name"\tAnatomical entity ID\t"Anatomical entity name"\tExpression\t'
            'Call quality\tFDR\tExpression score\tExpression rank\n')


def _bgee(tmp_path, name, rows):
    p = tmp_path / name
    with gzip.open(p, "wt") as fh:
        fh.write(BGEE_HDR)
        for g, e, call, q, s in rows:
            fh.write(f'{g}\t"x"\t{e}\t"n"\t{call}\t{q}\t0.01\t{s}\t1e3\n')
    return str(p)


def test_read_and_similarity(tmp_path):
    ents = [f"UBERON:{i:07d}" for i in range(6)]
    a = _bgee(tmp_path, "a.tsv.gz",
              [("GA1", e, "present", "gold quality", 10 + i) for i, e in enumerate(ents)]
              + [("GA2", e, "absent", "gold quality", 60 - i) for i, e in enumerate(ents)]
              + [("GA3", ents[0], "present", "bronze quality", 50)])          # low quality: dropped
    b = _bgee(tmp_path, "b.tsv.gz",
              [("GB1", e, "present", "silver quality", 5 + 2 * i) for i, e in enumerate(ents[:5])]
              + [("GB1", "CL:0000001", "present", "gold quality", 99)])      # entity only in B
    pa, pb = build_profiles(read_bgee(a)), build_profiles(read_bgee(b))
    assert list(pa.genes) == ["GA1", "GA2"]
    ca, cb, common = shared_columns(pa, pb)
    assert len(common) == 6 - 1 + 0 and "CL:0000001" not in common
    s = pair_similarity(pa, pb, pa.rows(["GA1", "GA2", "GAX"]), pb.rows(["GB1", "GB1", "GB1"]), min_shared=5)
    assert s.n_shared.tolist() == [5, 5, 0]
    assert s.spearman[0] == pytest.approx(1.0) and s.spearman[1] == pytest.approx(-1.0)
    assert np.isnan(s.spearman[2])
    assert s.jaccard[0] == pytest.approx(1.0) and s.jaccard[1] == pytest.approx(0.0)
    s6 = pair_similarity(pa, pb, pa.rows(["GA1"]), pb.rows(["GB1"]), min_shared=6)
    assert np.isnan(s6.spearman[0])


def test_spearman_ties_and_shuffle():
    assert spearman(np.array([1., 2, 2, 3]), np.array([1., 2, 2, 3])) == pytest.approx(1.0)
    assert np.isnan(spearman(np.array([1., 1, 1]), np.array([1., 2, 3])))
    df = pd.DataFrame({"gene": [f"G{i}" for i in range(50) for _ in range(3)],
                       "entity": ["E1", "E2", "E3"] * 50, "present": True,
                       "score": np.arange(150, dtype=np.float32)})
    p = build_profiles(df)
    q = shuffled_profiles(p, 11)
    assert list(q.genes) == list(p.genes)
    assert sorted(map(tuple, q.score.tolist())) == sorted(map(tuple, p.score.tolist()))
    assert not np.array_equal(q.score, p.score)


PEP = """>ENSP1.2 pep chromosome:X:1:1:10:1 gene:ENSG1.3 transcript:ENST1.1 gene_biotype:protein_coding
MKKLLA
>ENSP2.1 pep chromosome:X:1:1:10:1 gene:ENSG1.3 transcript:ENST2.1
MKKLLAWW
>FBpp9 pep chromosome:2L:1:10:1 gene:FBgn9 transcript:FBtr9
MAA
>Y1.1 pep chromosome:I:1:10:1 gene:WBGene1.1 transcript:Y1.1
MCC
"""


def test_parse_ensembl_pep(tmp_path):
    p = tmp_path / "x.pep.all.fa"
    p.write_text(PEP)
    d = parse_ensembl_pep(str(p))
    assert d.protein_id.tolist() == ["ENSP1", "ENSP2", "FBpp9", "Y1.1"]
    assert d.gene_id.tolist() == ["ENSG1", "ENSG1", "FBgn9", "WBGene1"]
    assert d.length.tolist() == [6, 8, 3, 3]
    assert d.transcript_id.tolist()[:2] == ["ENST1", "ENST2"]


def test_gene_pair_scores_and_decoys():
    pg = {"a0": "A", "a1": "A", "b0": "B", "b1": "B", "c0": "C"}
    q = ["a0", "a1", "b0", "a0", "a1"]
    t = ["b0", "b1", "a0", "c0", "a0"]          # last hit is within gene A: ignored
    bits = [100, 150, 110, 80, 999]
    s = gene_pair_scores(q, t, bits, pg, {"a0", "b0", "c0"}, [("A", "B"), ("B", "A"), ("A", "C"), ("B", "C")])
    assert s.canonical.tolist() == [110, 110, 80, 0]
    assert s.best_iso.tolist() == [150, 150, 80, 0]
    assert s.best_is_canonical.tolist() == [False, False, True, True]
    assert (s.best_prot_a[0], s.best_prot_b[0]) == ("a1", "b1")
    assert (s.best_prot_a[1], s.best_prot_b[1]) == ("b1", "a1")
    iso = {"A": ["a0", "a1"], "B": ["b0", "b1"], "C": ["c0"], "D": ["d0", "d1", "d2"]}
    canon = {"A": "a0", "B": "b0", "C": "c0", "D": "d0"}
    fam = {g: 1 for g in iso}
    sp = {"A": "s1", "B": "s2", "C": "s1", "D": "s1"}
    dec = decoy_isoform_sets(iso, canon, fam, sp, seed=1)
    assert dec["C"] == ["c0"]
    assert dec["A"][0] == "a0" and len(dec["A"]) == 2 and dec["A"][1] in ("d1", "d2")   # same species, other gene
    assert dec["B"][0] == "b0" and dec["B"][1] in ("a1", "d1", "d2")                   # fallback: any species
    assert set(dec["D"][1:]) <= {"a1"} and dec["D"][0] == "d0"
    h = pd.DataFrame({"q": q, "t": t, "bits": bits})
    assert set_pair_scores(h, {"A": ["a0"], "B": ["b0", "b1"]}, [("A", "B")]).tolist() == [110]
    assert set_pair_scores(h, iso, [("A", "B")]).tolist() == [150]
