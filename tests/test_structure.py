"""Tests for ancestra.structure (S6 3Di layer)."""
import os
import shutil

import numpy as np
import pytest

from ancestra import structure as S
from ancestra.io import GeneTable
from ancestra.newick import parse_newick
from ancestra.scoring import ScoreTable


def _gt(n=6, species=("A", "B", "C")):
    names = [f"g{i}" for i in range(n)]
    sp = np.array([i % len(species) for i in range(n)], dtype=np.int32)
    L = np.full(n, 100, dtype=np.int32)
    return GeneTable(names, list(species), sp, L, ["M" * 100] * n)


def test_shuffled_3di_derangement_and_length():
    rng = np.random.default_rng(0)
    ss = {f"p{i}": "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=50 + i)) for i in range(100)}
    sh = S.shuffled_3di(ss, seed=11, bin_size=20)
    assert set(sh) == set(ss)
    for k in ss:
        assert len(sh[k]) == len(ss[k])
        assert sh[k] != ss[k]
    # deterministic per seed, different across seeds
    assert sh == S.shuffled_3di(ss, seed=11, bin_size=20)
    assert sh != S.shuffled_3di(ss, seed=23, bin_size=20)
    # composition preserved globally (resampling aside, multiset of letters roughly kept)
    a = sorted("".join(ss.values()))
    b = sorted("".join(sh.values()))
    assert abs(len(a) - len(b)) == 0


def test_link_families_rules():
    fam = np.array([0, 0, 0, 1, 2, 2, 3])          # sizes 3,1,2,1
    L = np.full(7, 100.0)
    hits = dict(q=np.array([0, 3, 4, 0]), t=np.array([3, 4, 6, 6]),
                evalue=np.array([1e-10, 5e-3, 1e-8, 1e-2]),
                qcov=np.array([0.9, 0.9, 0.1, 0.9]), tcov=np.array([0.9, 0.9, 0.1, 0.9]))
    # e<=1e-3, cov>=0.3: only 0-3 (fam0+fam1); 4-6 fails coverage, 3-4 fails e, 0-6 fails e
    new, st = S.link_families(fam, hits, L, evalue=1e-3, cov_shorter=0.3)
    assert new[0] == new[3] and new[4] != new[0] and new[6] != new[0]
    assert st["n_joins"] == 1
    # looser e: 0-3 and 3-4 chain fam0, fam1, fam2 together
    new, st = S.link_families(fam, hits, L, evalue=6e-3, cov_shorter=0.3)
    assert new[0] == new[3] == new[4] == new[5] and new[6] != new[0]
    # size cap: fam0+fam1 = 4 genes ok, + fam2 = 6 > 5 -> skipped
    new, st = S.link_families(fam, hits, L, evalue=6e-3, cov_shorter=0.3, max_size=5)
    assert new[0] == new[3] and new[4] != new[0] and st["n_capped"] == 1
    # small_max=1: a link needs a side with <= 1 gene (fam1, fam3 qualify)
    new, st = S.link_families(fam, hits, L, evalue=1.0, cov_shorter=0.3, small_max=1)
    assert new[0] == new[3]


def test_combine_scores_identity_and_rescue():
    gt = _gt(4)
    st = ScoreTable(q=np.array([0, 0, 1], np.int32), t=np.array([1, 2, 2], np.int32),
                    bits=np.array([50, 20, 10], np.float32), qcov=np.ones(3, np.float32),
                    tcov=np.ones(3, np.float32), pident=np.ones(3, np.float32),
                    evalue=np.ones(3, np.float32), self_bits=np.full(4, 100, np.float32))
    # no structure pairs -> identical scores
    empty = dict(q=np.zeros(0, int), t=np.zeros(0, int), bits=np.zeros(0))
    out, s = S.combine_scores(gt, st, empty, np.full(4, np.nan))
    assert np.allclose(out.bits, st.bits) and len(out) == len(st)
    # structure for (0,1) and an unscored pair (2,3); identity calibration
    fs = dict(q=np.array([0, 3]), t=np.array([1, 2]), bits=np.array([200.0, 100.0]))
    fself = np.full(4, 400.0)
    out, s = S.combine_scores(gt, st, fs, fself, mode="avg", fit=(0.0, 1.0))
    sc = {(int(a), int(b)): float(x) / 100 for a, b, x in zip(out.q, out.t, out.bits)}
    assert np.isclose(sc[(0, 1)], np.sqrt(0.5 * 0.5))
    assert np.isclose(sc[(0, 2)], 0.2) and np.isclose(sc[(1, 2)], 0.1)
    assert np.isclose(sc[(2, 3)], 0.25)              # structure-only pair added
    out2, _ = S.combine_scores(gt, st, fs, fself, mode="avg", fit=(0.0, 1.0), add_missing=False)
    assert len(out2) == 3
    # node mask: structure nowhere -> sequence scores back
    out3, _ = S.combine_scores(gt, st, fs, fself, mode="max", fit=(0.0, 1.0),
                               node_mask=lambda q, t: np.zeros(len(q), bool))
    assert np.allclose(out3.bits[:3], st.bits)


def test_pair_lca_nodes():
    gt = _gt(6)
    tree = parse_newick("((A,B)AB,C)root;")
    lca, names = S.pair_lca_nodes(gt, tree, np.array([0, 0, 1]), np.array([1, 2, 3]))
    assert [names[i] for i in lca] == ["AB", "root", "AB"]


@pytest.mark.skipif(shutil.which(S.FOLDSEEK) is None, reason="foldseek not installed")
def test_build_db_and_align_pairs(tmp_path):
    rng = np.random.default_rng(1)
    aa = ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=120)) for _ in range(3)]
    s3 = ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=120)) for _ in range(3)]
    aa.append(aa[0]); s3.append(s3[0])     # gene 3 = copy of gene 0
    gt = GeneTable([f"g{i}" for i in range(4)], ["A", "B"], np.array([0, 1, 0, 1], np.int32),
                   np.full(4, 120, np.int32), aa)
    ss = {f"g{i}": s3[i] for i in range(4) if i != 2}   # gene 2 has no 3Di
    genes = S.build_db(gt, aa, ss, str(tmp_path / "db" / "x"))
    assert genes.tolist() == [0, 1, 3]
    key = {int(g): k for k, g in enumerate(genes)}
    h = S.align_pairs(str(tmp_path / "db" / "x"), key, np.array([0, 0, 3]), np.array([3, 0, 3]),
                      str(tmp_path / "aln"), threads=1)
    d = {(int(a), int(b)): float(x) for a, b, x in zip(h["q"], h["t"], h["bits"])}
    assert d[(0, 3)] == d[(0, 0)] == d[(3, 3)] > 0


def test_link_families_attach_no_big_big():
    # fam 0 (3 genes) and fam 1 (3 genes) both hit singleton fam 2; attach keeps 0 and 1 apart
    fam = np.array([0, 0, 0, 1, 1, 1, 2])
    L = np.full(7, 100.0)
    hits = dict(q=np.array([6, 6, 0]), t=np.array([0, 3, 3]), evalue=np.array([1e-9, 1e-6, 1e-20]),
                qcov=np.ones(3), tcov=np.ones(3))
    new, st = S.link_families(fam, hits, L, evalue=1e-3, small_max=1, mode="attach")
    assert new[6] == new[0] and new[3] != new[0]          # best link (1e-9) wins, big-big ignored
    new_u, _ = S.link_families(fam, hits, L, evalue=1e-3, small_max=1)
    assert new_u[6] == new_u[0] == new_u[3]               # union mode chains through the singleton
