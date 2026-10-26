"""Toy checks with known answers for ancestra/treebench.py (S4 tree benchmark)."""
import numpy as np

from ancestra.localtree import nj
from ancestra.treebench import (SpeciesTree, UTree, nrf, pair_labels, paired_bootstrap,
                                quartet_distance, root_overlap, root_reconcile)

ST = SpeciesTree("(((A,B)AB,C)ABC,D)root;")


def sp(name):
    return name[0]  # leaf "A1" -> species "A"


def test_rf_identity_and_one_nni():
    t1 = UTree.from_newick("((A1,B1),(C1,D1),E1);")
    t2 = UTree.from_newick("((A1,C1),(B1,D1),E1);")
    assert nrf(t1, t1)[0] == 0.0
    r, rf, mx, n = nrf(t1, t2)
    assert n == 5 and mx == 4 and rf == 4 and r == 1.0
    assert quartet_distance(t1, t1)[0] == 0.0
    assert quartet_distance(t1, t2)[0] > 0


def test_rf_ignores_rooting_and_prunes():
    t1 = UTree.from_newick("(((A1,B1),C1),(D1,E1));")
    t2 = UTree.from_newick("((A1,B1),(C1,(D1,E1)));")
    assert nrf(t1, t2)[0] == 0.0
    t3 = UTree.from_newick("(((A1,B1),X9),C1,(D1,E1));")
    assert nrf(t1, t3)[0] == 0.0  # X9 is not in t1 -> pruned


def test_root_reconcile_and_labels_duplication():
    # duplication before ABC: two copies (A1,B1,C1) and (A2,B2,C2), plus D1 outgroup
    t = UTree.from_newick("((((A1,B1),C1),((A2,B2),C2)),D1);")
    ch, r, info = root_reconcile(t, sp, ST)
    lab = pair_labels(t, ch, r, sp)
    assert lab[("A1", "B1")] == 1 and lab[("A1", "C1")] == 1
    assert lab[("A1", "B2")] == 0 and lab[("A2", "C1")] == 0
    assert lab[("A1", "D1")] == 1 and lab[("B2", "D1")] == 1
    assert info["dups"] == 1
    # same unrooted tree written with another root gives the same labels
    t2 = UTree.from_newick("(A1,B1,(C1,(((A2,B2),C2),D1)));")
    ch2, r2, _ = root_reconcile(t2, sp, ST)
    assert pair_labels(t2, ch2, r2, sp) == lab
    # species-overlap rooting alone cannot tell the root on D1 from the root between the two
    # copies (1 duplication each): why S4 roots by LCA reconciliation (dups, then losses)
    ch3, r3, info3 = root_overlap(t2, sp)
    assert info3["dups"] == 1


def test_nj_tree_recovers_additive_tree():
    # additive distances of ((A1,B1),(C1,D1)) with long C1 branch: NJ must recover it
    names = ["A1", "B1", "C1", "D1"]
    D = np.array([[0, 2, 7, 6], [2, 0, 7, 6], [7, 7, 0, 7], [6, 6, 7, 0]], float)
    edges, n = nj(D)
    t = UTree.from_edges(edges, n, names)
    assert nrf(t, UTree.from_newick("((A1,B1),(C1,D1));"))[0] == 0.0


def test_join_makes_root_duplication():
    t1 = UTree.from_newick("((A1,B1),C1);")
    t2 = UTree.from_newick("((A2,B2),C2);")
    j = UTree.join(t1, t1.edges()[0], t2, t2.edges()[0])
    assert j.leaf_names() == {"A1", "B1", "C1", "A2", "B2", "C2"}
    ch, r, _ = root_reconcile(j, sp, ST)
    lab = pair_labels(j, ch, r, sp)
    assert lab[("A1", "B2")] == 0 and lab[("A1", "B1")] == 1


def test_paired_bootstrap_sign():
    na = np.array([5, 5, 5, 5.0]); nb = na + 1; den = np.full(4, 10.0)
    d, lo, hi, p = paired_bootstrap(na, den, nb, den, B=200)
    assert abs(d - 0.1) < 1e-12 and lo > 0 and p == 0.0
    d, lo, hi, p = paired_bootstrap(na, den, na, den, B=200)
    assert d == 0 and lo == 0 and hi == 0


def test_3di_model_files_valid():
    import os
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = os.path.join(root, "work/S4/3DI_foldseek.paml")
    if not os.path.exists(p):
        import pytest
        pytest.skip("run scripts/S4_3di_model.py first")
    rows = [l.split() for l in open(p) if l.strip()]
    assert [len(r) for r in rows[:19]] == list(range(1, 20))
    pi = np.array(rows[19], float)
    assert len(pi) == 20 and abs(pi.sum() - 1) < 1e-4
    assert all(float(x) > 0 for r in rows[:19] for x in r)
    q = np.loadtxt(p.replace(".paml", ".fasttree.tsv"), skiprows=1, usecols=range(1, 22))
    assert np.allclose(q[:, :20].sum(0), 0, atol=1e-6)  # FastTree: columns sum to zero
