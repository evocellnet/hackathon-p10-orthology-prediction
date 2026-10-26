"""Tests for ancestra.speciestree (S1): toy data with known answers + controls."""
import os

import numpy as np
import pytest

from ancestra import speciestree as st
from ancestra.newick import parse_newick

T8 = "((((a:1,b:1)AB:1,(c:1,d:1)CD:1)ABCD:1,((e:1,f:1)EF:1,(g:1,h:1)GH:1)EFGH:1)root);"


def tree():
    return parse_newick(T8)


def test_compare_identity_is_zero():
    r = st.compare_trees(tree(), tree())
    assert r["rf"] == 0 and r["nrf"] == 0.0
    assert r["ref_split_recall"] == 1.0 and r["conflict_frac"] == 0.0
    assert r["quartet"]["qd_all"] == 0.0


def test_rooting_does_not_change_unrooted_splits():
    t = tree()
    rr = st.reroot(t, t.find("a"))
    assert sorted(rr.leaf_names()) == sorted(t.leaf_names())
    assert st.compare_trees(rr, t)["rf"] == 0


def test_known_rf_one_nni():
    # swap b and c: splits ab, cd removed; ac, bd added
    q = parse_newick("((((a,c),(b,d)),((e,f),(g,h))));")
    r = st.compare_trees(q, tree())
    assert r["rf"] == 4
    assert r["n_splits_query"] == 5 and r["n_splits_ref"] == 5
    assert abs(r["nrf"] - 0.4) < 1e-12


def test_polytomy_refinement_is_not_conflict():
    ref = parse_newick("((a,b,c,d),(e,f,g,h));")   # only one split
    r = st.compare_trees(tree(), ref)
    assert r["ref_split_recall"] == 1.0
    assert r["conflict_frac"] == 0.0           # tree refines the polytomies
    assert r["nrf"] > 0                        # but plain nRF counts refinements
    assert r["quartet"]["qd_resolved"] == 0.0


def test_quartets_brute_force_small():
    q = parse_newick("((a,c),(b,d),e);")
    r = parse_newick("((a,b),(c,d),e);")
    res = st.compare_trees(q, r)["quartet"]
    # 5 quartets; abcd differs; others: abce, abde, acde, bcde
    assert res["n_quartets"] == 5
    # brute force via four point on unit distances
    assert res["qd_all"] == 1.0   # all five quartets differ (worked by hand)


def test_nj_recovers_additive_tree():
    t = parse_newick("(((a:1,b:2):1,c:3):1,(d:2,e:1):2);")
    leaves = sorted(t.leaf_names())
    # path lengths
    D = np.zeros((5, 5))
    for i, x in enumerate(leaves):
        for j, y in enumerate(leaves):
            D[i, j] = _path_len(t, x, y)
    nt = st.nj(D, leaves)
    assert st.compare_trees(nt, t)["rf"] == 0


def _path_len(t, x, y):
    def anc(n):
        out, d = {}, 0.0
        while n is not None:
            out[id(n)] = d
            d += n.length or 0.0
            n = n.parent
        return out
    ax, ay = anc(t.find(x)), anc(t.find(y))
    return min(ax[k] + ay[k] for k in ax if k in ay)


def test_midpoint_root():
    t = parse_newick("((a:1,b:1):1,(c:1,d:9):1);")
    m = st.midpoint_root(t)
    # midpoint lies on d's branch: root children are d and the rest
    sizes = sorted(len(c.leaves()) for c in m.children)
    assert sizes == [1, 3]
    assert st.compare_trees(m, t)["rf"] == 0


def test_root_by_reference_recovers_root():
    t = tree()
    unrooted = st.reroot(t, t.find("a"))       # wrong root
    rr, info = st.root_by_reference(unrooted, t)
    names = st.clades(rr)
    assert frozenset("abcd") in names and frozenset("efgh") in names
    assert info["ref_clades_matched"] == info["ref_clades_total"]


def test_name_internal_nodes():
    q = parse_newick("((((a,c),(b,d)),((e,f),(g,h))));")
    st.name_internal_nodes(q, tree())
    names = {n.name for n in q.preorder() if n.children}
    assert {"EF", "GH", "EFGH", "ABCD", "root"} <= names
    assert any(n.startswith("ABCD~") for n in names)   # synthetic, named after LCA
    assert len(names) == sum(1 for n in q.preorder() if n.children)


def test_nni_perturb_deterministic_and_moves():
    t = tree()
    a1, d1 = st.nni_perturb(t, 2, 11)
    a2, _ = st.nni_perturb(t, 2, 11)
    assert a1.to_newick() == a2.to_newick()
    assert sorted(a1.leaf_names()) == sorted(t.leaf_names())
    assert d1 == 2
    # 0 NNIs -> identical splits
    a0, _ = st.nni_perturb(t, 0, 11)
    assert st.compare_trees(a0, t)["rf"] == 0


def test_nni_changes_rf_on_average():
    t = parse_newick(st.random_topology([f"x{i}" for i in range(30)], 5).to_newick())
    rfs = [st.compare_trees(st.nni_perturb(t, 3, s)[0], t, quartets=False)["rf"] for s in (11, 23, 47)]
    assert max(rfs) > 0 and all(r <= 2 * 3 for r in rfs)


def test_random_topology_negative_control():
    leaves = [f"x{i}" for i in range(40)]
    t1 = st.random_topology(leaves, 11)
    t2 = st.random_topology(leaves, 23)
    assert sorted(t1.leaf_names()) == sorted(leaves)
    assert all(len(n.children) == 2 for n in t1.preorder() if n.children)
    assert st.compare_trees(t1, t2, quartets=False)["nrf"] > 0.8   # near-maximal


def test_ncbi_tree_from_toy_taxdump(tmp_path):
    d = tmp_path
    # 1 root; 2 A(clade); 3 B(clade) under A; 10,11 under B; 12 under A; 20 under 1 (other)
    nodes = [(1, 1, "no rank"), (2, 1, "clade"), (3, 2, "clade"), (10, 3, "species"),
             (11, 3, "species"), (12, 2, "species"), (4, 1, "clade"), (20, 4, "species"),
             (5, 4, "clade"), (21, 5, "species")]
    with open(d / "nodes.dmp", "w") as fh:
        for t, p, r in nodes:
            fh.write(f"{t}\t|\t{p}\t|\t{r}\t|\tXX\t|\n")
    names = {1: "root", 2: "Alpha", 3: "Beta", 10: "sp10", 11: "sp11", 12: "sp12", 4: "Beta",
             20: "sp20", 5: "Gamma", 21: "sp21"}
    with open(d / "names.dmp", "w") as fh:
        for t, n in names.items():
            fh.write(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n")
    (d / "merged.dmp").write_text("99\t|\t21\t|\n")
    t = st.ncbi_tree({"L10": 10, "L11": 11, "L12": 12, "L20": 20, "L21": 99}, str(d))
    nw = t.to_newick(lengths=False)
    # Gamma (5) and node 4 are unary -> collapsed; duplicate 'Beta' gets taxid suffix
    internal = sorted(n.name for n in t.preorder() if n.children)
    assert internal == ["Alpha", "Beta_3", "Beta_4", "root"] or internal == ["Alpha", "Beta_3", "Beta_4", "root"]
    assert all(len(n.children) >= 2 for n in t.preorder() if n.children)
    assert st.clades(t)[frozenset({"L10", "L11"})] == "Beta_3"
    assert "L21" in nw


def test_families_from_s0(tmp_path):
    p = tmp_path / "fam.tsv"
    rows = [("g1", "F1"), ("g2", "F1"), ("g3", "F1"), ("g4", "F2"), ("g5", "F2"), ("g6", "F2"),
            ("g7", "F2"), ("g8", "F2")]
    p.write_text("#gene_id\tfamily_id\n" + "".join(f"{a}\t{b}\n" for a, b in rows))
    gs = {"g1": "A", "g2": "B", "g3": "C", "g4": "A", "g5": "A", "g6": "A", "g7": "B", "g8": "C"}
    f = st.families_from_s0(str(p), gs, 3)
    assert set(f) == {"F1"}       # F2 has 3 copies in A


def test_refine_polytomies_adds_only_compatible_clades():
    ref = parse_newick("(((a,b,c,d)X,(e,f,g,h)Y)Z,(i,j)W)R;")
    data = parse_newick("((((a,b),(c,d)),((e,g),(f,h))),(i,j));")
    t, added = st.refine_polytomies(ref, data)
    cl = st.clades(t)
    assert frozenset("ab") in cl and frozenset("cd") in cl and frozenset("eg") in cl
    assert added == 4
    # never removes a reference clade
    assert frozenset("abcd") in cl and frozenset("efgh") in cl
    # an incompatible data clade (cuts X) is not added
    bad = parse_newick("((((a,e),(b,c,d)),(f,g,h)),(i,j));")
    t2, _ = st.refine_polytomies(ref, bad)
    assert frozenset("ae") not in st.clades(t2)
    assert st.compare_trees(t2, ref)["conflict_frac"] == 0.0


def test_refine_polytomies_keeps_root_polytomy():
    ref = parse_newick("((a,b)A,(c,d)C,(e,f)E)R;")
    data = parse_newick("(((a,b),(c,d)),(e,f));")
    t, added = st.refine_polytomies(ref, data)
    assert added == 0 and len(t.children) == 3
    t2, added2 = st.refine_polytomies(ref, data, refine_root=True)
    assert added2 == 1 and frozenset("abcd") in st.clades(t2)
