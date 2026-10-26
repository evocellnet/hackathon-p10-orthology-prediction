"""Tests for ancestra.orthoxml: flattening, levels, round trip, schema, pyham, OrthoBench."""
import io
import os
import sys
from bisect import bisect_right  # noqa: F401  (used by the QfO reference code)

import pytest

from ancestra.newick import TreeIndex, parse_newick, read_newick
from ancestra.orthoxml import (GeneEntry, Group, OrthoXML, hogs_at_level, hogs_to_pairs,
                               normalize_groups, read_orthoxml, to_orthobench_ogs, validate_xsd,
                               write_orthoxml)

from conftest import ORTHOBENCH_DIR, ORTHOBENCH_POSITIVE, QFO_DIR

TREE = "((A:1,B:1)AB:1,C:2)ABC;"


def hand_example() -> OrthoXML:
    """ABC speciation; a duplication on the branch above AB (both copies kept in A and B).

        orthologGroup ABC
          paralogGroup (dup above AB)
            orthologGroup AB: A1 B1
            orthologGroup AB: A2 B2
          C1
    plus a second family with an in-paralog pair in C: orthologGroup ABC(A3, paralogGroup(C2, C3)).
    """
    genes = {g: GeneEntry(g, g[0]) for g in ["A1", "A2", "A3", "B1", "B2", "C1", "C2", "C3", "C4"]}
    fam1 = Group("ortholog", [
        Group("paralog", [Group("ortholog", ["A1", "B1"], level="AB"),
                          Group("ortholog", ["A2", "B2"], level="AB")], level="AB"),
        "C1"], level="ABC", id="F1")
    fam2 = Group("ortholog", ["A3", Group("paralog", ["C2", "C3"], level="C")], level="ABC", id="F2")
    return OrthoXML(genes=genes, groups=[fam1, fam2], species_tree=parse_newick(TREE))


def P(*pairs):
    return {tuple(sorted(p)) for p in pairs}


# ---------------------------------------------------------------- flattening
def test_flattening_hand_example():
    ox = hand_example()
    orth, para = hogs_to_pairs(ox, with_paralogs=True)
    assert orth == P(("A1", "B1"), ("A2", "B2"), ("A1", "C1"), ("B1", "C1"), ("A2", "C1"),
                     ("B2", "C1"), ("A3", "C2"), ("A3", "C3"))
    assert para == P(("A1", "A2"), ("A1", "B2"), ("B1", "A2"), ("B1", "B2"), ("C2", "C3"))
    lv = hogs_to_pairs(ox, with_level=True)
    assert lv[("A1", "B1")] == "AB" and lv[("A1", "C1")] == "ABC"
    # C4 is a singleton: in no relation
    assert not any("C4" in p for p in orth | para)


def test_same_species_pairs_under_orthologgroup_are_dropped():
    genes = {g: GeneEntry(g, g[0]) for g in ["A1", "A2", "B1"]}
    ox = OrthoXML(genes=genes, groups=[Group("ortholog", ["A1", "A2", "B1"])])
    assert hogs_to_pairs(ox) == P(("A1", "B1"), ("A2", "B1"))
    o, p = hogs_to_pairs(ox, with_paralogs=True)
    assert p == set()


def _qfo_pairs(xml_bytes: bytes, ox: OrthoXML):
    """Run the unmodified QfO map_relations.py extraction on an OrthoXML document."""
    sys.path.insert(0, QFO_DIR)
    try:
        import map_relations as mr
    finally:
        sys.path.pop(0)
    species = ox.species
    mapping, goff, k = {}, [0], 0
    for sp in species:
        for g in (g for g in ox.genes.values() if g.species == sp):
            k += 1
            mapping[g.prot_id] = k
        goff.append(k)
    inv = {v: kk for kk, v in mapping.items()}

    class DB:
        def __init__(self):
            self.pairs = set()

        def add_orthologs(self, a, b):
            self.pairs.add(tuple(sorted((inv[a], inv[b]))))

    db = DB()
    proc = mr.PairwiseOrthologRelationExtractor({"mapping": mapping, "Goff": goff[:-1],
                                                 "species": species}, db)
    mr.parse_orthoxml(io.BytesIO(xml_bytes), proc)
    return db.pairs


@pytest.mark.skipif(not os.path.exists(os.path.join(QFO_DIR, "map_relations.py")),
                    reason="QfO benchmark-webservice checkout not available")
@pytest.mark.parametrize("which", ["hand", "basic", "complex"])
def test_flattening_matches_qfo_reference(which, request):
    if which == "hand":
        ox = hand_example()
    else:
        _, out = request.getfixturevalue(f"sim_{which}")
        ox = read_orthoxml(os.path.join(out, "truth", "hogs.orthoxml"))
    buf = io.StringIO()
    write_orthoxml(ox, buf)
    assert _qfo_pairs(buf.getvalue().encode(), ox) == hogs_to_pairs(ox)


# ---------------------------------------------------------------- levels
def test_hogs_at_level_hand_example():
    ox = hand_example()
    S = lambda groups: {frozenset(g) for g in groups}  # noqa: E731
    assert S(hogs_at_level(ox, "ABC")) == S([{"A1", "B1", "A2", "B2", "C1"}, {"A3", "C2", "C3"}, {"C4"}])
    assert S(hogs_at_level(ox, "AB")) == S([{"A1", "B1"}, {"A2", "B2"}, {"A3"}])
    # dup in C is after the C "speciation" level? level C is a leaf: dup on branch above C -> separate
    assert S(hogs_at_level(ox, "C")) == S([{"C1"}, {"C2"}, {"C3"}, {"C4"}])
    assert S(hogs_at_level(ox, None)) == S([{"A1", "B1", "A2", "B2", "C1"}, {"A3", "C2", "C3"}, {"C4"}])


def test_hogs_at_level_mrca_fallback_for_paralog_without_level():
    ox = hand_example()
    for g in ox.groups[0].iter_groups():
        if g.kind == "paralog":
            g.level = None     # MRCA of {A,B} = AB -> dup not below AB -> still separate at AB
    S = lambda groups: {frozenset(g) for g in groups}  # noqa: E731
    assert S(hogs_at_level(ox, "AB", include_singletons=False)) == S([{"A1", "B1"}, {"A2", "B2"}])


def test_normalize_collapses_single_child_groups():
    g = Group("ortholog", [Group("paralog", ["x"]), Group("ortholog", ["y", Group("paralog", [])])])
    out = normalize_groups([g, Group("ortholog", ["z"])])
    assert len(out) == 1 and out[0].children == ["x", "y"]


# ---------------------------------------------------------------- round trip + schema
def test_round_trip_identical(sim_complex, tmp_path):
    _, out = sim_complex
    p = os.path.join(out, "truth", "hogs.orthoxml")
    ox = read_orthoxml(p)
    q = tmp_path / "rt.orthoxml"
    write_orthoxml(ox, str(q))
    assert open(p).read() == q.read_text()
    ox2 = read_orthoxml(str(q))
    assert hogs_to_pairs(ox2, with_paralogs=True) == hogs_to_pairs(ox, with_paralogs=True)
    assert [g.level for t in ox.groups for g in t.iter_groups()] == \
        [g.level for t in ox2.groups for g in t.iter_groups()]
    assert ox2.species_tree.to_newick(lengths=False) == ox.species_tree.to_newick(lengths=False)


def test_round_trip_hand_example():
    ox = hand_example()
    buf = io.StringIO()
    write_orthoxml(ox, buf)
    ox2 = read_orthoxml(io.BytesIO(buf.getvalue().encode()))
    assert set(ox2.genes) == set(ox.genes)
    assert hogs_to_pairs(ox2, with_paralogs=True) == hogs_to_pairs(ox, with_paralogs=True)
    lv = sorted((g.kind, g.level) for t in ox2.groups for g in t.iter_groups())
    assert lv == sorted((g.kind, g.level) for t in ox.groups for g in t.iter_groups())


@pytest.mark.parametrize("version", ["0.3", "0.4", "0.5"])
def test_schema_valid(sim_complex, tmp_path, version):
    _, out = sim_complex
    ox = read_orthoxml(os.path.join(out, "truth", "hogs.orthoxml"))
    q = str(tmp_path / f"v{version}.orthoxml")
    write_orthoxml(ox, q, version=version)
    ok, errs = validate_xsd(q)
    assert ok, errs[:5]
    # ortholog pairs do not depend on the version
    assert hogs_to_pairs(q) == hogs_to_pairs(ox)


def test_schema_negative_control():
    bad = b"""<?xml version="1.0"?><orthoXML xmlns="http://orthoXML.org/2011/" origin="x"
    originVersion="1" version="0.5"><species name="A"><database name="d" version="1"><genes>
    <gene id="1" protId="a"/></genes></database></species><groups><orthologGroup>
    <geneRef id="1"/></orthologGroup></groups></orthoXML>"""
    ok, errs = validate_xsd(bad)
    assert not ok and errs
    dangling = bad.replace(b'<geneRef id="1"/>', b'<geneRef id="1"/><geneRef id="7"/>')
    assert not validate_xsd(dangling)[0]


def test_truth_file_valid(sim_basic):
    _, out = sim_basic
    ok, errs = validate_xsd(os.path.join(out, "truth", "hogs.orthoxml"))
    assert ok, errs[:5]


# ---------------------------------------------------------------- pyham
@pytest.mark.parametrize("which", ["basic", "complex"])
def test_pyham_loads_and_agrees(which, request, tmp_path):
    pyham = pytest.importorskip("pyham")
    import logging
    logging.disable(logging.WARNING)
    try:
        _, out = request.getfixturevalue(f"sim_{which}")
        nwk = os.path.join(out, "species_tree.nwk")
        ox = read_orthoxml(os.path.join(out, "truth", "hogs.orthoxml"))
        q = str(tmp_path / "compat.orthoxml")
        write_orthoxml(ox, q, pyham_compat=True)
        assert validate_xsd(q)[0]
        assert hogs_to_pairs(q) == hogs_to_pairs(ox)
        ham = pyham.Ham(nwk, q, use_internal_name=True, tree_format="newick")
        # pyham places duplications at the MRCA of their genes and a top-level HOG at its own
        # level; apply the same conventions to our reader and compare HOGs at every level.
        oxc = read_orthoxml(q)
        for t in oxc.groups:
            for g in t.iter_groups():
                if g.kind == "paralog":
                    g.level = None
        tree = read_newick(nwk)
        tix = TreeIndex(tree)
        sp = {k: v.species for k, v in oxc.genes.items()}
        top_level = {}
        for t in oxc.groups:
            lvl = t.level or tix.mrca(sp[x] for x in t.genes())
            for x in t.genes():
                top_level[x] = lvl
        for lvl in tix.internal:
            ph = {frozenset(x.prot_id for x in h.get_all_descendant_genes())
                  for h in ham.get_ancestral_genome_by_name(lvl).genes}
            ours = {frozenset(h) for h in hogs_at_level(oxc, lvl, tree)
                    if next(iter(h)) in top_level and tix.is_ancestor(top_level[next(iter(h))], lvl)}
            assert ph == ours, lvl
    finally:
        logging.disable(logging.NOTSET)


# ---------------------------------------------------------------- OrthoBench export
def test_orthobench_export_format(tmp_path):
    ox = hand_example()
    f = tmp_path / "ogs.txt"
    ogs = to_orthobench_ogs(ox, out=str(f))
    lines = f.read_text().splitlines()
    assert lines[0] == "Orthogroup\tGenes"
    assert lines[1].startswith("OG0000000:\t")
    assert [sorted(x) for x in ogs] == [sorted(l.split("\t")[1:]) for l in lines[1:]]
    allg = [g for og in ogs for g in og]
    assert len(allg) == len(set(allg)) == len(ox.genes)
    ogs_ab = to_orthobench_ogs(ox, level="AB", include_singletons=False)
    assert sorted(map(sorted, ogs_ab)) == [["A1", "B1"], ["A2", "B2"]]


@pytest.mark.skipif(not (os.path.exists(ORTHOBENCH_POSITIVE) and os.path.isdir(ORTHOBENCH_DIR)),
                    reason="Open_Orthobench data not available")
def test_orthobench_official_reader_positive_control(tmp_path):
    """RefOGs -> OrthoXML -> to_orthobench_ogs -> official benchmark.py scoring gives F=100."""
    import contextlib
    import importlib.util
    import re
    spec = importlib.util.spec_from_file_location("ob", os.path.join(ORTHOBENCH_DIR, "benchmark.py"))
    ob = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ob)
    refogs = []
    for line in open(ORTHOBENCH_POSITIVE):
        genes = re.findall(r"[A-Za-z]+\d+(?:\.\d+)?", line.split(":", 1)[1]) if ":" in line else []
        if genes:
            refogs.append(genes)
    prefix = lambda g: re.match(r"[A-Za-z]+", g).group(0)  # noqa: E731
    genes = {g: GeneEntry(g, prefix(g)) for og in refogs for g in og}
    ox = OrthoXML(genes=genes, groups=[Group("ortholog", og) for og in refogs])
    buf = io.StringIO()
    write_orthoxml(ox, buf)
    ox = read_orthoxml(io.BytesIO(buf.getvalue().encode()))
    f = str(tmp_path / "ogs.txt")
    to_orthobench_ogs(ox, out=f)
    d_ref = os.path.join(ORTHOBENCH_DIR, "RefOGs") + os.sep
    ref = ob.read_refogs(d_ref)
    unc = ob.read_uncertain_refogs(d_ref + "low_certainty_assignments/")
    pred_smart = ob.read_orthogroups_smart(f)
    pred_basic = ob.read_orthogroups(f, set(genes), ob.get_n_col_skip(f))
    # the basic reader drops single-gene OGs by design; the regex reader keeps them
    assert {frozenset(x) for x in pred_smart if len(x) > 1} == {frozenset(x) for x in pred_basic}
    assert sum(len(x) for x in pred_smart) == len(genes)
    with contextlib.redirect_stdout(io.StringIO()):
        fscore, prec, rec = ob.calculate_benchmarks_pairwise(ref, unc, pred_basic)
    assert round(fscore, 6) == 100.0 and round(prec, 6) == 100.0 and round(rec, 6) == 100.0
