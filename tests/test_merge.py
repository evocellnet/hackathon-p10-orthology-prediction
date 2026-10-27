"""Tests for the core v0 engine: io, scoring cap, merge (exact hand-made HOGs), pipeline."""
import json
import os

import numpy as np
import pytest

from ancestra.io import GeneTable, load_families, load_guide_tree, load_proteomes
from ancestra.merge import ARMS, MergeConfig, decide_node, run_merge, to_orthoxml
from ancestra.newick import parse_newick
from ancestra.orthoxml import Group, hogs_to_pairs, read_orthoxml, validate_xsd, write_orthoxml
from ancestra.scoring import ScoreTable, ScoringConfig, apply_cap


# ----------------------------------------------------------------- helpers

def canon(x):
    """Order-free canonical form of a nested group (kind, level, children)."""
    if isinstance(x, Group):
        return (x.kind, x.level, frozenset(canon(c) for c in x.children))
    return x


def make_case(genes: dict, scores: dict, tree: str):
    """genes: species -> [gene names]; scores: {(g1, g2): bits}; self-score 1000."""
    names, sp_of, species = [], [], sorted(genes)
    for si, sp in enumerate(species):
        for g in genes[sp]:
            names.append(g)
            sp_of.append(si)
    idx = {g: i for i, g in enumerate(names)}
    gt = GeneTable(names, species, np.array(sp_of, dtype=np.int32),
                   np.full(len(names), 300, dtype=np.int32))
    q, t, b = [], [], []
    for (x, y), v in scores.items():
        a, c = sorted((idx[x], idx[y]))
        q.append(a)
        t.append(c)
        b.append(v)
    n = len(q)
    one = np.ones(n, dtype=np.float32)
    st = ScoreTable(np.array(q, np.int32), np.array(t, np.int32), np.array(b, np.float32), one, one,
                    one, one * 1e-30, np.full(len(names), 1000.0, dtype=np.float32))
    return gt, st, load_guide_tree(tree, species)


def predicted(gt, st, tree, cfg=None):
    res = run_merge(gt, tree, st, cfg or MergeConfig())
    ox = to_orthoxml(res, gt, tree)
    return ox, {canon(g) for g in ox.groups}


# ------------------------------------------------- hand-made exact cases

TREE = "((A,B)AB,C)root;"
GENES = {"A": ["a1", "a1p", "a2x", "a2y", "a3"], "B": ["b1", "b2x", "b2y", "b4"],
         "C": ["c1", "c2", "c3"]}
SCORES = {
    # family 1: duplication on the A branch (a1, a1p in-paralogs w.r.t. B)
    ("a1", "a1p"): 900, ("a1", "b1"): 700, ("a1p", "b1"): 640,
    ("a1", "c1"): 400, ("a1p", "c1"): 390, ("b1", "c1"): 400,
    # family 2: duplication between root and AB, both copies kept in A and B
    ("a2x", "b2x"): 700, ("a2y", "b2y"): 700, ("a2x", "b2y"): 500, ("a2y", "b2x"): 500,
    ("a2x", "a2y"): 500, ("b2x", "b2y"): 500,
    ("a2x", "c2"): 400, ("a2y", "c2"): 400, ("b2x", "c2"): 400, ("b2y", "c2"): 390,
    # family 3: lost in B
    ("a3", "c3"): 450,
    # family 4: b4 singleton (no scores)
}
EXPECTED = {
    ("ortholog", "root", frozenset({
        ("ortholog", "AB", frozenset({("paralog", "A", frozenset({"a1", "a1p"})), "b1"})),
        "c1"})),
    ("ortholog", "root", frozenset({
        ("paralog", "AB", frozenset({("ortholog", "AB", frozenset({"a2x", "b2x"})),
                                     ("ortholog", "AB", frozenset({"a2y", "b2y"}))})),
        "c2"})),
    ("ortholog", "root", frozenset({"a3", "c3"})),
}


@pytest.mark.parametrize("arm", ["max", "max+repair", "avg", "avg+repair", "max+repair/bits"])
def test_hand_made_exact_hogs(arm):
    gt, st, tree = make_case(GENES, SCORES, TREE)
    ox, got = predicted(gt, st, tree, ARMS[arm])
    assert got == EXPECTED
    # pairs: in-paralogs a1/a1p both orthologous to b1; a2x-b2y not orthologs
    pairs, paras = hogs_to_pairs(ox, with_paralogs=True)
    assert ("a1", "b1") in pairs and ("a1p", "b1") in pairs
    assert ("a1", "a1p") in paras
    assert ("a2x", "b2x") in pairs and ("a2x", "b2y") not in pairs and ("a2x", "b2y") in paras
    assert ("a2x", "c2") in pairs and ("a3", "c3") in pairs
    assert not any("b4" in p for p in pairs)


def test_rbh_only_misses_inparalog():
    gt, st, tree = make_case(GENES, SCORES, TREE)
    ox, got = predicted(gt, st, tree, ARMS["rbh"])
    pairs = hogs_to_pairs(ox)
    assert ("a1", "b1") in pairs and ("a1p", "b1") not in pairs


def test_unmatched_passes_up_and_xsd(tmp_path):
    gt, st, tree = make_case(GENES, SCORES, TREE)
    ox, _ = predicted(gt, st, tree)
    p = str(tmp_path / "x.orthoxml")
    write_orthoxml(ox, p)
    ok, errs = validate_xsd(p)
    assert ok, errs
    back = read_orthoxml(p)
    assert {canon(g) for g in back.groups} == EXPECTED


def test_polytomy_kway_with_duplication():
    genes = {"A": ["a"], "B": ["b"], "C": ["c", "cp"], "D": ["d"]}
    scores = {("a", "b"): 600, ("a", "c"): 600, ("b", "c"): 600, ("a", "cp"): 560,
              ("b", "cp"): 560, ("c", "cp"): 900, ("a", "d"): 300, ("b", "d"): 300,
              ("c", "d"): 300, ("cp", "d"): 290}
    gt, st, tree = make_case(genes, scores, "((A,B,C)ABC,D)root;")
    _, got = predicted(gt, st, tree)
    assert got == {("ortholog", "root", frozenset({
        ("ortholog", "ABC", frozenset({"a", "b", ("paralog", "C", frozenset({"c", "cp"}))})),
        "d"}))}


def test_old_dup_not_absorbed_without_inparalog_signal():
    # a_x/a_y paralogs older than AB, partner of a_y lost in B: a_y must stay apart at AB
    genes = {"A": ["ax", "ay"], "B": ["bx"], "C": ["c"]}
    scores = {("ax", "bx"): 700, ("ay", "bx"): 480, ("ax", "ay"): 470,
              ("ax", "c"): 400, ("bx", "c"): 400, ("ay", "c"): 380}
    gt, st, tree = make_case(genes, scores, TREE)
    _, got = predicted(gt, st, tree, ARMS["max"])
    # at AB: {ax,bx} and ay stay separate (ay is an out-paralog w.r.t. AB); at the root
    # both AB units match c and are in-paralogs relative to C -> duplication on root->AB
    assert got == {("ortholog", "root", frozenset({
        ("paralog", "AB", frozenset({("ortholog", "AB", frozenset({"ax", "bx"})), "ay"})),
        "c"}))}


def test_coverage_guard_blocks_merge():
    gt, st, tree = make_case({"A": ["a"], "B": ["b"]}, {("a", "b"): 500}, "(A,B)root;")
    st.qcov[:] = 0.3
    _, got = predicted(gt, st, tree, MergeConfig(min_cov=0.5))
    assert got == set()
    _, got = predicted(gt, st, tree, MergeConfig(min_cov=0.2))
    assert got == {("ortholog", "root", frozenset({"a", "b"}))}


def test_decide_node_repair_joins_leftovers():
    # child0: u0 (best partner taken), child1: v0, v1. u1 (child0) only weakly hits v1.
    child = np.array([0, 0, 1, 1])
    lu = np.array([0, 1, 1, 0])
    lv = np.array([2, 2, 3, 3])
    sc = np.array([1.0, 0.5, 0.4, 0.3])
    cfg = MergeConfig(repair=False)
    lab, _ = decide_node(4, child, lu, lv, sc, cfg)
    assert lab[0] == lab[2] and lab[1] != lab[3]
    lab, st = decide_node(4, child, lu, lv, sc, MergeConfig(repair=True, repair_ratio=0.5))
    assert lab[0] == lab[2] and lab[1] == lab[3] and st["r3_repaired"] == 1
    lab, _ = decide_node(4, child, lu, lv, sc, MergeConfig(repair=True, repair_ratio=0.9))
    assert lab[1] != lab[3]


# ----------------------------------------------------------------- io / scoring

def test_guide_tree_prune_collapse_and_names():
    t = load_guide_tree("(((A,X)AX,B),(C,D)CD)root;", ["A", "B", "C", "D"])
    names = {n.name for n in t.preorder()}
    assert "X" not in names and "AX" not in names       # pruned leaf, collapsed unary node
    assert set(t.leaf_names()) == {"A", "B", "C", "D"}
    assert all(n.name for n in t.preorder())             # unnamed internal node got a name
    with pytest.raises(ValueError):
        load_guide_tree("(A,B)r;", ["A", "Z"])


def test_load_proteomes_and_families(tmp_path):
    d = tmp_path / "p"
    d.mkdir()
    (d / "Sp2.fa").write_text(">sp|P1|X_HUMAN desc\nMKV*\n>tr|Q2|Y\nMAAA\n")
    (d / "Sp1.fasta").write_text(">g1\nMK\nLL\n")
    gt = load_proteomes(str(d))
    assert gt.species == ["Sp1", "Sp2"] and gt.names == ["g1", "P1", "Q2"]
    assert list(gt.lengths) == [4, 3, 4] and list(gt.species_of) == [0, 1, 1]
    f = tmp_path / "f.tsv"
    f.write_text("#gene_id\tfamily_id\ng1\tF1\nP1\tF1\n")
    fam = load_families(str(f), gt)
    assert fam[0] == fam[1] and fam[2] != fam[0]         # Q2 missing -> singleton


def test_apply_cap_top_per_species():
    # one family of 4 genes in 3 species; cap_all=3 forces the reduction
    sp = np.array([0, 1, 1, 2])
    fam = np.zeros(4, dtype=np.int64)
    q = np.array([0, 0, 0, 1, 1, 2], np.int32)
    t = np.array([1, 2, 3, 2, 3, 3], np.int32)
    bits = np.array([100, 90, 80, 70, 60, 50], np.float32)
    one = np.ones(6, np.float32)
    st = ScoreTable(q, t, bits, one, one, one, one, np.full(4, 200, np.float32))
    cfg = ScoringConfig(cap_all=3, top_per_species=1, top_same=0)
    stats = {}
    out = apply_cap(st, fam, sp, cfg, stats)
    kept = set(zip(out.q.tolist(), out.t.tolist()))
    # per query: best hit per other species; the same-species pair (1, 2) is dropped
    assert kept == {(0, 1), (0, 2), (0, 3), (1, 3), (2, 3)}
    assert stats["n_pairs_dropped_by_cap"] == 1 and stats["n_families_over_cap"] == 1


# ----------------------------------------------------------------- pipeline end-to-end

def test_pipeline_end_to_end_resume(sim_basic, tmp_path):
    from ancestra.evaluate_toy import evaluate
    from ancestra.pipeline import run_pipeline
    _, simdir = sim_basic
    out = str(tmp_path / "run")
    logs = []
    m = run_pipeline(os.path.join(simdir, "proteomes"), os.path.join(simdir, "species_tree.nwk"), out,
                     families=os.path.join(simdir, "truth", "gene_families.tsv"), fam_cols=(0, 2),
                     threads=4, log=logs.append)
    for f in ("hogs.orthoxml", "pairs.tsv", "orthobench_ogs.txt", "manifest.json"):
        assert os.path.exists(os.path.join(out, f))
    assert m["stages"]["output"]["stats"]["xsd_valid"]
    man = json.load(open(os.path.join(out, "manifest.json")))
    assert set(man["stages"]) == {"input", "families", "scores", "merge", "output"}
    assert len(man["nodes"]) == 7 and all("wall_s" in n for n in man["nodes"])
    # pairs.tsv == flattening of hogs.orthoxml
    ox = read_orthoxml(os.path.join(out, "hogs.orthoxml"))
    pairs = hogs_to_pairs(ox)
    with open(os.path.join(out, "pairs.tsv")) as fh:
        tsv = {tuple(l.split("\t")[:2]) for l in fh if not l.startswith("#")}
    assert tsv == pairs
    r = evaluate(ox, os.path.join(simdir, "truth", "hogs.orthoxml"))
    assert r["pairs"]["precision"] > 0.95 and r["pairs"]["recall"] > 0.85
    # resume: everything skipped, outputs identical
    logs.clear()
    m2 = run_pipeline(os.path.join(simdir, "proteomes"), os.path.join(simdir, "species_tree.nwk"), out,
                      families=os.path.join(simdir, "truth", "gene_families.tsv"), fam_cols=(0, 2),
                      threads=4, log=logs.append)
    assert all(v["skipped"] for v in m2["stages"].values())
    assert m2["outputs"] == m["outputs"]
    # changing the merge config reruns merge + output only
    m3 = run_pipeline(os.path.join(simdir, "proteomes"), os.path.join(simdir, "species_tree.nwk"), out,
                      families=os.path.join(simdir, "truth", "gene_families.tsv"), fam_cols=(0, 2),
                      threads=4, merge=MergeConfig(link="max"), log=logs.append)
    sk = {k: v["skipped"] for k, v in m3["stages"].items()}
    assert sk == {"input": True, "families": True, "scores": True, "merge": False, "output": False}


def test_join_families_guard():
    from ancestra.scoring import join_families
    fam = np.array([0, 0, 1, 1, 2, 3])
    cross = dict(q=np.array([0, 2, 4]), t=np.array([2, 4, 5]), qcov=np.array([0.9, 0.2, 0.9]),
                 tcov=np.array([0.8, 0.9, 0.9]), evalue=np.array([1e-20, 1e-20, 1.0]))
    new, n = join_families(fam, cross, cov=0.5, evalue=1e-5)
    assert n == 1                                  # 0-1 joined; 1-2 fails coverage; 2-3 fails e-value
    assert new[0] == new[2] and new[2] != new[4] and new[4] != new[5]
    assert sorted(set(new.tolist())) == [0, 1, 2]


# ----------------------------------------------------------------- synteny hook (S5 prep)

def test_gene_order_and_context_hook(sim_basic):
    from ancestra.genome import load_gene_order_gff3, neighbor_pairs
    from ancestra.merge import NodeContext
    from ancestra.orthoxml import hogs_to_pairs
    res_sim, simdir = sim_basic
    gt = load_proteomes(os.path.join(simdir, "proteomes"))
    tree = load_guide_tree(os.path.join(simdir, "species_tree.nwk"), gt.species)
    go = load_gene_order_gff3(os.path.join(simdir, "gff"), gt)
    assert all(v == 1.0 for v in go.coverage.values())
    # ranks follow the simulator's genomic order (protein IDs are numbered in that order)
    sp0 = gt.genes_of_species()[0]
    first_contig = sp0[go.contig[sp0] == go.contig[sp0[0]]]
    assert list(go.rank[first_contig]) == list(range(len(first_contig)))
    a, b, d = neighbor_pairs(go, 5)
    assert (go.contig[a] == go.contig[b]).all() and ((go.rank[b] - go.rank[a]) == d).all()
    assert d.min() == 1 and d.max() == 5
    fam = load_families(os.path.join(simdir, "truth", "gene_families.tsv"), gt, 0, 2)
    from ancestra.scoring import ScoringConfig, score_families
    st, _ = score_families(gt, fam, ScoringConfig(backend="mmseqs", threads=4),
                           os.path.join(simdir, "_hooktest"))
    seen = []

    def null_context(ctx: NodeContext):
        assert len(ctx.cand_a) == len(ctx.cand_b) == len(ctx.cand_score) == len(ctx.anchors)
        assert (ctx.child[ctx.cand_a] != ctx.child[ctx.cand_b]).all()
        ua, ub, nsp, dmin = ctx.adjacency
        assert (ctx.child[ua] == ctx.child[ub]).all()      # adjacency is within-child
        assert (nsp >= 1).all() and (dmin >= 1).all() and (dmin <= 5).all()
        assert len(ctx.members(0)) >= 1
        seen.append((ctx.node, int(ctx.anchors.sum()), len(ua)))
        return np.zeros(len(ctx.cand_a))
    base = run_merge(gt, tree, st, MergeConfig())
    hooked = run_merge(gt, tree, st, MergeConfig(), context=null_context, gene_order=go)
    assert len(seen) == 7 and all(n_adj > 0 for _, _, n_adj in seen)
    p0 = hogs_to_pairs(to_orthoxml(base, gt, tree))
    p1 = hogs_to_pairs(to_orthoxml(hooked, gt, tree))
    assert p0 == p1                                      # zero bonus = sequence-only result
    sh = go.shuffled(11)
    for s_ in range(len(gt.species)):                    # shuffle keeps the per-species gene set
        g = gt.genes_of_species()[s_]
        assert sorted(sh.contig[g]) == sorted(go.contig[g])


@pytest.mark.parametrize("arm", ["default", "rbh", "max+repair", "avg+repair/quartet",
                                 "avg+repair/either", "rd_both"])
def test_numba_engine_matches_python(sim_complex, arm, tmp_path):
    from dataclasses import replace
    from ancestra.orthoxml import hogs_to_pairs
    from ancestra.scoring import ScoringConfig, score_families
    _, simdir = sim_complex
    gt = load_proteomes(os.path.join(simdir, "proteomes"))
    tree = load_guide_tree(os.path.join(simdir, "species_tree.nwk"), gt.species)
    fam = load_families(os.path.join(simdir, "truth", "gene_families.tsv"), gt, 0, 2)
    st, _ = score_families(gt, fam, ScoringConfig(backend="mmseqs", threads=4), str(tmp_path))
    res = {e: hogs_to_pairs(to_orthoxml(run_merge(gt, tree, st, replace(ARMS[arm], engine=e)), gt, tree))
           for e in ("python", "numba")}
    assert res["python"] == res["numba"]


def test_localtree_labelling():
    from ancestra.localtree import tree_groups

    def D_from(pairs, m):
        D = np.zeros((m, m))
        for (i, j), v in pairs.items():
            D[i, j] = D[j, i] = v
        return D
    # fast-evolving in-paralog a2 of a1 (sisters), b1 in the other child, outgroup leaf 3
    D = D_from({(0, 1): 1.7, (0, 2): 1.0, (1, 2): 2.3, (0, 3): 2.5, (1, 3): 3.8, (2, 3): 2.5}, 4)
    assert sorted(map(sorted, tree_groups(D, [0, 0, 1, -1], 3))) == [[0, 1, 2]]
    # out-paralog a2 (duplication before X, partner lost): kept apart
    D = D_from({(0, 1): 2.0, (0, 2): 1.0, (1, 2): 2.0, (0, 3): 3.0, (1, 3): 3.0, (2, 3): 3.0}, 4)
    assert sorted(map(sorted, tree_groups(D, [0, 0, 1, -1], 3))) == [[0, 2], [1]]
    # old duplication with both copies kept, no outgroup: ((a1,b1),(a2,b2))
    D = D_from({(0, 1): 3, (0, 2): 1, (0, 3): 3, (1, 2): 3, (1, 3): 1, (2, 3): 3}, 4)
    assert sorted(map(sorted, tree_groups(D, [0, 0, 1, 1], None))) == [[0, 2], [1, 3]]


@pytest.mark.parametrize("mode", ["medoid", "all"])
def test_localtree_fasttree_backend_runs(sim_complex, tmp_path, mode):
    import shutil
    from dataclasses import replace
    if not (shutil.which("famsa") and shutil.which("FastTree")):
        pytest.skip("famsa/FastTree not available")
    from ancestra.scoring import ScoringConfig, score_families
    _, simdir = sim_complex
    gt = load_proteomes(os.path.join(simdir, "proteomes"))
    tree = load_guide_tree(os.path.join(simdir, "species_tree.nwk"), gt.species)
    fam = load_families(os.path.join(simdir, "truth", "gene_families.tsv"), gt, 0, 2)
    st, _ = score_families(gt, fam, ScoringConfig(backend="mmseqs", threads=4), str(tmp_path))
    cfg = replace(ARMS["default"], tree_require_outgroup=False, tree_min_d=0.0, tree_method="fasttree",
                  tree_ft_mode=mode, tree_ft_workers=2)
    res = run_merge(gt, tree, st, cfg)
    assert sum(s.get("tree_fasttree_jobs", 0) for s in res.node_stats) > 0
    ox = to_orthoxml(res, gt, tree)
    p = str(tmp_path / "x.orthoxml")
    write_orthoxml(ox, p)
    assert validate_xsd(p)[0]
    assert MergeConfig().key() == {k: v for k, v in MergeConfig().key().items()} and \
        "tree_method" not in MergeConfig().key()          # default hash unchanged


def test_localtree_fasttree_deterministic_across_workers(sim_complex, tmp_path):
    import shutil
    from dataclasses import replace
    from ancestra.orthoxml import hogs_to_pairs
    if not (shutil.which("famsa") and shutil.which("FastTree")):
        pytest.skip("famsa/FastTree not available")
    from ancestra.scoring import ScoringConfig, score_families
    _, simdir = sim_complex
    gt = load_proteomes(os.path.join(simdir, "proteomes"))
    tree = load_guide_tree(os.path.join(simdir, "species_tree.nwk"), gt.species)
    fam = load_families(os.path.join(simdir, "truth", "gene_families.tsv"), gt, 0, 2)
    st, _ = score_families(gt, fam, ScoringConfig(backend="mmseqs", threads=4), str(tmp_path))
    base = replace(ARMS["default"], tree_require_outgroup=False, tree_min_d=0.0, tree_method="fasttree")
    outs = []
    for w in (2, 4, 4):
        res = run_merge(gt, tree, st, replace(base, tree_ft_workers=w))
        assert sum(s.get("tree_fasttree_jobs", 0) for s in res.node_stats) > 0
        outs.append(([s.get("tree_changed") for s in res.node_stats],
                     hogs_to_pairs(to_orthoxml(res, gt, tree))))
    assert outs[0] == outs[1] == outs[2]


def test_pipeline_presets_do_not_share_cache(sim_basic, tmp_path):
    import shutil
    from ancestra.pipeline import run_pipeline
    if not (shutil.which("famsa") and shutil.which("FastTree")):
        pytest.skip("famsa/FastTree not available")
    _, simdir = sim_basic
    out = str(tmp_path / "run")
    kw = dict(families=os.path.join(simdir, "truth", "gene_families.tsv"), fam_cols=(0, 2), threads=2,
              log=lambda s: None)
    args = (os.path.join(simdir, "proteomes"), os.path.join(simdir, "species_tree.nwk"), out)
    m1 = run_pipeline(*args, preset="fast", **kw)
    m2 = run_pipeline(*args, preset="accurate", **kw)
    m3 = run_pipeline(*args, preset="accurate", **kw)
    assert m1["config"]["preset"] == "fast" and m2["config"]["preset"] == "accurate"
    assert m2["stages"]["merge"]["skipped"] is False            # preset change -> new merge
    assert m2["stages"]["scores"]["skipped"] is True
    assert m3["stages"]["merge"]["skipped"] is True              # same preset -> cached
    assert m2["stages"]["merge"]["key"] != m1["stages"]["merge"]["key"]
    with pytest.raises(ValueError):
        run_pipeline(*args, preset="nonsense", **kw)


def test_species_overlap_missed_speciation():
    from ancestra.orthoxml import hogs_to_pairs
    genes = {"A1": ["a1"], "A2": ["a2"], "B": ["b"]}
    scores = {("a1", "b"): 700, ("a2", "b"): 650}          # a1-a2 unscored: fission at A
    gt, st, tree = make_case(genes, scores, "((A1,A2)A,B)root;")
    base = MergeConfig(tree_mode=False, repair=False)
    ox0, _ = predicted(gt, st, tree, base)
    assert ("a2", "b") not in hogs_to_pairs(ox0)              # a2 left alone without the fix
    ox1, _ = predicted(gt, st, tree, MergeConfig(tree_mode=False, repair=False, so_absorb=True))
    o1, p1 = hogs_to_pairs(ox1, with_paralogs=True)
    assert ("a2", "b") in o1 and ("a1", "a2") in p1           # absorbed, but as paralogs
    ox2, got = predicted(gt, st, tree, MergeConfig(tree_mode=False, repair=False, so_absorb=True,
                                                   so_struct=True))
    o2 = hogs_to_pairs(ox2)
    assert {("a1", "a2"), ("a1", "b"), ("a2", "b")} <= o2     # species-disjoint -> speciation at A
    assert got == {("ortholog", "root", frozenset({("ortholog", "A", frozenset({"a1", "a2"})), "b"}))}


def test_deep_node_relaxed_tolerance():
    # child0: a ; child1: b1, b2. Scores (normalised) a-b1 0.15, a-b2 0.13: divergent node (d~1.9)
    child = np.array([0, 1, 1])
    lu, lv, sc = np.array([0, 0]), np.array([1, 2]), np.array([0.15, 0.13])
    _, st0 = decide_node(3, child, lu, lv, sc, MergeConfig(tree_mode=False, deep_min_d=2.5, deep_tol=0.8))
    assert "deep" not in st0 and st0["n_mb"] == 1                 # threshold not reached
    _, st1 = decide_node(3, child, lu, lv, sc, MergeConfig(tree_mode=False, deep_min_d=1.5, deep_tol=0.8))
    assert st1["deep"] == 1 and st1["n_mb_deep"] == 2              # a-b2 becomes mutual-best at 0.8
    assert "deep_min_d" not in MergeConfig().key()                 # default hash unchanged


def _nesting_violations(ox, tree):
    """(same-level OG in OG, paralogGroup level not strictly below enclosing OG or not containing
    members, OG TaxRange != LCA of its members)."""
    from ancestra.newick import TreeIndex
    tix = TreeIndex(tree)
    sp = {k: v.species for k, v in ox.genes.items()}
    bad = [0, 0, 0]
    for top in ox.groups:
        st = [(top, None)]
        while st:
            g, par = st.pop()
            mem = [sp[x] for x in g.genes()]
            if g.kind == "ortholog":
                if tix.mrca(mem) != g.level:
                    bad[2] += 1
                if par is not None and par.kind == "ortholog" and not tix.is_strict_descendant(g.level, par.level):
                    bad[0] += 1
            else:
                if not all(tix.is_ancestor(g.level, x) for x in mem) or \
                        (par is not None and not tix.is_strict_descendant(g.level, par.level)):
                    bad[1] += 1
            st.extend((c, g) for c in g.children if isinstance(c, Group))
    return bad


def test_so_struct_nesting_regression():
    # repro of results/qfo_run/repro_so_struct_nesting (Methanobacteriota polytomy with a
    # Stenosarchaea subclade): two separate ancestral genes {3,4} and {1,2} at "Metha" become
    # species-disjoint units at the root; the combined speciation must nest correctly.
    genes = {"S536": ["g1"], "S554": ["g2"], "S487": ["g3"], "S805": ["g4"], "O": ["o"]}
    scores = {("g1", "g2"): 600, ("g3", "g4"): 600, ("o", "g3"): 300, ("o", "g1"): 280}
    gt, st, tree = make_case(genes, scores, "((S536,(S554,S487)Steno,S805)Metha,O)root;")
    cfg = MergeConfig(tree_mode=False, so_struct=True, so_absorb=True)
    ox, got = predicted(gt, st, tree, cfg)
    assert got == {("ortholog", "root", frozenset({
        ("ortholog", "Metha", frozenset({"g1", "g4", ("ortholog", "Steno", frozenset({"g2", "g3"}))})),
        "o"}))}
    assert _nesting_violations(ox, tree) == [0, 0, 0]


def test_so_struct_rebucketed_paralog_level():
    # a paralogGroup created high up (level B, members all of species B1) that so_combine
    # re-buckets under a deeper orthologGroup must be relabelled strictly below it
    genes = {"B1": ["b1", "b1p"], "B2": ["b2"], "C": ["c"], "O": ["o"]}
    scores = {("b1", "b1p"): 900, ("b1", "c"): 500, ("b1p", "c"): 480, ("b2", "o"): 300,
              ("c", "o"): 250, ("b1", "o"): 240}
    gt, st, tree = make_case(genes, scores, "(((B1,B2)B,C)BC,O)root;")
    for cfg in (MergeConfig(tree_mode=False, so_struct=True, so_absorb=True),
                MergeConfig(tree_mode=False, so_struct=True)):
        ox, _ = predicted(gt, st, tree, cfg)
        assert _nesting_violations(ox, tree) == [0, 0, 0]


def test_fix_paralog_levels_flattens_and_relabels():
    from ancestra.merge import _fix_paralog_levels
    from ancestra.newick import TreeIndex, parse_newick
    from ancestra.orthoxml import GeneEntry, OrthoXML, hogs_to_pairs
    tree = parse_newick("((L,(A,B)Olf)Chordata,P)root;")
    tix = TreeIndex(tree)
    sp = {"l1": "L", "l2": "L", "l3": "L", "a": "A", "b": "B", "p": "P"}
    bad = Group("ortholog", [Group("paralog", [Group("paralog", ["l1", "l2"], level="Chordata"), "l3"],
                                   level="L"),
                             Group("ortholog", ["a", "b"], level="Olf")], level="Chordata")
    top = Group("ortholog", [bad, "p"], level="root")
    ox = OrthoXML(genes={g: GeneEntry(g, s) for g, s in sp.items()}, groups=[top], species_tree=tree)
    before = hogs_to_pairs(ox, with_paralogs=True)
    assert _fix_paralog_levels(top, sp, tix) >= 1
    par = top.children[0].children[0]
    assert par.kind == "paralog" and par.level == "L" and set(par.children) == {"l1", "l2", "l3"}
    assert hogs_to_pairs(ox, with_paralogs=True) == before          # pairs unchanged
