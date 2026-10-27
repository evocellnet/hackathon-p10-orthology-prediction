"""Tests for ancestra.rescue (S0 family rescue after S2; --family-rescue)."""
import os
import shutil

import numpy as np
import pytest

from ancestra import rescue as RS
from ancestra.io import GeneTable, load_guide_tree
from ancestra.merge import MergeResult
from ancestra.scoring import ScoreTable

AA = list("ACDEFGHIKLMNPQRSTVWY")


def _hits(q, t, e, cov=0.9, bits=None):
    q, t, e = np.asarray(q), np.asarray(t), np.asarray(e, float)
    return dict(q=q, t=t, evalue=e, bits=np.asarray(bits if bits is not None else -np.log10(e), float),
                qcov=np.full(len(q), cov), tcov=np.full(len(q), cov))


def test_attach_rules():
    # fam 0: 12 genes (target), fam 1: 11 genes (target), fam 2: singleton, fam 3: 2 genes, fam 4: singleton
    fam = np.array([0] * 12 + [1] * 11 + [2] + [3, 3] + [4])
    L = np.full(len(fam), 100.0)
    g2, g3a, g4 = 23, 24, 26
    cfg = RS.RescueConfig(mode="3di", small_max=10, max_size=7500)
    hits = [("a", _hits([g2, g2, g3a, 0, g4], [0, 12, g2, 12, 1], [1e-5, 1e-9, 1e-20, 1e-30, 1e-2]))]
    new, att, st = RS.attach_small(fam, hits, L, cfg)
    assert new[g2] == 1                          # best link (lowest e) wins: fam 1
    assert new[0] == 0 and new[12] == 1          # two targets are never joined (0-12 hit ignored)
    assert new[g3a] == 3 and new[g3a + 1] == 3   # small -> small links are not used
    assert new[g4] == 4                          # e > 1e-3: not attached
    assert st["n_attached"] == 1 and st["n_genes_attached"] == 1
    # coverage of the shorter protein < 0.3 -> no link
    h = _hits([g2], [0], [1e-10], cov=0.1)
    new, _, _ = RS.attach_small(fam, [("a", h)], L, cfg)
    assert new[g2] == 2
    # size cap: fam 1 (11) + singleton (1) = 12 ok; + 2-gene family -> 14 > 13 skipped
    cfg13 = RS.RescueConfig(mode="3di", small_max=10, max_size=13)
    h = _hits([g2, g3a], [12, 13], [1e-10, 1e-5])
    new, att, st = RS.attach_small(fam, [("a", h)], L, cfg13)
    assert new[g2] == 1 and new[g3a] == 3 and st["n_capped"] == 1
    # union of two sources: lower e-value wins whatever the source
    new, att, st = RS.attach_small(fam, [("a", _hits([g2], [0], [1e-6])), ("b", _hits([12], [g2], [1e-8]))], L, cfg)
    assert new[g2] == 1 and st["by_source"] == {"a": 0, "b": 1}
    # deterministic
    n1 = RS.attach_small(fam, hits, L, cfg)[0]
    assert (n1 == RS.attach_small(fam, hits, L, cfg)[0]).all()


def test_select_units():
    fam = np.array([0] * 12 + [1, 1])
    rep = np.array([0] * 6 + [6] * 6 + [12, 12])
    L = np.arange(14) + 100
    q, t, sizes = RS.select_units(fam, rep, L, small_max=10, reps_per_family=1)
    assert list(q) == [12]                       # one unit (stage-1 rep) for the small family
    assert len(t) == 1 and t[0] in (0, 6)        # one target rep: largest cluster (tie -> lower id)
    q, t, _ = RS.select_units(fam, rep, L, small_max=10, reps_per_family=20)
    assert sorted(t) == [0, 6]


def test_stage1_rep_of_canonical(tmp_path):
    seqs = ["MKV", "MKVLA", "MKVLA", "MA"]
    gt = GeneTable(["b", "a", "c", "d"], ["S"], np.zeros(4, np.int32), np.array([3, 5, 5, 2], np.int32), seqs)
    # canonical order: (len desc, seq, id): a(1), c(2), b(0), d(3) -> c0000000=a ...
    s0 = tmp_path / "s0"
    s0.mkdir()
    (s0 / "stage1_clusters.tsv").write_text("c0000000\tc0000000\nc0000000\tc0000001\nc0000000\tc0000002\n"
                                            "c0000003\tc0000003\n")
    rep = RS.stage1_rep_of(gt, str(s0))
    assert list(rep) == [1, 1, 1, 3]
    assert list(RS.stage1_rep_of(gt, None)) == [0, 1, 2, 3]


def test_revert_families_drops_cross_pairs():
    fam_old = np.array([0, 0, 0, 1, 2])
    fam_new = np.array([0, 0, 0, 0, 0])
    st = ScoreTable(q=np.array([0, 0, 3, 3], np.int32), t=np.array([1, 3, 4, 2], np.int32),
                    bits=np.ones(4, np.float32), qcov=np.ones(4, np.float32), tcov=np.ones(4, np.float32),
                    pident=np.ones(4, np.float32), evalue=np.ones(4, np.float32),
                    self_bits=np.ones(5, np.float32))
    fam, st2 = RS.revert_families(fam_old, fam_new, st, np.array([1]))
    assert list(fam) == [0, 0, 0, 1, 0]
    assert sorted(zip(st2.q.tolist(), st2.t.tolist())) == [(0, 1)]


def test_guard_detach():
    tree = load_guide_tree("((A,B)AB,C)root;", ["A", "B", "C"])
    names = ["a1", "b1", "c1", "b2", "c2"]
    sp = np.array([0, 1, 2, 1, 2], np.int32)
    gt = GeneTable(names, ["A", "B", "C"], sp, np.full(5, 100, np.int32), None)
    fam_old = np.array([0, 0, 0, 1, 2])                # target fam 0 (a1 b1 c1); attached b2 (fam1), c2 (fam2)
    fam_new = np.array([0, 0, 0, 0, 0])
    # scores: c2 is the RBH of a1 in C? no: a1-c1 0.9 > a1-c2 0.5 ; b2 RBH with c1 (b1-c1 lower)
    pairs = {(0, 1): 90, (0, 2): 90, (1, 2): 40, (0, 3): 50, (2, 3): 80, (0, 4): 50, (1, 4): 30}
    q, t = np.array([k[0] for k in pairs], np.int32), np.array([k[1] for k in pairs], np.int32)
    st = ScoreTable(q=q, t=t, bits=np.array(list(pairs.values()), np.float32), qcov=np.ones(len(q), np.float32),
                    tcov=np.ones(len(q), np.float32), pident=np.ones(len(q), np.float32),
                    evalue=np.ones(len(q), np.float32), self_bits=np.full(5, 100, np.float32))
    # merge result: b2 orthologous to c1 (O at root), c2 alone in a paralog group with c1
    res = MergeResult(roots=[("O", "root", [("O", "AB", [0, 1]), ("P", "C", [2, 4]), 3])],
                      node_stats=[], n_genes=5)
    det, info = RS.guard_detach(res, gt, tree, fam_old, fam_new, st)
    # b2-c1 predicted orthologs and RBH (c1's best in B is b2: 80 > 40) -> fam1 kept
    # c2: predicted ortholog of a1/b1/b2 (root O group) but a1's best in C is c1 -> not RBH; b1-c2 30 < b1-c1 40
    assert list(det) == [2]
    assert info["n_attached"] == 2 and info["n_confirmed"] == 1


def _rand(rng, n):
    return "".join(rng.choice(AA, size=n))


def _mutate(rng, s, p):
    return "".join(c if rng.random() > p else rng.choice(AA) for c in s)


@pytest.mark.skipif(shutil.which("mmseqs") is None or shutil.which("diamond") is None,
                    reason="mmseqs/diamond not installed")
def test_profile_rescue_end_to_end(tmp_path):
    """A remote homolog singleton of a 12-gene family is attached by the profile search and
    re-scored; an unrelated random singleton is not."""
    rng = np.random.default_rng(11)
    anc = _rand(rng, 250)
    seqs = [_mutate(rng, anc, 0.25) for _ in range(12)]
    seqs += [_mutate(rng, anc, 0.55), _rand(rng, 250)]
    seqs += [_rand(rng, 200)] * 0
    n = len(seqs)
    gt = GeneTable([f"g{i}" for i in range(n)], ["S1", "S2", "S3"], np.array([i % 3 for i in range(n)], np.int32),
                   np.array([len(s) for s in seqs], np.int32), seqs)
    fam = np.array([0] * 12 + [1, 2])
    st = ScoreTable(q=np.array([0], np.int32), t=np.array([1], np.int32), bits=np.array([300], np.float32),
                    qcov=np.ones(1, np.float32), tcov=np.ones(1, np.float32), pident=np.ones(1, np.float32),
                    evalue=np.ones(1, np.float32), self_bits=np.full(n, np.nan, np.float32))
    cfg = RS.RescueConfig(mode="profile", threads=2)
    f2, st2, att, info = RS.run_rescue(gt, fam, st, cfg, str(tmp_path / "w"), hits_cache=str(tmp_path / "c"),
                                       input_key="toy")
    assert f2[12] == 0 and f2[13] == 2
    assert info["attach"]["n_attached"] == 1
    new = (st2.q == 12) | (st2.t == 12)
    assert new.sum() >= 1 and np.isfinite(st2.self_bits[12])
    assert len(os.listdir(tmp_path / "c")) == 1
    # second call reuses the cached hits and gives the same result
    f3, _, _, info3 = RS.run_rescue(gt, fam, st, cfg, str(tmp_path / "w2"), hits_cache=str(tmp_path / "c"),
                                    input_key="toy")
    assert (f3 == f2).all() and "profile_from_cache" in info3
    # mode none is the identity
    f4, st4, _, _ = RS.run_rescue(gt, fam, st, RS.RescueConfig(), str(tmp_path / "w3"))
    assert f4 is fam and st4 is st


def test_pipeline_default_has_no_rescue_option_effect():
    from ancestra.pipeline import build_parser
    a = build_parser().parse_args(["--proteomes", "p", "--tree", "t", "--out", "o"])
    assert a.family_rescue == "none" and not a.family_rescue_guard and a.structure_3di is None


@pytest.mark.skipif(shutil.which("mmseqs") is None or shutil.which("diamond") is None,
                    reason="mmseqs/diamond not installed")
def test_pipeline_family_rescue_stage(sim_basic, tmp_path):
    """--family-rescue adds a 'rescue' stage after scores; 'none' leaves stages and keys as before;
    the guard reruns the merge only."""
    import json
    from ancestra.pipeline import run_pipeline
    _, simdir = sim_basic
    out = str(tmp_path / "run")
    kw = dict(families=os.path.join(simdir, "truth", "gene_families.tsv"), fam_cols=(0, 2), threads=2,
              log=lambda s: None, preset="fast")
    args = (os.path.join(simdir, "proteomes"), os.path.join(simdir, "species_tree.nwk"), out)
    m0 = run_pipeline(*args, **kw)
    assert set(m0["stages"]) == {"input", "families", "scores", "merge", "output"}
    assert "family_rescue" not in m0["config"]
    m1 = run_pipeline(*args, family_rescue="profile", **kw)
    assert set(m1["stages"]) == {"input", "families", "scores", "rescue", "merge", "output"}
    assert m1["stages"]["scores"]["skipped"] and not m1["stages"]["rescue"]["skipped"]
    assert m1["stages"]["merge"]["key"] != m0["stages"]["merge"]["key"]
    assert os.path.exists(os.path.join(out, "stages", "rescue", "families.tsv"))
    m2 = run_pipeline(*args, family_rescue="profile", rescue_guard=True, **kw)
    assert m2["stages"]["rescue"]["skipped"] and not m2["stages"]["merge"]["skipped"]
    assert "rescue_guard" in json.load(open(os.path.join(out, "stages", "merge", "stage.json")))["stats"]
    m3 = run_pipeline(*args, family_rescue="none", **kw)        # back to default: same merge key
    assert m3["stages"]["merge"]["key"] == m0["stages"]["merge"]["key"]
    with pytest.raises(ValueError):
        run_pipeline(*args, family_rescue="3di", **kw)            # 3di needs --structure-3di
