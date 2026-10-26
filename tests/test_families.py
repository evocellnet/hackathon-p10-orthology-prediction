"""Tests for ancestra.families: parsing, partitioning and stats on hand-made data."""
import json
import os

import pytest

from ancestra import families as F


def _write(p, text):
    p.write_text(text)
    return str(p)


def test_read_fasta_ids(tmp_path):
    fa = _write(tmp_path / "x.fa", ">a desc\nMK\nLL\n>b\nM\n>c extra words\nMKV\n")
    assert F.read_fasta_ids(fa) == ["a", "b", "c"]


def test_read_pairs_tsv_skips_comments(tmp_path):
    t = _write(tmp_path / "c.tsv", "#rep\tmember\na\ta\na\tb\n\nc\tc\textra\n")
    assert F.read_pairs_tsv(t) == [("a", "a"), ("a", "b"), ("c", "c")]


def test_partition_and_labels():
    genes = ["g1", "g2", "g3", "g4", "g5", "g6"]
    pairs = [("g1", "g2"), ("g2", "g3"), ("g4", "g5")]
    g2r = F.partition_from_pairs(genes, pairs)
    lab = F.label_families(g2r)
    assert lab["g1"] == lab["g2"] == lab["g3"] == "F0000001"   # largest first
    assert lab["g4"] == lab["g5"] == "F0000002"
    assert lab["g6"] == "F0000003"                              # singleton kept


def test_partition_strict_unknown_gene():
    with pytest.raises(KeyError):
        F.partition_from_pairs(["a"], [("a", "zzz")])


def test_label_tie_break_deterministic():
    g2r = {"b": "b", "a": "a", "c": "a", "d": "b"}  # two families of size 2
    lab = F.label_families(g2r)
    assert lab["a"] == lab["c"] == "F0000001"   # min gene id 'a' < 'b'
    assert lab["b"] == lab["d"] == "F0000002"


def test_family_stats():
    g2f = {"a": "F1", "b": "F1", "c": "F1", "d": "F2", "e": "F3"}
    st = F.family_stats(g2f)
    assert st["n_genes"] == 5 and st["n_families"] == 3
    assert st["n_singletons"] == 2
    assert st["sum_g2"] == 9 + 1 + 1
    assert st["sum_g2_over_N2"] == pytest.approx(11 / 25)
    assert st["max_family"] == 3 and st["top_sizes"] == [3, 1, 1]


def test_families_tsv_roundtrip(tmp_path):
    g2f = {"x": "F0000001", "y": "F0000001", "z": "F0000002"}
    p = str(tmp_path / "f.tsv")
    F.write_families_tsv(g2f, p)
    assert F.read_families_tsv(p) == g2f


def test_filter_hits_coverage_modes():
    # q s evalue bits qlen slen qstart qend sstart send
    rows = [
        ["r1", "r2", "1e-30", "200", "100", "400", "1", "90", "1", "90"],   # qcov .9, scov .225
        ["r2", "r1", "1e-30", "210", "400", "100", "1", "90", "1", "90"],   # reverse, higher bits
        ["r1", "r1", "0", "500", "100", "100", "1", "100", "1", "100"],     # self hit dropped
        ["r1", "r3", "1e-2", "50", "100", "100", "1", "100", "1", "100"],   # e-value fails
    ]
    cfg = F.FamilyConfig(s2_evalue=1e-5, s2_cov=0.5, s2_cov_mode="min")
    assert F.filter_hits(rows, cfg) == []
    cfg.s2_cov_mode = "max"
    assert F.filter_hits(rows, cfg) == [("r1", "r2", 210.0)]
    cfg.s2_cov_mode = "shorter"
    assert F.filter_hits(rows, cfg) == [("r1", "r2", 210.0)]
    cfg.s2_cov = 0.0
    assert len(F.filter_hits(rows, cfg)) == 1


def test_config_hash_ignores_resources():
    a = F.FamilyConfig(name="x", threads=1)
    b = F.FamilyConfig(name="y", threads=24, memory_gb=8)
    assert a.hash() == b.hash()
    assert a.hash() != F.FamilyConfig(approx_id=50).hash()


def test_leiden_split_two_cliques():
    a = [f"a{i}" for i in range(6)]
    b = [f"b{i}" for i in range(6)]
    edges = [(x, y, 100.0) for grp in (a, b) for i, x in enumerate(grp) for y in grp[i + 1:]]
    edges.append(("a0", "b0", 1.0))  # weak bridge
    w = {n: 10 for n in a + b}       # each rep stands for 10 genes -> component of 120 genes
    parts = F.leiden_split(a + b, edges, w, gmax=60, seed=11)
    assert sorted(sorted(p) for p in parts) == [sorted(a), sorted(b)]
    # under the cap: untouched
    assert F.leiden_split(a + b, edges, w, gmax=1000) == [a + b]


def test_best_hits_directed_and_filtered():
    rows = [
        ["g1", "r1", "1e-30", "100", "100", "100", "1", "100", "1", "100"],
        ["g1", "r2", "1e-40", "150", "100", "100", "1", "100", "1", "100"],
        ["g1", "r3", "1e-50", "300", "100", "1000", "1", "100", "1", "100"],  # scov 0.1 -> fails 'min'
        ["g2", "g2", "0", "500", "100", "100", "1", "100", "1", "100"],      # self
        ["g3", "r1", "1e-2", "40", "100", "100", "1", "100", "1", "100"],    # e-value fails
    ]
    cfg = F.FamilyConfig(s2_evalue=1e-5, s2_cov=0.5, s2_cov_mode="min")
    assert F.best_hits(rows, cfg) == {"g1": ("r2", 150.0)}
    cfg.s2_cov = 0.0
    assert F.best_hits(rows, cfg) == {"g1": ("r3", 300.0)}
