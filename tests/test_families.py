"""Tests for ancestra.families: parsing, partitioning, stats and RefOG metrics on hand-made data."""
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


def test_evaluate_refogs_hand_example():
    # families: F1={a1,a2,a3,b1}, F2={a4}, F3={b2,b3}, F4={x}; c1 missing from partition
    g2f = {"a1": "F1", "a2": "F1", "a3": "F1", "b1": "F1", "a4": "F2",
           "b2": "F3", "b3": "F3", "x": "F4"}
    refogs = {"A": {"a1", "a2", "a3", "a4"}, "B": {"b1", "b2", "b3"}, "C": {"x", "c1"}}
    r = F.evaluate_refogs(g2f, refogs)
    pa, pb, pc = r["per_refog"]["A"], r["per_refog"]["B"], r["per_refog"]["C"]
    assert pa["complete"] is False and pa["n_fam"] == 2 and pa["best_frac"] == 0.75
    assert pa["best_fam_size"] == 4 and pa["size_ratio"] == 1.0
    assert pa["foreign_refog_frac"] == pytest.approx(0.25) and pa["n_refogs_in_best"] == 2
    assert pb["n_fam"] == 2 and pb["best_frac"] == pytest.approx(2 / 3)
    assert pb["best_fam_size"] == 2 and pb["foreign_refog_frac"] == 0.0
    # c1 is absent from the partition -> its own singleton; tie (1,1) broken deterministically
    assert pc["n_fam"] == 2 and pc["complete"] is False
    agg = r["aggregate"]
    assert agg["n_refogs"] == 3 and agg["completeness"] == 0.0
    assert agg["mean_n_fam"] == 2.0
    assert agg["fusion_rate"] == pytest.approx(1 / 3)
    # OrthoBench q_even pair recall: A: TP=3 (3 of 6 pairs), FN=3 -> /3 ; B: TP=1, FN=2 -> /2; C: TP=0, FN=1 -> /1
    tp = 3 / 3 + 1 / 2 + 0
    fn = 3 / 3 + 2 / 2 + 1
    assert agg["pair_recall"] == pytest.approx(tp / (tp + fn))
    # FP: A: 3*1 (b1) + 1*0 = 3 -> /3; B: b1 in F1 with 3 others -> 1*3=3, F3: 0 -> /2; C: 0
    fp = 3 / 3 + 3 / 2
    assert agg["pair_precision"] == pytest.approx(tp / (tp + fp))


def test_evaluate_perfect_and_one_family_controls():
    refogs = {"A": {"a1", "a2"}, "B": {"b1", "b2", "b3"}}
    perfect = {"a1": "F1", "a2": "F1", "b1": "F2", "b2": "F2", "b3": "F2", "z": "F3"}
    r = F.evaluate_refogs(perfect, refogs)["aggregate"]
    assert r["completeness"] == 1.0 and r["pair_recall"] == 1.0 and r["pair_precision"] == 1.0
    assert r["fusion_rate"] == 0.0 and r["mean_n_fam"] == 1.0
    one = {g: "F1" for g in perfect}  # negative control: everything fused
    r = F.evaluate_refogs(one, refogs)["aggregate"]
    assert r["completeness"] == 1.0 and r["fusion_rate"] == 1.0 and r["pair_precision"] < 0.5
    singles = {g: g for g in perfect}  # negative control: everything split
    r = F.evaluate_refogs(singles, refogs)["aggregate"]
    assert r["completeness"] == 0.0 and r["pair_recall"] == 0.0


def test_read_refogs_removes_uncertain(tmp_path):
    d = tmp_path / "RefOGs"
    (d / "low_certainty_assignments").mkdir(parents=True)
    (d / "RefOG001.txt").write_text("a\nb\nc\n")
    (d / "RefOG002.txt").write_text("d\n")
    (d / "low_certainty_assignments" / "RefOG001.txt").write_text("c\n")
    r = F.read_refogs(str(d))
    assert r == {"RefOG001": {"a", "b"}, "RefOG002": {"d"}}
    assert F.read_refogs(str(d), ["RefOG001"], remove_uncertain=False) == {"RefOG001": {"a", "b", "c"}}


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
