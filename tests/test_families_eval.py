"""Tests for ancestra_eval.families_eval: RefOG reading and metrics on hand-made data."""
import pytest

from ancestra_eval import families_eval as E


def test_evaluate_refogs_hand_example():
    # families: F1={a1,a2,a3,b1}, F2={a4}, F3={b2,b3}, F4={x}; c1 missing from partition
    g2f = {"a1": "F1", "a2": "F1", "a3": "F1", "b1": "F1", "a4": "F2",
           "b2": "F3", "b3": "F3", "x": "F4"}
    refogs = {"A": {"a1", "a2", "a3", "a4"}, "B": {"b1", "b2", "b3"}, "C": {"x", "c1"}}
    r = E.evaluate_refogs(g2f, refogs)
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
    assert agg["gene_weighted_best_frac"] == pytest.approx((3 + 2 + 1) / 9)
    assert agg["frac_refogs_best_frac_ge95"] == 0.0
    assert agg["frac_refogs_le1_stray"] == pytest.approx(1.0)  # A: 1 stray, B: 1 stray, C: 1 stray
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
    r = E.evaluate_refogs(perfect, refogs)["aggregate"]
    assert r["completeness"] == 1.0 and r["pair_recall"] == 1.0 and r["pair_precision"] == 1.0
    assert r["fusion_rate"] == 0.0 and r["mean_n_fam"] == 1.0
    one = {g: "F1" for g in perfect}  # negative control: everything fused
    r = E.evaluate_refogs(one, refogs)["aggregate"]
    assert r["completeness"] == 1.0 and r["fusion_rate"] == 1.0 and r["pair_precision"] < 0.5
    singles = {g: g for g in perfect}  # negative control: everything split
    r = E.evaluate_refogs(singles, refogs)["aggregate"]
    assert r["completeness"] == 0.0 and r["pair_recall"] == 0.0


def test_read_refogs_removes_uncertain(tmp_path):
    d = tmp_path / "RefOGs"
    (d / "low_certainty_assignments").mkdir(parents=True)
    (d / "RefOG001.txt").write_text("a\nb\nc\n")
    (d / "RefOG002.txt").write_text("d\n")
    (d / "low_certainty_assignments" / "RefOG001.txt").write_text("c\n")
    r = E.read_refogs(str(d))
    assert r == {"RefOG001": {"a", "b"}, "RefOG002": {"d"}}
    assert E.read_refogs(str(d), ["RefOG001"], remove_uncertain=False) == {"RefOG001": {"a", "b", "c"}}
