"""Tests for ancestra.evaluate_toy: positive and negative controls, partition metrics."""
import os

import numpy as np

from ancestra.evaluate_toy import evaluate, main, partition_metrics
from ancestra.orthoxml import Group, OrthoXML, read_orthoxml


def _truth(out):
    return read_orthoxml(os.path.join(out, "truth", "hogs.orthoxml"))


def test_positive_control_truth_vs_truth(sim_complex, tmp_path):
    _, out = sim_complex
    t = os.path.join(out, "truth", "hogs.orthoxml")
    res = evaluate(t, t)
    pm = res["pairs"]
    assert pm["precision"] == pm["recall"] == pm["f1"] == 1.0
    for lvl, m in res["hogs"].items():
        assert m["exact"] == 1.0 and m["mean_best_jaccard"] == 1.0, lvl
        assert m["fissions"] == 0 and m["fusions"] == 0 and m["fusions_any"] == 0, lvl
        assert m["pair_f1"] == 1.0 or m["n_truth_hogs"] == 0, lvl
    js = tmp_path / "r.json"
    assert main(["--truth", t, "--pred", t, "--json", str(js)]) == 0
    assert js.exists()


def _relabel(group, perm):
    if isinstance(group, Group):
        return Group(group.kind, [_relabel(c, perm) for c in group.children], group.level, group.id)
    return perm[group]


def test_negative_control_random_groupings(sim_basic):
    """Same group structure, genes randomly permuted (within species, so the species
    composition of every group is preserved): pair precision and HOG recovery collapse."""
    _, out = sim_basic
    truth = _truth(out)
    rng = np.random.default_rng(11)
    perm = {}
    for sp in truth.species:
        ids = [g for g, e in truth.genes.items() if e.species == sp]
        for a, b in zip(ids, rng.permutation(ids)):
            perm[a] = b
    rand = OrthoXML(genes=truth.genes, groups=[_relabel(g, perm) for g in truth.groups],
                    species_tree=truth.species_tree)
    res = evaluate(rand, truth)
    assert res["pairs"]["precision"] < 0.1 and res["pairs"]["f1"] < 0.1
    assert res["hogs"]["families"]["exact"] < 0.1
    for lvl in truth.tree_index().internal:
        assert res["hogs"][lvl]["pair_f1"] < 0.2, lvl


def test_negative_control_all_singletons(sim_basic):
    _, out = sim_basic
    truth = _truth(out)
    empty = OrthoXML(genes=truth.genes, groups=[])
    res = evaluate(empty, truth)
    assert res["pairs"]["recall"] == 0.0 and res["pairs"]["n_pred"] == 0
    assert res["hogs"]["families"]["exact"] == 0.0


def test_flat_families_lose_precision_not_recall(sim_basic):
    """Families as flat orthologGroups: every true ortholog pair is found, but true
    paralogs become false orthologs (precision < 1 when duplications survived)."""
    _, out = sim_basic
    truth = _truth(out)
    flat = OrthoXML(genes=truth.genes, groups=[Group("ortholog", g.genes()) for g in truth.groups])
    res = evaluate(flat, truth)
    pm = res["pairs"]
    assert pm["recall"] == 1.0
    assert pm["precision"] < 1.0 and pm["fp_true_paralogs"] == pm["n_pred"] - pm["tp"]
    assert res["hogs"]["families"]["exact"] == 1.0


def test_partition_metrics_fusion_fission():
    truth = [{"a1", "a2", "a3"}, {"b1", "b2"}, {"c1"}]
    pred = [{"a1", "a2"}, {"a3", "b1", "b2"}, {"c1"}]
    m = partition_metrics(truth, pred)
    assert m["fissions"] == 1                     # A split
    assert m["fusions"] == 1                      # A and B merged in one predicted HOG
    assert m["n_exact"] == 0
    assert np.isclose(m["mean_best_jaccard"], (2 / 3 + 2 / 3) / 2)
    # pairs: truth {a1a2,a1a3,a2a3,b1b2}=4; pred {a1a2, a3b1,a3b2,b1b2}=4; tp=2
    assert np.isclose(m["pair_precision"], 0.5) and np.isclose(m["pair_recall"], 0.5)
    # genes left out by the prediction count as singletons
    m2 = partition_metrics(truth, [])
    assert m2["pair_recall"] == 0.0 and m2["fissions"] == 2


def test_unknown_and_missing_genes(sim_basic):
    _, out = sim_basic
    truth = _truth(out)
    from ancestra.orthoxml import GeneEntry
    genes = dict(truth.genes)
    genes["alien1"] = GeneEntry("alien1", "sp0")
    groups = list(truth.groups) + [Group("ortholog", ["alien1", truth.groups[0].genes()[0]])]
    # drop one family entirely (its genes are missing from the prediction)
    dropped = truth.groups[1].genes()
    for g in dropped:
        del genes[g]
    groups = [groups[0]] + groups[2:]
    res = evaluate(OrthoXML(genes=genes, groups=groups), truth)
    assert res["n_pred_unknown_genes"] == 1
    assert res["n_truth_genes_missing_in_pred"] == len(dropped)
    assert res["pairs"]["precision"] == 1.0 and res["pairs"]["recall"] < 1.0
