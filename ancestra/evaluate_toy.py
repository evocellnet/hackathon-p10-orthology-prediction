"""Compare predicted HOGs (OrthoXML) with simulated truth.

Usage::

    python -m ancestra.evaluate_toy --truth SIM/truth/hogs.orthoxml --pred pred.orthoxml \\
        [--species-tree SIM/species_tree.nwk] [--json out.json]

Metrics
-------
**Pairwise orthology** (QfO flattening, :func:`ancestra.orthoxml.hogs_to_pairs`):
precision, recall and F1 of predicted ortholog pairs against true pairs, and
recall broken down by the level of the true pair (the speciation node of its
LCA). Only genes present in the truth are scored; predicted genes unknown to
the truth are reported and ignored.

**HOGs per taxonomic level** ``L`` (every internal node of the species tree,
plus ``"families"`` = top-level groups): truth and predicted HOGs at ``L`` are
computed with :func:`ancestra.orthoxml.hogs_at_level` (each gene of the clade
below ``L`` belongs to exactly one HOG; genes absent from any predicted group
are singletons). For truth HOGs with >= 2 genes:

* ``exact``      -- fraction exactly reproduced by a predicted HOG;
* ``mean_best_jaccard`` -- mean over truth HOGs of max_P |T n P| / |T u P|;
* ``fissions``   -- number of truth HOGs whose genes fall in > 1 predicted HOG;
* ``fusions``    -- number of predicted HOGs containing genes of >= 2 truth
  HOGs, where at least two of them are multi-gene truth HOGs
  (``fusions_any`` also counts absorbed truth singletons);
* ``pair_precision/recall/f1`` -- co-membership pairs (two genes in the same
  HOG at ``L``), i.e. the OrthoBench-style pairwise measure at that level.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from typing import Union

from .newick import Node, TreeIndex, read_newick
from .orthoxml import OrthoXML, hogs_at_level, hogs_to_pairs, read_orthoxml

__all__ = ["pair_metrics", "partition_metrics", "evaluate", "main"]


def _prf(tp: int, n_pred: int, n_true: int) -> dict:
    p = tp / n_pred if n_pred else 0.0
    r = tp / n_true if n_true else 0.0
    f = 2 * p * r / (p + r) if p + r else 0.0
    return dict(tp=tp, n_pred=n_pred, n_true=n_true, precision=p, recall=r, f1=f)


def pair_metrics(pred_pairs, true_pairs) -> dict:
    """Precision/recall/F1 of two sets of unordered (sorted) gene pairs."""
    pred, true = set(pred_pairs), set(true_pairs)
    return _prf(len(pred & true), len(pred), len(true))


def _n_pairs(sizes) -> int:
    return sum(s * (s - 1) // 2 for s in sizes)


def partition_metrics(truth_groups: list[set], pred_groups: list[set]) -> dict:
    """Compare two partitions of the same gene set (see module docstring)."""
    t_of: dict[str, int] = {}
    for i, s in enumerate(truth_groups):
        for g in s:
            t_of[g] = i
    universe = set(t_of)
    p_of: dict[str, int] = {}
    pred = []
    for s in pred_groups:
        s = s & universe
        if s:
            pred.append(s)
    for j, s in enumerate(pred):
        for g in s:
            p_of[g] = j
    nxt = len(pred)
    for g in universe:            # genes the prediction left out -> singletons
        if g not in p_of:
            p_of[g] = nxt
            pred.append({g})
            nxt += 1
    # contingency
    cont: Counter = Counter((t_of[g], p_of[g]) for g in universe)
    tp_pairs = sum(n * (n - 1) // 2 for n in cont.values())
    n_true_pairs = _n_pairs(len(s) for s in truth_groups)
    n_pred_pairs = _n_pairs(len(s) for s in pred)
    by_t: dict[int, list[tuple[int, int]]] = defaultdict(list)
    by_p: dict[int, list[int]] = defaultdict(list)
    for (ti, pj), n in cont.items():
        by_t[ti].append((pj, n))
        by_p[pj].append(ti)
    multi = [i for i, s in enumerate(truth_groups) if len(s) >= 2]
    exact = 0
    jac_sum = 0.0
    fissions = 0
    for ti in multi:
        tsize = len(truth_groups[ti])
        best = 0.0
        for pj, n in by_t[ti]:
            j = n / (tsize + len(pred[pj]) - n)
            best = max(best, j)
        jac_sum += best
        if best == 1.0:
            exact += 1
        if len(by_t[ti]) > 1:
            fissions += 1
    fusions = fusions_any = 0
    for pj, tis in by_p.items():
        if len(tis) >= 2:
            fusions_any += 1
            if sum(1 for ti in tis if len(truth_groups[ti]) >= 2) >= 2:
                fusions += 1
    nm = len(multi)
    out = dict(n_genes=len(universe), n_truth_hogs=nm,
               n_pred_hogs=sum(1 for s in pred if len(s) >= 2),
               exact=exact / nm if nm else 1.0, n_exact=exact,
               mean_best_jaccard=jac_sum / nm if nm else 1.0,
               fissions=fissions, fusions=fusions, fusions_any=fusions_any)
    prf = _prf(tp_pairs, n_pred_pairs, n_true_pairs)
    out.update(pair_precision=prf["precision"], pair_recall=prf["recall"], pair_f1=prf["f1"])
    return out


def evaluate(pred: Union[OrthoXML, str], truth: Union[OrthoXML, str],
             species_tree: Union[Node, str, None] = None, levels: list[str] | None = None) -> dict:
    """Full evaluation of a predicted OrthoXML against the truth OrthoXML."""
    pred = pred if isinstance(pred, OrthoXML) else read_orthoxml(pred)
    truth = truth if isinstance(truth, OrthoXML) else read_orthoxml(truth)
    if isinstance(species_tree, str):
        species_tree = read_newick(species_tree)
    tree = species_tree or truth.species_tree
    if tree is None:
        raise ValueError("a species tree is required (argument or taxonomy in the truth file)")
    tix = TreeIndex(tree)
    pred = OrthoXML(genes=dict(pred.genes), groups=list(pred.groups), species_tree=pred.species_tree)
    universe = set(truth.genes)
    unknown = [g for g in pred.genes if g not in universe]
    missing = [g for g in universe if g not in pred.genes]
    # restrict prediction to known genes (species from the truth)
    for g in unknown:
        del pred.genes[g]
    for g in missing:
        pred.genes[g] = truth.genes[g]
    _prune_unknown(pred, set(unknown))

    true_pairs = hogs_to_pairs(truth, with_level=True)
    _, true_paralogs = hogs_to_pairs(truth, with_paralogs=True)
    pred_pairs = hogs_to_pairs(pred)
    res: dict = {"n_genes": len(universe), "n_pred_unknown_genes": len(unknown),
                 "n_truth_genes_missing_in_pred": len(missing)}
    pm = pair_metrics(pred_pairs, true_pairs.keys())
    # false-positive breakdown: true paralogs vs. pairs from different truth families
    fam_of = {g: i for i, s in enumerate(hogs_at_level(truth, None)) for g in s}
    fp = [q for q in pred_pairs if q not in true_pairs]
    pm["fp_true_paralogs"] = sum(1 for q in fp if q in true_paralogs)
    pm["fp_different_families"] = sum(1 for q in fp if fam_of[q[0]] != fam_of[q[1]])
    pm["fp_other"] = len(fp) - pm["fp_true_paralogs"] - pm["fp_different_families"]
    res["pairs"] = pm
    by_level: dict = defaultdict(lambda: [0, 0])
    for p, lvl in true_pairs.items():
        by_level[lvl][1] += 1
        if p in pred_pairs:
            by_level[lvl][0] += 1
    res["pair_recall_by_level"] = {lvl: dict(tp=v[0], n_true=v[1], recall=v[0] / v[1])
                                   for lvl, v in sorted(by_level.items(), key=lambda kv: str(kv[0]))}
    levels = levels or tix.internal
    hogs = {}
    fam_t = hogs_at_level(truth, None)
    fam_p = hogs_at_level(pred, None)
    hogs["families"] = partition_metrics(fam_t, fam_p)
    for lvl in levels:
        t = hogs_at_level(truth, lvl, species_tree=tree)
        p = hogs_at_level(pred, lvl, species_tree=tree)
        m = partition_metrics(t, p)
        m["depth"] = tix.depth[lvl]
        m["n_species"] = len(tix.clade(lvl))
        hogs[lvl] = m
    res["hogs"] = hogs
    return res


def _prune_unknown(ox: OrthoXML, unknown: set) -> None:
    if not unknown:
        return
    from .orthoxml import Group, normalize_groups

    def prune(x):  # builds new Group objects; the input is not modified
        if isinstance(x, Group):
            kids = [c for c in (prune(c) for c in x.children) if c is not None]
            return Group(x.kind, kids, x.level, x.id, dict(x.properties))
        return None if x in unknown else x

    ox.groups = normalize_groups(prune(g) for g in ox.groups)


def format_report(res: dict) -> str:
    pm = res["pairs"]
    lines = [f"genes={res['n_genes']}  pred_unknown={res['n_pred_unknown_genes']}  "
             f"missing_in_pred={res['n_truth_genes_missing_in_pred']}",
             f"ortholog pairs: P={pm['precision']:.4f} R={pm['recall']:.4f} F1={pm['f1']:.4f} "
             f"(tp={pm['tp']} pred={pm['n_pred']} true={pm['n_true']}; FP: paralogs="
             f"{pm['fp_true_paralogs']} cross-family={pm['fp_different_families']} other={pm['fp_other']})",
             "level\tdepth\tnsp\ttruthHOGs\tpredHOGs\texact\tmeanJacc\tfissions\tfusions\tpairP\tpairR\tpairF1"]
    for lvl, m in res["hogs"].items():
        lines.append(f"{lvl}\t{m.get('depth', '-')}\t{m.get('n_species', '-')}\t{m['n_truth_hogs']}\t"
                     f"{m['n_pred_hogs']}\t{m['exact']:.3f}\t{m['mean_best_jaccard']:.3f}\t"
                     f"{m['fissions']}\t{m['fusions']}\t{m['pair_precision']:.3f}\t"
                     f"{m['pair_recall']:.3f}\t{m['pair_f1']:.3f}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m ancestra.evaluate_toy", description=__doc__.split("\n")[0])
    ap.add_argument("--truth", required=True, help="truth OrthoXML (simulator truth/hogs.orthoxml)")
    ap.add_argument("--pred", required=True, help="predicted OrthoXML")
    ap.add_argument("--species-tree", help="Newick with internal node names (default: taxonomy in truth)")
    ap.add_argument("--json", help="write full results as JSON")
    a = ap.parse_args(argv)
    res = evaluate(a.pred, a.truth, a.species_tree)
    print(format_report(res))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(res, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
