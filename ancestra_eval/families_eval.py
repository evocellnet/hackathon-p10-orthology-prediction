"""Evaluation of S0 family partitions against OrthoBench RefOGs (split subset dev or test).

Moved out of ``ancestra/families.py`` so the method package holds no benchmark references.
See docs/S0_families.md for the metric definitions.

CLI
---
    python -m ancestra_eval.families_eval evaluate --families D/families.tsv --subset dev
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from typing import Iterable

from ancestra.families import read_families_tsv

__all__ = ["read_refogs", "evaluate_refogs"]

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_REFOGS = os.path.join(_ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs")
_SPLIT = os.path.join(_ROOT, "splits/orthobench_split.json")


def read_refogs(refog_dir: str, ids: Iterable[str] | None = None,
                remove_uncertain: bool = True) -> dict[str, set[str]]:
    """Read OrthoBench RefOG gene sets (optionally minus low_certainty_assignments)."""
    if ids is None:
        ids = sorted(f[:-4] for f in os.listdir(refog_dir) if f.startswith("RefOG") and f.endswith(".txt"))
    out = {}
    for rid in ids:
        genes = {l.strip() for l in open(os.path.join(refog_dir, rid + ".txt")) if l.strip()}
        unc = os.path.join(refog_dir, "low_certainty_assignments", rid + ".txt")
        if remove_uncertain and os.path.exists(unc):
            genes -= {l.strip() for l in open(unc) if l.strip()}
        out[rid] = genes
    return out


def evaluate_refogs(gene2fam: dict[str, str], refogs: dict[str, set[str]]) -> dict:
    """Score a family partition against reference groups.

    Per RefOG R (genes missing from ``gene2fam`` count as singleton families):
      complete     -- all of R lies in one family
      n_fam        -- number of distinct families R is spread over (fission; 1 = none)
                      (secondary: gene-weighted best_frac, share of RefOGs with best_frac >= 0.95,
                      share of RefOGs with at most one gene outside the best family)
      best_frac    -- largest fraction of R in a single family
      best_fam_size-- size of that family; size_ratio = best_fam_size / |R|
      foreign_refog_frac -- fraction of the best family's genes that belong to *other*
                       RefOGs in ``refogs`` (only the RefOGs passed in -- dev only during tuning)
      n_refogs_in_best -- number of RefOGs (of those passed) touching the best family
    Aggregate: completeness (= mean complete), mean / median fission, mean best_frac,
    median / max size_ratio, fusion rate (fraction of RefOGs whose best family also
    contains genes of another RefOG), and the OrthoBench-style pairwise
    recall / precision (per-RefOG normalised by |R|-1, as benchmark.py q_even).
    """
    fam_members = defaultdict(set)
    for g, f in gene2fam.items():
        fam_members[f].add(g)
    gene2ref = {}
    for rid, genes in refogs.items():
        for g in genes:
            gene2ref[g] = rid
    per = {}
    tp = fp = fn = 0.0
    for rid, R in refogs.items():
        if not R:
            continue
        cnt = Counter(gene2fam.get(g, "__single__" + g) for g in R)
        best_f, best_n = max(cnt.items(), key=lambda kv: (kv[1], kv[0]))
        fam = fam_members.get(best_f, {best_f[len("__single__"):]} if best_f.startswith("__single__") else set())
        others = Counter(gene2ref[g] for g in fam if g in gene2ref and gene2ref[g] != rid)
        n = len(R)
        # OrthoBench-style pair counts
        t_tp = t_fp = t_fn = 0.0
        for f, ov in cnt.items():
            fs = len(fam_members[f]) if not f.startswith("__single__") else 1
            t_tp += ov * (ov - 1) / 2
            t_fp += ov * (fs - ov)
            t_fn += (n - ov) * ov / 2
        if n > 1:
            tp += t_tp / (n - 1); fp += t_fp / (n - 1); fn += t_fn / (n - 1)
        per[rid] = {
            "size": n,
            "complete": len(cnt) == 1,
            "n_fam": len(cnt),
            "best_frac": best_n / n,
            "best_fam_size": len(fam),
            "size_ratio": len(fam) / n,
            "foreign_refog_frac": sum(others.values()) / len(fam) if fam else 0.0,
            "n_refogs_in_best": 1 + len(others),
        }
    vals = list(per.values())
    k = len(vals)

    def med(xs):
        s = sorted(xs)
        return s[len(s) // 2] if len(s) % 2 else 0.5 * (s[len(s) // 2 - 1] + s[len(s) // 2])

    agg = {
        "n_refogs": k,
        "completeness": sum(v["complete"] for v in vals) / k,
        "n_complete": sum(v["complete"] for v in vals),
        "mean_n_fam": sum(v["n_fam"] for v in vals) / k,
        "median_n_fam": med([v["n_fam"] for v in vals]),
        "max_n_fam": max(v["n_fam"] for v in vals),
        "mean_best_frac": sum(v["best_frac"] for v in vals) / k,
        "gene_weighted_best_frac": sum(v["best_frac"] * v["size"] for v in vals) / sum(v["size"] for v in vals),
        "frac_refogs_best_frac_ge95": sum(v["best_frac"] >= 0.95 for v in vals) / k,
        "frac_refogs_le1_stray": sum(v["size"] - round(v["best_frac"] * v["size"]) <= 1 for v in vals) / k,
        "median_size_ratio": med([v["size_ratio"] for v in vals]),
        "mean_size_ratio": sum(v["size_ratio"] for v in vals) / k,
        "max_size_ratio": max(v["size_ratio"] for v in vals),
        "fusion_rate": sum(v["n_refogs_in_best"] > 1 for v in vals) / k,
        "mean_foreign_refog_frac": sum(v["foreign_refog_frac"] for v in vals) / k,
        "pair_recall": tp / (tp + fn) if tp + fn else 0.0,
        "pair_precision": tp / (tp + fp) if tp + fp else 0.0,
    }
    return {"aggregate": agg, "per_refog": per}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ancestra_eval.families_eval", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("evaluate", help="score families.tsv vs OrthoBench RefOGs of one split subset")
    e.add_argument("--families", required=True)
    e.add_argument("--subset", choices=["dev", "test"], required=True)
    e.add_argument("--split", default=_SPLIT)
    e.add_argument("--refogs", default=_REFOGS)
    e.add_argument("--keep-uncertain", action="store_true")
    e.add_argument("--out")
    a = ap.parse_args(argv)
    split = json.load(open(a.split))
    refogs = read_refogs(a.refogs, split[a.subset], remove_uncertain=not a.keep_uncertain)
    res = evaluate_refogs(read_families_tsv(a.families), refogs)
    res["subset"] = a.subset
    res["remove_uncertain"] = not a.keep_uncertain
    s = json.dumps(res, indent=1)
    if a.out:
        open(a.out, "w").write(s)
    print(json.dumps(res["aggregate"], indent=1))


if __name__ == "__main__":
    sys.exit(main())
