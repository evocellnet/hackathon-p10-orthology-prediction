#!/usr/bin/env python3
"""TreeFam-A TPR/PPV restricted to a family subset (dev/test of splits/treefamA_split.json).

Uses the raw output of the official QfO RefPhyloTest (``raw/TreeFam-A/RefPhylo-TreeFam-A_*_raw.txt.gz``
written by bench_repo/scripts/run_qfo_local.sh): every reference pair of the test, labelled TP/FN
(reference speciation pair predicted / not predicted as orthologs) or FP/TN (reference duplication pair
predicted / not predicted). The reference pairs are therefore exactly the official ones.

Official formula (RefPhyloTest.drw; QfO 2022 TreeFam-A is ONE case): tp=fp=fn=tn=1 (uniform prior) plus
HALF the pair counts (each relation is stored under one protein only but added as length/2); TPR = tp/(tp+fn), PPV = tp/(tp+fp), CI half-width 1.96*sqrt(p(1-p)/n). The subset value
applies the same formula to the pairs of the subset's families.

Pooled check: ``all`` (= dev + test pooled) must equal the official assessment (assess/TreeFam-A.json).
Test-only values are refused unless ``--allow-test`` (DESIGN §4: test families are not inspected).

``--from-pairs``: RUN_DIR is a pipeline run dir; dev reference pairs are labelled directly from its
``pairs.tsv`` against the reference relations (no QfO benchmark run, dev families only). Validated: on arm A
(det_go) it gives the same dev TP/FP/FN/TN as the official raw output.

    python scripts/qfo_treefam_subset.py RUN_DIR [RUN_DIR ...] [--subset dev] [--json out.json]
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import math
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPLIT = os.path.join(ROOT, "splits/treefamA_split.json")


def load_split(path: str = SPLIT) -> tuple[dict, dict]:
    sp = json.load(open(path))
    gene2fam = {g: f for f, d in sp["families"].items() for g in d["members"]}
    fam2sub = {f: d["subset"] for f, d in sp["families"].items()}
    return gene2fam, fam2sub


def raw_pairs(run_dir: str):
    fs = glob.glob(os.path.join(run_dir, "raw", "TreeFam-A", "RefPhylo-TreeFam-A_*_raw.txt.gz"))
    if len(fs) != 1:
        raise SystemExit(f"{run_dir}: expected one TreeFam-A raw file, found {fs}")
    with gzip.open(fs[0], "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _, a, b, lab = line.rstrip("\n").split("\t")
            yield a, b, lab


def labels_from_pairs(pairs_tsv: str, gene2fam: dict, fam2sub: dict, subset: str = "dev"):
    """Label the reference pairs of ``subset`` families directly from a pipeline ``pairs.tsv`` (no QfO run):
    reference relation S/D (QfO 2022 ReconciledTrees_TreeFam-A.drw) x predicted ortholog yes/no ->
    TP/FN/FP/TN, exactly as RefPhyloTest labels them (validated against the official raw labels of
    arm A: identical dev counts). Yields (acc_a, acc_b, label)."""
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from make_treefamA_split import entry_to_acc, load_relations
    inv = entry_to_acc()
    ref = {}
    for x, y, e in load_relations():
        a, b = next(iter(inv[x])), next(iter(inv[y]))
        f = gene2fam.get(a)
        if f is None or fam2sub[f] != subset:
            continue
        ref[(a, b) if a < b else (b, a)] = e
    genes = {g for p in ref for g in p}
    pred = set()
    with open(pairs_tsv) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            a, b = line.split("\t", 2)[:2]
            if a in genes and b in genes:
                p = (a, b) if a < b else (b, a)
                if p in ref:
                    pred.add(p)
    for (a, b), e in ref.items():
        hit = (a, b) in pred
        yield a, b, ("TP" if hit else "FN") if e == "S" else ("FP" if hit else "TN")


def score(counts: Counter) -> dict:
    # RefPhyloTest stores each relation under its first protein only (specs[rel[1]]) but adds
    # length(...)/2 per protein: every reference pair therefore counts 1/2 (verified to 1e-8 vs official).
    tp, fp, fn, tn = (1 + counts[k] / 2 for k in ("TP", "FP", "FN", "TN"))
    tpr, ppv = tp / (tp + fn), tp / (tp + fp)
    return {"TPR": tpr, "TPR_ci": 1.96 * math.sqrt(tpr * (1 - tpr) / (tp + fn)),
            "PPV": ppv, "PPV_ci": 1.96 * math.sqrt(ppv * (1 - ppv) / (tp + fp)),
            **{f"n_{k}": counts[k] for k in ("TP", "FP", "FN", "TN")}}


def run(run_dir: str, subsets: list[str], gene2fam: dict, fam2sub: dict, from_pairs: bool = False) -> dict:
    c = {s: Counter() for s in ("all", "dev", "test")}
    unassigned = 0
    if from_pairs:
        if subsets != ["dev"]:
            raise SystemExit("--from-pairs scores dev families only")
        it = labels_from_pairs(os.path.join(run_dir, "pairs.tsv"), gene2fam, fam2sub, "dev")
    else:
        it = raw_pairs(run_dir)
    for a, b, lab in it:
        fa, fb = gene2fam.get(a), gene2fam.get(b)
        if fa is None or fa != fb:
            unassigned += 1
            continue
        c["all"][lab] += 1
        c[fam2sub[fa]][lab] += 1
    out = {"run_dir": run_dir, "n_pairs_not_in_one_family": unassigned}
    for s in subsets:
        out[s] = score(c[s])
    off = os.path.join(run_dir, "assess", "TreeFam-A.json")
    if os.path.exists(off) and "all" in subsets:
        o = {x["metrics"]["metric_id"]: x["metrics"]["value"] for x in json.load(open(off))
             if x.get("challenge_id") == "TreeFam-A"}
        out["official"] = o
        out["pooled_equals_official"] = all(abs(out["all"][m] - o[m]) < 1e-7 for m in ("TPR", "PPV"))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--subset", nargs="+", default=["dev", "all"], choices=["dev", "test", "all"])
    ap.add_argument("--allow-test", action="store_true")
    ap.add_argument("--split", default=SPLIT)
    ap.add_argument("--from-pairs", action="store_true",
                    help="RUN_DIR is a pipeline run: label dev reference pairs from its pairs.tsv (no QfO run)")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    if a.from_pairs:
        a.subset = ["dev"]
    if "test" in a.subset and not a.allow_test:
        raise SystemExit("test-family scores are not computed without --allow-test (DESIGN §4)")
    g2f, f2s = load_split(a.split)
    res = [run(r, a.subset, g2f, f2s, a.from_pairs) for r in a.run_dirs]
    for r in res:
        for s in a.subset:
            x = r[s]
            print(f"{os.path.basename(r['run_dir'].rstrip('/'))}\t{s}\tTPR {x['TPR']:.4f}±{x['TPR_ci']:.4f}\t"
                  f"PPV {x['PPV']:.4f}±{x['PPV_ci']:.4f}\tTP {x['n_TP']} FP {x['n_FP']} FN {x['n_FN']} TN {x['n_TN']}")
        if "official" in r:
            print(f"  pooled == official: {r['pooled_equals_official']} (official TPR {r['official']['TPR']:.8f}, "
                  f"PPV {r['official']['PPV']:.8f}; pooled {r['all']['TPR']:.8f}, {r['all']['PPV']:.8f}); "
                  f"pairs outside a single family: {r['n_pairs_not_in_one_family']}")
    if a.json:
        os.makedirs(os.path.dirname(os.path.abspath(a.json)), exist_ok=True)
        json.dump(res, open(a.json, "w"), indent=1)
    return 0 if all(r.get("pooled_equals_official", True) for r in res) else 1


if __name__ == "__main__":
    sys.exit(main())
