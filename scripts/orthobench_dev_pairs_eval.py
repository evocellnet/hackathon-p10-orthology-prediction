#!/usr/bin/env python
"""Pair-level OrthoBench DEV evaluation (DESIGN §3 S3 metric): ortholog-pair P/R/F1 inside DEV RefOGs.

Truth: every cross-species gene pair whose two genes are in the *same* DEV RefOG, labelled
ortholog / paralog by species overlap on the rooted RefOG reference tree
(scripts/S5_refog_labels.py: trees_tight, root minimising duplications, pruned to RefOG genes).
Prediction: a set of predicted ortholog pairs (pairs.tsv, an OrthoXML file flattened with
ancestra.orthoxml.hogs_to_pairs, or an in-memory OrthoXML / set).
Counts per truth pair: TP = ortholog predicted ortholog; FN = ortholog not predicted;
FP = paralog predicted ortholog; TN = paralog not predicted. Pairs outside DEV RefOGs and
pairs between RefOGs are ignored (not counted as FP).

Breakdown by species-pair class (mammal-mammal, tetrapod-tetrapod (with chicken), vert-teleost
(tetrapod x teleost), teleost-teleost, vert-invert (vertebrate x fly/worm/Ciona), invert-invert).
Paired bootstrap over RefOGs (B=2000, seed 11) for Delta between two runs.

    python scripts/orthobench_dev_pairs_eval.py PRED [PRED_B] [--json out.json]
      PRED = pairs.tsv | hogs.orthoxml | "truth" (self-check)
Benchmark-only code: never imported by the method (DESIGN §4 leakage rule).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

LABELS = os.path.join(ROOT, "work/S5/refog_pair_labels.parquet")
GROUP = {"Homo_sapiens": "mammal", "Pan_troglodytes": "mammal", "Mus_musculus": "mammal",
         "Rattus_norvegicus": "mammal", "Canis_familiaris": "mammal", "Monodelphis_domestica": "mammal",
         "Gallus_gallus": "bird", "Danio_rerio": "teleost", "Tetraodon_nigroviridis": "teleost",
         "Drosophila_melanogaster": "invert", "Caenorhabditis_elegans": "invert",
         "Ciona_intestinalis": "invert"}
CLASSES = ["mammal-mammal", "tetrapod-tetrapod", "vert-teleost", "teleost-teleost", "vert-invert",
           "invert-invert"]


def pair_class(sa: str, sb: str) -> str:
    ga, gb = sorted((GROUP[sa], GROUP[sb]))
    if ga == gb == "mammal":
        return "mammal-mammal"
    if {ga, gb} <= {"mammal", "bird"}:
        return "tetrapod-tetrapod"
    if ga == gb == "teleost":
        return "teleost-teleost"
    if "invert" in (ga, gb):
        return "invert-invert" if ga == gb else "vert-invert"
    return "vert-teleost"


def load_truth(subset: str = "dev") -> pd.DataFrame:
    """gene_a < gene_b, refog, label (1 ortholog / 0 paralog), cls. Within-RefOG pairs only."""
    if not os.path.exists(LABELS):
        import subprocess
        subprocess.run([sys.executable, os.path.join(ROOT, "scripts/S5_refog_labels.py"), "--subset", subset],
                       check=True)
    d = pd.read_parquet(LABELS)
    meta = json.load(open(LABELS.replace(".parquet", ".json")))
    assert meta["subset"] == subset, "label file was built for another subset"
    d = d[d.kind == "within"].copy()
    a = np.where(d.gene_a < d.gene_b, d.gene_a, d.gene_b)
    b = np.where(d.gene_a < d.gene_b, d.gene_b, d.gene_a)
    d["gene_a"], d["gene_b"] = a, b
    d["refog"] = d.refog_a
    d["cls"] = [pair_class(x, y) for x, y in zip(d.sp_a, d.sp_b)]
    d = d.drop_duplicates(["gene_a", "gene_b"])
    return d[["gene_a", "gene_b", "sp_a", "sp_b", "refog", "label", "cls"]].reset_index(drop=True)


def read_pred(src, genes: set | None = None) -> set:
    """Predicted ortholog pairs (sorted tuples). ``src``: path to pairs.tsv / *.orthoxml, an
    OrthoXML object, or a set. ``genes`` restricts the OrthoXML flattening to groups that
    contain at least one of these genes (speed)."""
    from ancestra.orthoxml import Group, OrthoXML, hogs_to_pairs, read_orthoxml
    if isinstance(src, set):
        return {tuple(sorted(p)) for p in src}
    if isinstance(src, str) and not src.endswith((".orthoxml", ".xml")):
        out = set()
        with open(src) as fh:
            for l in fh:
                if l.startswith("#") or not l.strip():
                    continue
                a, b = l.split("\t")[:2]
                out.add((a, b) if a < b else (b, a))
        return out
    ox = read_orthoxml(src) if isinstance(src, str) else src
    if genes is not None:
        keep = [g for g in ox.groups if isinstance(g, Group) and genes & set(g.genes())]
        ox = OrthoXML(genes=ox.genes, groups=keep, species_tree=ox.species_tree)
    return {tuple(sorted(p)) for p in hogs_to_pairs(ox)}


def counts(truth: pd.DataFrame, pred: set):
    """Per (RefOG, class) counts array [R, C, 4] = TP, FN, FP, TN, plus RefOG list."""
    key = list(zip(truth.gene_a, truth.gene_b))
    p = np.fromiter((k in pred for k in key), bool, len(key))
    y = truth.label.values == 1
    refogs = sorted(truth.refog.unique())
    ri = pd.Categorical(truth.refog, categories=refogs).codes
    ci = pd.Categorical(truth.cls, categories=CLASSES).codes
    cat = np.where(y & p, 0, np.where(y & ~p, 1, np.where(~y & p, 2, 3)))
    arr = np.zeros((len(refogs), len(CLASSES), 4), np.int64)
    np.add.at(arr, (ri, ci, cat), 1)
    return arr, refogs


def prf(c):
    """c[..., 4] -> precision, recall, F1 (nan when undefined)."""
    tp, fn, fp = c[..., 0].astype(float), c[..., 1].astype(float), c[..., 2].astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        P = tp / (tp + fp)
        R = tp / (tp + fn)
        F = 2 * P * R / (P + R)
    return P, R, F


def metrics(arr) -> dict:
    out = {}
    tot = arr.sum(axis=(0, 1))
    P, R, F = prf(tot)
    out["all"] = dict(P=_r(P), R=_r(R), F1=_r(F), tp=int(tot[0]), fn=int(tot[1]), fp=int(tot[2]), tn=int(tot[3]))
    for j, c in enumerate(CLASSES):
        t = arr[:, j].sum(axis=0)
        P, R, F = prf(t)
        out[c] = dict(P=_r(P), R=_r(R), F1=_r(F), tp=int(t[0]), fn=int(t[1]), fp=int(t[2]), tn=int(t[3]))
    return out


def _r(x):
    x = float(x)
    return None if np.isnan(x) else round(x, 4)


def paired_bootstrap(arr_a, arr_b, B: int = 2000, seed: int = 11) -> dict:
    """Delta (b - a) of P/R/F1 overall and per class; RefOGs resampled with replacement."""
    assert arr_a.shape == arr_b.shape
    R_ = arr_a.shape[0]
    rng = np.random.default_rng(seed)
    W = rng.multinomial(R_, np.full(R_, 1.0 / R_), size=B).astype(float)     # B x R weights
    out = {}
    groups = [("all", slice(None))] + [(c, j) for j, c in enumerate(CLASSES)]
    for name, sl in groups:
        ca = arr_a[:, sl].reshape(R_, -1, 4).sum(axis=1) if name == "all" else arr_a[:, sl]
        cb = arr_b[:, sl].reshape(R_, -1, 4).sum(axis=1) if name == "all" else arr_b[:, sl]
        if ca[:, :3].sum() == 0 and cb[:, :3].sum() == 0:
            continue
        ba, bb = W @ ca, W @ cb                                                    # B x 4
        res = {}
        for k, (va, vb) in zip(("P", "R", "F1"), zip(prf(ba), prf(bb))):
            d = vb - va
            d = d[~np.isnan(d)]
            pa, pb = prf(ca.sum(0)), prf(cb.sum(0))
            idx = ("P", "R", "F1").index(k)
            point = float(pb[idx] - pa[idx])
            res[k] = dict(delta=_r(point), ci95=[_r(np.percentile(d, 2.5)), _r(np.percentile(d, 97.5))]
                          if len(d) else None)
        out[name] = res
    return out


def evaluate(src, truth: pd.DataFrame | None = None) -> tuple[dict, np.ndarray, list]:
    truth = load_truth() if truth is None else truth
    if isinstance(src, str) and src == "truth":
        pred = set(zip(truth.gene_a[truth.label == 1], truth.gene_b[truth.label == 1]))
    else:
        genes = set(truth.gene_a) | set(truth.gene_b)
        pred = read_pred(src, genes)
    arr, refogs = counts(truth, pred)
    return metrics(arr), arr, refogs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pred")
    ap.add_argument("pred_b", nargs="?")
    ap.add_argument("--json")
    a = ap.parse_args()
    truth = load_truth()
    m, arr, refogs = evaluate(a.pred, truth)
    out = dict(pred=a.pred, n_truth_pairs=len(truth), n_refogs=len(refogs), metrics=m)
    if a.pred_b:
        mb, arrb, _ = evaluate(a.pred_b, truth)
        out.update(pred_b=a.pred_b, metrics_b=mb, delta_b_minus_a=paired_bootstrap(arr, arrb))
    print(json.dumps(out, indent=1))
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
