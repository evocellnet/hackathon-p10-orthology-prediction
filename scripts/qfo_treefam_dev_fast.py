"""Fast TreeFam-A DEV scoring of an in-memory / on-disk prediction (no Darwin run).

The official RefPhyloTest raw output of any earlier run lists every TreeFam-A reference
pair with its label; TP/FN pairs are reference *speciation* pairs and FP/TN pairs are
reference *duplication* pairs, independent of the prediction. We extract the DEV pairs
(splits/treefamA_split.json; test families are skipped at parse time) once, then score a
new prediction by checking which reference pairs it predicts as orthologs (QfO
flattening of the HOGs), with the official formula of scripts/qfo_treefam_subset.py
(pair counts / 2 + 1 pseudocount). Validated against the official dev values.

    python scripts/qfo_treefam_dev_fast.py PRED.orthoxml|pairs.tsv [--ref RAW.txt.gz]
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from qfo_treefam_subset import load_split  # noqa: E402

DEFAULT_RAW = os.path.join(ROOT, "work/qfo_peek1/det_go/raw/TreeFam-A/"
                           "RefPhylo-TreeFam-A_Ancestra-det-go_25e97643cce0_raw.txt.gz")
CACHE = os.path.join(ROOT, "work/qfo_recall/treefamA_dev_refpairs.json")


def dev_refpairs(raw: str = DEFAULT_RAW) -> list:
    """[(a, b, 'S'|'D', family)] for DEV families only."""
    if os.path.exists(CACHE):
        return [tuple(x) for x in json.load(open(CACHE))]
    gene2fam, fam2sub = load_split()
    out = []
    with gzip.open(raw, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _, a, b, lab = line.rstrip("\n").split("\t")
            f = gene2fam.get(a) or gene2fam.get(b)
            if f is None or fam2sub.get(f) != "dev":
                continue          # test families are never inspected
            out.append((a, b, "S" if lab in ("TP", "FN") else "D", f))
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    json.dump(out, open(CACHE, "w"))
    return out


def predicted_orthologs(src, genes: set) -> set:
    from ancestra.orthoxml import Group, OrthoXML, hogs_to_pairs, read_orthoxml
    if isinstance(src, str) and not src.endswith((".orthoxml", ".xml")):
        out = set()
        with open(src) as fh:
            for l in fh:
                if l.startswith("#"):
                    continue
                a, b = l.split("\t")[:2]
                if a in genes and b in genes:
                    out.add((a, b) if a < b else (b, a))
        return out
    ox = read_orthoxml(src) if isinstance(src, str) else src
    keep = [g for g in ox.groups if isinstance(g, Group) and genes & set(g.genes())]
    sub = OrthoXML(genes=ox.genes, groups=keep, species_tree=ox.species_tree)
    return {p for p in hogs_to_pairs(sub) if p[0] in genes and p[1] in genes}


def score(src, refs=None) -> dict:
    refs = refs if refs is not None else dev_refpairs()
    genes = {x for r in refs for x in r[:2]}
    pred = predicted_orthologs(src, genes)
    c = Counter()
    for a, b, t, _ in refs:
        p = ((a, b) if a < b else (b, a)) in pred
        c[("TP" if p else "FN") if t == "S" else ("FP" if p else "TN")] += 1
    tp, fp, fn, tn = (1 + c[k] / 2 for k in ("TP", "FP", "FN", "TN"))
    tpr, ppv = tp / (tp + fn), tp / (tp + fp)
    return {"TPR": tpr, "PPV": ppv, "TPR_ci": 1.96 * math.sqrt(tpr * (1 - tpr) / (tp + fn)),
            "PPV_ci": 1.96 * math.sqrt(ppv * (1 - ppv) / (tp + fp)), **{f"n_{k}": c[k] for k in ("TP", "FP", "FN", "TN")}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pred")
    ap.add_argument("--ref", default=DEFAULT_RAW)
    a = ap.parse_args()
    print(json.dumps(score(a.pred, dev_refpairs(a.ref)), indent=1))


if __name__ == "__main__":
    main()
