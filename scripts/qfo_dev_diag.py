#!/usr/bin/env python3
"""Diagnose missed (FN) and false (FP) TreeFam-A ortholog pairs of a run, on DEV families only.

Only pairs of dev families (splits/treefamA_split.json) are read from the official raw RefPhyloTest output;
test-family pairs are skipped at parse time (DESIGN §4).

FN categories (mutually exclusive, in this order):
  (i)   different families: the two genes are in different families of the partition the merge used
        (S0 after the S2 cross-family rescue, ``stages/scores/families_final.npy``); also reported: different
        already in raw S0 (``stages/families/families.npy``), and whether one of them is a singleton family.
  (iv)  same family, but at least one gene is a singleton in the output (not in any root HOG).
  (ii)  same family, both in root HOGs, different root HOGs (merge split).
  (iii) same root HOG, but the deepest common group is a paralogGroup (wrong duplication call).
For every pair: the deepest common taxonomic level = LCA of the two species in the guide tree, and whether
the pair was scored in S2 (``scores.npz``; with min(qcov,tcov) >= 0.5 = passes the merge coverage guard).
FP: by level (= LCA level = TaxRange of the orthologGroup that induced the pair).
Comparison: which of our FN pairs are TP for a comparator (OrthoFinder public predictions, same raw format).

    python scripts/qfo_dev_diag.py --run work/qfo2022_det/det_go --peek work/qfo_peek1/det_go \
        --comparator OrthoFinder=work/qfo_peek1/competitors/OrthoFinder --out results/qfo_dev_diag
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SPLIT = os.path.join(ROOT, "splits/treefamA_split.json")


def dev_raw(peek_dir: str, gene2fam: dict, fam2sub: dict) -> dict:
    fs = glob.glob(os.path.join(peek_dir, "raw", "TreeFam-A", "RefPhylo-TreeFam-A_*_raw.txt.gz"))
    assert len(fs) == 1, fs
    out = {}
    with gzip.open(fs[0], "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _, a, b, lab = line.rstrip("\n").split("\t")
            f = gene2fam.get(a)
            if f is None or fam2sub[f] != "dev":     # test pairs are never kept
                continue
            out[(a, b) if a < b else (b, a)] = (lab, f)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--peek", required=True)
    ap.add_argument("--comparator", action="append", default=[], help="NAME=peek_dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--examples", type=int, default=30)
    a = ap.parse_args(argv)
    from ancestra.newick import TreeIndex, read_newick
    from ancestra.orthoxml import Group, read_orthoxml

    sp = json.load(open(SPLIT))
    gene2fam = {g: f for f, d in sp["families"].items() for g in d["members"]}
    fam2sub = {f: d["subset"] for f, d in sp["families"].items()}
    ours = dev_raw(a.peek, gene2fam, fam2sub)
    comps = {n: dev_raw(p, gene2fam, fam2sub) for n, p in (c.split("=", 1) for c in a.comparator)}
    dev_genes = {g for p in ours for g in p}

    # run data
    gid, sp_of = {}, {}
    with open(os.path.join(a.run, "stages/input/genes.tsv")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            i, p, s, _ = line.rstrip("\n").split("\t")
            gid[p] = int(i)
            sp_of[p] = s
    fam0 = np.load(os.path.join(a.run, "stages/families/families.npy"))
    famF = np.load(os.path.join(a.run, "stages/scores/families_final.npy"))
    sz0, szF = np.bincount(fam0), np.bincount(famF)
    tree = read_newick(os.path.join(a.run, "stages/input/guide_tree.nwk"))
    tix = TreeIndex(tree)
    ox = read_orthoxml(os.path.join(a.run, "hogs.orthoxml"))
    root_of, top_size = {}, []
    for ti, top in enumerate(t for t in ox.groups if isinstance(t, Group)):
        gs = top.genes()
        top_size.append(len(gs))
        for g in gs:
            root_of[g] = ti
    tops = [t for t in ox.groups if isinstance(t, Group)]

    def deepest_common(top: Group, x: str, y: str):
        """(kind, level) of the deepest group containing both x and y."""
        best = (top.kind, top.level)
        node = top
        while True:
            nxt = None
            for c in node.children:
                if isinstance(c, Group):
                    gs = set(c.genes())
                    if x in gs and y in gs:
                        nxt = c
                        break
            if nxt is None:
                return best
            node = nxt
            best = (node.kind, node.level)

    # scores (S2): pair -> min coverage
    z = np.load(os.path.join(a.run, "stages/scores/scores.npz"))
    q, t = z["q"].astype(np.int64), z["t"].astype(np.int64)
    n = len(gid)
    want = set()
    for (x, y), (lab, _) in ours.items():
        if lab == "FN":
            i, j = sorted((gid[x], gid[y]))
            want.add(i * n + j)
    key = q * n + t
    sel = np.isin(key, np.fromiter(want, dtype=np.int64))
    mincov = dict(zip(key[sel].tolist(), np.minimum(z["qcov"][sel], z["tcov"][sel]).tolist()))
    del q, t, key, z

    # pairs.tsv levels for FP
    fp_pairs = {p for p, (lab, _) in ours.items() if lab == "FP"}
    lvl_of = {}
    with open(os.path.join(a.run, "pairs.tsv")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            x, y, lv = line.rstrip("\n").split("\t")
            if x in dev_genes and y in dev_genes:
                p = (x, y) if x < y else (y, x)
                if p in fp_pairs:
                    lvl_of[p] = lv

    rows = []
    for (x, y), (lab, f) in ours.items():
        if lab not in ("FN", "FP"):
            continue
        i, j = gid[x], gid[y]
        lca = tix.mrca([sp_of[x], sp_of[y]])
        r = {"pair": f"{x}|{y}", "label": lab, "family": f, "lca_level": lca, "lca_depth": tix.depth[lca],
             "sp": f"{sp_of[x]}|{sp_of[y]}"}
        for cn, cr in comps.items():
            r[f"{cn}_label"] = cr.get((x, y), ("NA",))[0]
        if lab == "FP":
            r["our_level"] = lvl_of.get((x, y), "?")
            rows.append(r)
            continue
        ii, jj = sorted((i, j))
        mc = mincov.get(ii * n + jj)
        r["scored"] = mc is not None
        r["cov_ok"] = mc is not None and mc >= 0.5
        r["diff_S0"] = bool(fam0[i] != fam0[j])
        r["diff_final"] = bool(famF[i] != famF[j])
        r["singleton_family"] = bool(szF[famF[i]] == 1 or szF[famF[j]] == 1)
        ri, rj = root_of.get(x), root_of.get(y)
        if r["diff_final"]:
            cat = "i_diff_family"
        elif ri is None or rj is None:
            cat = "iv_singleton_in_output"
        elif ri != rj:
            cat = "ii_diff_root_hog"
            r["root_sizes"] = f"{top_size[ri]}|{top_size[rj]}"
        else:
            kind, lvl = deepest_common(tops[ri], x, y)
            cat = "iii_paralog_call" if kind == "paralog" else "iii_other"
            r["dup_level"] = lvl
        r["category"] = cat
        rows.append(r)

    fn = [r for r in rows if r["label"] == "FN"]
    fp = [r for r in rows if r["label"] == "FP"]
    cats = Counter(r["category"] for r in fn)
    summ = {"n_dev_families": sum(1 for v in fam2sub.values() if v == "dev"),
            "dev_labels_ours": dict(Counter(l for l, _ in ours.values())),
            "FN_total": len(fn), "FN_by_category": dict(cats.most_common()),
            "FN_i_also_diff_in_raw_S0": sum(r["diff_S0"] for r in fn if r["category"] == "i_diff_family"),
            "FN_i_with_singleton_family": sum(r["singleton_family"] for r in fn if r["category"] == "i_diff_family"),
            "FN_rescue_fixed_family_split": sum(r["diff_S0"] and not r["diff_final"] for r in fn),
            "FN_scored_by_category": {c: sum(r["scored"] for r in fn if r["category"] == c) for c in cats},
            "FN_covok_by_category": {c: sum(r["cov_ok"] for r in fn if r["category"] == c) for c in cats},
            "FN_level_by_category": {c: dict(Counter(r["lca_level"] for r in fn if r["category"] == c).most_common(12))
                                     for c in cats},
            "FN_iii_dup_level": dict(Counter(r.get("dup_level") for r in fn
                                             if r["category"] == "iii_paralog_call").most_common(12)),
            "FN_families_top": dict(Counter(r["family"] for r in fn).most_common(10)),
            "FP_total": len(fp), "FP_by_level": dict(Counter(r["our_level"] for r in fp).most_common()),
            "FP_families_top": dict(Counter(r["family"] for r in fp).most_common(10))}
    for cn in comps:
        summ[f"FN_{cn}_TP_by_category"] = {c: sum(r[f"{cn}_label"] == "TP" for r in fn if r["category"] == c)
                                           for c in cats}
        summ[f"FN_{cn}_TP_total"] = sum(r[f"{cn}_label"] == "TP" for r in fn)
        summ[f"FP_{cn}_also_FP"] = sum(r[f"{cn}_label"] == "FP" for r in fp)
        oc = comps[cn]
        summ[f"{cn}_FN_that_we_get_TP"] = sum(1 for p, (lab, _) in oc.items() if lab == "FN"
                                             and ours.get(p, ("",))[0] == "TP")
        summ[f"{cn}_dev_labels"] = dict(Counter(l for l, _ in oc.values()))
    os.makedirs(a.out, exist_ok=True)
    json.dump(summ, open(os.path.join(a.out, "diag_summary.json"), "w"), indent=1)
    cols = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("pair", "label", "category"), k))
    with open(os.path.join(a.out, "diag_pairs.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
