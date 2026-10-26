"""Diagnose OrthoBench DEV errors of an `ancestra run` output directory.

    python scripts/diagnose_orthobench.py RUN_DIR [--out per_refog.tsv]

For every DEV RefOG (certain genes only; test RefOGs are never read):
* S0: number of final families (after cross-family rescue) holding its genes, largest share;
* merge: number of root HOGs (orthobench_ogs.txt groups) holding its genes within its main
  family -> a merge fission if > 1; for each piece, the MRCA clade of its species;
* FP: genes sharing a root HOG with the RefOG that are not in it -> split into genes of
  another DEV RefOG vs other genes, plus the MRCA of the FP genes' species.
Classification: 'exact', 'S0_split', 'merge_split', 'overmerge', or combinations.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import read_gene_table  # noqa: E402
from ancestra.newick import TreeIndex, read_newick  # noqa: E402

BENCH = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out")
    ap.add_argument("--ogs", help="alternative orthobench_ogs.txt (e.g. RUN/iter/<hash>/...)")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    split = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
    gt = read_gene_table(os.path.join(a.run, "stages/input/genes.tsv"))
    idx = gt.index()
    tree = read_newick(os.path.join(a.run, "stages/input/guide_tree.nwk"))
    tix = TreeIndex(tree)
    famf = os.path.join(a.run, "stages/scores/families_final.npy")
    fam = np.load(famf if os.path.exists(famf) else os.path.join(a.run, "stages/families/families.npy"))
    og_of = {}
    with open(a.ogs or os.path.join(a.run, "orthobench_ogs.txt")) as fh:
        next(fh)
        for k, line in enumerate(fh):
            for g in line.rstrip("\n").split("\t")[1:]:
                og_of[g] = k
    og_members = defaultdict(list)
    for g, k in og_of.items():
        og_members[k].append(g)
    dev_genes = {}
    refogs = {}
    for rid in split["dev"]:
        genes = {l.strip() for l in open(os.path.join(BENCH, "RefOGs", rid + ".txt")) if l.strip()}
        uf = os.path.join(BENCH, "RefOGs/low_certainty_assignments", rid + ".txt")
        unc = {l.strip() for l in open(uf)} if os.path.exists(uf) else set()
        refogs[rid] = genes - unc
        for g in genes:
            dev_genes[g] = rid
    sp = lambda g: gt.species[gt.species_of[idx[g]]]  # noqa: E731
    rows = []
    for rid, genes in refogs.items():
        genes = [g for g in genes if g in idx]
        fams = Counter(int(fam[idx[g]]) for g in genes)
        main_fam, main_n = fams.most_common(1)[0]
        ogs = Counter(og_of[g] for g in genes)
        in_main = [g for g in genes if fam[idx[g]] == main_fam]
        ogs_main = Counter(og_of[g] for g in in_main)
        pieces = []
        for k, c in ogs_main.most_common():
            pieces.append(f"{c}@{tix.mrca(sp(g) for g in in_main if og_of[g] == k)}")
        fp_dev = fp_other = 0
        fp_sp = []
        for k in ogs:
            for g in og_members[k]:
                if g in genes:
                    continue
                if g in dev_genes and dev_genes[g] != rid:
                    fp_dev += 1
                else:
                    fp_other += 1
                fp_sp.append(sp(g))
        cls = []
        if len(fams) > 1:
            cls.append("S0_split")
        if len(ogs_main) > 1:
            cls.append("merge_split")
        if fp_dev + fp_other > 0:
            cls.append("overmerge")
        rows.append(dict(refog=rid, n=len(genes), n_families=len(fams), main_family_share=main_n / len(genes),
                         n_root_hogs=len(ogs), n_root_hogs_in_main_family=len(ogs_main),
                         pieces=";".join(pieces), fp_dev_refog_genes=fp_dev, fp_other_genes=fp_other,
                         fp_mrca=tix.mrca(fp_sp) if fp_sp else "", cls="+".join(cls) or "exact"))
    c = Counter(r["cls"] for r in rows)
    print("DEV RefOG error classes:", dict(c))
    print("S0-split RefOGs:", sum(1 for r in rows if "S0_split" in r["cls"]),
          " merge-split:", sum(1 for r in rows if "merge_split" in r["cls"]),
          " overmerged:", sum(1 for r in rows if "overmerge" in r["cls"]))
    for r in sorted(rows, key=lambda r: r["refog"]):
        if r["cls"] != "exact" and not a.quiet:
            print(f"  {r['refog']} n={r['n']} {r['cls']}: fams={r['n_families']} share={r['main_family_share']:.2f} "
                  f"pieces={r['pieces']} FP(dev/other)={r['fp_dev_refog_genes']}/{r['fp_other_genes']} "
                  f"fpMRCA={r['fp_mrca']}")
    if a.out:
        with open(a.out, "w") as fh:
            keys = list(rows[0])
            fh.write("\t".join(keys) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[k]) for k in keys) + "\n")


if __name__ == "__main__":
    main()
