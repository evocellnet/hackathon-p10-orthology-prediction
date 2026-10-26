#!/usr/bin/env python3
"""Pivot the tidy OpenEBench TSV into one row per tool with the headline metrics.
Usage: make_reference_table.py reference_results/openebench_QfO_2022.tsv > reference_results/qfo2022_summary_wide.tsv"""
import csv, sys
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
tab = defaultdict(dict)
for r in rows:
    c, t = r["challenge"], r["tool"]
    for ax in ("x", "y"):
        m = r[ax + "_metric"].replace("STD_num_orthologs", "NR_ORTHOLOGS")
        tab[t][f"{c}:{m}"] = r[ax]
cols = ["SwissTree:TPR", "SwissTree:PPV", "TreeFam-A:TPR", "TreeFam-A:PPV", "VGNC:TPR", "VGNC:PPV",
        "GOtest:NR_ORTHOLOGS", "GOtest:avg Schlicker", "ECtest:NR_ORTHOLOGS", "ECtest:avg Schlicker"]
for c in ["G_STD2_Luca", "G_STD2_Eukaryota", "G_STD2_Fungi", "G_STD2_Vertebrata", "STD_Eukaryota", "STD_Fungi", "STD_Bacteria"]:
    cols += [f"{c}:NR_COMPLETED_TREE_SAMPLINGS", f"{c}:RF_DISTANCE", f"{c}:FRAC_INCORRECT_TREES"]
cols += ["FAS:NR_ORTHOLOGS", "FAS:FAS"]
w = csv.writer(sys.stdout, delimiter="\t")
w.writerow(["tool"] + cols)
for t in sorted(tab):
    w.writerow([t] + [("%.4g" % float(tab[t][c])) if tab[t].get(c) not in (None, "") else "" for c in cols])
