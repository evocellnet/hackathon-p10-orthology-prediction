#!/usr/bin/env python3
"""Compare a local QfO run (summary.tsv from run_qfo_local.sh) with the published OpenEBench
numbers of a participant (reference_results/openebench_QfO_<event>.tsv).
Usage: compare_to_published.py summary.tsv openebench_QfO_2022.tsv TOOL_ID
Prints benchmark, metric, local, published, published stderr, |diff|/stderr.
Also usable to compare our own method with any competitor (pass two summary files with --local2)."""
import csv, sys
summ, ref, tool = sys.argv[1:4]
ALIAS = {"GO": "GOtest", "EC": "ECtest", "SwissTrees": "SwissTree", "G_STD2_Luca": "G_STD2_LUCA"}
loc = {}
for r in csv.DictReader(open(summ), delimiter="\t"):
    loc[(r["benchmark"], r["metric"])] = float(r["value"])
pub = {}
for r in csv.DictReader(open(ref), delimiter="\t"):
    if r["tool"] != tool:
        continue
    for ax in ("x", "y"):
        m = r[ax + "_metric"].replace("STD_num_orthologs", "NR_ORTHOLOGS")
        pub[(r["challenge"], m)] = (float(r[ax]), float(r["stderr_" + ax] or 0))
print("benchmark\tmetric\tlocal\tpublished\tpub_stderr\tz")
for (b, m), v in sorted(loc.items()):
    key = None
    for c in (b, ALIAS.get(b, b), b.upper()):
        if (c, m) in pub:
            key = (c, m); break
    if key is None:
        print(f"{b}\t{m}\t{v:.4g}\tNA\tNA\tNA"); continue
    p, se = pub[key]
    z = abs(v - p) / se if se > 0 else (0.0 if abs(v - p) < 1e-9 else float("inf"))
    print(f"{b}\t{m}\t{v:.4g}\t{p:.4g}\t{se:.3g}\t{z:.2f}")
