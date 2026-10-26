#!/usr/bin/env python3
"""Thin wrapper around the official Open_Orthobench benchmark.py (unmodified).
Imports its functions, so scoring is identical; only adds (i) graceful handling of the
0/0 case (F reported as 0 when TP=0) and (ii) machine-readable TSV/JSON output.
Usage: orthobench_eval.py <BENCHMARKS dir> <orthogroups file> [--basic] [--json out.json] [--label NAME]
"""
import argparse, importlib.util, io, json, os, re, sys, contextlib

p = argparse.ArgumentParser()
p.add_argument("bench_dir"); p.add_argument("ogs"); p.add_argument("--basic", action="store_true")
p.add_argument("--json"); p.add_argument("--label", default=None)
a = p.parse_args()
spec = importlib.util.spec_from_file_location("ob", os.path.join(a.bench_dir, "benchmark.py"))
ob = importlib.util.module_from_spec(spec); spec.loader.exec_module(ob)
ob.__file__ = os.path.join(a.bench_dir, "benchmark.py")
d_ref = os.path.join(a.bench_dir, "RefOGs") + os.sep
ref = ob.read_refogs(d_ref)
unc = ob.read_uncertain_refogs(d_ref + "low_certainty_assignments/")
# get_expected_genes uses os.path.dirname(__file__) of the module
exp = ob.get_expected_genes()
if a.basic:
    pred = ob.read_orthogroups(a.ogs, exp, ob.get_n_col_skip(a.ogs))
else:
    pred = ob.read_orthogroups_smart(a.ogs)
ob.check_orthogroups(pred, exp)
buf = io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        f, prec, rec = ob.calculate_benchmarks_pairwise(ref, unc, pred)
    n_exact = int(re.search(r"(\d+) Orthogroups exactly correct", buf.getvalue()).group(1))
except ZeroDivisionError:
    f, prec, rec, n_exact = 0.0, 0.0, 0.0, 0
res = dict(label=a.label or os.path.basename(a.ogs), F=round(f, 2), precision=round(prec, 2),
           recall=round(rec, 2), exact_refogs=n_exact, n_pred_groups=len(pred),
           n_genes_assigned=sum(len(x) for x in pred))
print("label\tF\tprecision\trecall\texact_refogs")
print("%s\t%.1f\t%.1f\t%.1f\t%d" % (res["label"], f, prec, rec, n_exact))
if a.json:
    json.dump(res, open(a.json, "w"), indent=1)
