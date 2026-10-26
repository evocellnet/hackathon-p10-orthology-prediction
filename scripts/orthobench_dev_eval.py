"""OrthoBench scoring restricted to the DEV RefOGs (splits/orthobench_split.json).

Uses the official, unmodified ``BENCHMARKS/benchmark.py`` functions (read_refogs,
read_uncertain_refogs, read_orthogroups_smart, calculate_benchmarks_pairwise); the
only change is that the RefOG lists passed to calculate_benchmarks_pairwise contain
the dev RefOGs only. Test RefOGs are never scored by this script (hard guard).

Per dev RefOG it also reports (same uncertain-gene handling as benchmark.py):
TP/FP/FN pairs, #predicted groups it is split into, #genes missing (in no
multi-gene group), the largest predicted group's share, and the FP partner
genes' membership (another DEV RefOG vs. any other gene; test RefOGs are never looked at) -- used to
tell S0 splits from merge errors.

    python scripts/orthobench_dev_eval.py OGS_FILE [--json out.json] [--per-refog out.tsv]
"""
import argparse
import contextlib
import importlib.util
import io
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS")
SPLIT = os.path.join(ROOT, "splits/orthobench_split.json")


def load_bench():
    spec = importlib.util.spec_from_file_location("ob", os.path.join(BENCH, "benchmark.py"))
    ob = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ob)
    ob.__file__ = os.path.join(BENCH, "benchmark.py")
    return ob


def dev_eval(ogs_file: str, subset: str = "dev") -> dict:
    if subset != "dev":
        raise ValueError("only the dev subset may be scored before the freeze")
    ob = load_bench()
    d_ref = os.path.join(BENCH, "RefOGs") + os.sep
    with contextlib.redirect_stdout(io.StringIO()):
        ref = ob.read_refogs(d_ref)
        unc = ob.read_uncertain_refogs(d_ref + "low_certainty_assignments/")
    split = json.load(open(SPLIT))
    dev_ids = split["dev"]
    idx = [int(r.replace("RefOG", "")) - 1 for r in dev_ids]
    assert not set(dev_ids) & set(split["test"])
    ref_d = [ref[i] for i in idx]
    unc_d = [unc[i] for i in idx]
    pred = ob.read_orthogroups_smart(ogs_file)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        f, p, r = ob.calculate_benchmarks_pairwise(ref_d, unc_d, pred)
    n_exact = int(re.search(r"(\d+) Orthogroups exactly correct", buf.getvalue()).group(1))
    # per-RefOG diagnostics
    dev_ref_genes = {}            # only DEV RefOG membership is ever looked at
    for i in idx:
        for g in ref[i]:
            dev_ref_genes.setdefault(g, i)
    g2og = {}
    for k, og in enumerate(pred):
        for g in og:
            g2og[g] = k
    rows = []
    for rid, i in zip(dev_ids, idx):
        refog = ref[i] - unc[i]
        n = len(refog)
        tp = fp = fn = 0
        groups = {}
        for g in refog:
            if g in g2og:
                groups.setdefault(g2og[g], set()).add(g)
        missing = n - sum(len(v) for v in groups.values())
        fp_dev = fp_other = 0
        for k, mem in groups.items():
            og = pred[k] - unc[i]
            ov = len(mem)
            tp += ov * (ov - 1) // 2
            fp += ov * (len(og) - ov)
            fn += (n - ov) * ov
            for g in og - refog:
                if g in dev_ref_genes:
                    fp_dev += 1
                else:
                    fp_other += 1
        fn = (fn + missing * (n - 1)) // 2
        largest = max((len(v) for v in groups.values()), default=0)
        rows.append(dict(refog=rid, n=n, tp=tp, fp=fp, fn=fn, n_groups=len(groups), missing=missing,
                         largest_share=largest / n if n else 0, fp_genes_other_dev_refog=fp_dev,
                         fp_genes_other=fp_other,
                         exact=(fp == 0 and fn == 0)))
    n_fission = sum(1 for x in rows if x["n_groups"] + (1 if x["missing"] else 0) > 1)
    n_fusion = sum(1 for x in rows if x["fp"] > 0)
    return dict(subset="dev", n_refogs=len(dev_ids), F=f, precision=p, recall=r, exact=n_exact,
                n_refogs_split=n_fission, n_refogs_with_fp=n_fusion, n_pred_groups=len(pred),
                per_refog=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ogs")
    ap.add_argument("--json")
    ap.add_argument("--per-refog")
    a = ap.parse_args()
    res = dev_eval(a.ogs)
    print(f"DEV ({res['n_refogs']} RefOGs): F={res['F']:.2f} P={res['precision']:.2f} R={res['recall']:.2f} "
          f"exact={res['exact']} split={res['n_refogs_split']} with_FP={res['n_refogs_with_fp']}")
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)
    if a.per_refog:
        with open(a.per_refog, "w") as fh:
            keys = list(res["per_refog"][0])
            fh.write("\t".join(keys) + "\n")
            for r in res["per_refog"]:
                fh.write("\t".join(str(r[k]) for k in keys) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
