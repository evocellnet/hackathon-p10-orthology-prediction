#!/usr/bin/env python3
"""QfO peek 1 table: our arms next to the published QfO 2022 rows (OpenEBench OEBE0020000003).

Inputs
* ours: work/qfo_peek1/<arm>/{summary.tsv, raw/, orthologs.db} (bench_repo/scripts/run_qfo_local.sh)
* published: bench_repo/reference_results/openebench_QfO_2022.tsv
* competitors re-scored locally on SwissTrees/TreeFam-A/VGNC/SP (scripts/qfo_peek1_rescore_competitors.sh):
  work/qfo_peek1/competitors/<tool>/ -> evaluated-pair counts + reproduction check vs published.

Columns
* published metrics (competitors) / local metrics (ours): SwissTree, TreeFam-A, VGNC TPR/PPV with 95% CI half-width
  (QfO's reported "stderr" = 1.96*sqrt(var)); GO and EC avg Schlicker with n = NR_ORTHOLOGS (pairs with annotations
  used by the test).
* SP TPR/PPV (SwissProt-IDs test, not in the published event: local values for all rows).
* evaluated pairs = TP+FP from the raw files of SwissTree, TreeFam-A, VGNC, and SP "Nr orthologs among SwissProt
  entries"; n_pairs_mapped = distinct predicted ortholog pairs mapped by map_relations.py (orthologs.db).

Writes <out>/summary.tsv, summary.md, competitor_reproduction.tsv (default out: results/qfo_peek1).
Usage for later peeks: --ours "LABEL=work/qfo_peekN/<arm>=local (peek N)" ... --out results/qfo_peekN
"""
from __future__ import annotations

import csv
import glob
import gzip
import os
import sqlite3
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PEEK = os.path.join(ROOT, "work/qfo_peek1")
OUT = os.path.join(ROOT, "results/qfo_peek1")
PUB = os.path.join(ROOT, "bench_repo/reference_results/openebench_QfO_2022.tsv")
OURS = [("Ancestra A (gene order, det_go)", "det_go"), ("Ancestra B (no gene order, det_nogo)", "det_nogo")]
COMP = [("OrthoFinder", "OrthoFinder"), ("FastOMA_v0.1.4_nov", "FastOMA_v0.1.4_nov"),
        ("Hieranoid_2", "Hieranoid_2"), ("InParanoid5", "InParanoid5"),
        ("SonicParanoid2-sens-g", "SonicParanoid2-sens-g"), ("OMA_HOGs", "OMA_HOGs"),
        ("OMA_Pairs", "OMA_Pairs"), ("PANTHER_(v18)_all", "PANTHER_v18_all"),
        ("pbDOOD (eggnog v7)", "pbDOOD_eggnog_v7"), ("pbDOOD (midpoint)", "pbDOOD_midpoint")]
PUB_CH = {"SwissTree": "SwissTrees", "TreeFam-A": "TreeFam-A", "VGNC": "VGNC", "GOtest": "GO", "ECtest": "EC"}


def published() -> dict:
    out: dict = {}
    for r in csv.DictReader(open(PUB), delimiter="\t"):
        if r["challenge"] not in PUB_CH or "_agg_" not in r["dataset"]:
            continue
        b = PUB_CH[r["challenge"]]
        d = out.setdefault(r["tool"], {})
        d[(b, r["x_metric"])] = (float(r["x"]), float(r["stderr_x"] or 0))
        d[(b, r["y_metric"])] = (float(r["y"]), float(r["stderr_y"] or 0))
    return out


def local(run_dir: str) -> dict:
    d: dict = {}
    p = os.path.join(run_dir, "summary.tsv")
    if not os.path.exists(p):
        return d
    for r in csv.DictReader(open(p), delimiter="\t"):
        try:
            d[(r["benchmark"], r["metric"])] = (float(r["value"]), float(r["stderr"] or 0))
        except ValueError:
            pass
    for b in ("SwissTrees", "TreeFam-A", "VGNC"):
        fs = glob.glob(os.path.join(run_dir, "raw", b, "*raw.txt.gz"))
        if not fs:
            continue
        c: Counter = Counter()
        with gzip.open(fs[0], "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                lab = parts[3] if b != "VGNC" else parts[2]
                c[lab] += 1
        d[(b, "TP")], d[(b, "FP")], d[(b, "FN")] = c["TP"], c["FP"], c["FN"]
        d[(b, "evaluated")] = c["TP"] + c["FP"]
    db = os.path.join(run_dir, "orthologs.db")
    if os.path.exists(db):
        con = sqlite3.connect(db)
        n = con.execute("SELECT COUNT(*) FROM (SELECT DISTINCT prot_nr1, prot_nr2 FROM orthologs "
                        "WHERE prot_nr1 < prot_nr2)").fetchone()[0]
        con.close()
        d[("all", "n_pairs_mapped")] = n
    return d


def f(x, nd=3):
    if x is None:
        return ""
    if isinstance(x, tuple):
        v, se = x
        return f"{v:.{nd}f}±{se:.{nd}f}" if se else f"{v:.{nd}f}"
    return f"{x:,}" if isinstance(x, int) else f"{x:.{nd}f}"


def row(name: str, src: str, pub: dict | None, loc: dict) -> dict:
    g = (lambda k: pub.get(k)) if pub is not None else (lambda k: loc.get(k))
    gn = (lambda k: int(pub[k][0]) if pub and k in pub else (int(loc[k][0]) if k in loc else None))
    return {"method": name, "source": src,
            "SwissTree TPR": f(g(("SwissTrees", "TPR"))), "SwissTree PPV": f(g(("SwissTrees", "PPV"))),
            "TreeFam-A TPR": f(g(("TreeFam-A", "TPR"))), "TreeFam-A PPV": f(g(("TreeFam-A", "PPV"))),
            "VGNC TPR": f(g(("VGNC", "TPR"))), "VGNC PPV": f(g(("VGNC", "PPV"))),
            "GO avg Schlicker": f(g(("GO", "avg Schlicker")), 4), "GO n": f(gn(("GO", "NR_ORTHOLOGS"))),
            "EC avg Schlicker": f(g(("EC", "avg Schlicker")), 4), "EC n": f(gn(("EC", "NR_ORTHOLOGS"))),
            "SP TPR (local)": f(loc.get(("SP", "TPR"))), "SP PPV (local)": f(loc.get(("SP", "PPV"))),
            "SwissTree evaluated pairs (TP+FP)": f(loc.get(("SwissTrees", "evaluated"))),
            "TreeFam-A evaluated pairs (TP+FP)": f(loc.get(("TreeFam-A", "evaluated"))),
            "VGNC evaluated pairs (TP+FP)": f(loc.get(("VGNC", "evaluated"))),
            "SP pairs among SwissProt entries": f(int(loc[("SP", "Nr orthologs among SwissProt entries")][0])
                                                  if ("SP", "Nr orthologs among SwissProt entries") in loc else None),
            "predicted pairs mapped (all)": f(loc.get(("all", "n_pairs_mapped")))}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="QfO peek table (ours + published competitors)")
    ap.add_argument("--ours", action="append", default=[],
                    help="LABEL=QFO_RUN_DIR[=SOURCE] (default: the two peek-1 arms)")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    out_dir = a.out
    ours = [tuple(x.split("=", 2)) for x in a.ours] or \
        [(n, os.path.join(PEEK, d), "local (peek 1)") for n, d in OURS]
    os.makedirs(out_dir, exist_ok=True)
    pub = published()
    rows, repro = [], []
    for name, d, *src in ours:
        rows.append(row(name, src[0] if src else "local", None, local(d)))
    for name, loc_name in COMP:
        loc = local(os.path.join(PEEK, "competitors", loc_name))
        p = pub.get(name)
        if p is None:
            print(f"[table] no published row for {name}", file=sys.stderr)
            continue
        rows.append(row(name, "published (OpenEBench 2022); pair counts + SP local", p, loc))
        for b in ("SwissTrees", "TreeFam-A", "VGNC"):
            for m in ("TPR", "PPV"):
                if (b, m) in p and (b, m) in loc:
                    repro.append({"method": name, "benchmark": b, "metric": m, "published": p[(b, m)][0],
                                  "local": loc[(b, m)][0], "abs_diff": abs(p[(b, m)][0] - loc[(b, m)][0])})
    cols = list(rows[0])
    with open(os.path.join(out_dir, "summary.tsv"), "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(out_dir, "summary.md"), "w") as fh:
        fh.write("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n")
        for r in rows:
            fh.write("| " + " | ".join(r[c] for c in cols) + " |\n")
    if repro:
        with open(os.path.join(out_dir, "competitor_reproduction.tsv"), "w") as fh:
            w = csv.DictWriter(fh, list(repro[0]), delimiter="\t")
            w.writeheader()
            w.writerows(repro)
        print(f"[table] reproduction: max |local - published| = {max(r['abs_diff'] for r in repro):.2e} "
              f"over {len(repro)} values", file=sys.stderr)
    print(f"[table] wrote {out_dir}/summary.tsv ({len(rows)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
