"""Verifier: published competitor rows in results/qfo_final/summary.tsv vs bench_repo/reference_results/openebench_QfO_2022.tsv."""
import csv, json
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
ref = list(csv.DictReader(open(f"{ROOT}/bench_repo/reference_results/openebench_QfO_2022.tsv"), delimiter="\t"))
R = {}
for r in ref:
    R[(r["challenge"], r["x_metric"], r["y_metric"], r["tool"])] = r
summ = list(csv.DictReader(open(f"{ROOT}/results/qfo_final/summary.tsv"), delimiter="\t"))
def val(s): return float(s.split("±")[0].replace(",", "")) if s else None
checks = [("SwissTree TPR", ("SwissTree", "TPR", "PPV"), "x"), ("SwissTree PPV", ("SwissTree", "TPR", "PPV"), "y"),
          ("TreeFam-A TPR", ("TreeFam-A", "TPR", "PPV"), "x"), ("TreeFam-A PPV", ("TreeFam-A", "TPR", "PPV"), "y"),
          ("VGNC TPR", ("VGNC", "TPR", "PPV"), "x"), ("VGNC PPV", ("VGNC", "TPR", "PPV"), "y"),
          ("GO avg Schlicker", ("GOtest", "NR_ORTHOLOGS", "avg Schlicker"), "y"), ("GO n", ("GOtest", "NR_ORTHOLOGS", "avg Schlicker"), "x"),
          ("EC avg Schlicker", ("ECtest", "NR_ORTHOLOGS", "avg Schlicker"), "y"), ("EC n", ("ECtest", "NR_ORTHOLOGS", "avg Schlicker"), "x"),
          ("FAS", ("FAS", "NR_ORTHOLOGS", "FAS"), "y")]
for b in ("STD_Eukaryota", "STD_Fungi", "STD_Bacteria", "G_STD2_Luca", "G_STD2_Eukaryota", "G_STD2_Fungi", "G_STD2_Vertebrata"):
    checks += [(f"{b} frac incorrect", (b, "NR_COMPLETED_TREE_SAMPLINGS", "FRAC_INCORRECT_TREES"), "y"),
               (f"{b} RF", (b, "NR_COMPLETED_TREE_SAMPLINGS", "RF_DISTANCE"), "y"),
               (f"{b} n trees", (b, "NR_COMPLETED_TREE_SAMPLINGS", "FRAC_INCORRECT_TREES"), "x")]
n = 0; bad = []; missing = []; tools = []
for row in summ:
    if row["kind"].startswith("ours"): continue
    tools.append(row["method"])
    for col, key, ax in checks:
        r = R.get(key + (row["method"],))
        s = val(row.get(col, ""))
        if r is None:
            if s is not None: missing.append((row["method"], col))
            continue
        if s is None: missing.append((row["method"], col, "blank in summary")); continue
        x = float(r[ax]); n += 1
        tol = 0.5 if col.endswith(" n") or col.endswith("n trees") else 6e-4
        if abs(x - s) > tol: bad.append((row["method"], col, s, x))
ref_tools = sorted({k[3] for k in R})
out = dict(n_compared=n, mismatches=bad, missing=missing, tools_in_summary=tools,
           ref_tools_not_in_summary=[t for t in ref_tools if t not in tools])
print(json.dumps(out, indent=1))
json.dump(out, open(f"{ROOT}/results/verification/qfo_competitors.json", "w"), indent=1)
