"""Per-RefOG OG-level diff of two OrthoBench dev evals (dev_eval json or dev_eval_rescue json)."""
import json, sys
def per(p):
    d = json.load(open(p))
    rows = d["per_refog"] if "per_refog" in d else list(d["per_refog_og"].values())
    return {r["refog"]: r for r in rows}
a, b = per(sys.argv[1]), per(sys.argv[2])
for k in sorted(a):
    x, y = a[k], b[k]
    if (x["tp"], x["fp"], x["fn"]) != (y["tp"], y["fp"], y["fn"]):
        print(k, "n", x["n"], "| tp/fp/fn", x["tp"], x["fp"], x["fn"], "->", y["tp"], y["fp"], y["fn"],
              "| fp genes other", x["fp_genes_other"], "->", y["fp_genes_other"], "| groups", x["n_groups"], "->", y["n_groups"])
