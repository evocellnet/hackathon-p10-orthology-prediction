#!/usr/bin/env python3
"""Summarise results/S0_families/*.json into a markdown table + apply the pre-stated selection rule.

Rule (results/S0_families/configs.json "_rule"): among configs with sum_g2 <= 0.05 N^2 and
max_family <= 0.03 N, maximise DEV completeness; ties -> higher dev pair_recall, then lower wall.
Writes results/S0_families/summary.tsv and results/S0_families/selection.json; prints the markdown table.
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results/S0_families")


def load():
    rows = []
    for f in sorted(os.listdir(RES)):
        if not f.endswith(".json") or f.endswith(".TEST.json") or f in ("configs.json", "selection.json"):
            continue
        r = json.load(open(os.path.join(RES, f)))
        a, fs, rt = r["dev"]["aggregate"], r["family_stats"], r["runtime"]
        n = fs["n_genes"]
        rows.append({
            "name": r["name"],
            "dev_completeness": a["completeness"], "dev_n_complete": a["n_complete"], "dev_n": a["n_refogs"],
            "dev_complete_keep_uncertain": r["dev_keep_uncertain_aggregate"]["completeness"],
            "mean_fission": a["mean_n_fam"], "mean_best_frac": a["mean_best_frac"],
            "median_size_ratio": a["median_size_ratio"], "max_size_ratio": a["max_size_ratio"],
            "fusion_rate": a["fusion_rate"], "pair_recall": a["pair_recall"], "pair_precision": a["pair_precision"],
            "n_families": fs["n_families"], "n_singletons": fs["n_singletons"],
            "n_stage1": fs.get("n_stage1_clusters"), "n_reps": fs.get("n_reps"),
            "sum_g2": fs["sum_g2"], "sum_g2_over_N2": fs["sum_g2_over_N2"],
            "max_family": fs["max_family"], "max_over_N": fs["max_family"] / n, "top10": fs["top_sizes"],
            "n_split_families": fs.get("n_split_families", 0),
            "wall_s": rt.get("wall_incl_reused_stage1_s", rt.get("wall_s")), "peak_rss_gb": rt.get("peak_rss_gb"),
            "cpu_h": rt.get("cpu_h"), "load": rt.get("node_loadavg_at_start", ["?"])[0],
            "eligible": fs["sum_g2"] <= 0.05 * n * n and fs["max_family"] <= 0.03 * n,
        })
    return rows


def main():
    rows = load()
    elig = [r for r in rows if r["eligible"]]
    best = max(elig, key=lambda r: (r["dev_completeness"], r["pair_recall"], -r["wall_s"])) if elig else None
    cols = list(rows[0].keys())
    with open(os.path.join(RES, "summary.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(json.dumps(r[c]) if isinstance(r[c], list) else str(r[c]) for c in cols) + "\n")
    if best:
        json.dump({"rule": json.load(open(os.path.join(RES, "configs.json")))["_rule"], "selected": best["name"],
                   "selected_row": best}, open(os.path.join(RES, "selection.json"), "w"), indent=1)
    print("| config | dev complete | fission (mean #fam) | best frac | fusion rate | size ratio med/max | pair R/P | "
          "Σg²/N² | max fam (%N) | singletons | wall s (RSS GB) | eligible |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda r: -r["dev_completeness"]):
        mark = " **(selected)**" if best and r["name"] == best["name"] else ""
        print(f"| {r['name']}{mark} | {r['dev_completeness']:.3f} ({r['dev_n_complete']}/{r['dev_n']}) | "
              f"{r['mean_fission']:.2f} | {r['mean_best_frac']:.3f} | {r['fusion_rate']:.3f} | "
              f"{r['median_size_ratio']:.2f}/{r['max_size_ratio']:.1f} | {r['pair_recall']:.3f}/{r['pair_precision']:.3f} | "
              f"{r['sum_g2_over_N2']:.4f} | {r['max_family']} ({100 * r['max_over_N']:.1f}%) | {r['n_singletons']} | "
              f"{r['wall_s']:.0f} ({r['peak_rss_gb']:.1f}) | {'yes' if r['eligible'] else 'no'} |")
    if best:
        print(f"\nselected: {best['name']}")


if __name__ == "__main__":
    sys.exit(main())
