#!/usr/bin/env python3
"""Build results/orthobench_competitors/SUMMARY.md from <method>/score.json + runtime.json (+ paired bootstrap)."""
import json, os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orthobench_eval import bootstrap_delta_f
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "orthobench_competitors")
split = json.load(open(os.path.join(C, "..", "..", "..", "splits", "orthobench_split.json")))
order = ["orthofinder", "orthofinder_mcl", "sonicparanoid_default", "sonicparanoid_sensitive", "proteinortho", "fastoma_luca", "fastoma_luca_ogs"]
M = {}
for m in order + sorted(set(os.path.basename(d.rstrip("/")) for d in glob.glob(C + "/*/")) - set(order)):
    if os.path.exists(f"{C}/{m}/score.json"):
        M[m] = (json.load(open(f"{C}/{m}/score.json")), json.load(open(f"{C}/{m}/runtime.json")))
L = ["# OrthoBench competitor baselines (12 species, 251,378 proteins, 70 RefOGs)", "",
     "Scores: official Open_Orthobench `benchmark.py` (pairwise, RefOG-size normalised; percent). dev/test = sum of per-RefOG TP/FP/FN over the",
     "34/36 RefOGs of `splits/orthobench_split.json` (dev+test reproduces the all-70 numbers; verified identical to benchmark.py). Runs: 8 threads,",
     "`/usr/bin/time -v`, SHARED node (24 cores, load avg 20-55 during runs) -> wall times are inflated vs an idle node; CPU h = user+sys.", "",
     "| method | F all | F dev | F test | precision | recall | exact RefOGs | wall (h) | CPU (h) | peak RSS (GB) |", "|---|---|---|---|---|---|---|---|---|---|"]
for m, (s, r) in M.items():
    a = s["all"]
    L.append(f"| {m} | {a['F']:.1f} | {s['dev']['F']:.1f} | {s['test']['F']:.1f} | {a['precision']:.1f} | {a['recall']:.1f} | {a['exact_refogs']}/70 | "
             f"{r['wall_h']:.2f} | {r['cpu_h']:.1f} | {r['max_rss_gb']:.1f} |")
L += ["", "Precision/recall on dev / test:", "", "| method | P dev | R dev | P test | R test |", "|---|---|---|---|---|"]
for m, (s, r) in M.items():
    L.append(f"| {m} | {s['dev']['precision']:.1f} | {s['dev']['recall']:.1f} | {s['test']['precision']:.1f} | {s['test']['recall']:.1f} |")
if "orthofinder" in M:
    L += ["", "## Paired bootstrap vs OrthoFinder (dF = F(method) - F(OrthoFinder), percentage points; RefOGs resampled, B=2000, seed 11)", "",
          "| method | subset | dF | 95% CI | SE |", "|---|---|---|---|---|"]
    base = M["orthofinder"][0]["per_refog"]
    for m, (s, r) in M.items():
        if m == "orthofinder": continue
        for sub in ("all", "dev", "test"):
            ids = None if sub == "all" else split[sub]
            b = bootstrap_delta_f(s["per_refog"], base, ids, 2000, 11)
            L.append(f"| {m} | {sub} | {b['delta_point']:+.1f} | [{b['ci95'][0]:+.1f}, {b['ci95'][1]:+.1f}] | {b['se']:.1f} |")
L += ["", "## Notes", ""]
for m in M:
    v = f"{C}/{m}/version.txt"
    L.append(f"- **{m}**: " + (open(v).read().strip() if os.path.exists(v) else "") + " | " + M[m][1].get("note", ""))
nf = f"{C}/NOTES_extra.md"
if os.path.exists(nf): L += ["", open(nf).read().strip()]
open(C + "/SUMMARY.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
