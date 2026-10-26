"""Tables for docs/S6_structure.md from results/S6_structure/{s0_rescue.json, merge_runs.jsonl}.
    python scripts/S6_summarize.py  -> results/S6_structure/tables.md"""
import json, os
import sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results/S6_structure")
SH = ["shuf11", "shuf23", "shuf47"]
out = []
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from orthobench_dev_pairs_eval import paired_bootstrap  # noqa: E402
OFF_ARR = np.load(os.path.join(ROOT, "work/S6/merge/off_-/pair_counts.npy"))


def delta(r):
    arr = np.load(os.path.join(ROOT, "work/S6/merge", f"{r['tag']}{r['arm'].replace(':', '_')}_{r['variant']}", "pair_counts.npy"))
    return paired_bootstrap(OFF_ARR, arr)["all"]

# ---- exp 1
d = json.load(open(os.path.join(R, "s0_rescue.json")))
rows = d["rows"]
off = rows[0]
cfgs = []
for r in rows[1:]:
    if r["config"] not in cfgs:
        cfgs.append(r["config"])
out.append("## Exp 1 — S0 rescue (DEV RefOGs, n=34)\n")
out.append(f"off: complete {off['n_complete']}/34 = {off['completeness']:.3f}, fission {off['mean_n_fam']:.3f}, "
           f"fusion rate {off['fusion_rate']:.3f}, Σg²/N² {off['sum_g2_over_N2']:.4f}, families {off['n_families']}, singletons {off['n_singletons']}\n")
out.append("| config | pre-reg | real complete | shuffled complete (11/23/47) | Δ real−shuf mean | spread | keep | real fission | real fusion rate (shuf mean) | real Σg²/N² (shuf mean) | real joins (shuf mean) |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|")
for c in cfgs:
    by = {r["variant"]: r for r in rows[1:] if r["config"] == c}
    re_ = by["real"]
    sh = [by[s] for s in SH if s in by]
    sc = [x["completeness"] for x in sh]
    spread = max(sc) - min(sc)
    dsh = re_["completeness"] - np.mean(sc)
    doff = re_["completeness"] - off["completeness"]
    keep = "yes" if dsh > 0 and doff > spread else "no"
    pre = "no (post hoc)" if ("mode" in c or "cov_min" in c or '"small_max": 1,' in c) else "yes"
    out.append(f"| `{c}` | {pre} | {re_['n_complete']} ({re_['completeness']:.3f}) | {'/'.join(str(x['n_complete']) for x in sh)} | "
               f"{dsh:+.3f} | {spread:.3f} | {keep} | {re_['mean_n_fam']:.3f} | {re_['fusion_rate']:.3f} ({np.mean([x['fusion_rate'] for x in sh]):.3f}) | "
               f"{re_['sum_g2_over_N2']:.4f} ({np.mean([x['sum_g2_over_N2'] for x in sh]):.4f}) | {re_['link']['n_joins']} ({np.mean([x['link']['n_joins'] for x in sh]):.0f}) |")
out.append("")

# ---- exp 2/3
runs = [json.loads(l) for l in open(os.path.join(R, "merge_runs.jsonl"))]
base = [r for r in runs if r["arm"] == "off" and not r["tag"]][-1]
arms = []
for r in runs:
    if r["arm"] != "off" and r["arm"] not in arms:
        arms.append(r["arm"])
out.append("## Exp 2 — structure-aware scores in the merge (DEV; pair metrics on 34 RefOGs)\n")
out.append(f"off: OG-F {base['F']}, pair P/R/F1 {base['pairP']}/{base['pairR']}/{base['pairF1']}, deep pair F1 "
           f"{base['by_depth']['deep']['F1']}, shallow {base['by_depth']['shallow']['F1']}\n")
out.append("| arm | 3Di | OG-F | pair P | pair R | pair F1 | ΔF1 vs off [95% CI] | deep F1 | shallow F1 |")
out.append("|---|---|---|---|---|---|---|---|---|")
last = {}
for r in runs:
    if r["arm"] != "off":
        last[(r["arm"], r["variant"])] = r
for a in arms:
    for v in ["real"] + SH:
        r = last.get((a, v))
        if not r:
            continue
        dl = (r.get("delta_vs_off") or delta(r))["F1"]
        out.append(f"| {a} | {v} | {r['F']} | {r['pairP']} | {r['pairR']} | {r['pairF1']} | {dl['delta']:+.4f} [{dl['ci95'][0]:+.3f}, {dl['ci95'][1]:+.3f}] | "
                   f"{r['by_depth']['deep']['F1']} | {r['by_depth']['shallow']['F1']} |")
out.append("")
out.append("Keep rule on pair F1 (DESIGN §5): Δ(real − mean shuffled) > 0 and Δ(real − off) > seed spread (max − min of the 3 shuffled runs).\n")
out.append("| arm | real F1 | shuffled mean | spread | Δ real−shuf | Δ real−off | keep | deep: Δ real−off / Δ real−shuf | shallow: Δ real−off / Δ real−shuf |")
out.append("|---|---|---|---|---|---|---|---|---|")
for a in arms:
    re_ = last.get((a, "real"))
    sh = [last[(a, s)] for s in SH if (a, s) in last]
    if not re_ or not sh:
        continue
    f = [x["pairF1"] for x in sh]
    dsh, doff, spread = re_["pairF1"] - np.mean(f), re_["pairF1"] - base["pairF1"], max(f) - min(f)
    dd = lambda k: (re_["by_depth"][k]["F1"] - base["by_depth"][k]["F1"], re_["by_depth"][k]["F1"] - np.mean([x["by_depth"][k]["F1"] for x in sh]))
    out.append(f"| {a} | {re_['pairF1']} | {np.mean(f):.4f} | {spread:.4f} | {dsh:+.4f} | {doff:+.4f} | {'yes' if dsh > 0 and doff > spread else 'no'} | "
               f"{dd('deep')[0]:+.4f} / {dd('deep')[1]:+.4f} | {dd('shallow')[0]:+.4f} / {dd('shallow')[1]:+.4f} |")
out.append("")
out.append("## Exp 3 — pair F1 by 3Di coverage of the two genes (n = truth pairs per subset)\n")
nb = {k: base["by_coverage"][k]["tp"] + base["by_coverage"][k]["fn"] + base["by_coverage"][k]["fp"] + base["by_coverage"][k]["tn"] for k in base["by_coverage"]}
out.append(f"Truth pairs: both 3Di {nb['both_3di']}, one {nb['one_3di']}, none {nb['no_3di']}.\n")
out.append("| arm | 3Di | both F1 (Δoff) | one F1 (Δoff) | none F1 (Δoff) |")
out.append("|---|---|---|---|---|")
out.append(f"| off | - | {base['by_coverage']['both_3di']['F1']} | {base['by_coverage']['one_3di']['F1']} | {base['by_coverage']['no_3di']['F1']} |")
for a in arms:
    for v in ["real"] + SH:
        r = last.get((a, v))
        if not r:
            continue
        c = r["by_coverage"]; b = base["by_coverage"]
        out.append(f"| {a} | {v} | " + " | ".join(f"{c[k]['F1']} ({c[k]['F1'] - b[k]['F1']:+.4f})" for k in ("both_3di", "one_3di", "no_3di")) + " |")
out.append("")
down = [r for r in runs if r["tag"].startswith("s0down")]
if down:
    out.append("## Exp 1 downstream — merge on rescued families (attach e≤1e-3, ≤10 genes)\n")
    out.append("| families | OG-F | OG P | OG R | exact | pair P | pair R | pair F1 | deep F1 |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    out.append(f"| S0 + core rescue (off) | {base['F']} | {base['P_og']} | {base['R_og']} | {base['exact']} | {base['pairP']} | {base['pairR']} | {base['pairF1']} | {base['by_depth']['deep']['F1']} |")
    for r in down:
        out.append(f"| + 3Di attach, {r['tag'].split('_')[1]} | {r['F']} | {r['P_og']} | {r['R_og']} | {r['exact']} | {r['pairP']} | {r['pairR']} | {r['pairF1']} | {r['by_depth']['deep']['F1']} |")
open(os.path.join(R, "tables.md"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
