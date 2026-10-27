"""G-grid table + adoption rule -> results/giant/summary.{md,json} (from results/giant/runs.tsv, last row per arm)."""
import json
rows = [json.loads(l) for l in open("results/giant/runs.tsv")]
q = {r["label"]: r for r in rows if r["bench"] == "qfo"}
o = {r["label"]: r for r in rows if r["bench"] == "ob"}
L = ["| arm | QfO TF-A dev TPR / PPV (fast) | ΔTPR / ΔPPV | FN: diff family / diff root HOG / paralog call | "
     "OB OG-F (P/R) | OB macro pair F1 | ΔOG-F / Δmacro | rule |", "|---|---|---|---|---|---|---|---|"]
out = []
for a in ("G0", "G1", "G2", "G3"):
    Q, O = q[a], o[a]
    fc = Q["fn_classes"]
    dT, dP = Q["TPR"] - q["G0"]["TPR"], Q["PPV"] - q["G0"]["PPV"]
    dF, dm = O["OG_F"] - o["G0"]["OG_F"], O["pairF1_macro"] - o["G0"]["pairF1_macro"]
    ok = a != "G0" and dT > 0 and dP >= -0.01 and dF >= -1 and dm >= -0.005
    out.append(dict(arm=a, qfo=Q, ob=O, dTPR=round(dT, 4), dPPV=round(dP, 4), dOG_F=round(dF, 2), dmacro=round(dm, 4),
                    passes=None if a == "G0" else ok))
    dq = "–" if a == "G0" else f"{dT:+.4f} / {dP:+.4f}"
    do = "–" if a == "G0" else f"{dF:+.2f} / {dm:+.4f}"
    rule = "–" if a == "G0" else ("pass" if ok else "fail")
    L.append(f"| {a} | {Q['TPR']:.4f} / {Q['PPV']:.4f} | {dq} | {fc['i_diff_family']} / {fc['ii_diff_root_hog']} / "
             f"{fc['iii_paralog_call']} | {O['OG_F']:.2f} ({O['OG_P']:.1f}/{O['OG_R']:.1f}) | {O['pairF1_macro']:.4f} | {do} | {rule} |")
open("results/giant/summary.md", "w").write("\n".join(L) + "\n")
json.dump(out, open("results/giant/summary.json", "w"), indent=1)
print("\n".join(L))
