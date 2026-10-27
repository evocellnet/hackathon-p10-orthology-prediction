"""Arm table + adoption rule for the S0 rescue grid -> results/rescue/summary.{json,md}.
Rule (PI): vs R0, TreeFam-A dev TPR up with PPV drop <= 0.01 AND OrthoBench dev OG-F >= -1 AND
OrthoBench dev macro pair F1 >= -0.005."""
import csv, json, os, shutil
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results/rescue")


def last(path):
    out = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter="\t"):
            out[r["label"]] = r
    return out


ob, qf = last(os.path.join(RES, "ob_runs.tsv")), last(os.path.join(RES, "qfo_runs.tsv"))
arms = ["R0", "R1", "R2", "R3", "R4", "R2n", "R4n"]
rows = []
# merge.py changed at 02:08 (another agent's nesting fix) while the grid ran: QfO R1-R4 ran on the old
# merge.py and are compared with the old-code fast-preset baseline (results/core_v0_qfo/runs.tsv "B0 fast (bug
# repro)" = A8check, same code), R0/R2n/R4n on the new merge.py. OrthoBench R0 is identical on both codes
# (per-RefOG), and OB R*n reruns are reported too.
OLD_QFO = None
for r in csv.DictReader(open(os.path.join(ROOT, "results/core_v0_qfo/runs.tsv")), delimiter="\t"):
    if r["label"].startswith("B0 fast"):
        OLD_QFO = {"TPR": r["TF_dev_TPR"], "PPV": r["TF_dev_PPV"]}
b_ob = ob.get("R0")
for a in arms:
    b_qf = OLD_QFO if a in ("R1", "R2", "R3", "R4") else qf.get("R0")
    o, q = ob.get(a), qf.get(a)
    r = dict(arm=a)
    if q:
        r.update(TPR=float(q["TPR"]), PPV=float(q["PPV"]), qfo_rescue_s=q["rescue_s"], qfo_merge_s=q["merge_s"],
                 s0_fn=json.loads(q["fn_classes"]).get("i_diff_family"), fn_classes=json.loads(q["fn_classes"]),
                 qfo_attach=json.loads(q["rescue_attach"]) if q["rescue_attach"] not in ("None", "") else None,
                 qfo_guard=json.loads(q["guard"]) if q["guard"] not in ("None", "") else None)
    if o:
        r.update(OG_F=float(o["F"]), OG_P=float(o["P"]), OG_R=float(o["R"]), pairF1_macro=float(o["pairF1_macro"]),
                 pairF1=float(o["pairF1"]), ob_rescue_s=o["rescue_s"], ob_merge_s=o["merge_s"],
                 ob_attach=json.loads(o["rescue_stats"]) if o["rescue_stats"] not in ("None", "") else None,
                 ob_guard=json.loads(o["guard"]) if o["guard"] not in ("None", "") else None)
    r["qfo_code"] = "old merge.py (pre 02:08)" if a in ("R1", "R2", "R3", "R4") else "current merge.py"
    if a != "R0" and b_ob and b_qf and o and q:
        d = dict(dTPR=r["TPR"] - float(b_qf["TPR"]), dPPV=r["PPV"] - float(b_qf["PPV"]),
                 dOG_F=r["OG_F"] - float(b_ob["F"]), dmacro=r["pairF1_macro"] - float(b_ob["pairF1_macro"]))
        d = {k: round(v, 4) for k, v in d.items()}
        r.update(d, passes=bool(d["dTPR"] > 0 and d["dPPV"] >= -0.01 and d["dOG_F"] >= -1 and d["dmacro"] >= -0.005))
    rows.append(r)
json.dump(rows, open(os.path.join(RES, "summary.json"), "w"), indent=1)
L = ["| arm | QfO code | TF-A dev TPR | PPV | S0-split FN | OB OG-F (P/R) | OB macro pair F1 | ΔTPR | ΔPPV | ΔOG-F | Δmacro | rule |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    f = lambda k, fmt="{:.4f}": (fmt.format(r[k]) if k in r and r[k] is not None else "–")
    L.append(f"| {r['arm']} | {r['qfo_code']} | {f('TPR')} | {f('PPV')} | {r.get('s0_fn', '–')} | {f('OG_F', '{:.2f}')} ({f('OG_P', '{:.1f}')}/{f('OG_R', '{:.1f}')}) | "
             f"{f('pairF1_macro')} | {f('dTPR', '{:+.4f}')} | {f('dPPV', '{:+.4f}')} | {f('dOG_F', '{:+.2f}')} | {f('dmacro', '{:+.4f}')} | "
             f"{'pass' if r.get('passes') else ('–' if 'passes' not in r else 'fail')} |")
open(os.path.join(RES, "summary.md"), "w").write("\n".join(L) + "\n")
for bench in ("ob", "qfo"):
    for a in arms:
        for fn in ("dev_eval_rescue.json", "tf_dev_rescue.json", "manifest.json"):
            p = os.path.join(ROOT, "work/rescue", bench, a, fn)
            if os.path.exists(p):
                shutil.copyfile(p, os.path.join(RES, f"{bench}_{a}_{fn}"))
print("\n".join(L))
