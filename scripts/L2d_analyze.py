"""L2-domains: keep-rule analysis of the merge arms (docs/L2_domains.md §0.4).

    python scripts/L2d_analyze.py results/l2_domains/arms -> results/l2_domains/integration.json, tables.md

QfO: primary TreeFam-A dev pair F1 = 2TP/(2TP+FP+FN), paired bootstrap over dev families (B=2000, seed 11).
OrthoBench: OG-F guard; micro pair F1 with paired bootstrap over RefOGs; macro pair F1.
"""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys

import numpy as np

SEEDS = ("11", "23", "47")
WT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def f1(c):
    tp, fn, fp = c[..., 0], c[..., 1], c[..., 2]
    with np.errstate(invalid="ignore", divide="ignore"):
        return 2 * tp / (2 * tp + fp + fn)


def ppv(c):
    return c[..., 0] / (c[..., 0] + c[..., 2])


def tpr(c):
    return c[..., 0] / (c[..., 0] + c[..., 1])


def paired_boot(A, Bm, fn, B=2000, seed=11):
    """A, Bm: [units, 4] count arrays (same unit order). Delta = fn(sum B) - fn(sum A)."""
    rng = np.random.default_rng(seed)
    n = len(A)
    d = np.empty(B)
    for i in range(B):
        ii = rng.integers(0, n, n)
        d[i] = fn(Bm[ii].sum(0)) - fn(A[ii].sum(0))
    return [round(float(np.percentile(d, 2.5)), 4), round(float(np.percentile(d, 97.5)), 4)]


def macro(c):
    v = f1(c)
    return float(np.nanmean(v))


def main(d):
    runs = {}
    for p in sorted(glob.glob(os.path.join(d, "*.json"))):
        r = json.load(open(p))
        runs[(r["bench"], r["arm"], r["shuffle"])] = r
    out = {"qfo": {}, "ob": {}}
    rows = []
    for bench in ("qfo", "ob"):
        off = runs.get((bench, "off", "none"))
        if off is None:
            continue
        key = "per_family" if bench == "qfo" else "per_refog"
        units = sorted(off[key])
        C0 = np.array([off[key][u] for u in units], float)
        for arm in ("bonus", "gate", "gate_deep", "veto"):
            real = runs.get((bench, arm, "none"))
            if real is None:
                continue
            Cr = np.array([real[key][u] for u in units], float)
            sh = [runs.get((bench, arm, s)) for s in SEEDS]
            sh = [s for s in sh if s is not None]
            Cs = [np.array([s[key][u] for u in units], float) for s in sh]
            r = dict(n_shuffled=len(Cs), stats=real["stats"],
                     shuffled_stats=[s["stats"] for s in sh])
            r["F1_off"], r["F1_real"] = round(float(f1(C0.sum(0))), 4), round(float(f1(Cr.sum(0))), 4)
            r["dF1"] = round(r["F1_real"] - r["F1_off"], 4)
            r["dF1_ci"] = paired_boot(C0, Cr, f1)
            shv = [float(f1(c.sum(0))) for c in Cs]
            r["F1_shuf"] = [round(x, 4) for x in shv]
            r["seed_spread"] = round(max(shv) - min(shv), 4) if shv else None
            r["real_minus_shuf"] = round(r["F1_real"] - float(np.mean(shv)), 4) if shv else None
            r["dF1_shuf_vs_off_ci"] = [paired_boot(C0, c, f1) for c in Cs]
            r["PPV_off"], r["PPV_real"] = round(float(ppv(C0.sum(0))), 4), round(float(ppv(Cr.sum(0))), 4)
            r["TPR_off"], r["TPR_real"] = round(float(tpr(C0.sum(0))), 4), round(float(tpr(Cr.sum(0))), 4)
            r["dPPV"], r["dTPR"] = round(r["PPV_real"] - r["PPV_off"], 4), round(r["TPR_real"] - r["TPR_off"], 4)
            r["dPPV_ci"] = paired_boot(C0, Cr, ppv)
            r["dTPR_ci"] = paired_boot(C0, Cr, tpr)
            if bench == "qfo":
                r["TPR_official"] = round(real["TF_TPR_official"], 4)
                r["PPV_official"] = round(real["TF_PPV_official"], 4)
                r["TPR_official_off"] = round(off["TF_TPR_official"], 4)
                r["PPV_official_off"] = round(off["TF_PPV_official"], 4)
                r["primary_ci_pass"] = r["dF1_ci"][0] > 0
                r["shuffle_pass"] = (r["real_minus_shuf"] is not None and len(Cs) == 3
                                     and r["real_minus_shuf"] > r["seed_spread"])
                r["ppv_guard_pass"] = r["dPPV"] >= -0.01
            else:
                r["OG_F_off"], r["OG_F_real"] = round(off["OG_F"], 2), round(real["OG_F"], 2)
                r["dOG_F"] = round(real["OG_F"] - off["OG_F"], 2)
                r["OG_F_shuf"] = [round(s["OG_F"], 2) for s in sh]
                r["macroF1_off"], r["macroF1_real"] = round(macro(C0), 4), round(macro(Cr), 4)
                r["macroF1_shuf"] = [round(macro(c), 4) for c in Cs]
                r["og_guard_pass"] = r["dOG_F"] >= -1.0
            out[bench][arm] = r
            rows.append((bench, arm, r))
    # keep verdict per arm (needs both benches)
    verdict = {}
    for arm in ("bonus", "gate", "gate_deep", "veto"):
        q = out["qfo"].get(arm)
        o = out["ob"].get(arm)
        if q is None:
            continue
        ok = q["primary_ci_pass"] and q["shuffle_pass"] and q["ppv_guard_pass"]
        ob_ok = True if arm == "gate_deep" and o is None else (o is not None and o["og_guard_pass"])
        verdict[arm] = dict(primary_ci=q["primary_ci_pass"], shuffle=q["shuffle_pass"],
                            ppv_guard=q["ppv_guard_pass"], ob_og_guard=ob_ok,
                            keep_before_orthoxml_check=bool(ok and ob_ok))
    out["verdict"] = verdict
    git = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    out["manifest"] = dict(git_commit=git, script="scripts/L2d_analyze.py", inputs=sorted(
        os.path.basename(p) for p in glob.glob(os.path.join(d, "*.json"))), boot=dict(B=2000, seed=11))
    od = os.path.dirname(os.path.abspath(d))
    json.dump(out, open(os.path.join(od, "integration.json"), "w"), indent=1)
    # markdown table
    L = ["| bench | arm | F1 off | F1 real | ΔF1 [95% CI] | F1 shuffled (11/23/47) | real−shuf | seed spread | "
         "ΔTPR | ΔPPV | OG-F off→real | macro F1 off→real | events (real) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for bench, arm, r in rows:
        ev = {k: v for k, v in r["stats"].items() if k.startswith("arch_")}
        og = f"{r['OG_F_off']}→{r['OG_F_real']}" if bench == "ob" else "–"
        mc = f"{r['macroF1_off']}→{r['macroF1_real']}" if bench == "ob" else "–"
        L.append(f"| {bench} | {arm} | {r['F1_off']} | {r['F1_real']} | {r['dF1']:+.4f} [{r['dF1_ci'][0]:+.4f}, "
                 f"{r['dF1_ci'][1]:+.4f}] | {', '.join(str(x) for x in r['F1_shuf'])} | {r['real_minus_shuf']} | "
                 f"{r['seed_spread']} | {r['dTPR']:+.4f} | {r['dPPV']:+.4f} | {og} | {mc} | {ev} |")
    open(os.path.join(od, "integration_table.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    print(json.dumps(verdict, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
