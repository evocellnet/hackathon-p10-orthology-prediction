#!/usr/bin/env python
"""L2-gpu task 1 summary: evaluates every arm against the CPU frozen path (OrthoBench accurate+go,
QfO fast+go) and applies the pre-registered equivalence rule E1-E3 (docs/L2_gpu.md §0.1).
Writes results/l2_gpu/task1_summary.json and prints a markdown table.
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W = os.path.join(ROOT, "..")          # work/
REF_OB = os.path.join(W, "ob_final/main_accurate_go")
REF_QFO = os.path.join(W, "gpuv2/qfo_cpu_fast")
ARMS = [  # name, OB run, QfO run, kind
    ("CPU no rescue (A)", "gpuv2/diagA_cpu_norescue", "l2gpu/qfo/A", "reference"),
    ("GPU split (gpu_max_family 1000)", "gpuv2/e2e_gpu_acc", "gpuv2/qfo_gpusw_fast", "existing"),
    ("all-GPU, unsparsified", "gpuv2/ob_gpusw_acc", "gpuv2/qfo_allgpu_fast", "existing"),
    ("K39 top-39/species (either)", "l2gpu/ob/K39", "l2gpu/qfo/K39", "arm"),
    ("K10 top-10/species (either)", "l2gpu/ob/K10", "l2gpu/qfo/K10", "arm"),
    ("R207 bits >= 0.207 best/species", "l2gpu/ob/R207", "l2gpu/qfo/R207", "arm"),
    ("F0528 norm_max >= 0.0528", "l2gpu/ob/F0528", "l2gpu/qfo/F0528", "arm"),
    ("Crand random 57.5% (control)", "l2gpu/ob/Crand", "l2gpu/qfo/Crand", "control"),
    ("MK39 mutual top-39 (D1)", "l2gpu/ob/MK39", "l2gpu/qfo/MK39", "arm"),
    ("SP3000 splice (D1)", "l2gpu/ob/SP3000", "l2gpu/qfo/SP3000", "arm"),
    ("SP5000 splice (D1)", "l2gpu/ob/SP5000", "l2gpu/qfo/SP5000", "arm"),
    ("J11 CPU A + 1% jitter (D2)", "l2gpu/ob/J11", "l2gpu/qfo/J11", "noise"),
    ("J23 CPU A + 1% jitter (D2)", "l2gpu/ob/J23", None, "noise"),
    ("J47 CPU A + 1% jitter (D2)", "l2gpu/ob/J47", None, "noise"),
    ("X1 GPU scores on DIAMOND pair set", "l2gpu/ob/X1", None, "diagnostic"),
    ("X2 GPU table + CPU self-scores", "l2gpu/ob/X2", None, "diagnostic"),
    ("X3 CPU table + GPU self-scores", "l2gpu/ob/X3", None, "diagnostic"),
    ("Dmax_cpuA max-link, CPU A", "l2gpu/ob/Dmax_cpuA", None, "diagnostic"),
    ("Dmax_gpu max-link, all-GPU", "l2gpu/ob/Dmax_gpu", None, "diagnostic"),
]


def done(p):
    return p and os.path.exists(os.path.join(W, p, "hogs.orthoxml"))


def main():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from gpu_dev_eval import boot, one
    from orthobench_dev_pairs_eval import load_truth
    truth = load_truth()
    a = one(REF_OB, truth)
    a0 = one(os.path.join(W, "gpuv2/diagA_cpu_norescue"), truth)   # same-rescue (none) CPU reference
    ob, obA = {}, {}
    for name, o, q, kind in ARMS:
        if done(o):
            b = one(os.path.join(W, o), truth)
            d = {k: v for k, v in b.items() if not k.startswith("_")}
            d.update(boot(a, b))
            ob[name] = d
            obA[name] = boot(a0, b)
    qruns = [(name, os.path.join(W, q)) for name, o, q, kind in ARMS if done(q)]
    tmp = os.path.join(ROOT, "results/l2_gpu/.tf_tmp.json")
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts/l2gpu_treefam_paired.py"), REF_QFO,
                    *[p for _, p in qruns], "--json", tmp], check=True, capture_output=True)
    tf = json.load(open(tmp))
    os.remove(tmp)
    qfo = {name: r for (name, _), r in zip(qruns, tf[1:])}
    rows = []
    print("| arm | kind | OB OG-F | ΔOG-F [95% CI] | OB macro pair F1 | Δmacro [95% CI] | TF-A dev TPR / PPV | ΔTPR [CI] | ΔPPV [CI] | E1 | E2 | E3 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for name, o, q, kind in ARMS:
        x, y = ob.get(name), qfo.get(name)
        e1 = e2 = e3 = None
        if x:
            e1 = bool(x["dOG_F"] >= -1.0 and x["dOG_F_ci"][0] >= -4.0)
            e2 = bool(x["dmacro"] >= -0.005 and x["dmacro_ci"][0] >= -0.02)
        if y:
            e3 = bool(y["dTPR"] >= -0.01 and y["dPPV"] >= -0.01 and y["dTPR_ci"][0] >= -0.02 and y["dPPV_ci"][0] >= -0.02)
        rows.append(dict(arm=name, kind=kind, ob_run=o, qfo_run=q, ob=x, ob_vs_cpuA=obA.get(name), qfo=y,
                         E1=e1, E2=e2, E3=e3))
        f = lambda v: "–" if v is None else ("pass" if v else "fail")  # noqa: E731
        obs = (f"{x['F']:.2f} | {x['dOG_F']:+.2f} [{x['dOG_F_ci'][0]:+.1f}, {x['dOG_F_ci'][1]:+.1f}] | "
               f"{x['pairF1_macro']:.4f} | {x['dmacro']:+.4f} [{x['dmacro_ci'][0]:+.3f}, {x['dmacro_ci'][1]:+.3f}]"
               if x else "– | – | – | –")
        qs = (f"{y['TPR']:.4f} / {y['PPV']:.4f} | {y['dTPR']:+.4f} [{y['dTPR_ci'][0]:+.3f}, {y['dTPR_ci'][1]:+.3f}] | "
              f"{y['dPPV']:+.4f} [{y['dPPV_ci'][0]:+.3f}, {y['dPPV_ci'][1]:+.3f}]" if y else "– | – | –")
        print(f"| {name} | {kind} | {obs} | {qs} | {f(e1)} | {f(e2)} | {f(e3)} |")
    print()
    print("| arm | ΔOG-F vs CPU A (no rescue) [CI] | Δmacro vs CPU A [CI] |")
    print("|---|---|---|")
    for name, o, q, kind in ARMS:
        z = obA.get(name)
        if z:
            print(f"| {name} | {z['dOG_F']:+.2f} [{z['dOG_F_ci'][0]:+.1f}, {z['dOG_F_ci'][1]:+.1f}] | "
                  f"{z['dmacro']:+.4f} [{z['dmacro_ci'][0]:+.3f}, {z['dmacro_ci'][1]:+.3f}] |")
    ref = {"ob": {k: v for k, v in a.items() if not k.startswith("_")}, "qfo": tf[0]}
    git = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    json.dump({"reference": ref, "rows": rows,
               "rule": "E1 dOG-F>=-1 & CI lo>=-4; E2 dmacro>=-0.005 & CI lo>=-0.02; E3 dTPR,dPPV>=-0.01 & CI lo>=-0.02",
               "manifest": {"git_commit": git, "ref_ob": REF_OB, "ref_qfo": REF_QFO,
                            "bootstrap": "paired, B=2000, seed 11 (OB over 34 dev RefOGs; TF-A over 359 dev families)"}},
              open(os.path.join(ROOT, "results/l2_gpu/task1_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
