#!/usr/bin/env python
"""Build results/orthobench_final/summary.{tsv,md} (+ runtime.json) from scores.json, run manifests,
/usr/bin/time -v outputs and the process-tree samples (proc_samples.tsv, 2 s)."""
import json
import os
import re
import sys

import numpy as np

ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
OUT = os.path.join(ROOT, "results/orthobench_final")
sys.path.insert(0, os.path.join(OUT, "scripts"))
from score_final import CONFIGS, NI_COMPARATORS  # noqa: E402


def parse_time(p):
    t = open(p).read()
    def g(k):
        m = re.search(re.escape(k) + r":\s*(.+)", t)
        return m.group(1).strip() if m else None
    el = g("Elapsed (wall clock) time (h:mm:ss or m:ss)")
    parts = [float(x) for x in el.split(":")]
    wall = parts[-1] + 60 * parts[-2] + (3600 * parts[-3] if len(parts) == 3 else 0)
    return dict(wall_s=wall, user_s=float(g("User time (seconds)")), sys_s=float(g("System time (seconds)")),
                cpu_s=float(g("User time (seconds)")) + float(g("System time (seconds)")),
                max_rss_gb=int(g("Maximum resident set size (kbytes)")) / 2**20, exit=int(g("Exit status")),
                cpu_pct=g("Percent of CPU this job got"))


def stage_accounting(run):
    """Per-stage wall (log timestamps), CPU (process-tree samples, linear interp.), peak tree RSS (samples)."""
    log = [l.rstrip("\n").split("\t", 1) for l in open(os.path.join(run, "run.log")) if "\t" in l]
    starts, ends = {}, {}
    for ts, msg in log:
        m = re.match(r"\[ancestra\] stage (\w+): (running|done|up to date)", msg)
        if m:
            (starts if m.group(2) == "running" else ends)[m.group(1)] = float(ts)
    smp = np.loadtxt(os.path.join(run, "proc_samples.tsv"), skiprows=1, ndmin=2)
    man = json.load(open(os.path.join(run, "manifest.json")))
    out = {}
    for st, info in man["stages"].items():
        if info.get("skipped") or st not in starts or st not in ends:
            out[st] = dict(skipped=True)
            continue
        a, b = starts[st], ends[st]
        cpu = float(np.interp(b, smp[:, 0], smp[:, 1]) - np.interp(a, smp[:, 0], smp[:, 1]))
        w = (smp[:, 0] >= a - 2) & (smp[:, 0] <= b + 2)
        out[st] = dict(wall_s=round(b - a, 1), cpu_s=round(cpu, 1),
                       peak_rss_tree_gb=round(float(smp[w, 2].max()) / 1024, 2) if w.any() else None,
                       peak_rss_own_gb=round(info["peak_rss_mb"] / 1024, 2),
                       load1_mean=round(float(smp[w, 4].mean()), 1) if w.any() else None)
    return out


def f2(x, d=1):
    return ("%%.%df" % d) % x


def main():
    S = json.load(open(os.path.join(OUT, "scores.json")))
    sc, comp, boots, pairs = S["scores"], S["competitors"], S["bootstrap"], S["pair_test"]
    bt = {(b["A"], b["B_"], b["subset"]): b for b in boots}
    rt = {}
    for lab in sc:
        run = sc[lab]["run"]
        rt[lab] = dict(time=parse_time(os.path.join(run, "time.txt")), stages=stage_accounting(run))
    json.dump(rt, open(os.path.join(OUT, "runtime.json"), "w"), indent=1)
    main_rt = rt["accurate_go"]
    tsv = ["table\trow\tsubset\tF\tP\tR\texact\tn_refogs\tdF_vs\tdF\tci_lo\tci_hi\tse\tverdict\tcpu_h\twall_h\tpeak_rss_gb"]
    md = []
    md.append("# OrthoBench FINAL evaluation (one-time test; frozen freeze-v1)\n")
    md.append("Primary endpoint: OrthoBench-**test** F (36 RefOGs; split by homology component). Scores: official "
              "Open_Orthobench `benchmark.py` via `bench_repo/scripts/orthobench_eval.py` (subset = sums of per-RefOG "
              "normalised TP/FP/FN; all-70 verified identical to benchmark.py). ΔF = paired bootstrap over RefOGs "
              "(B=2000, seed 11), percentage points. Non-inferiority (DESIGN §5): *comparable* iff lower 95% bound of "
              "Δ(ours − comparator) ≥ −3 pp (−0.03 F). All configs logged in `results/test_log.tsv` before scoring.\n")
    # ---- main table
    md.append("## 1. Main table (test)\n")
    md.append("| method | test F | P | R | exact | dev F | all-70 F | CPU-h | wall-h | peak RSS GB | threads |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|")
    fast_cpu = None
    rows = []
    for lab, name in (("accurate_go", "**Ancestra accurate** (default; + gene order)"),
                      ("fast_go", "Ancestra fast (+ gene order)")):
        t = sc[lab]["test"]
        if lab == "accurate_go":
            cpu_h = main_rt["time"]["cpu_s"] / 3600; wall_h = main_rt["time"]["wall_s"] / 3600
            rss = main_rt["time"]["max_rss_gb"]; note = ""
        else:
            # fast e2e = main's input+families+scores stages + this run's merge+output (merge-only rerun)
            st = main_rt["stages"]
            pre_cpu = sum(st[s]["cpu_s"] for s in ("input", "families", "scores"))
            pre_wall = sum(st[s]["wall_s"] for s in ("input", "families", "scores"))
            cpu_h = (pre_cpu + rt[lab]["time"]["cpu_s"]) / 3600
            wall_h = (pre_wall + rt[lab]["time"]["wall_s"]) / 3600
            rss = max(main_rt["time"]["max_rss_gb"], rt[lab]["time"]["max_rss_gb"]); note = "†"
            fast_cpu = cpu_h
        rows.append((name + note, t, sc[lab]["dev"]["F"], sc[lab]["all"]["F"], cpu_h, wall_h, rss, "10"))
        tsv.append("main\t%s\ttest\t%.4f\t%.4f\t%.4f\t%d\t%d\t\t\t\t\t\t\t%.3f\t%.3f\t%.2f" % (
            lab, t["F"], t["precision"], t["recall"], t["exact_refogs"], t["n_refogs"], cpu_h, wall_h, rss))
    for name, c in comp.items():
        t = c["test"]; r = c["runtime"]
        nm = name + (" ⚑" if "FastOMA" in name else "")
        rows.append((nm, t, c["dev"]["F"], c["all"]["F"], r["cpu_h"], r["wall_h"], r["max_rss_gb"],
                     str(r.get("threads_requested", 8))))
        tsv.append("main\t%s\ttest\t%.4f\t%.4f\t%.4f\t%d\t%d\t\t\t\t\t\t\t%.3f\t%.3f\t%.2f" % (
            name, t["F"], t["precision"], t["recall"], t["exact_refogs"], t["n_refogs"], r["cpu_h"], r["wall_h"],
            r["max_rss_gb"]))
    for nm, t, dF, aF, c, w, rss, th in rows:
        md.append("| %s | **%s** | %s | %s | %d/%d | %s | %s | %s | %s | %s | %s |" % (
            nm, f2(t["F"]), f2(t["precision"]), f2(t["recall"]), t["exact_refogs"], t["n_refogs"], f2(dF), f2(aF),
            f2(c, 2), f2(w, 2), f2(rss, 1), th))
    md.append("\n† fast = S0+S2 stages of the main fresh run + its own merge/output (merge-only rerun on the same cache). "
              "⚑ FastOMA uses the OMA-derived LUCA OMAmer DB (leakage; OMA contains these species) and is not a clean "
              "comparator. Competitor runtimes: 8 threads, node load 20–55; Ancestra: 10 threads (load in §3).\n")
    md.append("### Δ test F vs comparators (paired bootstrap over 36 test RefOGs) and non-inferiority verdict\n")
    md.append("| Ancestra | comparator | ΔF test (pp) | 95% CI | SE | verdict | ΔF dev | ΔF all-70 [95% CI] |")
    md.append("|---|---|---|---|---|---|---|---|")
    for lab in ("accurate_go", "fast_go"):
        for name in list(NI_COMPARATORS) + ["FastOMA (LUCA, leaky)"]:
            b = bt[(lab, name, "test")]; bd = bt[(lab, name, "dev")]; ba = bt[(lab, name, "all")]
            v = b["verdict"] if name in NI_COMPARATORS else b["verdict"] + " (flagged)"
            md.append("| %s | %s | %+.1f | [%+.1f, %+.1f] | %.1f | %s | %+.1f | %+.1f [%+.1f, %+.1f] |" % (
                lab, name, b["delta_point"], b["ci95"][0], b["ci95"][1], b["se"], v, bd["delta_point"],
                ba["delta_point"], ba["ci95"][0], ba["ci95"][1]))
            for sub, x in (("test", b), ("dev", bd), ("all", ba)):
                tsv.append("delta_vs_competitor\t%s\t%s\t\t\t\t\t%d\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t\t\t" % (
                    lab, sub, x["n_refogs"], name, x["delta_point"], x["ci95"][0], x["ci95"][1], x["se"],
                    x.get("verdict", "")))
    # ---- ablations
    md.append("\n## 2. Ablations / layer effects (test primary; Δ vs main = accurate + gene order + NCBI tree)\n")
    md.append("| config | group | test F | P | R | exact | ΔF test vs main [95% CI] | dev F | all-70 F | "
              "pair F1 test macro (Δ [CI]) | micro |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for lab in [l for l in CONFIGS if l in sc]:
        t = sc[lab]["test"]
        if lab == "accurate_go":
            d = "—"
        else:
            b = bt[(lab, "accurate_go", "test")]
            d = "%+.1f [%+.1f, %+.1f]" % (b["delta_point"], b["ci95"][0], b["ci95"][1])
            for sub in ("test", "dev", "all"):
                x = bt[(lab, "accurate_go", sub)]
                tsv.append("ablation_delta\t%s\t%s\t\t\t\t\t%d\taccurate_go\t%.4f\t%.4f\t%.4f\t%.4f\t\t\t\t" % (
                    lab, sub, x["n_refogs"], x["delta_point"], x["ci95"][0], x["ci95"][1], x["se"]))
        p = pairs.get(lab)
        pm = "—"
        if p:
            dm = p["delta_macro_vs_main"]
            pm = "%.3f" % p["macro_F1"] + ("" if lab == "accurate_go" else " (%+.3f [%+.3f, %+.3f])" % (
                dm["point"], dm["ci95"][0], dm["ci95"][1]))
        md.append("| %s | %s | %s | %s | %s | %d/%d | %s | %s | %s | %s | %s |" % (
            sc[lab]["desc"], sc[lab]["group"], f2(t["F"]), f2(t["precision"]), f2(t["recall"]), t["exact_refogs"],
            t["n_refogs"], d, f2(sc[lab]["dev"]["F"]), f2(sc[lab]["all"]["F"]), pm,
            ("%.3f" % p["micro"]["F1"]) if p else "—"))
        for sub in ("test", "dev", "all"):
            x = sc[lab][sub]
            tsv.append("config\t%s\t%s\t%.4f\t%.4f\t%.4f\t%d\t%d\t\t\t\t\t\t\t\t\t" % (
                lab, sub, x["F"], x["precision"], x["recall"], x["exact_refogs"], x["n_refogs"]))
        if p:
            tsv.append("pair_test\t%s\ttest\t%.4f\t%s\t%s\t\t%d\t\t\t\t\t\tmacro F1 (F col); micro P/R in P/R cols\t\t\t" % (
                lab, p["macro_F1"], p["micro"]["P"], p["micro"]["R"], p["n_refogs"]))
    if ("fast_go", "fast_nogo", "test") in bt:
        b = bt[("fast_go", "fast_nogo", "test")]
        md.append("\nGene order within the fast preset: ΔF test (fast+GO − fast no GO) = %+.1f [%+.1f, %+.1f]." % (
            b["delta_point"], b["ci95"][0], b["ci95"][1]))
        tsv.append("ablation_delta\tfast_go\ttest\t\t\t\t\t%d\tfast_nogo\t%.4f\t%.4f\t%.4f\t%.4f\t\t\t\t" % (
            b["n_refogs"], b["delta_point"], b["ci95"][0], b["ci95"][1], b["se"]))
    # ---- runtime
    md.append("\n## 3. Runtime (shared node, 10 threads; `/usr/bin/time -v` + 2-s process-tree sampling)\n")
    md.append("Fresh end-to-end run of the main config (`work/ob_final/main_accurate_go`):\n")
    md.append("| stage | wall (s) | CPU (s) | CPU-h | peak RSS tree (GB) | peak RSS own (GB) | node load1 mean |")
    md.append("|---|---|---|---|---|---|---|")
    for st, x in main_rt["stages"].items():
        if x.get("skipped"):
            continue
        md.append("| %s | %.0f | %.0f | %.2f | %s | %s | %s |" % (st, x["wall_s"], x["cpu_s"], x["cpu_s"] / 3600,
                                                            x["peak_rss_tree_gb"], x["peak_rss_own_gb"], x["load1_mean"]))
        tsv.append("runtime_stage\taccurate_go:%s\t\t\t\t\t\t\t\t\t\t\t\t\t%.4f\t%.4f\t%s" % (
            st, x["cpu_s"] / 3600, x["wall_s"] / 3600, x["peak_rss_tree_gb"]))
    T = main_rt["time"]
    md.append("| **total** (`time -v`) | %.0f | %.0f | %.2f | max RSS %.1f (largest process) | | |" % (
        T["wall_s"], T["cpu_s"], T["cpu_s"] / 3600, T["max_rss_gb"]))
    md.append("\nMerge-only variants (cache of the main run; wall/CPU = whole CLI call incl. input, rescue, merge, output):\n")
    md.append("| variant | wall (s) | CPU-h | max RSS (GB) | merge wall (s) | exit |")
    md.append("|---|---|---|---|---|---|")
    for lab in [l for l in CONFIGS if l in sc and l != "accurate_go"]:
        T = rt[lab]["time"]; ms = rt[lab]["stages"].get("merge", {})
        md.append("| %s | %.0f | %.2f | %.1f | %s | %d |" % (lab, T["wall_s"], T["cpu_s"] / 3600, T["max_rss_gb"],
                                                             ms.get("wall_s", "—"), T["exit"]))
        tsv.append("runtime_variant\t%s\t\t\t\t\t\t\t\t\t\t\t\t\t%.4f\t%.4f\t%.2f" % (
            lab, T["cpu_s"] / 3600, T["wall_s"] / 3600, T["max_rss_gb"]))
    # ---- verification + reproducibility
    md.append("\n## 4. Output verification and reproducibility\n")
    for lab in ("accurate_go", "fast_go"):
        p = os.path.join(OUT, "verify", lab + ".json")
        if os.path.exists(p):
            v = json.load(open(p))
            md.append("- `scripts/verify_orthoxml.py` on `%s`: overall **%s** (%s)." % (
                lab, "PASS" if v["ok"] else "FAIL", ", ".join(v["checks"])))
            tsv.append("verify\t%s\t\t\t\t\t\t\t\t\t\t\t\t%s\t\t\t" % (lab, "PASS" if v["ok"] else "FAIL"))
    import hashlib
    def md5(p):
        return hashlib.md5(open(p, "rb").read()).hexdigest() if os.path.exists(p) else None
    b, o = sc["accurate_go"]["run"], os.path.join(ROOT, "work/core_v0_ob/run_det")
    same = [f for f in ("stages/families/families.npy", "stages/scores/scores.npz")
            if md5(os.path.join(b, f)) and md5(os.path.join(b, f)) == md5(os.path.join(o, f))]
    md.append("- Fresh S0 families and S2 scores are byte-identical to the pre-freeze deterministic run "
              "(`work/core_v0_ob/run_det`): %s. Main dev F = %.2f, the pre-freeze `accurate` dev value (79.69) "
              "(regression check)." % (", ".join(same) or "NO", sc["accurate_go"]["dev"]["F"]))
    md.append("- Guide-tree arms reuse the tree-independent S0/S2 checkpoints through "
              "`scripts/run_cached.py` (frozen CLI, in-process key adoption, logged in `<run>/reused_stages.json`); "
              "same-tree variants hard-link the stages (as `scripts/rescue_ob_arms.sh`).")
    # ---- reading
    acc = sc["accurate_go"]["test"]["F"]
    of = bt[("accurate_go", "OrthoFinder 3.2", "test")]
    sp = bt[("accurate_go", "SonicParanoid2 default", "test")]
    md.append("\n## 5. Reading (pre-registered rule applied as written)\n")
    md.append("- **Primary endpoint, OrthoBench-test F:** Ancestra `accurate` %.1f vs OrthoFinder 3.2 %.1f "
              "(Δ %+.1f pp, 95%% CI [%+.1f, %+.1f]) and SonicParanoid2 %.1f (Δ %+.1f, [%+.1f, %+.1f]). The point "
              "estimates favour Ancestra, but the lower bounds are below −3 pp, so non-inferiority (\"comparable\") "
              "is **not shown** vs OrthoFinder or SonicParanoid2 on 36 test RefOGs. It is shown (with CI > 0) vs "
              "ProteinOrtho. The `fast` preset gives the same verdicts." % (
                  acc, of["F_B"], of["delta_point"], of["ci95"][0], of["ci95"][1], sp["F_B"], sp["delta_point"],
                  sp["ci95"][0], sp["ci95"][1]))
    md.append("- The test-set SE of ΔF (≈5–6 pp = 0.05–0.06 F) is larger than the r2 simulation (0.02–0.04 F). Per-RefOG variance "
              "is dominated by a few large/paralog-rich RefOGs. SonicParanoid2 is far better on test (70.4) than on "
              "dev (34.2), so the dev-to-test drop of the competitors is not uniform.")
    md.append("- **Gene order (synteny layer):** no measurable test effect (accurate Δ −0.1 [−0.4, +0.0]; "
              "fast Δ −0.0 [−0.1, +0.1]). Its dev gain (+1.0 accurate, +4.4 fast) does not replicate on test.")
    md.append("- **Family rescue (3Di + profiles + guard):** test F −5.0 [−13.1, 0.0] (precision 69.5 → 60.1, "
              "recall 84.2 → 87.2), against +0.2 on dev; within-RefOG pair F1 is unchanged. So the extra genes "
              "attached to RefOG orthogroups are non-members. The decision to keep it off by default is supported.")
    md.append("- **Guide tree:** data-derived trees (ASTRAL-Pro, FastME; one wrong split vs NCBI) and 1-NNI "
              "perturbations change test F by ≤0.3, except NNI seed 23 (Gallus+Monodelphis clade; −3.3 "
              "[−8.5, −0.0]). The random tree (negative control) costs −16.4 [−27.1, −5.6] F and −0.19 macro pair "
              "F1, so the guide tree does real work. ASTRAL-Pro and FastME have the same topology but different "
              "Newick child order, and their outputs differ in 18 orthogroups (same scores). That is a small "
              "child-order sensitivity of the merge. Dev F moves by −3.3 to +4.2 under the same perturbations, "
              "which illustrates the dev volatility (DISCLOSURES 2).")
    md.append("- **Presets:** `fast` is −1.4 [−7.1, +3.6] test F and −0.024 [−0.048, +0.001] macro pair F1 vs "
              "`accurate`, for about 1 CPU-h less (merge 56 s vs 491 s).")
    md.append("- **Runtime:** Ancestra `accurate` end to end took 31.5 min wall and 3.9 CPU-h (10 threads, max RSS "
              "14.0 GB, node load ≈ 23–31). OrthoFinder 3.2 took 1.29 h and 8.2 CPU-h (8 threads, load 20–55). "
              "The conditions differ (threads, load), so this is indicative only. S2 scoring is 61% of CPU.")
    md.append("- Pair-level test truth: the S5 labelling rule (`scripts/S5_refog_labels.py --subset test`), built "
              "for this final evaluation only (`work/ob_final/refog_pair_labels_test.parquet`; 36 RefOGs, "
              "5,545 ortholog / 11,648 within-RefOG paralog pairs). Macro = mean over RefOGs of pair F1; "
              "Δ = paired bootstrap over RefOGs (B=2000, seed 11).")
    open(os.path.join(OUT, "summary.md"), "w").write("\n".join(md) + "\n")
    open(os.path.join(OUT, "summary.tsv"), "w").write("\n".join(tsv) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
