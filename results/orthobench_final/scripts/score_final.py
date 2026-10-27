#!/usr/bin/env python
"""OrthoBench FINAL scoring (one-time test evaluation; configs logged in results/test_log.tsv first).

For every config in CONFIGS (run dir must hold orthobench_ogs.txt + hogs.orthoxml + manifest.json):
  * official benchmark.py scores via bench_repo/scripts/orthobench_eval.py (all-70 / dev / test, per-RefOG TP/FP/FN)
    -> results/orthobench_final/scores/<label>.json + <label>.per_refog.tsv
  * paired bootstrap of dF over TEST RefOGs (B=2000, seed 11; also dev/all for information):
    Ancestra main configs vs OrthoFinder 3.2 / SonicParanoid2 default / ProteinOrtho (competitor JSONs in
    bench_repo/results/orthobench_competitors/), and every ablation vs the main config.
  * non-inferiority (DESIGN §5): comparable iff lower 95% bound of dF >= -0.03 (fraction) = -3 pp.
  * macro / micro ortholog-pair F1 on TEST RefOGs (scripts/orthobench_dev_pairs_eval.py logic with test labels
    built by scripts/S5_refog_labels.py --subset test).
Writes results/orthobench_final/scores.json (everything) for make_summary.py.
"""
import json
import os
import sys

import numpy as np

ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
sys.path.insert(0, os.path.join(ROOT, "bench_repo/scripts"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, ROOT)
import orthobench_eval as OE  # noqa: E402

OUT = os.path.join(ROOT, "results/orthobench_final")
BENCH = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS")
COMP = os.path.join(ROOT, "bench_repo/results/orthobench_competitors")
W = os.path.join(ROOT, "work/ob_final")
TEST_LABELS = os.path.join(W, "refog_pair_labels_test.parquet")
MARGIN_PP = 3.0

# label -> (run dir, description, group)
CONFIGS = {
    "accurate_go": (f"{W}/main_accurate_go", "Ancestra accurate + gene order (MAIN, fresh e2e)", "main"),
    "fast_go": (f"{W}/fast_go", "Ancestra fast + gene order", "main"),
    "accurate_nogo": (f"{W}/accurate_nogo", "accurate, no gene order", "layer"),
    "fast_nogo": (f"{W}/fast_nogo", "fast, no gene order", "layer"),
    "accurate_go_rescue": (f"{W}/accurate_go_rescue", "accurate + gene order + family rescue (3Di+profile) + guard",
                           "layer"),
    "tree_astralpro": (f"{W}/tree_astralpro", "guide tree ASTRAL-Pro (armB, NCBI-rooted)", "tree"),
    "tree_fastme": (f"{W}/tree_fastme", "guide tree FastME (armA, NCBI-rooted)", "tree"),
    "tree_nni1_s11": (f"{W}/tree_nni1_s11", "guide tree NCBI + 1 NNI (seed 11)", "tree"),
    "tree_nni1_s23": (f"{W}/tree_nni1_s23", "guide tree NCBI + 1 NNI (seed 23)", "tree"),
    "tree_nni1_s47": (f"{W}/tree_nni1_s47", "guide tree NCBI + 1 NNI (seed 47)", "tree"),
    "tree_random_s11": (f"{W}/tree_random_s11", "guide tree random (seed 11; negative control)", "tree"),
}
COMPETITORS = {"OrthoFinder 3.2": "orthofinder", "SonicParanoid2 default": "sonicparanoid_default",
               "ProteinOrtho": "proteinortho", "FastOMA (LUCA, leaky)": "fastoma_luca"}
NI_COMPARATORS = ["OrthoFinder 3.2", "SonicParanoid2 default", "ProteinOrtho"]


def split_ids():
    sp = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
    return {"all": None, "dev": sp["dev"], "test": sp["test"]}


def verdict(lo):
    return OE.ni_verdict(lo, MARGIN_PP)


def pair_eval_test(labels):
    import orthobench_dev_pairs_eval as PE
    if not os.path.exists(TEST_LABELS):
        import subprocess
        subprocess.run([sys.executable, os.path.join(ROOT, "scripts/S5_refog_labels.py"), "--subset", "test",
                        "--out", TEST_LABELS], check=True)
    PE.LABELS = TEST_LABELS
    truth = PE.load_truth("test")
    res = {}
    for lab in labels:
        m, arr, refogs = PE.evaluate(os.path.join(CONFIGS[lab][0], "hogs.orthoxml"), truth)
        per = arr.sum(axis=1).astype(float)
        tp, fn, fp = per[:, 0], per[:, 1], per[:, 2]
        with np.errstate(invalid="ignore", divide="ignore"):
            f1 = 2 * tp / (2 * tp + fp + fn)
        res[lab] = dict(micro=m["all"], macro_F1=float(np.nanmean(f1)), n_refogs=len(refogs),
                        n_refogs_defined=int(np.isfinite(f1).sum()),
                        per_refog_F1={r: (None if not np.isfinite(x) else round(float(x), 4))
                                      for r, x in zip(refogs, f1)}, arr=arr)
    # paired bootstrap of macro F1 vs main (RefOGs resampled, B=2000 seed 11)
    base = res["accurate_go"]
    rng = np.random.default_rng(11)
    R = base["arr"].shape[0]
    S = rng.integers(0, R, size=(2000, R))
    def f1v(arr):
        per = arr.sum(axis=1).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            return 2 * per[:, 0] / (2 * per[:, 0] + per[:, 2] + per[:, 1])
    fb = f1v(base["arr"])
    for lab in labels:
        fa = f1v(res[lab]["arr"])
        d = np.nanmean(fa[S], axis=1) - np.nanmean(fb[S], axis=1)
        res[lab]["delta_macro_vs_main"] = dict(point=float(np.nanmean(fa) - np.nanmean(fb)),
                                               ci95=[float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))])
        del res[lab]["arr"]
    del base
    return res


def main():
    os.makedirs(os.path.join(OUT, "scores"), exist_ok=True)
    ids = split_ids()
    scores, per = {}, {}
    for lab, (run, desc, grp) in CONFIGS.items():
        if not os.path.exists(os.path.join(run, "orthobench_ogs.txt")):
            print("missing", lab, file=sys.stderr)
            continue
        js = os.path.join(OUT, "scores", lab + ".json")
        r = OE.score_orthogroups(BENCH, os.path.join(run, "orthobench_ogs.txt"), quiet=True)
        r["label"] = lab
        json.dump(r, open(js, "w"), indent=1)
        with open(os.path.join(OUT, "scores", lab + ".per_refog.tsv"), "w") as fh:
            fh.write("refog\tsubset\tn_genes\tTP\tFP\tFN\tTP_raw\tFP_raw\tFN_raw\tn_pred_fragments\texact\n")
            for x in r["per_refog"]:
                s = "test" if x["refog"] in ids["test"] else "dev"
                fh.write("%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.1f\t%.1f\t%.1f\t%d\t%d\n" % (
                    x["refog"], s, x["n_genes"], x["TP"], x["FP"], x["FN"], x["TP_raw"], x["FP_raw"], x["FN_raw"],
                    x["n_pred_fragments"], x["exact"]))
        assert r["verify"]["identical_to_official"] or r["verify"]["max_abs_diff"] < 1e-9
        scores[lab] = {k: r[k] for k in ("all", "dev", "test", "verify", "n_pred_groups", "n_genes_assigned")}
        scores[lab].update(desc=desc, group=grp, run=run)
        per[lab] = r["per_refog"]
        print(lab, "test F %.2f dev %.2f all %.2f" % (r["test"]["F"], r["dev"]["F"], r["all"]["F"]))
    comp = {}
    for name, d in COMPETITORS.items():
        r = json.load(open(os.path.join(COMP, d, "score.json")))
        comp[name] = {k: r[k] for k in ("all", "dev", "test")}
        comp[name]["runtime"] = json.load(open(os.path.join(COMP, d, "runtime.json")))
        per[name] = r["per_refog"]
    boots = []
    for lab in ("accurate_go", "fast_go"):
        if lab not in per:
            continue
        for name in COMPETITORS:
            for sub in ("test", "dev", "all"):
                b = OE.bootstrap_delta_f(per[lab], per[name], ids[sub], B=2000, seed=11)
                b.update(A=lab, B_=name, subset=sub, verdict=verdict(b["ci95"][0]))
                boots.append(b)
    for lab in per:
        if lab in COMPETITORS or lab == "accurate_go" or "accurate_go" not in per:
            continue
        base = "accurate_go"
        for sub in ("test", "dev", "all"):
            b = OE.bootstrap_delta_f(per[lab], per[base], ids[sub], B=2000, seed=11)
            b.update(A=lab, B_=base, subset=sub)
            boots.append(b)
    # gene-order effect within the fast preset as well
    if "fast_nogo" in per and "fast_go" in per:
        for sub in ("test", "dev", "all"):
            b = OE.bootstrap_delta_f(per["fast_go"], per["fast_nogo"], ids[sub], B=2000, seed=11)
            b.update(A="fast_go", B_="fast_nogo", subset=sub)
            boots.append(b)
    labs = [l for l in CONFIGS if l in scores]
    pairs = pair_eval_test(labs) if "accurate_go" in scores else {}
    json.dump(dict(scores=scores, competitors=comp, bootstrap=boots, pair_test=pairs, margin_pp=MARGIN_PP),
              open(os.path.join(OUT, "scores.json"), "w"), indent=1, default=str)
    for b in boots:
        if b["subset"] == "test":
            print("%-20s - %-24s dF %+6.2f [%+6.2f, %+6.2f] %s" % (b["A"], b["B_"], b["delta_point"], b["ci95"][0],
                                                                  b["ci95"][1], b.get("verdict", "")))
    for lab, p in pairs.items():
        print("pair-test %-20s macro %.4f micro F1 %s" % (lab, p["macro_F1"], p["micro"]["F1"]))


if __name__ == "__main__":
    main()
