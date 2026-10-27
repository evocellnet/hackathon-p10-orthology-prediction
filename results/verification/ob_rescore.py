#!/usr/bin/env python
"""Independent OrthoBench re-score (verifier). Imports ONLY the official upstream benchmark.py
(readers, read_refogs, read_uncertain_refogs, calculate_benchmarks_pairwise). Subset scores are
obtained by calling the official calculate_benchmarks_pairwise on the subset of RefOG lists
(the official function pools normalised TP/FP/FN over whatever RefOGs it is given).
Per-RefOG normalised TP/FP/FN for the bootstrap come from a line-by-line copy of the official
per-RefOG loop; their pooled sums are checked against the official function.
Does not import bench_repo/scripts/orthobench_eval.py."""
import contextlib, io, json, os, re, sys
import numpy as np
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
BENCH = f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS"
sys.path.insert(0, BENCH)
import benchmark as ob  # official

def per_refog(refogs, uncerts, preds):
    out = []
    for refOg, uncert in zip(refogs, uncerts):
        refOg = refOg.difference(uncert)
        n = len(refOg); notp = set(refOg); tp = fp = fn = 0.
        for p in preds:
            ov = len(refOg & p)
            if ov > 0:
                p = p.difference(uncert); ov = len(refOg & p)
            if ov > 0:
                notp -= p; tp += ov*(ov-1)/2; fp += ov*(len(p)-ov); fn += (n-ov)*ov
        fn += len(notp)*(n-1); fn /= 2; N = float(n-1)
        out.append((tp/N, fp/N, fn/N, int(tp>0 and fp == 0 and fn == 0) if False else int(fp == 0 and fn == 0)))
    return np.array(out)

def official(refogs, uncerts, preds):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        f, p, r = ob.calculate_benchmarks_pairwise(refogs, uncerts, preds)
    ex = int(re.search(r"(\d+) Orthogroups exactly correct", buf.getvalue()).group(1))
    return dict(F=f, P=p, R=r, exact=ex)

def fpr(a):
    tp, fp, fn = a[..., 0].sum(-1), a[..., 1].sum(-1), a[..., 2].sum(-1)
    P = tp/(tp+fp); R = tp/(tp+fn); return 100*2*P*R/(P+R)

def main():
    split = json.load(open(f"{ROOT}/splits/orthobench_split.json"))
    names = ["RefOG%03d" % i for i in range(1, 71)]
    d = BENCH + "/RefOGs/"
    with contextlib.redirect_stdout(io.StringIO()):
        refogs = ob.read_refogs(d); unc = ob.read_uncertain_refogs(d + "low_certainty_assignments/")
        exp = ob.get_expected_genes()
    files = {
        "accurate_go": f"{ROOT}/work/ob_final/main_accurate_go/orthobench_ogs.txt",
        "fast_go": f"{ROOT}/work/ob_final/fast_go/orthobench_ogs.txt",
        "OrthoFinder": f"{ROOT}/data/orthobench/runs/orthofinder/Orthogroups.og",
        "SonicParanoid2": f"{ROOT}/data/orthobench/runs/sonicparanoid_default.og",
        "ProteinOrtho": f"{ROOT}/data/orthobench/runs/proteinortho.og",
    }
    res, per = {}, {}
    for lab, fn in files.items():
        preds = ob.read_orthogroups_smart(fn)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ob.check_orthogroups(preds, exp)
        a = per_refog(refogs, unc, preds); per[lab] = a
        res[lab] = {"check": buf.getvalue().strip().splitlines()[-1], "n_groups": len(preds)}
        for sub in ("all", "dev", "test"):
            ids = names if sub == "all" else split[sub]
            idx = [names.index(x) for x in ids]
            o = official([refogs[i] for i in idx], [unc[i] for i in idx], preds)
            o["F_from_perrefog"] = float(fpr(a[idx][None])[0]); o["n"] = len(idx)
            assert abs(o["F_from_perrefog"] - o["F"]) < 1e-9, (lab, sub, o)
            res[lab][sub] = o
        print(lab, {s: round(res[lab][s]["F"], 2) for s in ("test", "dev", "all")}, res[lab]["test"])
    # paired bootstrap over RefOGs, B=2000, seed 11 (own implementation)
    boots = []
    for A in ("accurate_go", "fast_go"):
        for Bn in ("OrthoFinder", "SonicParanoid2", "ProteinOrtho"):
            for sub in ("test", "dev", "all"):
                ids = names if sub == "all" else split[sub]
                idx = np.array([names.index(x) for x in ids])
                rng = np.random.default_rng(11)
                S = idx[rng.integers(0, len(idx), size=(2000, len(idx)))]
                dd = fpr(per[A][S]) - fpr(per[Bn][S])
                pt = fpr(per[A][idx][None])[0] - fpr(per[Bn][idx][None])[0]
                lo, hi = np.percentile(dd, [2.5, 97.5])
                boots.append(dict(A=A, B=Bn, subset=sub, delta=float(pt), lo=float(lo), hi=float(hi),
                                  se=float(dd.std(ddof=1)), NI=bool(lo >= -3.0)))
                if sub == "test":
                    print("%s - %s %s dF %+.2f [%+.2f, %+.2f] SE %.2f NI=%s" % (A, Bn, sub, pt, lo, hi, dd.std(ddof=1), lo >= -3))
    json.dump(dict(scores=res, bootstrap=boots, split_test=split["test"]), open(f"{ROOT}/results/verification/ob_rescore.json", "w"), indent=1)

if __name__ == "__main__":
    main()
