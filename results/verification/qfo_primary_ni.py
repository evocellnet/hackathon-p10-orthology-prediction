"""Verifier: pre-registered QfO primary endpoints with CIs from official raw RefPhylo outputs.
TreeFam-A TEST families (splits/treefamA_split.json) and SwissTrees (18 cases), Ancestra vs comparators,
paired bootstrap over families / cases (B=2000, seed 11), official RefPhylo formula (prior 1, counts/2).
Non-inferiority: lower 95% bound of delta >= -0.03."""
import collections, glob, gzip, json
import numpy as np
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
sp = json.load(open(f"{ROOT}/splits/treefamA_split.json"))
g2f = {g: f for f, d in sp["families"].items() for g in d["members"]}
f2s = {f: d["subset"] for f, d in sp["families"].items()}
runs = {"Ancestra_accurate": f"{ROOT}/work/qfo_final/cand2_accurate/raw",
        "Ancestra_fast": f"{ROOT}/work/qfo_final/cand2_fast/raw"}
for c in ("OrthoFinder", "FastOMA_v0.1.4_nov", "OMA_HOGs", "OMA_Pairs", "Hieranoid_2", "InParanoid5", "SonicParanoid2-sens-g"):
    runs[c] = f"{ROOT}/work/qfo_peek1/competitors/{c}/raw"
LAB = ["TP", "FP", "FN", "TN"]
def load(raw, test, unit):
    f = glob.glob(f"{raw}/{test}/RefPhylo-{test}_*_raw.txt.gz")[0]
    cnt = collections.defaultdict(lambda: np.zeros(4))
    for l in gzip.open(f, "rt"):
        if l.startswith("#"): continue
        case, a, b, lab = l.rstrip("\n").split("\t")
        if unit == "family":
            fa = g2f.get(a)
            if fa is None or f2s[fa] != "test": continue
            cnt[fa][LAB.index(lab)] += 1
        else:
            cnt[case][LAB.index(lab)] += 1
    return cnt
def rates(c):  # c: (..., 4) counts; official prior + halves
    tp, fp, fn = 1 + c[..., 0] / 2, 1 + c[..., 1] / 2, 1 + c[..., 2] / 2
    return tp / (tp + fn), tp / (tp + fp)
out = {}
for test, unit in (("TreeFam-A", "family"), ("SwissTrees", "case")):
    data = {k: load(v, test, unit) for k, v in runs.items()}
    units = sorted(set().union(*[set(d) for d in data.values()]))
    M = {k: np.array([d.get(u, np.zeros(4)) for u in units]) for k, d in data.items()}
    rng = np.random.default_rng(11); S = rng.integers(0, len(units), size=(2000, len(units)))
    res = {"n_units": len(units)}
    B_ = {}
    for k, m in list(M.items()):
        if unit == "family":   # one pooled case (as official TreeFam-A)
            r, p = rates(m.sum(0)); rb, pb = rates(m[S].sum(1))
        else:                  # macro over cases with >5 proteins: approximate by cases present in raw
            r, p = [x.mean() for x in rates(m)]; rb, pb = [x.mean(1) for x in rates(m[S])]
        res[k] = dict(TPR=round(float(r), 4), PPV=round(float(p), 4))
        B_[k] = (rb, pb)
    for a in ("Ancestra_accurate", "Ancestra_fast"):
        for c in [k for k in runs if not k.startswith("Ancestra")]:
            for i, m in enumerate(("TPR", "PPV")):
                d = B_[a][i] - B_[c][i]
                lo, hi = np.percentile(d, [2.5, 97.5])
                res[f"{a} - {c} {m}"] = dict(delta=round(res[a][m] - res[c][m], 4), ci95=[round(lo, 4), round(hi, 4)],
                                             NI=bool(lo >= -0.03))
    out[test + (" TEST families" if unit == "family" else "")] = res
print(json.dumps(out, indent=1))
json.dump(out, open(f"{ROOT}/results/verification/qfo_primary_ni.json", "w"), indent=1)
