"""Verifier: S4 pooled pair-label accuracy (single sets, c3di scope, matched truth) a_nj vs b_ft, own paired bootstrap."""
import json, re
import numpy as np, pandas as pd
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
d = pd.read_csv(f"{ROOT}/results/S4_trees/per_set.tsv", sep="\t")
split = json.load(open(f"{ROOT}/splits/orthobench_split.json"))
refs = {r for s in d["set"] for r in re.findall(r"RefOG\d+", s)}
s = d[(d.kind == "single") & (d.scope == "c3di")]
out = {"all_sets_dev": refs <= set(split["dev"])}
piv = s.pivot_table(index="set", columns="arm", values=["m_all_n", "m_all_ok"])
for scope in ("c3di", "full"):
    s = d[(d.kind == "single") & (d.scope == scope)]
    n = s.pivot(index="set", columns="arm", values="m_all_n"); ok = s.pivot(index="set", columns="arm", values="m_all_ok")
    res = {}
    for arm in ("a_nj", "b_ft", "b_iq", "d_part"):
        if arm not in n: continue
        res[arm] = round(ok[arm].sum() / n[arm].sum(), 4)
    rng = np.random.default_rng(11); k = len(n); S = rng.integers(0, k, size=(2000, k))
    A = (ok["b_ft"].values[S].sum(1) / n["b_ft"].values[S].sum(1)) - (ok["a_nj"].values[S].sum(1) / n["a_nj"].values[S].sum(1))
    res["n_sets"] = k
    res["delta_bft_minus_anj"] = round(res["b_ft"] - res["a_nj"], 4)
    res["ci95"] = [round(float(x), 4) for x in np.percentile(A, [2.5, 97.5])]
    out[scope] = res
print(out)
json.dump(out, open(f"{ROOT}/results/verification/s4_acc.json", "w"), indent=1)
