#!/usr/bin/env python
"""S5 phase 2 analysis: pair-level Delta(real-off), Delta(real-shuffled) with paired bootstrap over
DEV RefOGs (B=2000, seed 11), per species-pair class, and the DESIGN §7 keep rule:
Delta(real-shuffled) > 0 (reported with CI) AND Delta(real-off) > seed spread
(seed spread = max-min pair F1 over the 3 shuffled seeds of the same config).
Writes results/S5_synteny/phase2/summary.json and phase2/tables.md.
"""
import glob, json, os, sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from orthobench_dev_pairs_eval import CLASSES, metrics, paired_bootstrap, prf  # noqa: E402

D = os.path.join(ROOT, "results/S5_synteny/phase2")
A = os.path.join(D, "arms")
SEEDS = (11, 23, 47)
FOCUS = ["all", "mammal-mammal", "tetrapod-tetrapod", "vert-teleost", "teleost-teleost", "vert-invert"]


def load(name):
    return np.load(os.path.join(A, f"{name}.npy")).astype(float), json.load(open(os.path.join(A, f"{name}.json")))


def f1(arr, cls="all"):
    c = arr.sum(axis=(0, 1)) if cls == "all" else arr[:, CLASSES.index(cls)].sum(axis=0)
    return float(prf(c)[2])


def main():
    off, offj = load("off")
    reals = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(A, "real_*.npy")))
    rows, summ = [], {"off": offj["pair"], "arms": {}}
    for r in reals:
        arr, j = load(r)
        sh = [f"{r.replace('real_', 'shuf_')}_s{s}" for s in SEEDS]
        if not all(os.path.exists(os.path.join(A, f"{x}.npy")) for x in sh):
            continue
        sarr = [load(x)[0] for x in sh]
        smean = np.mean(sarr, axis=0)
        d_off = paired_bootstrap(off, arr)
        d_shuf = paired_bootstrap(smean, arr)
        d_shuf_off = paired_bootstrap(off, smean)
        ent = dict(pair=j["pair"], og=j["og"], n_ctx_bonus=j["n_ctx_bonus"], delta_vs_off=d_off,
                   delta_vs_shuffled=d_shuf, delta_shuffled_vs_off=d_shuf_off, per_class={})
        for cls in FOCUS:
            sF = [f1(x, cls) for x in sarr]
            spread = float(np.nanmax(sF) - np.nanmin(sF))
            do = d_off.get(cls, {}).get("F1", {})
            ds = d_shuf.get(cls, {}).get("F1", {})
            keep = (ds.get("delta") is not None and ds["delta"] > 0 and do.get("delta") is not None
                    and do["delta"] > spread)
            keep_strict = keep and ds.get("ci95") and ds["ci95"][0] > 0
            ent["per_class"][cls] = dict(F1=round(f1(arr, cls), 4), F1_off=round(f1(off, cls), 4),
                                         F1_shuf=[round(x, 4) for x in sF], seed_spread=round(spread, 4),
                                         d_off=do, d_shuf=ds, keep=bool(keep), keep_ci=bool(keep_strict))
            rows.append(dict(arm=r, cls=cls, F1=round(f1(arr, cls), 4), F1_off=round(f1(off, cls), 4),
                             F1_shuf_mean=round(float(np.mean(sF)), 4), seed_spread=round(spread, 4),
                             d_off=do.get("delta"), d_off_ci=do.get("ci95"), d_shuf=ds.get("delta"),
                             d_shuf_ci=ds.get("ci95"), dP_off=d_off.get(cls, {}).get("P", {}).get("delta"),
                             dR_off=d_off.get(cls, {}).get("R", {}).get("delta"), keep=keep, keep_ci=keep_strict,
                             OG_F=round(j["og"]["F"], 2)))
        summ["arms"][r] = ent
    df = pd.DataFrame(rows)
    json.dump(summ, open(os.path.join(D, "summary.json"), "w"), indent=1, default=str)
    out = ["### Phase 2: pair-level F1 on OrthoBench DEV (all classes), real arms vs off and shuffled\n"]
    a = df[df.cls == "all"].sort_values("d_off", ascending=False)
    out.append(a.drop(columns=["cls"]).to_markdown(index=False))
    best = a.arm.iloc[0]
    for arm in dict.fromkeys([best, "real_w10_l0.3", "real_w10_l0.5", "real_w5_l0.1", "real_w10_l0.1"]):
        if arm not in set(df.arm):
            continue
        out.append(f"\n\n### Per species-pair class: `{arm}`\n")
        out.append(df[df.arm == arm].drop(columns=["arm"]).to_markdown(index=False))
    # off metrics per class
    om = summ["off"]
    out.append("\n\n### Off arm (sequence-only, R1-only base, min_cov=0) per class\n")
    out.append(pd.DataFrame([dict(cls=k, **v) for k, v in om.items()]).to_markdown(index=False))
    open(os.path.join(D, "tables.md"), "w").write("\n".join(out))
    print("\n".join(out))


if __name__ == "__main__":
    main()
