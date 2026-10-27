#!/usr/bin/env python
"""L2-genome keep-rule table: each real arm vs `off` and vs its shuffled seeds (docs/L2_genome.md §1.6).

    python scripts/L2_compare.py ob  F_c F_u ...     (OrthoBench dev, accurate base; primary macro pair F1)
    python scripts/L2_compare.py qfo fast_R_c ...    (TreeFam-A dev; primary TPR)

Paired bootstrap B=2000, seed 11, over dev RefOGs (ob) or dev families (qfo). Shuffled arms are
<arm>_s11/_s23/_s47. Writes results/l2_genome/keep_<ob|qfo>.json and prints a markdown table.
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B, SEED = 2000, 11


def ob_load(arm):
    return json.load(open(os.path.join(ROOT, "results/l2_genome/ob_arms", arm + ".json")))


def qfo_load(arm):
    return json.load(open(os.path.join(ROOT, "results/l2_genome/qfo_arms", arm + ".json")))


def macro(per, w):
    tp, fn, fp = per[:, 0], per[:, 1], per[:, 2]
    den = 2 * tp + fp + fn
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(den > 0, 2 * tp / den, np.nan)
    ok = ~np.isnan(f)
    return (w[:, ok] * f[ok]).sum(1) / w[:, ok].sum(1)


def tf_rates(c):          # c: (..., 4) TP FN FP TN; official pseudocount formula
    tp, fn, fp = 1 + c[..., 0] / 2, 1 + c[..., 1] / 2, 1 + c[..., 2] / 2
    return tp / (tp + fn), tp / (tp + fp)


def ci(x):
    return [round(float(np.percentile(x, 2.5)), 4), round(float(np.percentile(x, 97.5)), 4)]


def main():
    kind, arms = sys.argv[1], sys.argv[2:]
    out = {}
    if kind == "ob":
        off = ob_load("off")
        keys = sorted(off["per_refog"])
        po = np.array([off["per_refog"][k] for k in keys], float)
        R = len(keys)
        rng = np.random.default_rng(SEED)
        W = rng.multinomial(R, np.full(R, 1 / R), size=B).astype(float)
        m_off = macro(po, np.ones((1, R)))[0]
        print("| arm | OG-F (Δ) | pair P micro (Δ) | macro pair F1 (Δ [95% CI]) | shuffled macro (3 seeds) | Δ real−mean shuffled vs spread | S8 cand | GSR joins | keep |")
        print("|---|---|---|---|---|---|---|---|---|")
        for arm in arms:
            r = ob_load(arm)
            pr = np.array([r["per_refog"][k] for k in keys], float)
            m = macro(pr, np.ones((1, R)))[0]
            d = macro(pr, W) - macro(po, W)
            sh = []
            for s in (11, 23, 47):
                p = os.path.join(ROOT, "results/l2_genome/ob_arms", f"{arm}_s{s}.json")
                if os.path.exists(p):
                    sh.append(ob_load(f"{arm}_s{s}"))
            shm = [x["pair_macro_F1"] for x in sh]
            spread = (max(shm) - min(shm)) if shm else None
            dsh = (m - np.mean(shm)) if shm else None
            dP = r["pair_micro"]["P"] - off["pair_micro"]["P"]
            dF = r["OG_F"] - off["OG_F"]
            lo = ci(d)[0]
            keep = bool(lo > 0 and shm and dsh > spread and dP >= -0.01 and dF >= -1.0)
            e = dict(OG_F=r["OG_F"], dOG_F=round(dF, 2), pairP=r["pair_micro"]["P"], dpairP=round(dP, 4),
                     macro=round(m, 4), dmacro=round(m - m_off, 4), dmacro_ci=ci(d), shuffled_macro=shm,
                     d_real_minus_shuf=None if dsh is None else round(dsh, 4),
                     shuffle_spread=None if spread is None else round(spread, 4),
                     s8=r["s8_fission_cand"], s8_shuf=[x["s8_fission_cand"] for x in sh],
                     gsr_joins=r["stats"].get("gsr_repaired"), gsr_joins_shuf=[x["stats"].get("gsr_repaired") for x in sh],
                     dOG_F_shuf=[round(x["OG_F"] - off["OG_F"], 2) for x in sh], keep_rule_ob=keep)
            out[arm] = e
            print(f"| {arm} | {r['OG_F']:.2f} ({dF:+.2f}) | {r['pair_micro']['P']:.4f} ({dP:+.4f}) | {m:.4f} "
                  f"({m - m_off:+.4f} [{e['dmacro_ci'][0]:+.4f}, {e['dmacro_ci'][1]:+.4f}]) | "
                  f"{', '.join(f'{x:.4f}' for x in shm)} | {e['d_real_minus_shuf']} vs {e['shuffle_spread']} | "
                  f"{e['s8']} (shuf {e['s8_shuf']}) | {e['gsr_joins']} (shuf {e['gsr_joins_shuf']}) | {keep} |")
        out["off"] = dict(OG_F=off["OG_F"], pairP=off["pair_micro"]["P"], macro=round(m_off, 4),
                          s8=off["s8_fission_cand"])
    else:
        base = arms[0].split("_")[0] + "_off"
        off = qfo_load(base)
        fams = sorted(off["per_family"])
        co = np.array([off["per_family"][f] for f in fams], float)
        R = len(fams)
        rng = np.random.default_rng(SEED)
        W = rng.multinomial(R, np.full(R, 1 / R), size=B).astype(float)
        tpr0, ppv0 = tf_rates(co.sum(0))
        bo = tf_rates(W @ co)
        print("| arm | TPR (Δ [95% CI]) | PPV (Δ [95% CI]) | shuffled TPR | Δ real−mean shuffled vs spread | GSR joins | keep |")
        print("|---|---|---|---|---|---|---|")
        for arm in arms:
            r = qfo_load(arm)
            c = np.array([r["per_family"].get(f, [0, 0, 0, 0]) for f in fams], float)
            tpr, ppv = tf_rates(c.sum(0))
            br = tf_rates(W @ c)
            dt, dp = br[0] - bo[0], br[1] - bo[1]
            sh = [qfo_load(f"{arm}_s{s}") for s in (11, 23, 47)
                  if os.path.exists(os.path.join(ROOT, "results/l2_genome/qfo_arms", f"{arm}_s{s}.json"))]
            sht = [x["TF_dev"]["TPR"] for x in sh]
            spread = (max(sht) - min(sht)) if sht else None
            dsh = (tpr - np.mean(sht)) if sht else None
            keep = bool(ci(dt)[0] > 0 and sht and dsh > spread and (ppv - ppv0) >= -0.01)
            e = dict(TPR=round(tpr, 4), dTPR=round(tpr - tpr0, 4), dTPR_ci=ci(dt), PPV=round(ppv, 4),
                     dPPV=round(ppv - ppv0, 4), dPPV_ci=ci(dp), shuffled_TPR=[round(x, 4) for x in sht],
                     shuffled_PPV=[round(x["TF_dev"]["PPV"], 4) for x in sh],
                     d_real_minus_shuf=None if dsh is None else round(dsh, 4),
                     shuffle_spread=None if spread is None else round(spread, 4),
                     gsr_joins=r["stats"].get("gsr_repaired"), gsr_joins_shuf=[x["stats"].get("gsr_repaired") for x in sh],
                     keep_rule_tf=keep)
            out[arm] = e
            print(f"| {arm} | {tpr:.4f} ({tpr - tpr0:+.4f} [{e['dTPR_ci'][0]:+.4f}, {e['dTPR_ci'][1]:+.4f}]) | "
                  f"{ppv:.4f} ({ppv - ppv0:+.4f} [{e['dPPV_ci'][0]:+.4f}, {e['dPPV_ci'][1]:+.4f}]) | "
                  f"{', '.join(f'{x:.4f}' for x in sht)} | {e['d_real_minus_shuf']} vs {e['shuffle_spread']} | "
                  f"{e['gsr_joins']} (shuf {e['gsr_joins_shuf']}) | {keep} |")
        out[base] = dict(TPR=round(tpr0, 4), PPV=round(ppv0, 4))
    p = os.path.join(ROOT, f"results/l2_genome/keep_{kind}.json")
    old = json.load(open(p)) if os.path.exists(p) else {}
    old.update(out)
    json.dump(old, open(p, "w"), indent=1)


if __name__ == "__main__":
    main()
