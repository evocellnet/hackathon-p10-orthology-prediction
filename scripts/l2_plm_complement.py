#!/usr/bin/env python
"""L2-plm diagnostic: is a pLM arm complementary to sequence? On the triplet dump of l2_plm_discrim.py --tag TAG:
accuracy of each arm where sequence is wrong vs right (ties excluded), real vs shuffled mean, per group.
    python scripts/l2_plm_complement.py --tag final  -> results/l2_plm/complement_final.json"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2_plm_common import WORK, RES, SEEDS


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="final"); a = ap.parse_args()
    z = np.load(os.path.join(WORK, f"triplets_{a.tag}.npz"))
    names = list(z["names"]); C = dict(zip(names, z["cols"])); T = z["T"]
    seq = C["seq"]
    groups = {"all": T[:, 2], "deep": T[:, 2] & T[:, 0], "hard": T[:, 2] & T[:, 1]}
    arms = sorted({n.split("|")[0] for n in names if "|" in n and not n.startswith("comb_")})
    out = {}
    for g, m in groups.items():
        w, r = m & (seq == 0), m & (seq == 1)
        out[g] = dict(n_seq_wrong=int(w.sum()), n_seq_right=int(r.sum()))
        for arm in arms:
            real = C[f"{arm}|real"]; sh = np.mean([C[f"{arm}|shuf{s}"] for s in SEEDS], 0)
            out[g][arm] = dict(acc_when_seq_wrong=round(float(real[w].mean()), 4), shuf_when_seq_wrong=round(float(sh[w].mean()), 4),
                               acc_when_seq_right=round(float(real[r].mean()), 4), shuf_when_seq_right=round(float(sh[r].mean()), 4))
    json.dump(out, open(os.path.join(RES, f"complement_{a.tag}.json"), "w"), indent=1)
    for g in out:
        print(g, out[g]["n_seq_wrong"], out[g]["n_seq_right"])
        for arm in arms:
            x = out[g][arm]
            print(f"  {arm:22s} seq-wrong: {x['acc_when_seq_wrong']:.3f} (shuf {x['shuf_when_seq_wrong']:.3f})  "
                  f"seq-right: {x['acc_when_seq_right']:.3f} (shuf {x['shuf_when_seq_right']:.3f})")


if __name__ == "__main__":
    main()
