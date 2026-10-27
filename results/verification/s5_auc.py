"""Verifier: synteny ROC-AUC (Mann-Whitney, ties=0.5) per species pair from the stored per-pair score tables."""
import glob, json
import numpy as np, pandas as pd
from scipy.stats import rankdata
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
split = json.load(open(f"{ROOT}/splits/orthobench_split.json"))
def auc(s, y):
    r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
out = {}
for f in sorted(glob.glob(f"{ROOT}/work/S5/pair_scores_*.parquet")):
    d = pd.read_parquet(f); y = d["label"].values.astype(int)
    refs = set(d.refog_a) | set(d.refog_b)
    k = f.split("pair_scores_")[1][:-8]
    out[k] = dict(n_orth=int(y.sum()), n_par=int((1 - y).sum()), refogs_all_dev=refs <= set(split["dev"]),
                  seq=round(auc(d.seq.values, y), 4), w20=round(auc(d.syn_w20.values, y), 4),
                  w20_shuf11=round(auc(d.shuf11_w20.values, y), 4), w10=round(auc(d.syn_w10.values, y), 4))
    print(k, out[k])
json.dump(out, open(f"{ROOT}/results/verification/s5_auc.json", "w"), indent=1)
