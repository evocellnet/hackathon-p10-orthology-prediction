"""S6 diagnostic: how many DEV truth pairs (orthologs / paralogs, deep / shallow) have a sequence score,
and how many of the unscored ones get a structure-only score (Foldseek e <= 1e-3, 'miss' set) with real
vs shuffled 3Di.  -> results/S6_structure/pair_coverage.json"""
import json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import read_gene_table
from ancestra.scoring import ScoreTable
from orthobench_dev_pairs_eval import load_truth
RUN = os.path.join(ROOT, "work/core_v0_ob/run_v2")
DEEP = {"vert-teleost", "vert-invert", "invert-invert"}

gt = read_gene_table(os.path.join(RUN, "stages/input/genes.tsv")); n = gt.n
idx = {x: i for i, x in enumerate(gt.names)}
st = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
kseq = set((st.q.astype(np.int64) * n + st.t.astype(np.int64)).tolist())
has = np.zeros(n, bool); has[np.load(os.path.join(ROOT, "work/S6/db/real/genes.npy"))] = True
truth = load_truth()
a = np.array([idx[x] for x in truth.gene_a]); b = np.array([idx[x] for x in truth.gene_b])
k = np.minimum(a, b) * n + np.maximum(a, b)
inseq = np.fromiter((x in kseq for x in k.tolist()), bool, len(k))
out = {}
miss = {}
for v in ["real", "shuf11", "shuf23", "shuf47"]:
    z = np.load(os.path.join(ROOT, "work/S6/pairs/miss", v, "fs.npz"))
    q, t, e = z["q"].astype(np.int64), z["t"].astype(np.int64), z["evalue"]
    ok = e <= 1e-3
    miss[v] = set((np.minimum(q[ok], t[ok]) * n + np.maximum(q[ok], t[ok])).tolist())
deep = np.isin(truth.cls.values, list(DEEP))
for lab, ln in ((1, "ortholog"), (0, "paralog")):
    for dn, dm in (("deep", deep), ("shallow", ~deep)):
        m = (truth.label.values == lab) & dm
        both3 = has[a] & has[b]
        r = dict(n=int(m.sum()), with_seq_score=int((m & inseq).sum()), unscored=int((m & ~inseq).sum()),
                 unscored_both_3di=int((m & ~inseq & both3).sum()))
        for v, S in miss.items():
            r[f"struct_rescued_{v}"] = int(sum(1 for x, mm in zip(k.tolist(), (m & ~inseq).tolist()) if mm and x in S))
        out[f"{ln}_{dn}"] = r
json.dump(out, open(os.path.join(ROOT, "results/S6_structure/pair_coverage.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
