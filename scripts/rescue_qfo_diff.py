"""Where does a rescue arm's TreeFam-A DEV change come from? Pair status changes (dev pairs only) between
R0 and an arm, split by whether the pair touches an attached gene, a gene of a family that received
attachments (receiving family), or neither (merge perturbation elsewhere, e.g. via deep-node medians).

    python scripts/rescue_qfo_diff.py work/rescue/qfo/R0 work/rescue/qfo/R2
"""
import json, os, sys
from collections import Counter
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import read_gene_table  # noqa
from qfo_treefam_dev_fast import dev_refpairs, predicted_orthologs  # noqa

r0, r1 = sys.argv[1], sys.argv[2]
refs = dev_refpairs()
genes = {x for r in refs for x in r[:2]}
p0 = predicted_orthologs(os.path.join(r0, "hogs.orthoxml"), genes)
p1 = predicted_orthologs(os.path.join(r1, "hogs.orthoxml"), genes)
gt = read_gene_table(os.path.join(r1, "stages/input/genes.tsv")); idx = gt.index()
fb = np.load(os.path.join(r1, "stages/rescue/families_before.npy"))
fa = np.load(os.path.join(r1, "stages/rescue/families_final.npy"))
moved = fb != fa
recv = np.zeros(int(fa.max()) + 1, bool); recv[np.unique(fa[moved])] = True
c = Counter()
for a, b, t, _ in refs:
    k = (a, b) if a < b else (b, a)
    x, y = k in p0, k in p1
    if x == y:
        continue
    ia, ib = idx.get(a), idx.get(b)
    cls = "attached" if (moved[ia] or moved[ib]) else ("receiving_family" if (recv[fa[ia]] or recv[fa[ib]]) else "elsewhere")
    lab = ("TP" if y else "FN") if t == "S" else ("FP" if y else "TN")
    c[(lab, cls)] += 1
out = {f"{l}_{k}": v for (l, k), v in sorted(c.items())}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(ROOT, "results/rescue", f"qfo_diff_{os.path.basename(r1)}_vs_{os.path.basename(r0)}.json"), "w"), indent=1)
