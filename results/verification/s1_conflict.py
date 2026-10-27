"""Verifier: NCBI guide tree vs published QfO SwissTree species tree, own bipartition code via ete3 parsing only."""
import json
from ete3 import Tree
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
def load(p):
    t = Tree(open(p).read().strip(), format=1); return t
q = load(f"{ROOT}/results/S1_speciestree/qfo/ncbi.nwk")
r = load(f"{ROOT}/results/S1_speciestree/eval_only/qfo_swisstree_lineage_tree.nwk")
common = sorted(set(q.get_leaf_names()) & set(r.get_leaf_names()))
print("leaves", len(q), len(r), "common", len(common))
q.prune(common, preserve_branch_length=False); r.prune(common, preserve_branch_length=False)
L = frozenset(common)
def splits(t):
    out = set()
    for n in t.traverse():
        s = frozenset(n.get_leaf_names())
        if 1 < len(s) < len(L) - 1:
            a = s if min(s) == min(L) else L - s  # canonical side containing first leaf
            out.add(frozenset(a))
    return out
sq, sr = splits(q), splits(r)
def compat(a, b):
    A, B = a, b; return A <= B or B <= A or not (A & B) or (A | B) == L
conf = [a for a in sq if any(not compat(a, b) for b in sr)]
res = dict(n=len(common), n_splits_ncbi=len(sq), n_splits_pub=len(sr), conflict_frac=len(conf) / len(sq),
           nrf=(len(sq - sr) + len(sr - sq)) / (len(sq) + len(sr)))
print(res)
json.dump(res, open(f"{ROOT}/results/verification/s1_conflict.json", "w"), indent=1)
