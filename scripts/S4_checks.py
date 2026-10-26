#!/usr/bin/env python
"""S4 verification (single DEV RefOGs, c3di scope):
1. nRF cross-check with ete3 (independent implementation) for every arm tree.
2. Negative control: random topologies (NJ on uniform random distances, seeds 11/23/47),
   same rooting/labelling -> label accuracy floor; and the all-ortholog baseline.
3. Primary accuracies recomputed from the saved Newick files (not from per_set.tsv).
    python scripts/S4_checks.py -> results/S4_trees/checks.json
"""
import json, os, sys
import numpy as np
import pandas as pd
from ete3 import Tree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import read_gene_table  # noqa: E402
from ancestra.localtree import nj  # noqa: E402
from ancestra.treebench import SpeciesTree, UTree, nrf, pair_labels, root_reconcile  # noqa: E402

RUN = f"{ROOT}/work/core_v0_ob/run_v2"
OB = f"{ROOT}/data/orthobench/Open_Orthobench"
W = f"{ROOT}/work/S4/sets"
ARMS = ["a_nj", "b_ft", "b_iq", "c_ft3di", "c_iq3di", "d_part", "d_ftcat", "e_fident", "e_fsbits"]


def main():
    gt = read_gene_table(f"{RUN}/stages/input/genes.tsv")
    idx = gt.index()
    sp_of = {n: gt.species[s].split(".")[0] for n, s in zip(gt.names, gt.species_of)}
    stree = SpeciesTree(open(f"{RUN}/stages/input/guide_tree.nwk").read(), rename=lambda s: s.split(".")[0])
    df = pd.read_csv(f"{ROOT}/results/S4_trees/per_set.tsv", sep="\t")
    prim = df[(df.kind == "single") & (df.scope == "c3di")]

    def leaf_name(nm):
        parts = nm.split("_")
        for L in range(len(parts), 0, -1):
            for i in range(0, len(parts) - L + 1):
                g = "_".join(parts[i:i + L])
                if g in idx:
                    return g
        return nm

    spf = lambda nm: sp_of.get(nm) or "_".join(nm.split("_")[:2])  # noqa: E731
    ete_bad, n_cmp = [], 0
    acc = {a: [0, 0] for a in ARMS + ["random11", "random23", "random47", "all_ortholog"]}
    for s in sorted(prim.set.unique()):
        key = f"{s}__c3di"
        t_ref_full = UTree.from_newick(open(f"{OB}/Supporting_Data/Data_for_RefOGs/trees_tight/{s}.tre").read(),
                                       rename=leaf_name)
        trees = {a: UTree.from_newick(open(f"{W}/{key}/{a}.nwk").read()) for a in ARMS}
        leaves = trees["a_nj"].leaf_names()
        ref = t_ref_full.prune(leaves)
        rch, rr, _ = root_reconcile(ref, spf, stree)
        truth = pair_labels(ref, rch, rr, spf)
        acc["all_ortholog"][0] += sum(truth.values())
        acc["all_ortholog"][1] += len(truth)
        for seed in (11, 23, 47):
            rng = np.random.default_rng(seed + len(leaves))
            names = sorted(leaves)
            D = rng.uniform(0.1, 1.0, (len(names), len(names)))
            D = (D + D.T) / 2
            np.fill_diagonal(D, 0)
            e, nn = nj(D)
            trees[f"random{seed}"] = UTree.from_edges(e, nn, names)
        for a, t in trees.items():
            ch, r, _ = root_reconcile(t, spf, stree)
            lab = pair_labels(t, ch, r, spf)
            acc[a][0] += sum(int(lab[p] == y) for p, y in truth.items())
            acc[a][1] += len(truth)
            if a in ARMS:
                ours_n, ours_rf, ours_mx, _ = nrf(t, ref)
                # ete3 independent RF
                e1 = Tree(open(f"{W}/{key}/{a}.nwk").read())
                from ancestra.treebench import UTree as U  # noqa
                e2 = Tree(_to_newick(ref))
                rf, mx = e1.robinson_foulds(e2, unrooted_trees=True)[:2]
                n_cmp += 1
                # compare raw RF counts (ete3's max_rf for multifurcating trees is the number of
                # splits of both trees, ours is the binary 2(n-3) bound)
                if rf != ours_rf:
                    ete_bad.append((key, a, ours_rf, rf, ours_mx, mx))
    out = dict(ete3_rf_checked=n_cmp, ete3_rf_mismatches=ete_bad,
               label_acc={a: v[0] / v[1] for a, v in acc.items()}, n_pairs=acc["a_nj"][1])
    json.dump(out, open(f"{ROOT}/results/S4_trees/checks.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


def _to_newick(t):
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from S4_trees import to_newick
    return to_newick(t)


if __name__ == "__main__":
    main()
