#!/usr/bin/env python3
"""Nested, domain-stratified species subsets of QfO 2022 for the runtime-scaling benchmark (DESIGN §5).

Sizes 8/16/32/64/78. Per domain a random permutation (random.Random(11), over the sorted species list);
subset n takes the first k_d(n) species of each domain, with k_d(n) from largest-remainder proportional
allocation (checked monotone in n, so the subsets are nested). Writes results/scaling/subsets/sub{n}.{txt,nwk,json}
and symlink dirs work/scaling/inputs/sub{n}/ (the 78 = data/inputs/qfo2022/proteomes). Guide tree = the NCBI tree
pruned to the subset (ancestra.speciestree.prune_to).
"""
import json, os, random, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.newick import read_newick
from ancestra.speciestree import prune_to, write_tree

SRC = os.path.join(ROOT, "data/qfo/proteomes_2022_02_updZEBRAFISH")
DOMS = ["Eukaryota", "Bacteria", "Archaea"]
SIZES = [8, 16, 32, 64, 78]
SEED = 11
OUT = os.path.join(ROOT, "results/scaling/subsets")
INP = os.path.join(ROOT, "work/scaling/inputs")


def canon():
    out = {}
    for d in DOMS:
        for f in sorted(os.listdir(os.path.join(SRC, d))):
            if f.startswith("UP") and f.endswith(".fasta") and "_additional" not in f and "_DNA" not in f:
                out[f[:-6]] = (d, os.path.join(SRC, d, f))
    return out


def nprot(p):
    n = 0
    with open(p) as fh:
        for l in fh:
            n += l.startswith(">")
    return n


def alloc(n, counts):
    tot = sum(counts.values())
    q = {d: n * c / tot for d, c in counts.items()}
    k = {d: int(q[d]) for d in q}
    for d in sorted(q, key=lambda d: -(q[d] - k[d]))[: n - sum(k.values())]:
        k[d] += 1
    return k


def main():
    cf = canon()
    assert len(cf) == 78
    by = {d: sorted(s for s, (dd, _) in cf.items() if dd == d) for d in DOMS}
    counts = {d: len(v) for d, v in by.items()}
    rng = random.Random(SEED)
    perm = {}
    for d in DOMS:
        p = list(by[d]); rng.shuffle(p); perm[d] = p
    ks = {n: alloc(n, counts) for n in SIZES}
    for a, b in zip(SIZES, SIZES[1:]):
        assert all(ks[a][d] <= ks[b][d] for d in DOMS), (a, b, ks)
    tree = read_newick(os.path.join(ROOT, "data/inputs/qfo2022/guide_tree_ncbi.nwk"))
    os.makedirs(OUT, exist_ok=True)
    info = {}
    prev = set()
    for n in SIZES:
        sub = sorted(s for d in DOMS for s in perm[d][: ks[n][d]])
        assert len(sub) == n and prev <= set(sub)
        prev = set(sub)
        npr = {s: nprot(cf[s][1]) for s in sub}
        open(f"{OUT}/sub{n}.txt", "w").write("\n".join(sub) + "\n")
        write_tree(prune_to(tree, sub), f"{OUT}/sub{n}.nwk", lengths=False)
        d = f"{INP}/sub{n}"
        os.makedirs(d, exist_ok=True)
        for s in sub:
            dst = f"{d}/{s}.fasta"
            if not os.path.lexists(dst):
                os.symlink(cf[s][1], dst)
        info[n] = {"n_species": n, "domains": ks[n], "n_proteins": sum(npr.values()), "seed": SEED,
                   "species": {s: {"domain": cf[s][0], "n_proteins": npr[s]} for s in sub}}
        json.dump(info[n], open(f"{OUT}/sub{n}.json", "w"), indent=1)
        print(n, ks[n], sum(npr.values()))


main()
