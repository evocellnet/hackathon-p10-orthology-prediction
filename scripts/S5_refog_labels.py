#!/usr/bin/env python
"""Ortholog / paralog labels for gene pairs inside OrthoBench RefOGs (S5 benchmark truth).

Benchmark-only code (reads RefOGs; never imported by the method).

For each RefOG (dev split by default):
1. Read the curated RefOG tree ``Data_for_RefOGs/trees_tight/RefOGxxx.tre`` (unrooted
   in meaning; includes some outgroup / non-RefOG genes). Leaf = ``Genus_species_ID``.
2. Root it: among all edges, choose the root minimising the number of
   species-overlap duplication nodes (all leaves, outgroups included);
   ties -> most balanced root (|max root-to-tip left - right| smallest).
3. Prune to the RefOG genes of the 12 OrthoBench species (root position kept).
4. Species-overlap labelling: an internal node is a duplication iff its
   children's species sets overlap. A cross-species pair is an ortholog iff its
   LCA is a speciation node.
Pairs of genes in *different* RefOGs of the same homology component
(splits/orthobench_refog_components.tsv) are labelled paralog ("deep").

Output: work/S5/refog_pair_labels.parquet (gene_a, gene_b, sp_a, sp_b, refog_a, refog_b,
label 1=ortholog 0=paralog, kind=within|deep, lca_n_species, lca_dup_depth) + stats json.
"""
import argparse, json, os, re, sys
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.newick import parse_newick  # noqa: E402

OB = f"{ROOT}/data/orthobench/Open_Orthobench"
SPECIES = ["Caenorhabditis_elegans", "Canis_familiaris", "Ciona_intestinalis", "Danio_rerio",
           "Drosophila_melanogaster", "Gallus_gallus", "Homo_sapiens", "Monodelphis_domestica",
           "Mus_musculus", "Pan_troglodytes", "Rattus_norvegicus", "Tetraodon_nigroviridis"]


G2SP = {}


def leaf_species_id(name):
    """(species, gene id): the gene id is the longest '_'-token substring that is an
    OrthoBench protein ID (handles 'Genus_species_ID', bare 'ID', 'ID_2' suffixes);
    other leaves (outgroups) -> ('Genus_species', rest)."""
    parts = name.split("_")
    for L in range(len(parts), 0, -1):
        for i in range(0, len(parts) - L + 1):
            g = "_".join(parts[i:i + L])
            if g in G2SP:
                return G2SP[g], g
    return "_".join(parts[:2]), "_".join(parts[2:])


class UTree:
    """Undirected tree from a parsed Newick (multifurcating root allowed)."""

    def __init__(self, root):
        self.adj = defaultdict(dict)
        self.name = {}
        ids = {}
        for i, n in enumerate(root.preorder()):
            ids[id(n)] = i
            self.name[i] = n.name if not n.children else ""
        for n in root.preorder():
            for c in n.children:
                L = c.length if c.length is not None else 0.0
                self.adj[ids[id(n)]][ids[id(c)]] = L
                self.adj[ids[id(c)]][ids[id(n)]] = L
        # suppress a degree-2 root
        r = 0
        if len(self.adj[r]) == 2:
            (a, la), (b, lb) = self.adj[r].items()
            del self.adj[a][r], self.adj[b][r], self.adj[r]
            self.adj[a][b] = la + lb
            self.adj[b][a] = la + lb
        self.leaves = [i for i in self.adj if len(self.adj[i]) == 1]

    def edges(self):
        return [(a, b) for a in self.adj for b in self.adj[a] if a < b]

    def rooted(self, a, b):
        """Children lists of a tree rooted on edge (a, b): returns (children dict, root id)."""
        root = -1
        ch = {root: [a, b]}
        stack = [(a, b), (b, a)]
        while stack:
            u, par = stack.pop()
            ch[u] = [v for v in self.adj[u] if v != par]
            stack.extend((v, u) for v in ch[u])
        return ch, root


def species_sets(ch, root, sp_of):
    out = {}
    order = []
    stack = [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack.extend(ch.get(u, []))
    for u in reversed(order):
        kids = ch.get(u, [])
        out[u] = {sp_of[u]} if not kids else set().union(*(out[k] for k in kids))
    return out, order


def n_dups(ch, root, sp_of):
    ss, order = species_sets(ch, root, sp_of)
    d = 0
    for u in order:
        kids = ch.get(u, [])
        if len(kids) >= 2:
            seen = set()
            for k in kids:
                if seen & ss[k]:
                    d += 1
                    break
                seen |= ss[k]
    return d


def tip_depth(ch, u, ut, par):
    best = 0.0
    stack = [(u, 0.0)]
    while stack:
        x, dd = stack.pop()
        kids = ch.get(x, [])
        if not kids:
            best = max(best, dd)
        for k in kids:
            stack.append((k, dd + ut.adj[x][k]))
    return best


def label_refog(tree_path, refog_genes):
    root = parse_newick(open(tree_path).read().strip())
    ut = UTree(root)
    sp_of = {}
    gid = {}
    for l in ut.leaves:
        sp, g = leaf_species_id(ut.name[l])
        sp_of[l] = sp
        gid[l] = g
    best = None
    for a, b in ut.edges():
        ch, r = ut.rooted(a, b)
        nd = n_dups(ch, r, sp_of)
        bal = abs(tip_depth(ch, a, ut, b) + ut.adj[a][b] / 2 - tip_depth(ch, b, ut, a) - ut.adj[a][b] / 2)
        key = (nd, bal)
        if best is None or key < best[0]:
            best = (key, a, b)
    (ndup_full, _), a, b = best
    ch, r = ut.rooted(a, b)
    # prune to RefOG genes of the 12 species
    keep_leaf = {l for l in ut.leaves if gid[l] in refog_genes and sp_of[l] in SPECIES}
    ss_full, order = species_sets(ch, r, sp_of)
    # pruned children: compress unary
    pr = {}
    for u in reversed(order):
        kids = ch.get(u, [])
        if not kids:
            pr[u] = u if u in keep_leaf else None
            continue
        ks = [pr[k] for k in kids if pr[k] is not None]
        if not ks:
            pr[u] = None
        elif len(ks) == 1:
            pr[u] = ks[0]
        else:
            ch[u] = ks
            pr[u] = u
    proot = pr[r]
    chp = {}
    stack = [proot] if proot is not None else []
    while stack:
        u = stack.pop()
        kids = ch.get(u, []) if u not in keep_leaf else []
        chp[u] = kids
        stack.extend(kids)
    ss, order_p = species_sets(chp, proot, sp_of) if proot is not None else ({}, [])
    # label each internal node; depth of duplication
    is_dup = {}
    for u in order_p:
        kids = chp.get(u, [])
        if len(kids) >= 2:
            seen = set(); dup = False
            for k in kids:
                if seen & ss[k]:
                    dup = True
                seen |= ss[k]
            is_dup[u] = dup
    # pair labels via leaf sets per node
    leafs_under = {}
    pairs = []
    for u in reversed(order_p):
        kids = chp.get(u, [])
        if not kids:
            leafs_under[u] = [u]
            continue
        for i in range(len(kids)):
            for j in range(i + 1, len(kids)):
                for x in leafs_under[kids[i]]:
                    for y in leafs_under[kids[j]]:
                        if sp_of[x] != sp_of[y]:
                            pairs.append((gid[x], gid[y], sp_of[x], sp_of[y], 0 if is_dup[u] else 1,
                                          len(ss[u])))
        leafs_under[u] = [x for k in kids for x in leafs_under[k]]
        for k in kids:
            del leafs_under[k]
    found = {gid[l] for l in keep_leaf}
    return pairs, dict(n_leaves_tree=len(ut.leaves), n_refog_genes=len(refog_genes),
                       n_refog_in_tree=len(found), n_dup_full=ndup_full,
                       n_dup_pruned=int(sum(is_dup.values())), n_internal=len(is_dup))


def read_refog(path):
    return [l.strip() for l in open(path) if l.strip()]


def gene_species_map():
    m = {}
    for sp in SPECIES:
        fa = [f for f in os.listdir(f"{OB}/BENCHMARKS/Input") if f.startswith(sp + ".")][0]
        for l in open(f"{OB}/BENCHMARKS/Input/{fa}"):
            if l.startswith(">"):
                m[l[1:].split()[0]] = sp
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", default="dev")
    ap.add_argument("--out", default=f"{ROOT}/work/S5/refog_pair_labels.parquet")
    a = ap.parse_args()
    split = json.load(open(f"{ROOT}/splits/orthobench_split.json"))
    refogs = split[a.subset]
    g2sp = gene_species_map()
    G2SP.update(g2sp)
    rows, stats = [], {}
    members = {}
    for r in refogs:
        genes = set(read_refog(f"{OB}/BENCHMARKS/RefOGs/{r}.txt"))
        members[r] = genes
        tp = f"{OB}/Supporting_Data/Data_for_RefOGs/trees_tight/{r}.tre"
        pairs, st = label_refog(tp, genes)
        stats[r] = st
        for ga, gb, sa, sb, lab, nsp in pairs:
            rows.append((ga, gb, sa, sb, r, r, lab, "within", nsp))
    # deep paralogs: different RefOGs of the same component (both dev)
    comp = pd.read_csv(f"{ROOT}/splits/orthobench_refog_components.tsv", sep="\t")
    ndeep = 0
    for _, c in comp.iterrows():
        rs = [x for x in c["refogs"].split(",") if x in members]
        for r1, r2 in combinations(rs, 2):
            for x in members[r1]:
                for y in members[r2]:
                    if x in g2sp and y in g2sp and g2sp[x] != g2sp[y]:
                        rows.append((x, y, g2sp[x], g2sp[y], r1, r2, 0, "deep", -1))
                        ndeep += 1
    df = pd.DataFrame(rows, columns=["gene_a", "gene_b", "sp_a", "sp_b", "refog_a", "refog_b",
                                     "label", "kind", "lca_n_species"])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    df.to_parquet(a.out, index=False)
    summ = dict(subset=a.subset, n_refogs=len(refogs), n_pairs=len(df),
                n_ortholog=int((df.label == 1).sum()), n_paralog_within=int(((df.label == 0) & (df.kind == "within")).sum()),
                n_paralog_deep=ndeep,
                refog_genes_in_tree=sum(s["n_refog_in_tree"] for s in stats.values()),
                refog_genes=sum(s["n_refog_genes"] for s in stats.values()),
                per_refog=stats)
    json.dump(summ, open(a.out.replace(".parquet", ".json"), "w"), indent=1)
    print({k: v for k, v in summ.items() if k != "per_refog"})


if __name__ == "__main__":
    main()
