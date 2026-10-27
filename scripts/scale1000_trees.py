#!/usr/bin/env python
"""NCBI-taxdump guide trees for each scale1000 subset + loader check (ancestra.io, no pipeline run)."""
import os, sys, time, shutil
sys.path.insert(0, os.getcwd())
from ancestra.speciestree import ncbi_tree, load_taxdump, write_tree
from ancestra.io import load_guide_tree, TreeIndex

ROOT = "data/scale1000"; RES = "results/scaling_large"
os.makedirs(f"{ROOT}/trees", exist_ok=True); os.makedirs(f"{RES}/subsets", exist_ok=True); os.makedirs(f"{RES}/trees", exist_ok=True)
tax = load_taxdump("data/taxonomy")
for n in (125, 250, 500, 1000):
    stems = [l.strip() for l in open(f"{ROOT}/subsets/n{n}.txt")]
    leaf_taxid = {s: int(s.split("_")[1]) for s in stems}
    t0 = time.time()
    tree = ncbi_tree(leaf_taxid, tax=tax)
    p = f"{ROOT}/trees/ncbi_n{n}.nwk"
    write_tree(tree, p, lengths=False)
    t1 = time.time()
    # loader check: file list = proteome stems (as io.species_name would give from file names)
    files = [f"{s}.fasta.gz" for s in stems]
    from ancestra.io import species_name
    species = sorted(species_name(f) for f in files)
    g = load_guide_tree(p, species)
    TreeIndex(g)
    nl = sum(1 for x in g.postorder() if not x.children); ni = sum(1 for x in g.postorder() if x.children)
    maxk = max(len(x.children) for x in g.postorder())
    print(f"n{n}: leaves {nl} (expect {n}), internal {ni}, max children {maxk}, build {t1-t0:.1f}s, load {time.time()-t1:.1f}s")
    assert nl == n and set(g.leaf_names()) == set(species)
    shutil.copy(p, f"{RES}/trees/ncbi_n{n}.nwk")
    shutil.copy(f"{ROOT}/subsets/n{n}.txt", f"{RES}/subsets/n{n}.txt")
for f in ("selection.tsv", "subset_totals.tsv", "md5sums.txt"):
    shutil.copy(f"{ROOT}/{f}", f"{RES}/subsets/{f}")
