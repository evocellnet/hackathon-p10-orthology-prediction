#!/usr/bin/env python3
"""FastOMA (patched run) -> OrthoBench OG files.
Each TOP-LEVEL orthologGroup of FastOMA_HOGs.orthoxml (= root-level HOG of the species tree, here Bilateria; identical to
RootHOGs.tsv) becomes one orthogroup; genes in no group are singletons (absent from the file).
Two variants from ONE patched run:
  <out>_patched.og   : all top-level HOGs (incl. subHOGs rescued by the gene-loss patch, logged in rescued_subhogs.tsv)
  <out>_released.og  : top-level HOGs minus those whose gene set equals a rescued subHOG = exactly what unpatched
                       FastOMA 0.5.1 outputs (rescued subHOGs are frozen and never merged, so they do not change any other HOG).
Usage: fastoma_v2_to_og.py FastOMA_HOGs.orthoxml rescued_subhogs.tsv <out_prefix> [<all_input_fasta_dir>]
"""
import sys, glob, os
from lxml import etree
src, resc, outp = sys.argv[1:4]
fa_dir = sys.argv[4] if len(sys.argv) > 4 else None
id2prot = {}; tops = []
for ev, el in etree.iterparse(src, events=("end",)):
    tag = etree.QName(el).localname
    if tag == "gene" and el.get("protId") is not None:
        id2prot[el.get("id")] = el.get("protId")
    elif tag == "groups":
        for top in el:
            genes = [id2prot[g.get("id")] for g in top.iter() if etree.QName(g).localname == "geneRef"]
            if genes: tops.append(genes)
rescued = set()
n_resc_genes = 0
for l in open(resc):
    a = l.rstrip("\n").split("\t")
    mem = frozenset(m.split("||")[0] for m in a[3].split(","))
    n_resc_genes += len(mem)
    if len(mem) >= 2: rescued.add(mem)
n_rel = 0; genes_p = genes_r = 0
with open(outp + "_patched.og", "w") as fp, open(outp + "_released.og", "w") as fr:
    for i, g in enumerate(tops):
        fp.write("HOG%07d: %s\n" % (i + 1, " ".join(g))); genes_p += len(g)
        if frozenset(g) not in rescued:
            n_rel += 1; genes_r += len(g); fr.write("HOG%07d: %s\n" % (i + 1, " ".join(g)))
n_in = None
if fa_dir:
    n_in = sum(1 for f in glob.glob(os.path.join(fa_dir, "*.fa")) for l in open(f) if l.startswith(">"))
print("top-level HOGs patched=%d (genes %d) released=%d (genes %d); rescued subHOGs logged=%d (size>=2: %d, genes %d); input proteins=%s"
      % (len(tops), genes_p, n_rel, genes_r, sum(1 for _ in open(resc)), len(rescued), n_resc_genes, n_in))
