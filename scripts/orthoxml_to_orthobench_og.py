#!/usr/bin/env python3
"""OrthoXML (HOG) -> OrthoBench OG file: each TOP-LEVEL orthologGroup (root HOG of the species tree given to the method)
becomes one orthogroup. Usage: orthoxml_to_orthobench_og.py in.orthoxml out.og"""
import sys
from lxml import etree
src, dst = sys.argv[1], sys.argv[2]
id2prot = {}
n = 0
with open(dst, "w") as fo:
    for ev, el in etree.iterparse(src, events=("end",)):
        tag = etree.QName(el).localname
        if tag == "gene" and el.get("protId") is not None:
            id2prot[el.get("id")] = el.get("protId")
        elif tag == "groups":
            for top in el:
                genes = [id2prot[g.get("id")] for g in top.iter() if etree.QName(g).localname == "geneRef"]
                if genes:
                    n += 1; fo.write("HOG%07d: %s\n" % (n, " ".join(genes)))
print("top-level HOGs:", n)
