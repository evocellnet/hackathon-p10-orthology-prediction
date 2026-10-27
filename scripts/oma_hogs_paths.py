#!/usr/bin/env python3
"""Stream the OMA browser HOG OrthoXML (oma-hogs.orthoXML.gz) and write, for every gene of the selected species,
its path of nested orthologGroups (root -> deepest): one TSV row per gene:
  omaid  species  root_hog_id  path   where path = '|'.join(f"{groupid}@{TaxRange}") from the root HOG down
a paralogGroup (duplication) on the path is written as the token 'PARALOG@-'. Usage: oma_hogs_paths.py oma-hogs.orthoXML.gz out.tsv "Homo sapiens" ...
"""
import sys, gzip
from lxml import etree
src, dst, species = sys.argv[1], sys.argv[2], set(sys.argv[3:])
NS = "{http://orthoXML.org/2011/}"
gid2prot = {}; gid2sp = {}
cur_sp = None
out = open(dst, "w"); n_top = n_genes = 0

def walk(el, path):
    global n_genes
    tag = el.tag[len(NS):] if el.tag.startswith(NS) else el.tag
    if tag == "orthologGroup":
        tr = None
        for c in el:
            if c.tag == NS + "property" and c.get("name") == "TaxRange":
                tr = c.get("value"); break
        path = path + ["%s@%s" % (el.get("id"), tr)]
    for c in el:
        ct = c.tag[len(NS):]
        if ct == "geneRef":
            g = c.get("id")
            if g in gid2prot:
                out.write("%s\t%s\t%s\t%s\n" % (gid2prot[g], gid2sp[g], path[0].split("@")[0], "|".join(path))); n_genes += 1
        elif ct == "orthologGroup":
            walk(c, path)
        elif ct == "paralogGroup":
            walk(c, path + ["PARALOG@-"])

ctx = etree.iterparse(gzip.open(src, "rb"), events=("start", "end"), huge_tree=True)
depth_groups = False
for ev, el in ctx:
    t = el.tag
    if ev == "start":
        if t == NS + "species":
            cur_sp = el.get("name") if el.get("name") in species else None
        elif t == NS + "groups":
            depth_groups = True
        continue
    # end events
    if t == NS + "gene":
        if cur_sp is not None:
            gid2prot[el.get("id")] = el.get("protId"); gid2sp[el.get("id")] = cur_sp
        el.clear()
    elif t == NS + "species":
        el.clear()
    elif depth_groups and t == NS + "orthologGroup" and el.getparent() is not None and el.getparent().tag == NS + "groups":
        n_top += 1
        walk(el, [])
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
print("genes of selected species in file:", len(gid2prot), "top-level HOGs:", n_top, "gene rows written:", n_genes, file=sys.stderr)
