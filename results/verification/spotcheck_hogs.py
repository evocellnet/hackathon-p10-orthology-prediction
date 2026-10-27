#!/usr/bin/env python
"""Verifier spot-check of 20 random root HOGs in the OrthoBench accurate output (seed 11),
independent re-derivation of root groups vs orthobench_ogs.txt, and pyham load of the 20-HOG subset."""
import collections, json, random, re, sys, os, copy, logging
import xml.etree.ElementTree as ET
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
RUN = f"{ROOT}/work/ob_final/main_accurate_go"
NS = "{http://orthoXML.org/2011/}"
OUT = f"{ROOT}/results/verification"
ET.register_namespace("", "http://orthoXML.org/2011/")
tree = ET.parse(f"{RUN}/hogs.orthoxml"); root = tree.getroot()
gid2prot, gid2sp = {}, {}
for sp in root.findall(NS + "species"):
    name = sp.get("name").split(".")[0]
    for g in sp.iter(NS + "gene"):
        gid2prot[g.get("id")] = g.get("protId"); gid2sp[g.get("id")] = name
groups = root.find(NS + "groups")
tops = [g for g in groups if g.tag == NS + "orthologGroup"]
def genes(el): return [r.get("id") for r in el.iter(NS + "geneRef")]
# 1. root groups vs orthobench_ogs.txt
ogs = []
for l in open(f"{RUN}/orthobench_ogs.txt"):
    if l.startswith("OG"):
        g = l.split()[1:]
        if len(g) > 1: ogs.append(frozenset(g))
mine = [frozenset(gid2prot[x] for x in genes(t)) for t in groups]  # top-level orthologGroups AND paralogGroups
mine_multi = [m for m in mine if len(m) > 1]
same = set(ogs) == set(mine_multi)
print("top-level orthologGroups", len(tops), "multi-gene", len(mine_multi), "OB export OGs(>1)", len(ogs), "identical sets:", same)
# RefOG membership for annotation
refog = {}
for i in range(1, 71):
    for gname in open(f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs/RefOG%03d.txt" % i):
        refog[gname.strip()] = "RefOG%03d" % i
sizes = sorted((len(m) for m in mine), reverse=True)
big = sorted(tops, key=lambda t: -len(genes(t)))[:6]
for t in big:
    gs = genes(t); rf = collections.Counter(refog.get(gid2prot[x], "-") for x in gs)
    print("BIG", t.get("id"), "tax", t.get("taxonId"), len(gs), dict(collections.Counter(gid2sp[x] for x in gs)), {k: v for k, v in rf.items() if k != "-"})
print("largest 10 root HOGs:", sizes[:10], "n>=500:", sum(s >= 500 for s in sizes), "genes in groups>=500:", sum(s for s in sizes if s >= 500))
random.seed(11)
pick = random.sample(range(len(tops)), 20)
rows = []
for i in pick:
    t = tops[i]; gs = genes(t)
    spc = collections.Counter(gid2sp[x] for x in gs)
    n_par = sum(1 for _ in t.iter(NS + "paralogGroup"))
    n_ortho = sum(1 for _ in t.iter(NS + "orthologGroup"))
    rf = collections.Counter(refog.get(gid2prot[x], "-") for x in gs)
    rows.append(dict(id=t.get("id"), taxonId=t.get("taxonId"), size=len(gs), n_species=len(spc),
                     max_copies=max(spc.values()), species=dict(sorted(spc.items())),
                     n_paralogGroups=n_par, n_orthologGroups=n_ortho,
                     refogs={k: v for k, v in rf.items() if k != "-"}))
for r in rows:
    print(r["id"], "tax", r["taxonId"], "size", r["size"], "nsp", r["n_species"], "maxcp", r["max_copies"],
          "par", r["n_paralogGroups"], "og", r["n_orthologGroups"], r["refogs"], r["species"])
# level name for taxonId: find from a nested element? print level names
# 2. pyham on the 20-HOG subset
sub = copy.deepcopy(root); sg = sub.find(NS + "groups")
keep = {tops[i].get("id") for i in pick}
for g in list(sg):
    if g.get("id") not in keep: sg.remove(g)
used = {x for g in sg for x in genes(g)}
for sp in sub.findall(NS + "species"):
    gg = sp.find(NS + "database").find(NS + "genes")
    for g in list(gg):
        if g.get("id") not in used: gg.remove(g)
    if len(gg) == 0: sub.remove(sp)
p = f"{OUT}/spotcheck_20hogs.orthoxml"
ET.ElementTree(sub).write(p, xml_declaration=True, encoding="UTF-8")
ok_ham = None
try:
    import pyham
    logging.getLogger("pyham").setLevel(logging.ERROR)
    ham = pyham.Ham(f"{RUN}/stages/input/guide_tree.nwk", p, use_internal_name=True, tree_format="newick")
    tl = ham.get_list_top_level_hogs()
    ok_ham = dict(n_top=len(tl), expected=20, sizes=sorted(len(h.get_all_descendant_genes()) for h in tl))
except Exception as e:
    ok_ham = dict(error=repr(e)[:500])
print("pyham (raw subset, no compat rewrite):", ok_ham)
json.dump(dict(root_groups_identical_to_ob_export=same, n_top=len(tops), n_multi=len(mine_multi),
               largest=sizes[:10], sample=rows, pyham=ok_ham), open(f"{OUT}/spotcheck_hogs.json", "w"), indent=1)
