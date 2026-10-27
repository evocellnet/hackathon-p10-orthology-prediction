#!/usr/bin/env python3
"""OMA browser HOGs -> OrthoBench OG files, restricted to OrthoBench genes (mapping from oma_browser_map_*.py).
  bilateria: the HOG at the Bilateria level (see bilateria_label; ancestors of Bilateria in OMA: Eumetazoa, Metazoa,
             Opisthokonta, Eukaryota, LUCA); families younger than Bilateria -> their root HOG.
  root:      the root-level HOG (top-level orthologGroup).
Unmapped genes / genes whose OMA entry is in no HOG are singletons (absent from the file). If two OrthoBench proteins map to
the same OMA entry, both are placed in that entry's group.
Usage: oma_browser_to_orthobench_og.py mapping_final.tsv hog_paths.tsv out_prefix
"""
import sys, collections
mp, paths, outp = sys.argv[1:4]
ANC = {"Bilateria", "Eumetazoa", "Metazoa", "Opisthokonta", "Eukaryota", "LUCA"}
bil, root = {}, {}

def bilateria_label(oid, toks):
    """toks: root->leaf list of (groupid, TaxRange) incl. ('PARALOG', '-') for paralogGroups.
    Deepest orthologGroup with TaxRange in ANC = A (Bilateria or older). If A is Bilateria: label A. If A is older than
    Bilateria, OMA has no Bilateria-level group on this path, i.e. the lineage below A was not resolved at Bilateria:
    without a duplication between A and the gene, the gene is in A's single Bilateria lineage (label A); if a paralogGroup
    follows A, each paralog child is a separate lineage (label = the first orthologGroup below that paralogGroup, or the
    gene itself if the paralog child is a bare gene). If no group is Bilateria-or-older: the root HOG (younger family)."""
    ia = None
    for i, (g, tr) in enumerate(toks):
        if g != "PARALOG" and tr in ANC: ia = i
    if ia is None: return toks[0][0]
    if toks[ia][1] == "Bilateria": return toks[ia][0]
    rest = toks[ia + 1:]
    if not rest or rest[0][0] != "PARALOG": return toks[ia][0]
    for g, tr in rest:
        if g != "PARALOG": return g
    return "gene:" + oid
for l in open(paths):
    oid, sp, r, path = l.rstrip("\n").split("\t")
    toks = [g.rsplit("@", 1) for g in path.split("|")]
    bil[oid] = bilateria_label(oid, toks); root[oid] = r
for name, lab in (("bilateria", bil), ("root", root)):
    G = collections.defaultdict(list); n_map = 0
    for l in list(open(mp))[1:]:
        ob, sp, oid, meth, ih = l.rstrip("\n").split("\t")
        if oid in lab: G[lab[oid]].append(ob); n_map += 1
    with open(f"{outp}_{name}.og", "w") as fo:
        k = 0
        for g, mem in sorted(G.items()):
            if len(mem) >= 2: k += 1; fo.write(f"{g}: {' '.join(mem)}\n")
    sz = sorted((len(m) for m in G.values()), reverse=True)
    print(name, "genes in HOGs:", n_map, "groups (size>=2):", k, "largest:", sz[:5])
