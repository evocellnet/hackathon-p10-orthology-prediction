#!/usr/bin/env python3
"""Map OrthoBench (Ensembl 99) proteins to OMA browser entries (All.May2026).
Order per protein: (1) Ensembl protein ID (unversioned) in oma-ensembl xrefs, (2) Ensembl 99 gene ID (from the Ensembl 99
pep.all.fa headers; C. elegans OrthoBench IDs are WBGene IDs), (3) Ensembl 99 transcript ID, (4) exact sequence identity
against oma-seqs.fa.gz (same species), (5) exact sequence identity ignoring a leading 'M'/terminal '*' (not used if (4) hits).
The first route (in this order) that reaches an entry present in the HOG file (OMA main isoforms) is used; otherwise the
first route that maps at all. Among several candidates, entries in the HOG file are preferred; ties -> smallest OMA ID.
Writes mapping TSV: orthobench_id species oma_id method  (oma_id '-' = unmapped).
Usage: oma_browser_map_orthobench.py <OB Input dir> <ensembl99 pep dir> <ens12.tsv> <hog_paths.tsv> <oma-seqs.fa.gz|-> <out.tsv>
"""
import sys, os, gzip, glob, re, collections
ob_dir, pep_dir, ens_tsv, paths_tsv, seqs_fa, out = sys.argv[1:7]
SP = {"Caenorhabditis_elegans": "CAEEL", "Canis_familiaris": "CANLF", "Ciona_intestinalis": "CIOIN", "Danio_rerio": "DANRE",
      "Drosophila_melanogaster": "DROME", "Gallus_gallus": "CHICK", "Homo_sapiens": "HUMAN", "Monodelphis_domestica": "MONDO",
      "Mus_musculus": "MOUSE", "Pan_troglodytes": "PANTR", "Rattus_norvegicus": "RATNO", "Tetraodon_nigroviridis": "TETNG"}
unv = lambda x: re.sub(r"\.\d+$", "", x)

def read_fa(fh):
    h = None; s = []
    for l in fh:
        if l.startswith(">"):
            if h is not None: yield h, "".join(s)
            h = l[1:].rstrip("\n"); s = []
        else: s.append(l.strip())
    if h is not None: yield h, "".join(s)

# OrthoBench proteins
ob = {}  # id -> (species, seq)
for f in sorted(glob.glob(os.path.join(ob_dir, "*.fa"))):
    sp = os.path.basename(f).split(".")[0]
    for h, s in read_fa(open(f)):
        ob[h.split()[0]] = (sp, s.rstrip("*"))
# Ensembl 99 protein -> gene/transcript
p2g, p2t = {}, {}
for f in glob.glob(os.path.join(pep_dir, "*.pep.all.fa.gz")):
    for h, s in read_fa(gzip.open(f, "rt")):
        a = h.split(); pid = unv(a[0])
        d = dict(x.split(":", 1) for x in a[1:] if ":" in x and x.split(":", 1)[0] in ("gene", "transcript"))
        p2g[pid] = unv(d.get("gene", "")); p2t[pid] = unv(d.get("transcript", ""))
# OMA xrefs
x2oma = collections.defaultdict(set)
for l in open(ens_tsv):
    o, x = l.rstrip("\n").split("\t")[:2]
    x2oma[(o[:5], unv(x))].add(o)
in_hog = set(l.split("\t", 1)[0] for l in open(paths_tsv))
seq2oma = collections.defaultdict(set)
if seqs_fa != "-":
    codes = set(SP.values())
    for h, s in read_fa(gzip.open(seqs_fa, "rt")):
        o = h.split()[0]
        if o[:5] in codes and o[5:].isdigit():
            seq2oma[(o[:5], s.rstrip("*"))].add(o)
def pick(c):
    c = sorted(c); h = [x for x in c if x in in_hog]
    return (h or c)[0]
res = []; stat = collections.Counter(); spstat = collections.defaultdict(collections.Counter)
for pid, (sp, seq) in ob.items():
    code = SP[sp]; key = unv(pid); m = None; how = "unmapped"
    if sp == "Caenorhabditis_elegans":
        cands = [("gene", key)]
    else:
        cands = [("protein", key), ("gene", p2g.get(key, "")), ("transcript", p2t.get(key, ""))]
    routes = []
    for meth, k in cands:
        if k and (code, k) in x2oma: routes.append((meth, x2oma[(code, k)]))
    if (code, seq) in seq2oma: routes.append(("sequence", seq2oma[(code, seq)]))
    if seq.startswith("M") and (code, seq[1:]) in seq2oma: routes.append(("sequence_noM", seq2oma[(code, seq[1:])]))
    # first route that reaches an entry present in the HOG file (main isoform); else first route at all
    for meth, c in routes:
        if any(x in in_hog for x in c):
            m = pick(c); how = meth; break
    else:
        if routes: m = pick(routes[0][1]); how = routes[0][0]
    res.append((pid, sp, m or "-", how, "inhog" if m in in_hog else "nothog"))
    stat[how] += 1; spstat[sp][how] += 1; spstat[sp]["inhog" if m in in_hog else "nothog"] += 1
with open(out, "w") as fo:
    fo.write("orthobench_id\tspecies\toma_id\tmethod\tin_hog_file\n")
    for r in res: fo.write("\t".join(r) + "\n")
print("total", len(res), dict(stat))
for sp in sorted(spstat): print(sp, dict(spstat[sp]))
