"""Verifier: GO/EC avg Schlicker + SEM from the official raw per-pair files, and check that a random sample of
the scored pairs are orthologs in our own flattening of hogs.orthoxml (provenance of the raw file)."""
import gzip, glob, json, math, random, sys
import xml.etree.ElementTree as ET
label, hogs, rawdir = sys.argv[1:4]
OUT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/results/verification"
res = {}; sample = set(); pairs = []
random.seed(11)
for b in ("GO", "EC"):
    f = glob.glob(f"{rawdir}/{b}/{b}_*_raw.txt.gz")[0]
    v = []; rows = []
    for l in gzip.open(f, "rt"):
        if l.startswith("#"): continue
        a, c, s = l.split("\t"); v.append(float(s)); rows.append((a, c))
    m = sum(v) / len(v); sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))
    res[b] = dict(n=len(v), avg=m, sem=sd / math.sqrt(len(v)))
    for p in random.sample(rows, 300): pairs.append((b, p)); sample |= set(p)
NS = "{http://orthoXML.org/2011/}"
gid = {}; paths = {}; stack = []; uid = 0
for ev, el in ET.iterparse(hogs, events=("start", "end")):
    t = el.tag
    if ev == "start":
        if t == NS + "gene" and el.get("protId") in sample: gid[el.get("id")] = el.get("protId")
        elif t in (NS + "orthologGroup", NS + "paralogGroup"): uid += 1; stack.append((uid, t == NS + "orthologGroup"))
        elif t == NS + "geneRef" and el.get("id") in gid: paths[gid[el.get("id")]] = tuple(stack)
    else:
        if t in (NS + "orthologGroup", NS + "paralogGroup"):
            stack.pop()
            if not stack: el.clear()
def orth(a, b):
    pa, pb = paths.get(a, ()), paths.get(b, ()); last = None
    for x, y in zip(pa, pb):
        if x != y: break
        last = x
    return bool(last and last[1])
for b in ("GO", "EC"):
    ok = sum(orth(*p) for bb, p in pairs if bb == b); res[b]["sample_pairs_orthologs_in_hogs"] = f"{ok}/300"
print(label, res)
json.dump(res, open(f"{OUT}/qfo_goec_{label}.json", "w"), indent=1)
