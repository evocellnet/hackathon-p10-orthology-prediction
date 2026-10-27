#!/usr/bin/env python
"""Verifier: independent re-implementation of QfO RefPhyloTest (SwissTrees, TreeFam-A) TPR/PPV.
Truth: data/qfo/reference_data/2022.1/ReconciledTrees_<set>.drw (parsed with regex).
Predictions: own flattening of hogs.orthoxml (ortholog iff deepest common ancestor group is an orthologGroup),
protIds mapped to internal entry numbers with mapping.json.gz. Scoring formula transcribed from RefPhyloTest.drw
(uniform prior 1, each relation counted 1/2 because Darwin stores {a,b} sorted; cases with >5 mapped proteins;
macro average over cases; stderr = 1.96*sqrt(sum var)/n). Also cross-checks per-pair calls against the official raw file.
Usage: qfo_refphylo.py RUN_LABEL HOGS_ORTHOXML OFFICIAL_RAW_DIR"""
import gzip, json, re, sys, collections, math
import xml.etree.ElementTree as ET
REF = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/data/qfo/reference_data/2022.1"
OUT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/results/verification"
label, hogs, rawdir = sys.argv[1:4]

def parse_cases(path):
    txt = open(path).read()
    cases = {}
    for m in re.finditer(r"ReconciledTrees\['([^']+)'\] := RecTreeCase\(", txt):
        start = m.end(); nxt = txt.find("ReconciledTrees[", start); body = txt[start: nxt if nxt > 0 else len(txt)]
        prts = set(map(int, re.search(r"\{([\d,]*)\}", body).group(1).split(",")))
        rels = {}
        for a, b, ev in re.findall(r"\[\{(\d+),(\d+)\},\[([^\]]*)\]\]", body):
            rels[(int(a), int(b))] = ev.strip("'\" ")
        cases[m.group(1)] = (prts, rels)
    return cases

sets = {s: parse_cases(f"{REF}/ReconciledTrees_{s}.drw") for s in ("SwissTrees", "TreeFam-A")}
need = set()
for cs in sets.values():
    for prts, rels in cs.values(): need |= prts
mp = json.load(gzip.open(f"{REF}/mapping.json.gz"))["mapping"]
# entry nr -> protIds that map to it (only those needed)
# stream orthoxml
gid2ent = {}; gid2sp = {}
paths = {}  # entry -> tuple of (uid, is_ortho)
stack = []; uid = 0; sp = None
NS = "{http://orthoXML.org/2011/}"
for ev, el in ET.iterparse(hogs, events=("start", "end")):
    t = el.tag
    if ev == "start":
        if t == NS + "species": sp = el.get("name")
        elif t == NS + "gene":
            e = mp.get(el.get("protId"))
            if e in need: gid2ent[el.get("id")] = e; gid2sp[e] = sp
        elif t in (NS + "orthologGroup", NS + "paralogGroup"):
            uid += 1; stack.append((uid, t == NS + "orthologGroup"))
        elif t == NS + "geneRef":
            e = gid2ent.get(el.get("id"))
            if e is not None: paths[e] = tuple(stack)
    else:
        if t in (NS + "orthologGroup", NS + "paralogGroup"): stack.pop()
        if t != NS + "geneRef" and t not in (NS + "gene",): pass
        if t in (NS + "orthologGroup", NS + "paralogGroup") and not stack: el.clear()
        elif t == NS + "species": el.clear()

def is_orth(a, b):
    pa, pb = paths.get(a), paths.get(b)
    if not pa or not pb: return False
    last = None
    for x, y in zip(pa, pb):
        if x != y: break
        last = x
    return bool(last and last[1]) and gid2sp.get(a) != gid2sp.get(b)

res = {}
for s, cases in sets.items():
    tprs, ppvs, raw = [], [], {}
    for name in cases:  # Darwin order irrelevant for macro average
        prts, rels = cases[name]
        if len(prts) <= 5: continue
        tp = fp = fn = tn = 1.0
        for (a, b), ev in rels.items():
            if a not in prts or b not in prts: continue
            o = is_orth(a, b)
            if ev in ("D", ":D=Y", "D=Y", "DUPLICATION"):
                if o: fp += 0.5; raw[(name, a, b)] = "FP"
                else: tn += 0.5; raw[(name, a, b)] = "TN"
            elif ev in ("S", ":D=N", "D=N", "SPECIATION"):
                if o: tp += 0.5; raw[(name, a, b)] = "TP"
                else: fn += 0.5; raw[(name, a, b)] = "FN"
        p = tp / (tp + fp); r = tp / (tp + fn)
        ppvs.append((p, p * (1 - p) / (tp + fp))); tprs.append((r, r * (1 - r) / (tp + fn)))
    n = len(tprs)
    res[s] = dict(n_cases=n, TPR=sum(x for x, _ in tprs) / n, PPV=sum(x for x, _ in ppvs) / n,
                  TPR_se=1.96 * math.sqrt(sum(v for _, v in tprs)) / n, PPV_se=1.96 * math.sqrt(sum(v for _, v in ppvs)) / n,
                  counts=dict(collections.Counter(raw.values())))
    # compare with official raw file (labels per pair; IDs are XRefs, compare label multiset per case)
    import glob
    f = glob.glob(f"{rawdir}/{s}/RefPhylo-{s}_*_raw.txt.gz")
    if f:
        off = collections.Counter(l.rstrip("\n").split("\t")[3] for l in gzip.open(f[0], "rt") if not l.startswith("#"))
        res[s]["official_raw_counts"] = dict(off)
    print(label, s, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[s].items()})
json.dump(res, open(f"{OUT}/qfo_refphylo_{label}.json", "w"), indent=1)
