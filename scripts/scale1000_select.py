#!/usr/bin/env python
"""Select 1,000 UniProt reference proteomes for the scalability set (fixed rule, random.Random(11)).

Input : data/scale1000/README_refprot (UniProt reference_proteomes README, table of proteome ID, taxid, size)
        data/taxonomy/{nodes,names,merged}.dmp
Output: data/scale1000/selection.tsv (ordered; `rank_in_domain` defines the nested subsets),
        data/scale1000/subsets/n{125,250,500,1000}.txt (stems, one per line, sorted)
Rule  : per domain, eligible = 'canonical entries' (#(1)) in [lo, hi] (Bacteria 1000-40000, Archaea 700-40000,
        Eukaryota 3000-30000); one proteome per genus (random.Random(11) picks among a genus' members);
        genera are grouped by phylum (fallbacks: class, kingdom, order, ...) and the phyla are visited
        round-robin in a seeded shuffled order, one genus per phylum per cycle, until the domain quota
        (800/100/100) is filled. rank_in_domain = position in that order, so any prefix is phylum-balanced;
        subset n takes the first q_n(domain) of each domain (nested): 125 = 100/12/13, 250 = 200/25/25,
        500 = 400/50/50, 1000 = 800/100/100 (Bacteria/Archaea/Eukaryota).
"""
import os, random, sys
from collections import defaultdict
sys.path.insert(0, os.getcwd())
from ancestra.speciestree import load_taxdump

ROOT = "data/scale1000"
LIM = {"bacteria": (1000, 40000), "archaea": (700, 40000), "eukaryota": (3000, 30000)}
QUOTA = {"bacteria": 800, "archaea": 100, "eukaryota": 100}
SUB = {125: {"bacteria": 100, "archaea": 12, "eukaryota": 13},
       250: {"bacteria": 200, "archaea": 25, "eukaryota": 25},
       500: {"bacteria": 400, "archaea": 50, "eukaryota": 50},
       1000: QUOTA}
rng = random.Random(11)

parent, rank, name, merged = load_taxdump("data/taxonomy")
rows = []
with open(f"{ROOT}/README_refprot") as fh:
    started = False
    for line in fh:
        if line.startswith("Proteome_ID\t"):
            started = True
            continue
        if not started:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or not f[0].startswith("UP"):
            continue
        rows.append(dict(upid=f[0], taxid=int(f[1]), dom=f[3], n=int(f[4]), sp=f[7]))
print(len(rows), "proteomes in README")

def lineage(t):
    t = merged.get(t, t)
    out = {}
    while t in parent:
        r = rank[t]
        if r != "no rank" and r not in out:
            out[r] = t
        if parent[t] == t: break
        t = parent[t]
    return out

PH = ["phylum", "class", "kingdom", "order", "subkingdom", "superphylum", "clade"]
sel_all = []
for dom, (lo, hi) in LIM.items():
    el = [r for r in rows if r["dom"] == dom and lo <= r["n"] <= hi]
    byg = defaultdict(list)
    for r in sorted(el, key=lambda r: r["upid"]):
        lin = lineage(r["taxid"])
        g = lin.get("genus")
        r["genus"] = name.get(g, "") if g else r["sp"].split()[0]
        key = None
        for k in PH:
            if k in lin: key = name.get(lin[k], str(lin[k])); break
        r["phylum"] = key or dom
        byg[(r["phylum"], r["genus"])].append(r)
    byp = defaultdict(list)
    for (ph, g), members in sorted(byg.items()):
        byp[ph].append(rng.choice(members))
    for ph in sorted(byp): rng.shuffle(byp[ph])
    phyla = sorted(byp); rng.shuffle(phyla)
    order = []
    while len(order) < QUOTA[dom] and any(byp[p] for p in phyla):
        for p in phyla:
            if byp[p] and len(order) < QUOTA[dom]:
                order.append(byp[p].pop())
    print(f"{dom}: eligible {len(el)}, genera {len(byg)}, phyla {len(byp)}, selected {len(order)}, "
          f"phyla covered {len({r['phylum'] for r in order})}")
    for i, r in enumerate(order):
        r["rank_in_domain"] = i; sel_all.append(r)

os.makedirs(f"{ROOT}/subsets", exist_ok=True)
with open(f"{ROOT}/selection.tsv", "w") as fo:
    fo.write("stem\tproteome\ttaxid\tdomain\tphylum\tgenus\tn_canonical\trank_in_domain\tspecies\n")
    for r in sel_all:
        fo.write(f"{r['upid']}_{r['taxid']}\t{r['upid']}\t{r['taxid']}\t{r['dom']}\t{r['phylum']}\t{r['genus']}\t{r['n']}\t{r['rank_in_domain']}\t{r['sp']}\n")
for n, q in SUB.items():
    stems = [f"{r['upid']}_{r['taxid']}" for r in sel_all if r["rank_in_domain"] < q[r["dom"]]]
    assert len(stems) == n, (n, len(stems))
    open(f"{ROOT}/subsets/n{n}.txt", "w").write("\n".join(sorted(stems)) + "\n")
    print("subset", n, len(stems))
