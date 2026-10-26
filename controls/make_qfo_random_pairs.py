#!/usr/bin/env python3
"""Negative control for the QfO benchmark: uniformly random inter-species protein pairs.

Reads the reference-set mapping.json.gz (from the QfO Zenodo reference data), picks one
UniProt accession per internal protein, and samples N random pairs of proteins from two
different species. Output: TSV with two UniProt accessions per line.
Usage: make_qfo_random_pairs.py mapping.json.gz out.tsv [--n 2000000] [--seed 42]
"""
import argparse, bisect, gzip, json, random, re

RE_UP = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
p = argparse.ArgumentParser()
p.add_argument("mapping"); p.add_argument("out")
p.add_argument("--n", type=int, default=2000000); p.add_argument("--seed", type=int, default=42)
a = p.parse_args()
d = json.load(gzip.open(a.mapping, "rt"))
acc = {}
for k, v in d["mapping"].items():
    if RE_UP.match(k) and v not in acc:
        acc[v] = k
ids = sorted(acc)
goff = d["Goff"]
sp = [bisect.bisect_right(goff, i - 1) - 1 for i in ids]
rng = random.Random(a.seed)
n = 0; seen = set()
with open(a.out, "w") as f:
    while n < a.n:
        i, j = rng.randrange(len(ids)), rng.randrange(len(ids))
        if sp[i] == sp[j] or (i, j) in seen or (j, i) in seen:
            continue
        seen.add((i, j)); n += 1
        f.write("%s\t%s\n" % (acc[ids[i]], acc[ids[j]]))
print("proteins with accession:", len(ids), "species:", len(d["species"]), "pairs:", n)
