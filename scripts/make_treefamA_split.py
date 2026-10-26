#!/usr/bin/env python3
"""TreeFam-A dev/test split for QfO 2022 (DESIGN.md §4: 30% hash split of TreeFam-A as QfO dev).

The QfO 2022 reference (``ReconciledTrees_TreeFam-A.drw``, Zenodo v2022.1) stores TreeFam-A as ONE case
(``TreeFamA``): a set of mapped proteins (Darwin entry numbers) and a relation table {pair -> S|D}; the
TreeFam family IDs are not included. Families are therefore reconstructed as the **connected components of
the relation graph** (each reconciled tree relates its mapped proteins pairwise; 1,013 of 1,108 components
are complete cliques, the rest are connected but miss some relations). Entry numbers are mapped to UniProt
accessions through ``mapping.json.gz`` (1:1 for all 11,182 proteins).

Family ID = ``TFA:<lexicographically smallest UniProt accession in the component>`` (stable, deterministic).
dev  <=> int(sha256(family_id), 16) % 10 < 3        (pre-registered rule, PI 2026-10-26)
test <=> otherwise.

Output: ``splits/treefamA_split.json`` (rule, counts, per family: subset, members, #relations S/D) and the
file's sha256 appended/updated in ``splits/SHA256SUMS``. Reads evaluation data: scripts/ only.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data/qfo/reference_data/2022.1/ReconciledTrees_TreeFam-A.drw")
MAPPING = os.path.join(ROOT, "data/qfo/reference_data/2022.1/mapping.json.gz")
OUT = os.path.join(ROOT, "splits/treefamA_split.json")
SUMS = os.path.join(ROOT, "splits/SHA256SUMS")
ACC = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def load_relations(path: str = REF):
    s = open(path).read()
    i = s.index("table([")
    return [(int(a), int(b), e) for a, b, e in re.findall(r"\[\{(\d+),(\d+)\},\['(\w+)'\]\]", s[i:])]


def entry_to_acc(mapping_path: str = MAPPING) -> dict:
    with gzip.open(mapping_path, "rt") as fh:
        m = json.load(fh)["mapping"]
    inv: dict = {}
    for k, v in m.items():
        if ACC.match(k):
            inv.setdefault(v, set()).add(k)
    return inv


def families(rels, inv) -> dict:
    par: dict = {}

    def f(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b, _ in rels:
        par[f(a)] = f(b)
    comp = defaultdict(set)
    for x in list(par):
        comp[f(x)].add(x)
    nrel = defaultdict(lambda: {"S": 0, "D": 0})
    for a, b, e in rels:
        nrel[f(a)][e] += 1
    out = {}
    for r, members in comp.items():
        accs = []
        for x in members:
            a = inv.get(x, set())
            if len(a) != 1:
                raise SystemExit(f"entry {x}: {len(a)} accessions")
            accs.append(next(iter(a)))
        accs.sort()
        fid = "TFA:" + accs[0]
        n = len(accs)
        out[fid] = {"members": accs, "n_prot": n, "n_rel_S": nrel[r]["S"], "n_rel_D": nrel[r]["D"],
                    "clique": nrel[r]["S"] + nrel[r]["D"] == n * (n - 1) // 2}
    return out


def subset_of(fid: str) -> str:
    return "dev" if int(hashlib.sha256(fid.encode()).hexdigest(), 16) % 10 < 3 else "test"


def main() -> int:
    rels = load_relations()
    fams = families(rels, entry_to_acc())
    for fid, d in fams.items():
        d["subset"] = subset_of(fid)
    dev = sorted(f for f, d in fams.items() if d["subset"] == "dev")
    test = sorted(f for f, d in fams.items() if d["subset"] == "test")
    pairs = lambda fs, e: sum(fams[f][f"n_rel_{e}"] for f in fs)
    doc = {"rule": "family_id = 'TFA:' + min UniProt accession of a connected component of the TreeFam-A "
                   "relation graph (QfO 2022 ReconciledTrees_TreeFam-A.drw); dev iff int(sha256(family_id),16) % 10 < 3",
           "reference": os.path.relpath(REF, ROOT),
           "reference_md5": hashlib.md5(open(REF, "rb").read()).hexdigest(),
           "n_families": len(fams), "n_dev": len(dev), "n_test": len(test),
           "n_relations": len(rels),
           "dev_relations": {"S": pairs(dev, "S"), "D": pairs(dev, "D")},
           "test_relations": {"S": pairs(test, "S"), "D": pairs(test, "D")},
           "n_non_clique_families": sum(not d["clique"] for d in fams.values()),
           "dev": dev, "test": test, "families": dict(sorted(fams.items()))}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(doc, fh, indent=1)
    h = hashlib.sha256(open(OUT, "rb").read()).hexdigest()
    lines = [l for l in open(SUMS).read().splitlines() if not l.endswith("treefamA_split.json")] \
        if os.path.exists(SUMS) else []
    lines.append(f"{h}  treefamA_split.json")
    with open(SUMS, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[split] {len(fams)} families: dev {len(dev)} ({100 * len(dev) / len(fams):.1f}%), test {len(test)}; "
          f"dev S/D relations {doc['dev_relations']}; sha256 {h[:16]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
