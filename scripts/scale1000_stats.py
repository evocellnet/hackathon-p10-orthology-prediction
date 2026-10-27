#!/usr/bin/env python
"""Per-proteome protein/residue counts (+ duplicate-ID check) -> data/scale1000/proteome_stats.tsv and per-subset totals."""
import gzip, sys, os
ROOT = "data/scale1000"
sel = [l.rstrip("\n").split("\t") for l in open(f"{ROOT}/selection.tsv")][1:]
stats, seen, dup = {}, set(), 0
for f in sel:
    stem = f[0]
    n = res = 0
    for line in gzip.open(f"{ROOT}/proteomes/{stem}.fasta.gz", "rt"):
        if line[0] == ">":
            n += 1
            pid = line[1:].split()[0].split("|")[1]
            if pid in seen: dup += 1
            seen.add(pid)
        else:
            res += len(line.strip())
    stats[stem] = (n, res)
with open(f"{ROOT}/proteome_stats.tsv", "w") as fo:
    fo.write("stem\tn_proteins\tn_residues\n")
    for s, (n, r) in stats.items(): fo.write(f"{s}\t{n}\t{r}\n")
print("duplicate protein IDs across proteomes:", dup)
out = ["subset\tn_proteomes\tn_proteins\tn_residues\tmax_proteome_proteins"]
for n in (125, 250, 500, 1000):
    st = [l.strip() for l in open(f"{ROOT}/subsets/n{n}.txt")]
    out.append(f"{n}\t{len(st)}\t{sum(stats[s][0] for s in st)}\t{sum(stats[s][1] for s in st)}\t{max(stats[s][0] for s in st)}")
open(f"{ROOT}/subset_totals.tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
