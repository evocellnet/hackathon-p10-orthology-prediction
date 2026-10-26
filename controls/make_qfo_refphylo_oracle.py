#!/usr/bin/env python3
"""Oracle / anti-oracle controls for the QfO reference-gene-tree benchmarks (SwissTrees, TreeFam-A).

Input: RefPhylo raw files written by any run (raw/<testset>/RefPhylo-*_raw.txt.gz), which label
every evaluated pair (both proteins in a reference tree) as TP/FN (speciation = true orthologs) or
FP/TN (duplication = paralogs). Writes
  oracle_orthologs.tsv : all speciation pairs  -> expect TPR = PPV ~ 1 (minus the +1 pseudocounts)
  anti_oracle.tsv      : all duplication pairs -> expect TPR ~ 0, PPV ~ 0
NB: circular by construction; this only validates ID mapping, parsing and metric direction.
Usage: make_qfo_refphylo_oracle.py outdir raw1.txt.gz [raw2.txt.gz ...]
"""
import gzip, os, sys
outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
orth, para = set(), set()
for fn in sys.argv[2:]:
    for l in gzip.open(fn, "rt"):
        if l.startswith("#"):
            continue
        fam, a, b, lab = l.rstrip("\n").split("\t")
        pair = tuple(sorted((a, b)))
        (orth if lab in ("TP", "FN") else para).add(pair)
for name, s in (("oracle_orthologs.tsv", orth), ("anti_oracle.tsv", para)):
    with open(os.path.join(outdir, name), "w") as f:
        for a, b in sorted(s):
            f.write(f"{a}\t{b}\n")
print("orthologous pairs:", len(orth), "paralogous pairs:", len(para))
