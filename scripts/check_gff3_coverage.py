#!/usr/bin/env python3
"""Check that every OrthoBench input protein ID is found (as protein_id of a CDS) in the GFF3
files, and report the fraction per species. Usage: check_gff3_coverage.py <BENCHMARKS/Input> <gff3 dir>"""
import glob, gzip, os, re, sys
inp, gdir = sys.argv[1:3]
pid_re = re.compile(r"protein_id=([^;\s]+)")
for fa in sorted(glob.glob(os.path.join(inp, "*.fa"))):
    sp = os.path.basename(fa).split(".")[0]
    ids = {l[1:].strip() for l in open(fa) if l.startswith(">")}
    g = [x for x in glob.glob(os.path.join(gdir, sp.split("_")[0] + "*.gff3.gz")) if not os.path.islink(x) or True]
    g = sorted([x for x in g if sp.split("_")[1] in x], key=lambda x: "chr_patch_hapl_scaff" not in x)
    found = set()
    for fn in g[:1]:
        with gzip.open(fn, "rt") as f:
            for l in f:
                if "\tCDS\t" in l:
                    m = pid_re.search(l)
                    if m:
                        found.add(m.group(1))
    # WormBase OrthoBench IDs are gene IDs with '.1' suffix
    hit = len(ids & found)
    if sp.startswith("Caenorhabditis"):
        gid = set()
        with gzip.open(g[0], "rt") as f:
            for l in f:
                if "\tgene\t" in l:
                    m = re.search(r"gene_id=([^;\s]+)", l)
                    if m: gid.add(m.group(1) + ".1")
        hit = len(ids & gid)
    print(f"{sp}\t{len(ids)}\t{hit}\t{100*hit/len(ids):.2f}%\t{os.path.basename(g[0]) if g else 'NA'}")
