#!/usr/bin/env python3
"""Convert a table of orthogroups (OrthoFinder Orthogroups.tsv / MCL_Orthogroups.tsv, SonicParanoid
multi_species_clusters / ortholog_groups.tsv, ProteinOrtho .proteinortho.tsv, ...) to the OrthoBench format:
one orthogroup per line, "<OGid>: gene1 gene2 ...", header lines dropped. Genes are extracted with the same regex as
Open_Orthobench benchmark.py::read_orthogroups_smart, so only the 12 OrthoBench proteomes' IDs are kept.
Usage: to_orthobench_og.py in.tsv out.og [--first-col-id]
Reports #groups, #genes, #genes duplicated across groups (must be 0 for the official check)."""
import re, sys, collections
GENE = re.compile(r"WBGene00\d+\.1|ENSCAFP\d+|ENSCINP\d+|ENSDARP\d+|FBpp0\d+|ENSGALP\d+|ENSP000\d+|ENSMODP\d+|ENSMUSP\d+|ENSPTRP\d+|ENSRNOP\d+|ENSTNIP\d+")
def main():
    src, dst = sys.argv[1], sys.argv[2]
    n = 0; cnt = collections.Counter()
    with open(src) as fi:
        first = fi.readline()
    if first.rstrip("\n").split("\t")[:2] in (["Group", "Protein"], ["RootHOG", "Protein"]):
        # long format (FastOMA OrthologousGroups.tsv / RootHOGs.tsv): one "group<TAB>protein" per line
        groups = collections.OrderedDict()
        with open(src) as fi:
            next(fi)
            for line in fi:
                f = line.rstrip("\n").split("\t")
                groups.setdefault(f[0], []).extend(GENE.findall(f[1]))
        with open(dst, "w") as fo:
            for k, g in groups.items():
                if g:
                    n += 1; cnt.update(g); fo.write("%s: %s\n" % (k, " ".join(g)))
        print("groups=%d genes=%d unique=%d duplicated_genes=%d" % (n, sum(cnt.values()), len(cnt), sum(1 for v in cnt.values() if v > 1)))
        return
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith("#"): continue
            g = GENE.findall(line)
            if not g: continue
            n += 1
            cnt.update(g)
            fo.write("OG%07d: %s\n" % (n, " ".join(g)))
    dup = sum(1 for v in cnt.values() if v > 1)
    print("groups=%d genes=%d unique=%d duplicated_genes=%d" % (n, sum(cnt.values()), len(cnt), dup))
if __name__ == "__main__":
    main()
