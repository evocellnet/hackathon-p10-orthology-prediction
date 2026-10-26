#!/usr/bin/env python
"""S8 follow-up: is a human 'loss' shared with chimp (chimp gene = pseudogene modelled as coding)?

For every human candidate called B_pseudo / C_deleted (+ GULO/UOX), search the chimp locus of the HOG's chimp
gene with the *non-hominid* proteins only (mouse, rat, dog, opossum), ignoring the chimp annotation
(no CDS index -> unannotated rule with the majority vote). B_pseudo in chimp too => the loss is at least
Homininae-level (taxonomic placement of the loss moves up), and the chimp protein is a pseudogene model.
Writes results/S8_geneloss/homininae_check.tsv.
"""
import os, sys, pickle
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ancestra.geneloss import (FastaIndex, GeneOrderIndex, LocusSearcher, MAMMALS, SOURCE_ORDER, read_fasta,
                               syntenic_locus)
from run_S8_geneloss import GENOME, PROTEOME, KNOWN, gff_overlap, parse_blocks

R, W = "results/S8_geneloss", "work/S8"
hi = pickle.load(open(f"{W}/hogindex.pkl", "rb"))
go = GeneOrderIndex.from_dir("data/genomes/orthobench", MAMMALS)
prot = {}
for sp in MAMMALS:
    prot.update(read_fasta(os.path.join("data/orthobench/Open_Orthobench/BENCHMARKS/Input", PROTEOME[sp])))
T = "Pan_troglodytes"
g = FastaIndex(os.path.join("data/genomes_dna", GENOME[T]))
gff = pd.read_parquet(f"{W}/gff_genes_{T}.parquet")
S = LocusSearcher(hi, go, T, g, prot, f"{W}/tmp/homininae", cds=None)
c = pd.read_csv(f"{R}/final_calls.tsv", sep="\t")
c = c[(c.target == "Homo_sapiens") & c.final_cls.isin(["B_pseudo", "C_deleted"])]
jobs = [(int(h), cl, "") for h, cl in zip(c.hog, c.final_cls)]
for nm, (mp, _) in KNOWN.items():
    jobs.append((hi.hog_of[mp], "known", nm))
recs = []
for h, cl, nm in jobs:
    chimp = hi.in_species(h, T)
    hide = frozenset(chimp) | frozenset(hi.in_species(h, "Homo_sapiens"))
    loc = syntenic_locus(hi, go, h, T, [s for s in SOURCE_ORDER[T] if s != "Homo_sapiens"], contig_len=g.lengths())
    r = S.search(f"hominine{h}", h, loc, hide=hide)
    if r.get("hit_blocks"):
        r["gff_names"], r["gff_biotypes"] = gff_overlap(gff, r["contig"], parse_blocks(r["hit_blocks"]))
    # does the best hit fall on the annotated chimp gene of this HOG?
    on_gene = False
    for cg in chimp:
        rr = go.row(T, cg)
        if rr is not None and r.get("hit_start") and rr["contig"] == r["contig"] and rr["start"] <= r["hit_end"] and rr["end"] >= r["hit_start"]:
            on_gene = True
    r.update({"human_cls": cl, "known": nm, "chimp_genes": ",".join(chimp), "hit_on_chimp_gene": on_gene})
    recs.append(r)
d = pd.DataFrame(recs)
d.to_csv(f"{R}/homininae_check.tsv", sep="\t", index=False)
print(pd.crosstab(d.human_cls, d.cls))
print(d[["hog", "known", "human_cls", "cls", "hit_on_chimp_gene", "fs", "st", "nq_intact", "nq_disrupted", "gff_names"]].to_string())
