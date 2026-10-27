#!/usr/bin/env python
"""L2-genome idea 3: intron positions (protein coordinates) for OrthoBench DEV proteins.

Ensembl 99 GFF3 CDS segments -> per protein, the list of CDS-segment junctions as
(residue index, phase): junction after L coding nt (5'->3', strand-aware, after removing the
first segment's GFF phase) -> residue L // 3, phase L % 3 (phase 0 = between codons).

Protein -> CDS: by `protein_id`; otherwise (C. elegans WBGene...1 IDs, or a miss) gene ->
transcripts (mRNA/transcript lines' Parent) and the transcript whose CDS length matches
3 x len(protein) (or + 3 with the stop codon), else the closest length (flagged).

Output: work/l2_genome/introns.parquet (protein, species, transcript, how, cds_nt, prot_len,
len_match, n_introns, pos[list], phase[list]) and work/l2_genome/introns_qc.json.
"""
import gzip, json, os, sys, time
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GFF = f"{ROOT}/data/orthobench/ensembl_99/gff3"
INP = f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/Input"
OUT = f"{ROOT}/work/l2_genome"
FILES = {   # alt-haplotype proteins are in OrthoBench -> full GFF for human/mouse/zebrafish
    "Caenorhabditis_elegans": "Caenorhabditis_elegans.WBcel235.99.gff3.gz",
    "Canis_familiaris": "Canis_familiaris.CanFam3.1.99.gff3.gz",
    "Ciona_intestinalis": "Ciona_intestinalis.KH.99.gff3.gz",
    "Danio_rerio": "Danio_rerio.GRCz11.99.chr_patch_hapl_scaff.gff3.gz",
    "Drosophila_melanogaster": "Drosophila_melanogaster.BDGP6.28.99.gff3.gz",
    "Gallus_gallus": "Gallus_gallus.GRCg6a.99.gff3.gz",
    "Homo_sapiens": "Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gff3.gz",
    "Monodelphis_domestica": "Monodelphis_domestica.ASM229v1.99.gff3.gz",
    "Mus_musculus": "Mus_musculus.GRCm38.99.chr_patch_hapl_scaff.gff3.gz",
    "Pan_troglodytes": "Pan_troglodytes.Pan_tro_3.0.99.gff3.gz",
    "Rattus_norvegicus": "Rattus_norvegicus.Rnor_6.0.99.gff3.gz",
    "Tetraodon_nigroviridis": "Tetraodon_nigroviridis.TETRAODON8.99.gff3.gz",
}


def attr(s, key):
    i = s.find(key + "=")
    if i < 0:
        return None
    return s[i + len(key) + 1:].split(";")[0].strip()


def read_fasta(path, want):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if name in want:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name in want:
        seqs[name] = "".join(buf)
    return seqs


def parse_gff(path):
    cds = defaultdict(list)       # transcript -> [(start, end, strand, phase)]
    t_prot, t_gene, g_ts = {}, {}, defaultdict(list)
    for line in gzip.open(path, "rt"):
        if line[0] == "#":
            continue
        x = line.rstrip("\n").split("\t")
        if len(x) < 9:
            continue
        a = x[8]
        if x[2] == "CDS":
            par = attr(a, "Parent")
            if not par:
                continue
            t = par.split(",")[0].replace("transcript:", "")
            ph = int(x[7]) if x[7] in "012" else 0
            cds[t].append((int(x[3]), int(x[4]), x[6], ph))
            p = attr(a, "protein_id")
            if p:
                t_prot[t] = p
        else:
            i, par = attr(a, "ID"), attr(a, "Parent")
            if i and i.startswith("transcript:") and par and par.startswith("gene:"):
                t = i[len("transcript:"):]
                gid = par[len("gene:"):]
                t_gene[t] = gid
                g_ts[gid].append(t)
    return cds, t_prot, t_gene, g_ts


def introns_of(segs):
    strand = segs[0][2]
    segs = sorted(segs, key=lambda s: s[0], reverse=(strand == "-"))
    lens = [e - s + 1 for s, e, _, _ in segs]
    ph0 = segs[0][3]
    cum = np.cumsum(lens)[:-1] - ph0
    total = int(sum(lens) - ph0)
    pos = [int(L // 3) for L in cum]
    phase = [int(L % 3) for L in cum]
    return pos, phase, total


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    lab = pd.read_parquet(f"{ROOT}/work/S5/refog_pair_labels.parquet")
    genes = pd.concat([lab[["gene_a", "sp_a"]].set_axis(["g", "sp"], axis=1),
                       lab[["gene_b", "sp_b"]].set_axis(["g", "sp"], axis=1)]).drop_duplicates()
    rows, qc = [], {}
    for sp, f in FILES.items():
        want = set(genes[genes.sp == sp].g)
        fa = [x for x in os.listdir(INP) if x.startswith(sp)][0]
        seqs = read_fasta(f"{INP}/{fa}", want)
        cds, t_prot, t_gene, g_ts = parse_gff(f"{GFF}/{f}")
        p_t = {p: t for t, p in t_prot.items()}
        go = f"{ROOT}/data/genomes/orthobench/{sp}.parquet"
        p2g = dict(pd.read_parquet(go, columns=["protein_id", "gene_id"]).values) if os.path.exists(go) else {}
        c = dict(n=len(want), by_protein_id=0, by_gene_len=0, by_gene_closest=0, unmapped=0, len_match=0)
        for p in sorted(want):
            L = len(seqs.get(p, "").rstrip("*"))
            t, how = p_t.get(p), "protein_id"
            if t is None:
                gid = p2g.get(p) or (p[:-2] if p.endswith(".1") else p)
                ts = [x for x in g_ts.get(gid, []) if x in cds]
                best = None
                for x in ts:
                    tot = introns_of(cds[x])[2]
                    d = min(abs(tot - 3 * L), abs(tot - 3 * L - 3))
                    if best is None or d < best[0]:
                        best = (d, x)
                if best is None:
                    c["unmapped"] += 1
                    continue
                t, how = best[1], ("gene_len" if best[0] == 0 else "gene_closest")
            pos, phase, tot = introns_of(cds[t])
            match = tot in (3 * L, 3 * L + 3)
            c["by_" + how] += 1
            c["len_match"] += int(match)
            rows.append(dict(protein=p, species=sp, transcript=t, how=how, cds_nt=tot, prot_len=L,
                             len_match=bool(match), n_introns=len(pos), pos=pos, phase=phase))
        c["mapped_share"] = round(1 - c["unmapped"] / max(1, c["n"]), 4)
        c["len_match_share"] = round(c["len_match"] / max(1, c["n"] - c["unmapped"]), 4)
        qc[sp] = c
        print(sp, c, flush=True)
    d = pd.DataFrame(rows)
    d.to_parquet(f"{OUT}/introns.parquet", index=False)
    qc["_total"] = dict(n=int(len(genes)), mapped=int(len(d)), len_match=int(d.len_match.sum()),
                        mean_introns=round(float(d.n_introns.mean()), 2), wall_s=round(time.time() - t0, 1))
    json.dump(qc, open(f"{OUT}/introns_qc.json", "w"), indent=1)
    print(qc["_total"])


if __name__ == "__main__":
    main()
