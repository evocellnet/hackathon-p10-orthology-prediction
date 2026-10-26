#!/usr/bin/env python
"""S7 prep: all Ensembl 99 isoforms of the OrthoBench DEV-RefOG genes + all-vs-all DIAMOND.

* Genes = every gene in work/S5/refog_pair_labels.parquet (DEV RefOGs, scripts/S5_refog_labels.py).
* OrthoBench input protein -> Ensembl 99 gene (pep.all header; C. elegans inputs are
  WBGene IDs + '.1'); the canonical protein is the isoform whose sequence equals the
  OrthoBench input sequence (checked for every species).
* Genome-wide isoform statistics per species (how often the OrthoBench protein is the longest).
* DIAMOND blastp --sensitive -k 0 --max-hsps 1 -e 1e-3 (all-vs-all over the isoform set;
  identical sequences inside a gene are collapsed onto one representative).
Outputs: work/S7/{isoforms.parquet, isoforms.fa, hits.tsv, prep_stats.json}
"""
import glob, json, os, subprocess, sys, time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.isoforms import parse_ensembl_pep  # noqa: E402

IN = f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/Input"
PEP = f"{ROOT}/data/orthobench/ensembl_99/pep"
W = f"{ROOT}/work/S7"


def read_fasta(path):
    out, k, s = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if k:
                out[k] = "".join(s)
            k, s = l[1:].split()[0], []
        else:
            s.append(l.strip())
    if k:
        out[k] = "".join(s)
    return out


def main(threads=6):
    os.makedirs(W, exist_ok=True)
    t0 = time.time()
    lab = pd.read_parquet(f"{ROOT}/work/S5/refog_pair_labels.parquet")
    gsp = pd.concat([lab[["gene_a", "sp_a"]].set_axis(["ob", "sp"], axis=1),
                     lab[["gene_b", "sp_b"]].set_axis(["ob", "sp"], axis=1)]).drop_duplicates()
    stats, rows = {}, []
    for f in sorted(glob.glob(f"{IN}/*.pep.all.fa")):
        sp = os.path.basename(f).split(".")[0]
        ob = read_fasta(f)
        pep = parse_ensembl_pep(glob.glob(f"{PEP}/{sp}.*.pep.all.fa.gz")[0])
        if sp == "Caenorhabditis_elegans":
            ob_gene = {k: k.rsplit(".", 1)[0] for k in ob}
        else:
            p2g = dict(zip(pep.protein_id, pep.gene_id))
            ob_gene = {k: p2g.get(k) for k in ob}
        # canonical isoform = same sequence within the gene (prefer same ID)
        pep["ob"] = None
        by_gene = {g: d for g, d in pep.groupby("gene_id")}
        n_seq_match = n_longest = n_multi = 0
        canon = {}
        for k, s in ob.items():
            g = ob_gene.get(k)
            if g is None or g not in by_gene:
                continue
            d = by_gene[g]
            m = d[d.seq == s]
            if len(m):
                n_seq_match += 1
                pick = m[m.protein_id == k].index if (m.protein_id == k).any() else m.index[:1]
                canon[k] = pep.loc[pick[0], "protein_id"]
            n_multi += int(d.seq.nunique() > 1)
            n_longest += int(len(s) >= d.length.max())
        stats[sp] = dict(ob_proteins=len(ob), mapped_to_gene=sum(v is not None and v in by_gene for v in ob_gene.values()),
                         canonical_seq_matched=n_seq_match, genes_with_multiple_distinct_isoforms=n_multi,
                         ob_is_longest=n_longest, pep_all_proteins=len(pep), pep_all_genes=int(pep.gene_id.nunique()))
        print(sp, stats[sp], flush=True)
        for k in gsp[gsp.sp == sp].ob:
            g = ob_gene.get(k)
            if g is None or g not in by_gene:
                rows.append((sp, k, g, k, True, len(ob.get(k, "")), ob.get(k, "")))
                continue
            d = by_gene[g].drop_duplicates("seq")
            c = canon.get(k)
            cs = ob[k]
            rows.append((sp, k, g, k, True, len(cs), cs))           # canonical = OrthoBench ID + seq
            for pid, s in zip(d.protein_id, d.seq):
                if s != cs:
                    rows.append((sp, k, g, pid, False, len(s), s))
    iso = pd.DataFrame(rows, columns=["species", "ob_id", "gene_id", "protein_id", "canonical", "length", "seq"])
    iso = iso.drop_duplicates(["ob_id", "protein_id"])
    iso["uid"] = [f"{o}|{p}" for o, p in zip(iso.ob_id, iso.protein_id)]
    iso.to_parquet(f"{W}/isoforms.parquet", index=False)
    with open(f"{W}/isoforms.fa", "w") as fh:
        for u, s in zip(iso.uid, iso.seq):
            if s:
                fh.write(f">{u}\n{s}\n")
    per = iso.groupby("ob_id").size()
    stats["_label_genes"] = dict(n_genes=int(len(per)), n_proteins=int(len(iso)),
                                 genes_with_alt_isoforms=int((per > 1).sum()),
                                 mean_isoforms=round(float(per.mean()), 3), max_isoforms=int(per.max()))
    print(stats["_label_genes"], flush=True)
    cmd = (f"diamond makedb --in {W}/isoforms.fa -d {W}/isoforms --quiet && "
           f"/usr/bin/time -v diamond blastp -q {W}/isoforms.fa -d {W}/isoforms -o {W}/hits.tsv --sensitive -k 0 "
           f"--max-hsps 1 -e 1e-3 -p {threads} --quiet --outfmt 6 qseqid sseqid bitscore pident length qlen slen "
           f"2> {W}/diamond.time")
    subprocess.run(cmd, shell=True, check=True)
    stats["_diamond"] = dict(cmd=cmd, wall_s=round(time.time() - t0, 1))
    json.dump(stats, open(f"{W}/prep_stats.json", "w"), indent=1)


if __name__ == "__main__":
    main()
