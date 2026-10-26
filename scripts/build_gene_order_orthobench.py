#!/usr/bin/env python
"""Build OrthoBench gene-order tables (data/genomes/orthobench/<species>.parquet).

Ensembl 99 GFF3: the chr_patch_hapl_scaff file where present (alt-haplotype genes
are in OrthoBench), primary flag = contig in the primary-assembly GFF3.
Writes coverage.tsv + coverage.json next to the tables.
"""
import glob, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.genomes import build_ensembl_species, write_gene_order  # noqa: E402

INP = f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/Input"
GFF = f"{ROOT}/data/orthobench/ensembl_99/gff3"
OUT = f"{ROOT}/data/genomes/orthobench"


def one(fa):
    sp = os.path.basename(fa).split(".")[0]
    full = glob.glob(f"{GFF}/{sp}.*.99.chr_patch_hapl_scaff.gff3.gz")
    prim = [p for p in glob.glob(f"{GFF}/{sp}.*.99.gff3.gz") if "chr_patch" not in p]
    gfull = full[0] if full else prim[0]
    gprim = prim[0] if full else None
    t0 = time.time()
    df, cov = build_ensembl_species(fa, gfull, gprim)
    write_gene_order(df, f"{OUT}/{sp}.parquet")
    cov.update(species=sp, gff=os.path.basename(gfull),
               gff_primary=os.path.basename(gprim) if gprim else None, wall_s=round(time.time() - t0, 1))
    return cov


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fas = sorted(glob.glob(f"{INP}/*.pep.all.fa"))
    with ProcessPoolExecutor(int(os.environ.get("THREADS", 6))) as ex:
        covs = list(ex.map(one, fas))
    with open(f"{OUT}/coverage.json", "w") as fh:
        json.dump(covs, fh, indent=1)
    cols = ["species", "n_proteins", "n_mapped", "coverage", "n_mapped_via_gene_id", "n_non_primary", "n_projected", "n_non_primary_unprojected",
            "n_contigs", "n_genes_on_contigs_ge10", "gff"]
    with open(f"{OUT}/coverage.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in covs:
            fh.write("\t".join(str(c[k]) for k in cols) + "\n")
    for c in covs:
        print(c["species"], c["n_mapped"], "/", c["n_proteins"], c["coverage"], "alt", c["n_non_primary"], "proj", c["n_projected"],
              "missing", c["missing_examples"][:3])
