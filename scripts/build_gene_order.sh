#!/usr/bin/env bash
# Build gene-order tables (S5 synteny input) under data/genomes/.
#   data/genomes/orthobench/<species>.parquet  (Ensembl 99 GFF3; OrthoBench 12 species)
#   data/genomes/qfo2022/<UPid_taxid>.parquet  (QfO 2022 gcoord.xml, then 2026 GFF where seq MD5 equal)
# Schema: ancestra/genomes.py (COLUMNS). Coverage reports: coverage.tsv in each directory.
set -euo pipefail
cd "$(dirname "$0")/.."
source /opt/conda/etc/profile.d/conda.sh
conda activate "$PWD/conda_envs/orth"
THREADS=${THREADS:-6} python scripts/build_gene_order_orthobench.py

# QfO 2022_02 (updated zebrafish), 78 species: 2022 _gcoord.xml first, then 2026_02 GFF/GBFF
# (via 2026 canonical FASTA MD5) and 2026 _gcoord.xml (via its sequence MD5), never mixing
# assemblies. Needs data/qfo/proteomes_current/{Archaea,Bacteria,Eukaryota} extracted from
# QfO_release_2026_02.tar.gz and gbff/ (see README_local_2026_tarball.txt there). ~15 min, 4 procs.
python -m ancestra.genomes_qfo --out data/genomes/qfo2022 --threads "${QFO_THREADS:-4}"
