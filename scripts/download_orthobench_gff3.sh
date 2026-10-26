#!/usr/bin/env bash
# Download Ensembl GFF3 (gene coordinates, for synteny) + DNA-independent annotation for the
# 12 OrthoBench species, from the Ensembl release whose assemblies match the OrthoBench
# input proteomes (default: release 99, Jan 2020 — OrthoBench proteomes were downloaded March 2020;
# release 99 gives 99.99% protein-ID coverage vs 99.8% for release 100; assemblies GRCh38, GRCm38, Rnor_6.0,
# CanFam3.1, Pan_tro_3.0, ASM229v1, GRCg6a, GRCz11, TETRAODON8, KH, BDGP6.28, WBcel235).
# Usage: download_orthobench_gff3.sh [RELEASE]   -> $DATA_DIR/orthobench/ensembl_<rel>/gff3/
set -euo pipefail
source "$(dirname "$0")/config.sh"
REL="${1:-99}"
OUT="$DATA_DIR/orthobench/ensembl_$REL/gff3"; mkdir -p "$OUT"; cd "$OUT"
declare -A ASM=( [homo_sapiens]=Homo_sapiens.GRCh38 [mus_musculus]=Mus_musculus.GRCm38
  [rattus_norvegicus]=Rattus_norvegicus.Rnor_6.0 [canis_familiaris]=Canis_familiaris.CanFam3.1
  [pan_troglodytes]=Pan_troglodytes.Pan_tro_3.0 [monodelphis_domestica]=Monodelphis_domestica.ASM229v1
  [gallus_gallus]=Gallus_gallus.GRCg6a [danio_rerio]=Danio_rerio.GRCz11
  [tetraodon_nigroviridis]=Tetraodon_nigroviridis.TETRAODON8 [ciona_intestinalis]=Ciona_intestinalis.KH
  [drosophila_melanogaster]=Drosophila_melanogaster.BDGP6.28 [caenorhabditis_elegans]=Caenorhabditis_elegans.WBcel235 )
# Human, mouse and zebrafish OrthoBench proteomes include genes on alternative haplotypes /
# patches (e.g. HSCHR6_MHC_*), which are only in the *.chr_patch_hapl_scaff.gff3.gz files.
for sp in homo_sapiens mus_musculus danio_rerio; do
  f="${ASM[$sp]}.$REL.chr_patch_hapl_scaff.gff3.gz"
  [[ -s "$f" ]] || wget -q "https://ftp.ensembl.org/pub/release-$REL/gff3/$sp/$f" || echo "FAILED $sp $f"
done
for sp in "${!ASM[@]}"; do
  f="${ASM[$sp]}.$REL.gff3.gz"
  [[ -s "$f" ]] && continue
  wget -q "https://ftp.ensembl.org/pub/release-$REL/gff3/$sp/$f" && continue
  # dog moved to canis_lupus_familiaris in Ensembl >= 99
  if [[ $sp == canis_familiaris ]]; then
    g="Canis_lupus_familiaris.CanFam3.1.$REL.gff3.gz"
    wget -q "https://ftp.ensembl.org/pub/release-$REL/gff3/canis_lupus_familiaris/$g" && ln -sf "$g" "$f" && continue
  fi
  echo "FAILED $sp $f"
done
wget -q -O CHECKSUMS_note.txt "https://ftp.ensembl.org/pub/release-$REL/gff3/homo_sapiens/CHECKSUMS" || true
md5sum *.gff3.gz > md5sums.txt
ls -la
