#!/usr/bin/env bash
# Extra data layers (expression, structure) — small parts downloaded, large parts only sized.
# Usage: download_extra_layers.sh [bgee|afdb_sizes|afdb <UPID_TAX_CODE>...]
#   bgee        : Bgee v15.2 "expr_simple" present/absent expression calls for the 10 OrthoBench
#                 species covered by Bgee (Tetraodon and Ciona are not in Bgee). ~460 MB.
#   afdb_sizes  : list AlphaFold DB v6 per-proteome tar sizes for QfO model organisms (no download)
#   afdb X      : download one AFDB proteome tar, e.g. UP000002311_559292_YEAST_v6 (~1 GB)
set -euo pipefail
source "$(dirname "$0")/config.sh"
what="${1:-bgee}"; shift || true
case "$what" in
  bgee)
    OUT="$DATA_DIR/extra/bgee_v15_2"; mkdir -p "$OUT"; cd "$OUT"
    for sp in Homo_sapiens Mus_musculus Rattus_norvegicus Canis_lupus_familiaris Pan_troglodytes \
              Monodelphis_domestica Gallus_gallus Danio_rerio Drosophila_melanogaster Caenorhabditis_elegans; do
      f=${sp}_expr_simple.tsv.gz
      [[ -s $f ]] || wget -q "https://www.bgee.org/ftp/bgee_v15_2/download/calls/expr_calls/$f" \
                  || wget -q "https://www.bgee.org/ftp/current/download/calls/expr_calls/$f" || echo "FAILED $f"
    done
    md5sum *.tsv.gz > md5sums.txt; ls -la ;;
  afdb_sizes)
    curl -s https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/ | sed 's/<[^>]*>/ /g' | grep -E '_v[0-9]+\.tar' \
      | awk '{print $1"\t"$4}' ;;
  afdb)
    OUT="$DATA_DIR/extra/afdb"; mkdir -p "$OUT"; cd "$OUT"
    for x in "$@"; do wget -c -q "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/$x.tar"; done; ls -la ;;
esac
