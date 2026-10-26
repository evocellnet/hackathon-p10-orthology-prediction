#!/usr/bin/env bash
# Quest for Orthologs benchmark data, everything needed to run the benchmarks locally:
#  1. QfO reference proteomes 2022_02 (the release used by the public benchmark service, with the
#     updated zebrafish proteome UP000000437 — verified: the 2022 benchmark mapping has 26,355
#     DANRE proteins = updated release) incl. gene2acc, _DNA.fasta, _additional.fasta (isoforms),
#     .idmapping, .xml, _gcoord.xml (genome coordinates, 20 species)
#  2. benchmark-webservice code (pinned commit) and Darwin (static binary + library)
#  3. reference data for the benchmark (Zenodo 10.5281/zenodo.10517603, v2022.1; optional 2018)
# Usage: download_qfo.sh [all|proteomes|code|refdata|refdata2018|current|gff_current]
set -euo pipefail
source "$(dirname "$0")/config.sh"
Q="$DATA_DIR/qfo"; mkdir -p "$Q"
WS_COMMIT=c0854a96c1a0fd7f2a891d971af0863002fabc90   # qfo/benchmark-webservice master, 2025-03-26
FTP=https://ftp.ebi.ac.uk/pub/databases/reference_proteomes
what="${1:-all}"
if [[ $what == all || $what == proteomes ]]; then
  R=qfo_release-2022_02_with_updated_UP000000437
  mkdir -p "$Q/proteomes_2022_02_updZEBRAFISH"; cd "$Q/proteomes_2022_02_updZEBRAFISH"
  for f in README STATS QfO_release_2022_02_with_updated_UP000000437.tar.gz; do wget -c -q "$FTP/previous_releases/$R/$f"; done
  [[ -d Eukaryota ]] || tar xzf QfO_release_2022_02_with_updated_UP000000437.tar.gz
fi
if [[ $what == current ]]; then   # newest QfO release (2026_02, 81 species) — not used by the 2022 benchmark
  mkdir -p "$Q/proteomes_current"; cd "$Q/proteomes_current"
  for f in README STATS PROTEOME_ASSEMBLY_MAPPING; do wget -c -q "$FTP/QfO/$f"; done
  echo "tarballs: $FTP/QfO/QfO_release_2026_02.tar.gz (2.2G), QfO_release_2026_02_additional.tar.gz (20G)"
fi
if [[ $what == all || $what == code ]]; then
  [[ -d "$Q/benchmark-webservice/.git" ]] || git clone https://github.com/qfo/benchmark-webservice "$Q/benchmark-webservice"
  git -C "$Q/benchmark-webservice" checkout -q "$WS_COMMIT"
  T="$DATA_DIR/tools/darwin"; mkdir -p "$T"; cd "$T"
  [[ -s darwin.linux64 ]] || { curl -sSL -O https://biorecipes.com/darwin/darwin.linux64; chmod +x darwin.linux64; }
  [[ -d lib ]] || { curl -sSL -O https://biorecipes.com/darwin/darwin-lib.tgz; tar xzf darwin-lib.tgz; \
      sed -i -e "s/\(datadirname\s*:=\s*\)'[^']*'[:;]/\1'\/refset':/" lib/darwinit; }   # as in Dockerfile_qfo_darwin
  printf '#!/bin/bash\nulimit -s unlimited\nD=$(dirname $(readlink -f $0))\nexec $D/darwin.linux64 -l $D/lib "$@"\n' > darwin
  chmod +x darwin
  md5sum darwin.linux64 darwin-lib.tgz > md5sums.txt
fi
if [[ $what == all || $what == refdata ]]; then
  cd "$Q/benchmark-webservice"
  $PYTHON fetch_reference_data.py 2022.1 --out-dir "$Q/reference_data/2022.1"
  ( cd "$Q/reference_data/2022.1" && md5sum $(ls | grep -v fas_annotations$) > ../md5_2022.1.txt )
fi
if [[ $what == refdata2018 ]]; then   # used for the positive control with the shipped OMA example file
  cd "$Q/benchmark-webservice"; $PYTHON fetch_reference_data.py 2018 --out-dir "$Q/reference_data/2018"
fi
if [[ $what == gff_current ]]; then  # GFF/GFF3 of the current (2026_02) QfO release, 74 files, ~490 MB
  # The 2022_02 release ships no GFF (only _gcoord.xml for 19 Ensembl eukaryotes); these GFFs belong to
  # the 2026 proteome versions: map to 2022 proteins via UniProt accession / Ensembl IDs (.idmapping).
  mkdir -p "$Q/proteomes_current/gff"; cd "$Q/proteomes_current/gff"
  for d in Archaea Bacteria Eukaryota; do
    for f in $(curl -s "$FTP/QfO/$d/" | grep -o 'UP[0-9]*\.gff3\?\.gz' | sort -u); do wget -c -q "$FTP/QfO/$d/$f"; done
  done
  md5sum *.gz > md5sums.txt
fi
echo "done: $what"
