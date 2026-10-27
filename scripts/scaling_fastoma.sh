#!/usr/bin/env bash
# FastOMA 0.5.1 (nextflow, local executor, 8 cpus) on nested subset N; full workflow incl. OMAmer placement.
# Same call as the OrthoBench competitor run (data/orthobench/runs/fastoma/logs/time.txt). Usage: scaling_fastoma.sh N RUNDIR
set -eu
N="$1"; RUN="$(realpath -m "$2")"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FO="$ROOT/data/orthobench/runs/fastoma"
IN="$ROOT/work/scaling/inputs/fo_$N"
mkdir -p "$RUN"; cd "$RUN"
export NXF_HOME="$FO/nxf_home"
exec nextflow run "$FO/src/FastOMA.nf" -c "$FO/local.config" --input "$IN" --proteome_folder "$IN/proteome" \
  --species_tree "$IN/species_tree.nwk" --omamer_db "$FO/db/LUCA.h5" --output_folder "$RUN/Output" -with-trace -with-report
