#!/usr/bin/env bash
# Score orthogroups against the 70 Open_Orthobench RefOGs with the unmodified upstream benchmark.py
# Usage: run_orthobench.sh <orthogroups file> [label] [--basic]
#   orthogroups file: one orthogroup per line (OrthoFinder Orthogroups.txt/.tsv, or any format whose
#   lines contain the Ensembl/FlyBase/WormBase protein IDs of BENCHMARKS/Input). Output JSON is
#   written to $BENCH_ROOT/results/orthobench/<label>.json
set -euo pipefail
source "$(dirname "$0")/config.sh"
OGS="$1"; LABEL="${2:-$(basename "$OGS")}"; shift; shift || true
B="${ORTHOBENCH_DIR:-$DATA_DIR/orthobench/Open_Orthobench/BENCHMARKS}"
mkdir -p "$BENCH_ROOT/results/orthobench"
$PYTHON -W ignore "$BENCH_ROOT/scripts/orthobench_eval.py" "$B" "$OGS" --label "$LABEL" \
    --json "$BENCH_ROOT/results/orthobench/$LABEL.json" "$@"
