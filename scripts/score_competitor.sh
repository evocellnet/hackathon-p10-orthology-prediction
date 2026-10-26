#!/usr/bin/env bash
# Usage: score_competitor.sh <method> <orthogroup-table> <runs-logs-dir> <threads> "<note>" "<version string>" [extra logs dirs...]
# Converts to OrthoBench format (data/orthobench/runs/<method>.og), scores with run_orthobench.sh (official benchmark.py)
# and orthobench_eval.py (per-RefOG, dev/test), writes results/orthobench_competitors/<method>/{score.json,per_refog.tsv,runtime.json,version.txt,time.txt,load.log}
set -euo pipefail
M="$1"; TAB="$2"; LOGS="$3"; THR="$4"; NOTE="$5"; VER="$6"; shift 6
cd "$(dirname "$0")/.."
source scripts/config.sh
C=results/orthobench_competitors/$M; mkdir -p "$C"
OG="$DATA_DIR/orthobench/runs/$M.og"
$PYTHON scripts/to_orthobench_og.py "$TAB" "$OG"
bash scripts/run_orthobench.sh "$OG" "competitor_$M" | tail -2
$PYTHON -W ignore scripts/orthobench_eval.py "$DATA_DIR/orthobench/Open_Orthobench/BENCHMARKS" "$OG" --label "$M" --json "$C/score.json" --per-refog-tsv "$C/per_refog.tsv" | grep -v "^Examples\|^[A-Z]*[0-9]*$"
INC=(); [ $# -gt 0 ] && INC=(--include "$@")
$PYTHON scripts/make_runtime_json.py "$LOGS" "$C/runtime.json" --threads "$THR" --note "$NOTE" "${INC[@]}" > /dev/null
echo "$VER" > "$C/version.txt"
cp "$LOGS/time.txt" "$LOGS/load.log" "$C/"
