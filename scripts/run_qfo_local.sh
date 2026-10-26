#!/usr/bin/env bash
# Run the Quest-for-Orthologs benchmarks LOCALLY, without docker or nextflow.
#
# Re-implements the orchestration of qfo/benchmark-webservice main.nf (DSL1, docker) as a
# plain bash script calling the *unmodified* upstream code natively:
#   validate.py -> map_relations.py -> {GoTest,EcTest,RefPhyloTest,SpeciesTreeDiscordanceTest}.sh
#   (Darwin, static linux64 binary from biorecipes.com) + vgnc_benchmark.py / swissprot_benchmark.py
#
# Usage: run_qfo_local.sh -i predictions.(orthoxml[.gz]|tsv[.gz]) -n METHOD_NAME [-o OUTDIR]
#                         [-y 2022.1] [-b "GO EC SwissTrees ..."] [-j NPAR]
# Env:   QFO_WS     path to benchmark-webservice clone (default $DATA_DIR/qfo/benchmark-webservice)
#        QFO_REF    path to reference data dir           (default $DATA_DIR/qfo/reference_data/<year>)
#        DARWIN_DIR path containing the 'darwin' wrapper (default $DATA_DIR/tools/darwin)
#        PYTHON     python with lxml, dendropy, biopython, jsonschema, pandas (default conda_envs/data)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"

YEAR=2022.1
BENCHES="GO EC VGNC SwissTrees TreeFam-A STD_Eukaryota STD_Fungi STD_Bacteria G_STD2_Luca G_STD2_Fungi G_STD2_Eukaryota G_STD2_Vertebrata"
NPAR=8
OUT=""
INPUT=""; NAME=""
while getopts "i:n:o:y:b:j:h" opt; do
  case $opt in
    i) INPUT="$(readlink -f "$OPTARG")";; n) NAME="$OPTARG";; o) OUT="$OPTARG";;
    y) YEAR="$OPTARG";; b) BENCHES="$OPTARG";; j) NPAR="$OPTARG";;
    h|*) sed -n 2,20p "$0"; exit 0;;
  esac
done
[[ -z "$INPUT" || -z "$NAME" ]] && { sed -n 2,20p "$0"; exit 1; }
QFO_WS="${QFO_WS:-$DATA_DIR/qfo/benchmark-webservice}"
QFO_REF="$(readlink -f "${QFO_REF:-$DATA_DIR/qfo/reference_data/$YEAR}")"
DARWIN_DIR="${DARWIN_DIR:-$DATA_DIR/tools/darwin}"
# darwin on PATH; mafft (needed by G_STD2_* tests) is expected next to $PYTHON (conda_envs/data)
export PATH="$DARWIN_DIR:$(dirname "$PYTHON"):$PATH"
export DARWIN_ORTHOLOG_BENCHMARK_REPO_PATH="$(readlink -f "$QFO_WS")"
METHOD="${NAME// /_}"
OUT="$(readlink -f -m "${OUT:-$BENCH_ROOT/results/qfo/$METHOD}")"
mkdir -p "$OUT/logs" "$OUT/raw" "$OUT/assess" "$OUT/tmp"
# sqlite index creation for large inputs overflows small /tmp on cluster nodes
export TMPDIR="$OUT/tmp" SQLITE_TMPDIR="$OUT/tmp"
cd "$OUT"
ulimit -s unlimited || true
command -v mafft >/dev/null || echo "[qfo] WARNING: mafft not on PATH; G_STD2_* tests will fail"
command -v darwin >/dev/null || { echo "darwin not found in $DARWIN_DIR (run scripts/download_qfo.sh)"; exit 1; }
echo "[qfo] input=$INPUT method=$METHOD refset=$QFO_REF out=$OUT"

t0=$(date +%s)
# 1. validate
$PYTHON "$QFO_WS/validate.py" --com QfO --challenges_ids "$BENCHES" --participant "$METHOD" \
    --out participant.json "$QFO_REF/mapping.json.gz" "$INPUT" > logs/validate.log 2>&1 \
  || { echo "[qfo] VALIDATION FAILED, see $OUT/logs/validate.log"; tail -5 logs/validate.log; exit 2; }
# 2. convert predictions to darwin db + sqlite
( cd "$QFO_WS" && $PYTHON map_relations.py --out "$OUT/predictions.db" --db "$OUT/orthologs.db" \
    "$QFO_REF/mapping.json.gz" "$INPUT" ) > logs/map_relations.log 2>&1
echo "[qfo] converted in $(( $(date +%s)-t0 ))s"

# 3. benchmarks (each writes <bench>.json assessment stub in assess/)
run_one() {
  local b="$1" db="$OUT/predictions.db" ws="$QFO_WS" ref="$QFO_REF"
  case "$b" in
    GO)  bash "$ws/GoTest.sh" -o "$OUT/raw/GO" -a "assess/GO.json" -c QfO -e exp "$db" "$METHOD" "$ref";;
    EC)  bash "$ws/EcTest.sh" -o "$OUT/raw/EC" -a "assess/EC.json" -c QfO "$db" "$METHOD" "$ref";;
    SwissTrees|TreeFam-A)
         bash "$ws/RefPhyloTest.sh" -o "$OUT/raw/$b" -a "assess/$b.json" -t "$b" "$db" "$METHOD" "$ref";;
    STD_*)    bash "$ws/SpeciesTreeDiscordanceTest.sh" -o "$OUT/raw/$b" -a "assess/$b.json" -c QfO -p "${b#STD_}" -m 0 "$db" "$METHOD" "$ref";;
    G_STD2_*) bash "$ws/SpeciesTreeDiscordanceTest.sh" -o "$OUT/raw/$b" -a "assess/$b.json" -c QfO -p "${b#G_STD2_}" -m 2 "$db" "$METHOD" "$ref";;
    G_STD_*)  bash "$ws/SpeciesTreeDiscordanceTest.sh" -o "$OUT/raw/$b" -a "assess/$b.json" -c QfO -p "${b#G_STD_}" -m 1 "$db" "$METHOD" "$ref";;
    VGNC) ( cd "$ws" && $PYTHON vgnc_benchmark.py --com QfO --participant "$METHOD" --assessment-out "$OUT/assess/VGNC.json" \
             --outdir "$OUT/raw/VGNC" --vgnc-orthologs "$ref/vgnc-orthologs.txt.gz" --db "$OUT/orthologs.db" );;
    SP)   ( cd "$ws" && $PYTHON swissprot_benchmark.py --com QfO --participant "$METHOD" --assessment-out "$OUT/assess/SP.json" \
             --outdir "$OUT/raw/SwissProtIDs" --mapping "$ref/mapping.json.gz" --sp-entries "$ref/swissprot.txt.gz" \
             --strategy ids_exist_in_both --lineage-tree "$ref/lineage_tree.phyloxml" --db "$OUT/orthologs.db" );;
    *) echo "unknown benchmark $b"; return 1;;
  esac
}
export -f run_one; export OUT QFO_WS QFO_REF METHOD PYTHON PATH DARWIN_ORTHOLOG_BENCHMARK_REPO_PATH
for b in $BENCHES; do
  rm -f "assess/$b.json"
  ( s=$(date +%s); run_one "$b" > "logs/$b.log" 2>&1 || true   # darwin scripts may exit non-zero on success (RefPhyloTest uses exit(1))
      if [[ -s "assess/$b.json" ]]; then
      echo "[qfo] $b done ($(( $(date +%s)-s ))s)"; else echo "[qfo] $b FAILED (see $OUT/logs/$b.log)"; fi ) &
  while (( $(jobs -rp | wc -l) >= NPAR )); do sleep 2; done
done
wait
# 4. summary table
$PYTHON "$HERE/summarize_qfo.py" "$OUT/assess" > summary.tsv
cat summary.tsv
echo "[qfo] total $(( $(date +%s)-t0 ))s; results in $OUT"
