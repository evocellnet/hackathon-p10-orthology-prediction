#!/usr/bin/env bash
# Re-run all benchmark controls end-to-end and store small summaries in results/controls/.
# Usage: run_controls.sh [orthobench|qfo_fast|qfo_full|all]
#   qfo_fast: oracle / anti-oracle (SwissTrees, TreeFam-A) + random pairs on fast tests (~10 min)
#   qfo_full: also all STD/G_STD2 tests on random pairs, 2018 OMA example (exact reproduction),
#             FastOMA & OrthoFinder 2022 submissions vs published (several hours)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; source "$HERE/config.sh"
RC="$BENCH_ROOT/results/controls"; mkdir -p "$RC"
R="$DATA_DIR/qfo/runs"; C="$DATA_DIR/qfo/controls"; mkdir -p "$C"
what="${1:-all}"
FAST="GO EC VGNC SP SwissTrees TreeFam-A"
ALL="$FAST STD_Eukaryota STD_Fungi STD_Bacteria G_STD2_Luca G_STD2_Fungi G_STD2_Eukaryota G_STD2_Vertebrata"
if [[ $what == orthobench || $what == all ]]; then
  B="$DATA_DIR/orthobench/Open_Orthobench/BENCHMARKS"
  $PYTHON "$BENCH_ROOT/controls/make_orthobench_controls.py" "$B" "$DATA_DIR/orthobench/controls"
  for f in "$DATA_DIR"/orthobench/controls/*.txt; do
    bash "$HERE/run_orthobench.sh" "$f" "control_$(basename "$f" .txt)" | tail -1
  done | sort | { echo -e "label\tF\tprecision\trecall\texact_refogs"; cat; } > "$RC/orthobench_controls.tsv"; cat "$RC/orthobench_controls.tsv"
fi
if [[ $what == qfo_fast || $what == qfo_full || $what == all ]]; then
  REF="$DATA_DIR/qfo/reference_data/2022.1"
  $PYTHON "$BENCH_ROOT/controls/make_qfo_random_pairs.py" "$REF/mapping.json.gz" "$C/random_pairs_2022.tsv" --n 2000000 --seed 42
  BEN="$FAST"
  [[ $what == qfo_full || $what == all ]] && BEN="$ALL"
  bash "$HERE/run_qfo_local.sh" -i "$C/random_pairs_2022.tsv" -n random_pairs -o "$R/random_pairs_2022" -b "$BEN" -j 13
  cp "$R/random_pairs_2022/summary.tsv" "$RC/qfo2022_random_pairs_summary.tsv"
  $PYTHON "$BENCH_ROOT/controls/make_qfo_refphylo_oracle.py" "$C" "$R"/random_pairs_2022/raw/SwissTrees/*.gz "$R"/random_pairs_2022/raw/TreeFam-A/*.gz
  for c in oracle_orthologs anti_oracle; do
    bash "$HERE/run_qfo_local.sh" -i "$C/$c.tsv" -n $c -o "$R/${c}_2022" -b "SwissTrees TreeFam-A" -j 2
    cp "$R/${c}_2022/summary.tsv" "$RC/qfo2022_${c}_summary.tsv"
  done
fi
if [[ $what == qfo_full || $what == all ]]; then
  bash "$HERE/run_qfo_local.sh" -y 2018 -i "$DATA_DIR/qfo/benchmark-webservice/example/oma-groups.orthoxml.gz" \
       -n OMA_Groups_example2018 -o "$R/oma_groups_2018" -b "GO EC SwissTrees TreeFam-A STD_Eukaryota STD_Fungi STD_Bacteria G_STD2_Luca G_STD2_Fungi G_STD2_Eukaryota G_STD2_Vertebrata" -j 11
  $PYTHON "$HERE/compare_to_published.py" "$R/oma_groups_2018/summary.tsv" \
       "$BENCH_ROOT/reference_results/openebench_QfO_2018-07-07.tsv" OMA_Groups_2.0 > "$RC/qfo2018_OMA_Groups_example_vs_published.tsv"
  $PYTHON "$HERE/fetch_qfo_participant_predictions.py" --download FastOMA_v0.1.4_nov OrthoFinder --outdir "$DATA_DIR/qfo/participants_2022"
  P="$DATA_DIR/qfo/participants_2022"
  bash "$HERE/run_qfo_local.sh" -i "$P/FastOMA_v0.1.4_qfo2022_17Nov3.orthoxml" -n FastOMA_v0.1.4_nov -o "$R/FastOMA_2022" -b "$ALL" -j 7
  bash "$HERE/run_qfo_local.sh" -i "$P/OrthoFinder.txt" -n OrthoFinder -o "$R/OrthoFinder_2022" -b "$ALL" -j 7
  for t in FastOMA_v0.1.4_nov:FastOMA_2022 OrthoFinder:OrthoFinder_2022; do
    $PYTHON "$HERE/compare_to_published.py" "$R/${t#*:}/summary.tsv" "$BENCH_ROOT/reference_results/openebench_QfO_2022.tsv" "${t%%:*}" \
       > "$RC/qfo2022_${t%%:*}_vs_published.tsv"
  done
fi
