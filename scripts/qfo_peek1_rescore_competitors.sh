#!/usr/bin/env bash
# Re-score published QfO 2022 participant prediction files (data/qfo/participants_2022, B2SHARE, CC BY 4.0)
# locally on the fast pair-level tests that report raw TP/FP/FN (SwissTrees, TreeFam-A, VGNC, SP), so that
# the number of *evaluated* predicted pairs (TP+FP) is available next to ours. Also a reproduction check
# of the published TPR/PPV. 3 tools at a time x 4 tests = <=12 processes.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
P=data/qfo/participants_2022; OUT=work/qfo_peek1/competitors; mkdir -p "$OUT"
declare -A F=(
  [OrthoFinder]=OrthoFinder.txt [FastOMA_v0.1.4_nov]=FastOMA_v0.1.4_qfo2022_17Nov3.orthoxml
  [Hieranoid_2]=pairs_hieranoid-diamond.tsv [InParanoid5]=InParanoid_QFO22.pairs_1.gz
  [SonicParanoid2-sens-g]=qfo22-challenge.sp2-sens-g.tsv.gz [OMA_HOGs]=HierarchicalGroups.orthoxml.gz
  [OMA_Pairs]=OMA.2.5.0-VPairs.txt.gz [PANTHER_v18_all]=PANTHER18_all_1
  [pbDOOD_eggnog_v7]=pbdood_pairs.tsv [pbDOOD_midpoint]=pbdood_pairs_midpoint.tsv)
for t in "${!F[@]}"; do
  [[ -s "$OUT/$t/summary.tsv" ]] && { echo "[rescore] $t cached"; continue; }
  ( /usr/bin/time -v bash bench_repo/scripts/run_qfo_local.sh -i "$P/${F[$t]}" -n "$t" -o "$OUT/$t" -j 4 \
      -b "SwissTrees TreeFam-A VGNC SP" > "$OUT/$t.log" 2>&1; echo "[rescore] $t exit $?" ) &
  while (( $(jobs -rp | wc -l) >= 3 )); do sleep 5; done
done
wait
echo "[rescore] all done"
