#!/usr/bin/env bash
# Final QfO layer / guide-tree ablations (PI 2026-10-27, freeze-v1; no code changes, no tuning).
# Each arm: merge-only rerun on the cached deterministic S0/S2 (work/qfo2022_det/det_nogo) with the frozen
# preset, verify_orthoxml (must PASS), test-log line (peek=final) BEFORE scoring, then the fast QfO tests.
#
#   bash scripts/qfo_final_ablations.sh LANE_THREADS ARM...
# ARMs: a_nogo  b_rescue  c_astral  d_fastme  e_nni5  f_random  g_accurate_nogo
# Baseline = cand2_fast (fast preset, gene order on, NCBI tree): work/qfo2022_det/cand2_fast.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
source /opt/conda/etc/profile.d/conda.sh; conda activate "$ROOT/conda_envs/orth"
TH="$1"; shift
P=data/inputs/qfo2022/proteomes; NCBI=data/inputs/qfo2022/guide_tree_ncbi.nwk; GO=data/inputs/qfo2022/gene_order
S1=results/S1_speciestree/qfo; SRC=work/qfo2022_det/det_nogo
W=work/qfo2022_det/ablation; Q=work/qfo_final/ablation; L=$W/logs
mkdir -p "$W" "$Q" "$L"
FAST="GO EC VGNC SwissTrees TreeFam-A"

for arm in "$@"; do
  tree=$NCBI; preset=fast; go=$GO; extra=(--reuse-stages-from "$SRC")
  case "$arm" in
    a_nogo)          go=none;;
    b_rescue)        extra+=(--family-rescue both --family-rescue-guard --structure-3di data/structure/qfo_3di_best.fa);;
    c_astral)        tree=$S1/armB_astralpro_ncbiroot.nwk; extra=(--transplant-from "$SRC");;
    d_fastme)        tree=$S1/armA_fastme_ncbiroot.nwk;    extra=(--transplant-from "$SRC");;
    e_nni5)          tree=$S1/perturbed_ncbi_nni5_s11.nwk; extra=(--transplant-from "$SRC");;
    f_random)        tree=$S1/random_s11.nwk;              extra=(--transplant-from "$SRC");;
    g_accurate_nogo) go=none; preset=accurate;;
    *) echo "unknown arm $arm"; continue;;
  esac
  out=$W/$arm
  echo "[abl] $arm start $(date -Iseconds) load $(cut -d' ' -f1-3 /proc/loadavg) threads=$TH tree=$tree preset=$preset go=$go"
  /usr/bin/time -v python scripts/qfo_pipeline_driver.py --proteomes $P --tree "$tree" --out "$out" --threads "$TH" \
      --gene-order "$go" --preset "$preset" "${extra[@]}" > "$L/$arm.time.log" 2>&1 \
    || { echo "[abl] $arm PIPELINE FAILED"; tail -5 "$L/$arm.time.log"; continue; }
  grep -E "\[driver\]|stage (rescue|merge|output): done|Elapsed|Maximum res" "$L/$arm.time.log"
  if ! python scripts/verify_orthoxml.py "$out" --json "results/qfo_final/ablation_${arm}_verify_orthoxml.json" > "$L/$arm.verify.log" 2>&1; then
    echo "[abl] $arm verify_orthoxml FAIL - not scored"; grep FAIL "$L/$arm.verify.log" | cut -c1-200; continue
  fi
  echo "[abl] $arm verify_orthoxml PASS"
  H=$(python -c "
import json,hashlib; m=json.load(open('$out/manifest.json')); d=m['driver']
print(hashlib.sha256(json.dumps({'config':m['config'],'families_config_hash':d['families_config_hash'],'code_keys':d['stage_code_keys']},sort_keys=True,default=str).encode()).hexdigest()[:12])")
  printf "%s\tQfO peek=final\tablation %s (preset %s, gene order %s, tree %s)\t%s\tQfO2022 fast: GO EC VGNC SwissTrees TreeFam-A\tfinal ablation vs cand2_fast, freeze-v1; git=%s; out=%s/hogs.orthoxml\n" \
    "$(date -Iseconds)" "$arm" "$preset" "$([ "$go" = none ] && echo off || echo on)" "$(basename "$tree")" "$H" \
    "$(git rev-parse --short HEAD)" "$out" >> results/test_log.tsv
  bash bench_repo/scripts/run_qfo_local.sh -i "$out/hogs.orthoxml" -n "Ancestra_abl_$arm" -o "$Q/$arm" -j 1 -b "$FAST" \
      > "$Q/$arm.log" 2>&1
  grep -E "\[qfo\] (total|.*FAILED)" "$Q/$arm.log"
  echo "[abl] $arm done $(date -Iseconds)"
done
