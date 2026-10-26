#!/usr/bin/env bash
# End-to-end Ancestra run on QfO 2022 (78 species), resumable. Does NOT run the QfO benchmark
# (peeks are PI-authorised only and logged in results/test_log.tsv; see docs/qfo_run.md).
#
#   bash scripts/run_qfo_full.sh [full|sub10] [--gene-order both|on|off] [--threads N]
#        [--stop-after input|families|scores|merge|output] [--config cfg.json] [--name PREFIX]
#        [--adopt-unrecorded] [--workdir work/qfo2022]
#
# Run dirs: <WORKDIR>/<PREFIX>_nogo (merge without gene order) and <PREFIX>_go (merge with
# the S5 gene-order tables, constrained synteny = current default merge). PREFIX defaults to
# "full" / "sub10_s11". With --gene-order both (default) the nogo run computes S0+S2 and the go
# run hard-links them (--reuse-stages-from), so S0 and scoring run once.
#
# Steps
#   1. scripts/prepare_qfo_inputs.py   symlinks of the 78 canonical FASTAs, the S5 gene-order
#                                      tables and the NCBI guide tree (data/inputs/qfo2022/)
#   2. scripts/verify_ids.py           input check (DESIGN §4, + gene-order tables); aborts on failure
#   3. scripts/qfo_pipeline_driver.py  ancestra.pipeline under /usr/bin/time -v. Stage cache in
#                                      <run>/stages/ (pipeline keys) + code keys (driver: sha256 of
#                                      the ancestra/ files each stage runs + effective families
#                                      config): a cached stage is recomputed when its code changes.
#                                      Rerunning this script resumes after the last finished stage.
#                                      Git commit, code hashes and families-config hash -> manifest.json.
#   4. scripts/verify_orthoxml.py      DESIGN §6 output checks -> results/qfo_run/<run>_verify_orthoxml.json
#   5. scripts/summarize_qfo_run.py    families / HOGs per level / runtime -> results/qfo_run/<run>_summary.{json,md}
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

MODE="full"; STOP=""; THREADS="${THREADS:-12}"; CONFIG=""; PREFIX=""; GO="both"; ADOPT=""; WD="work/qfo2022"
while [[ $# -gt 0 ]]; do
  case "$1" in
    full|sub10) MODE="$1"; shift ;;
    --gene-order) GO="$2"; shift 2 ;;
    --stop-after) STOP="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --name) PREFIX="$2"; shift 2 ;;
    --adopt-unrecorded) ADOPT="--adopt-unrecorded"; shift ;;
    --workdir) WD="$2"; shift 2 ;;
    -h|--help) sed -n '2,27p' "$0"; exit 0 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
case "$GO" in both|on|off) ;; *) echo "--gene-order must be both|on|off" >&2; exit 2 ;; esac

# shellcheck disable=SC1091
source /opt/conda/etc/profile.d/conda.sh
conda activate "$ROOT/conda_envs/orth"

IN="data/inputs/qfo2022"
TREE="$IN/guide_tree_ncbi.nwk"
GODIR="$IN/gene_order"
if [[ "$MODE" == "full" ]]; then PROT="$IN/proteomes"; PREFIX="${PREFIX:-full}";
else PROT="$IN/subsets/sub10_s11"; PREFIX="${PREFIX:-sub10_s11}"; fi
LOGD="$WD/logs"; RES="results/qfo_run"
mkdir -p "$LOGD" "$RES"
TS="$(date +%Y%m%dT%H%M%S)"
LOG="$LOGD/${PREFIX}.${TS}.log"

echo "[run_qfo] $(bash scripts/clock.sh)" | tee -a "$LOG"
echo "[run_qfo] workdir=$WD mode=$MODE prefix=$PREFIX gene_order=$GO threads=$THREADS stop_after=${STOP:-none} config=${CONFIG:-defaults} $ADOPT" | tee -a "$LOG"
echo "[run_qfo] git $(git rev-parse HEAD 2>/dev/null || echo none); load: $(cat /proc/loadavg)" | tee -a "$LOG"

# 1-2. inputs + verification (always on the full 78-species set: DESIGN §4)
python scripts/prepare_qfo_inputs.py 2>&1 | tee -a "$LOG"
if ! python scripts/verify_ids.py --tree "$TREE" --gene-order "$GODIR" --json "$RES/verify_ids.json" 2>&1 | tee -a "$LOG"; then
  echo "[run_qfo] verify_ids FAILED - not running the pipeline" | tee -a "$LOG"; exit 1
fi

run_one() {   # $1 = go|nogo, $2 = run dir to reuse S0+S2 from (optional)
  local v="$1" reuse="${2:-}" name="${PREFIX}_$1" out="$WD/${PREFIX}_$1"
  local args=(--proteomes "$PROT" --tree "$TREE" --out "$out" --threads "$THREADS")
  [[ "$v" == "go" ]] && args+=(--gene-order "$GODIR") || args+=(--gene-order none)
  [[ -n "$reuse" ]] && args+=(--reuse-stages-from "$reuse")
  [[ -n "$STOP" ]] && args+=(--stop-after "$STOP")
  [[ -n "$CONFIG" ]] && args+=(--config "$CONFIG")
  [[ -n "$ADOPT" ]] && args+=("$ADOPT")
  echo "[run_qfo] === $name: $(date -Iseconds) load $(cut -d' ' -f1-3 /proc/loadavg)" | tee -a "$LOG"
  local tlog="$LOGD/${name}.${TS}.time.log"
  /usr/bin/time -v python scripts/qfo_pipeline_driver.py "${args[@]}" > "$tlog" 2>&1 \
    || { echo "[run_qfo] $name pipeline FAILED (see $tlog)" | tee -a "$LOG"; tail -20 "$tlog" >> "$LOG"; return 1; }
  grep -E "\[driver\]|\[ancestra\]" "$tlog" >> "$LOG" || true
  echo "[run_qfo] $name pipeline: $(grep -E 'Elapsed|Maximum resident' "$tlog" | tr -s ' ' | tr '\n' ';')" | tee -a "$LOG"
  if [[ -f "$out/stages/output/DONE" && -f "$out/hogs.orthoxml" ]]; then
    if python scripts/verify_orthoxml.py "$out" --json "$RES/${name}_verify_orthoxml.json" >> "$LOG" 2>&1; then
      echo "[run_qfo] $name verify_orthoxml PASS" | tee -a "$LOG"
    else
      echo "[run_qfo] $name verify_orthoxml FAIL (see $RES/${name}_verify_orthoxml.json)" | tee -a "$LOG"
    fi
    python scripts/summarize_qfo_run.py "$out" --json "$RES/${name}_summary.json" --md "$RES/${name}_summary.md" \
      --timelog "$tlog" 2>&1 | tee -a "$LOG"
  fi
}

case "$GO" in
  off)  run_one nogo ;;
  on)   run_one go ;;
  both) run_one nogo && run_one go "$WD/${PREFIX}_nogo" ;;
esac
echo "[run_qfo] done $(bash scripts/clock.sh); log $LOG" | tee -a "$LOG"
