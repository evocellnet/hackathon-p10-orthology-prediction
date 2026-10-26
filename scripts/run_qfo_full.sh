#!/usr/bin/env bash
# End-to-end Ancestra run on QfO 2022 (78 species), resumable. Does NOT run the QfO benchmark
# (peeks are PI-authorised only and logged in results/test_log.tsv; see docs/qfo_run.md).
#
#   bash scripts/run_qfo_full.sh [full|sub10] [--stop-after input|families|scores|merge|output]
#        [--threads N] [--config cfg.json] [--name RUN_NAME]
#
# Steps
#   1. scripts/prepare_qfo_inputs.py   symlinks of the 78 canonical FASTAs + NCBI guide tree
#                                      (data/inputs/qfo2022/, md5 vs release md5sums.txt)
#   2. scripts/verify_ids.py           input check (DESIGN §4); aborts the run on failure
#   3. scripts/qfo_pipeline_driver.py  ancestra.pipeline (input, families=S0, scores=S2, merge=S3,
#                                      output) under /usr/bin/time -v; per-stage wall/CPU/RSS/load in
#                                      work/qfo2022/<name>/stage_accounting.json. Stages are cached in
#                                      work/qfo2022/<name>/stages/ (hash-keyed DONE files): rerunning
#                                      this script resumes after the last finished stage.
#   4. scripts/verify_orthoxml.py      DESIGN §6 output checks -> results/qfo_run/<name>_verify_orthoxml.json
#   5. scripts/summarize_qfo_run.py    families / HOGs per level / runtime -> results/qfo_run/<name>_summary.{json,md}
# Steps 4-5 run only once the output stage exists.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

MODE="full"; STOP=""; THREADS="${THREADS:-12}"; CONFIG=""; NAME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    full|sub10) MODE="$1"; shift ;;
    --stop-after) STOP="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done

# shellcheck disable=SC1091
source /opt/conda/etc/profile.d/conda.sh
conda activate "$ROOT/conda_envs/orth"

IN="data/inputs/qfo2022"
TREE="$IN/guide_tree_ncbi.nwk"
if [[ "$MODE" == "full" ]]; then PROT="$IN/proteomes"; NAME="${NAME:-full}";
else PROT="$IN/subsets/sub10_s11"; NAME="${NAME:-sub10_s11}"; fi
OUT="work/qfo2022/$NAME"
LOGD="work/qfo2022/logs"; RES="results/qfo_run"
mkdir -p "$OUT" "$LOGD" "$RES"
TS="$(date +%Y%m%dT%H%M%S)"
LOG="$LOGD/${NAME}.${TS}.log"

echo "[run_qfo] $(bash scripts/clock.sh)" | tee -a "$LOG"
echo "[run_qfo] mode=$MODE name=$NAME threads=$THREADS stop_after=${STOP:-none} config=${CONFIG:-defaults}" | tee -a "$LOG"
echo "[run_qfo] load: $(cat /proc/loadavg)" | tee -a "$LOG"

# 1-2. inputs + verification (always on the full 78-species set: DESIGN §4)
python scripts/prepare_qfo_inputs.py 2>&1 | tee -a "$LOG"
if ! python scripts/verify_ids.py --tree "$TREE" --json "$RES/verify_ids.json" 2>&1 | tee -a "$LOG"; then
  echo "[run_qfo] verify_ids FAILED - not running the pipeline" | tee -a "$LOG"; exit 1
fi
grep -q "overall: PASS" "$LOG" || { echo "[run_qfo] verify_ids did not pass" | tee -a "$LOG"; exit 1; }

# 3. pipeline (resumable)
ARGS=(--proteomes "$PROT" --tree "$TREE" --out "$OUT" --threads "$THREADS")
[[ -n "$STOP" ]] && ARGS+=(--stop-after "$STOP")
[[ -n "$CONFIG" ]] && ARGS+=(--config "$CONFIG")
/usr/bin/time -v python scripts/qfo_pipeline_driver.py "${ARGS[@]}" >> "$LOG" 2>&1
echo "[run_qfo] pipeline step finished: $(grep -E 'Elapsed|Maximum resident' "$LOG" | tail -2 | tr -s ' ' | tr '\n' ';')" | tee -a "$LOG"
echo "[run_qfo] load: $(cat /proc/loadavg)" | tee -a "$LOG"

# 4-5. verification + summary once the output exists
if [[ -f "$OUT/stages/output/DONE" && -f "$OUT/hogs.orthoxml" ]]; then
  /usr/bin/time -v python scripts/verify_orthoxml.py "$OUT" --json "$RES/${NAME}_verify_orthoxml.json" >> "$LOG" 2>&1 \
    && echo "[run_qfo] verify_orthoxml PASS" | tee -a "$LOG" \
    || echo "[run_qfo] verify_orthoxml FAIL (see $RES/${NAME}_verify_orthoxml.json)" | tee -a "$LOG"
  python scripts/summarize_qfo_run.py "$OUT" --json "$RES/${NAME}_summary.json" --md "$RES/${NAME}_summary.md" \
    --timelog "$LOG" 2>&1 | tee -a "$LOG"
fi
echo "[run_qfo] done $(bash scripts/clock.sh); log $LOG" | tee -a "$LOG"
