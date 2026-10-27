#!/usr/bin/env bash
# QfO FAS benchmark for the final runs (not wired in bench_repo/scripts/run_qfo_local.sh).
# Exactly the main.nf process `fas_benchmark` (upstream fas_benchmark.py, unmodified) with greedyFAS 1.18.14
# (pip wheel, separate venv conda_envs/fas_venv on top of conda_envs/data) providing fas.runMultiTaxa for the
# <= 9,000 pairs missing from fas_precomputed.json.gz. NB: upstream samples precomputed pairs with an unseeded
# random.shuffle -> FAS is stochastic (like STD); --cpus ${FAS_CPUS:-1} (shared-node budget).
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
WS=$ROOT/data/qfo/benchmark-webservice; REF=$ROOT/data/qfo/reference_data/2022.1
export PATH="$ROOT/conda_envs/fas_venv/bin:$PATH"
for arm in "$@"; do
  O=$ROOT/work/qfo_final/$arm
  [[ -s "$O/orthologs.db" ]] || { echo "[fas] $arm: no orthologs.db"; continue; }
  mkdir -p "$O/raw/FAS" "$O/assess"
  ( cd "$WS" && /usr/bin/time -v "$ROOT/conda_envs/fas_venv/bin/python" fas_benchmark.py --com QfO \
      --participant "Ancestra_$arm" --assessment-out "$O/assess/FAS.json" --outdir "$O/raw/FAS" \
      --fas-precomputed-scores "$REF/fas_precomputed.json.gz" --fas-data "$REF/fas_annotations/" \
      --db "$O/orthologs.db" --cpus "${FAS_CPUS:-1}" ) > "$O/logs/FAS.log" 2>&1
  echo "[fas] $arm exit $? $(grep -E 'FAS_mean|Elapsed' "$O/logs/FAS.log" | tr '\n' ' ' | cut -c1-200)"
done
