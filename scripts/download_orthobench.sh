#!/usr/bin/env bash
# Open_Orthobench (Emms & Kelly 2020, GBE): 12 input proteomes, 70 RefOGs, benchmark.py
set -euo pipefail
source "$(dirname "$0")/config.sh"
COMMIT=872d6f30592ab5ff837224db16a514b3f2bb916a     # master on 2020-10-14 (pinned)
D="$DATA_DIR/orthobench/Open_Orthobench"
mkdir -p "$DATA_DIR/orthobench"
[[ -d "$D/.git" ]] || git clone https://github.com/davidemms/Open_Orthobench "$D"
git -C "$D" checkout -q "$COMMIT"
B="$D/BENCHMARKS"
echo "proteins per input proteome:"; for f in "$B"/Input/*.fa; do echo "  $(basename $f) $(grep -c '>' $f)"; done
( cd "$B" && md5sum Input/*.fa benchmark.py RefOGs/*.txt RefOGs/low_certainty_assignments/*.txt > "$DATA_DIR/orthobench/md5_inputs.txt" )
# controls (positive = RefOGs themselves, negatives = random partitions)
$PYTHON "$BENCH_ROOT/controls/make_orthobench_controls.py" "$B" "$DATA_DIR/orthobench/controls"
