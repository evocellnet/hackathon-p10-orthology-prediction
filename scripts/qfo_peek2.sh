#!/usr/bin/env bash
# QfO PEEK 2 (PI-authorised 2026-10-27): frozen candidate cand2_{accurate,fast}.
# Gate: verify_orthoxml must PASS (0 structural violations, pyham loads) for both arms; then the test log is
# appended BEFORE any scoring; then the fast QfO tests; then the table and the TreeFam-A dev/test subsets.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
source /opt/conda/etc/profile.d/conda.sh; conda activate "$ROOT/conda_envs/orth"
ARMS="cand2_accurate cand2_fast"
for a in $ARMS; do
  python - "$a" <<'PY'
import json, sys
a = sys.argv[1]
v = json.load(open(f"results/qfo_run/{a}_verify_orthoxml.json"))
t = v["checks"]["taxrange_lca"]
ok = v["ok"] and t["n_nesting_violations"] == 0 and t["n_paralog_level_violations"] == 0 and v["checks"]["pyham"]["ok"]
print(f"[peek2] gate {a}: verify ok={v['ok']} nesting={t['n_nesting_violations']} paralog={t['n_paralog_level_violations']} pyham={v['checks']['pyham']['ok']}")
sys.exit(0 if ok else 3)
PY
done
python scripts/verify_ids.py --tree data/inputs/qfo2022/guide_tree_ncbi.nwk --gene-order data/inputs/qfo2022/gene_order \
  --json results/qfo_run/verify_ids_peek2.json | tail -1 | grep -q "overall: PASS"
TS=$(date -Iseconds); GIT=$(git rev-parse --short HEAD)
for a in $ARMS; do
  H=$(python -c "
import json,hashlib; m=json.load(open('work/qfo2022_det/$a/manifest.json')); d=m['driver']
print(hashlib.sha256(json.dumps({'config':m['config'],'families_config_hash':d['families_config_hash'],'code_keys':d['stage_code_keys']},sort_keys=True,default=str).encode()).hexdigest()[:12])")
  G=$(python -c "import json; print(json.load(open('work/qfo2022_det/$a/manifest.json'))['driver']['git']['commit'][:7])")
  printf "%s\tQfO peek=2\t%s (preset %s, gene order)\t%s\tQfO2022 fast: GO EC VGNC SwissTrees TreeFam-A SP\tfrozen candidate accurate/fast after FastTree + species-overlap fixes; git=%s; out=work/qfo2022_det/%s/hogs.orthoxml\n" \
    "$TS" "$a" "${a#cand2_}" "$H" "$G" "$a" >> results/test_log.tsv
done
tail -2 results/test_log.tsv
mkdir -p work/qfo_peek2
for a in $ARMS; do
  /usr/bin/time -v bash bench_repo/scripts/run_qfo_local.sh -i "work/qfo2022_det/$a/hogs.orthoxml" -n "Ancestra_$a" \
     -o "work/qfo_peek2/$a" -j 6 -b "GO EC VGNC SwissTrees TreeFam-A SP" > "work/qfo_peek2/$a.log" 2>&1 &
done
wait
grep -hE "\[qfo\] (total|.*FAILED)|Elapsed" work/qfo_peek2/*.log
python scripts/make_qfo_peek1_table.py --out results/qfo_peek2 \
  --ours "Ancestra cand2_accurate (peek 2)=work/qfo_peek2/cand2_accurate=local (peek 2)" \
  --ours "Ancestra cand2_fast (peek 2)=work/qfo_peek2/cand2_fast=local (peek 2)" \
  --ours "Ancestra A det_go (peek 1)=work/qfo_peek1/det_go=local (peek 1)" \
  --ours "Ancestra B det_nogo (peek 1)=work/qfo_peek1/det_nogo=local (peek 1)"
python scripts/qfo_treefam_subset.py work/qfo_peek2/cand2_accurate work/qfo_peek2/cand2_fast work/qfo_peek1/det_go \
  --subset dev test all --allow-test --json results/qfo_peek2/treefamA_dev_test.json
echo "[peek2] done"
