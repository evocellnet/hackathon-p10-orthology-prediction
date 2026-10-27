#!/bin/bash
# Re-run merges of rescue arms on the current merge.py (nesting fix of 02:08 landed mid-grid), reusing the
# arm's rescue stage (rescue does not depend on merge.py). usage: BENCH THREADS ARM... (ARM in R2 R4 R1 R3)
set -e
cd "$(dirname "$0")/.."
source /opt/conda/etc/profile.d/conda.sh && conda activate ./conda_envs/orth
B=$1; T=$2; shift 2
if [ $B = qfo ]; then
  base=work/rescue/qfo/base/stages; ss=data/structure/qfo_3di_best.fa
  common="--proteomes data/inputs/qfo2022/proteomes --tree data/inputs/qfo2022/guide_tree_ncbi.nwk --gene-order data/inputs/qfo2022/gene_order --preset fast"
else
  base=work/core_v0_ob/run_det/stages; ss=data/structure/orthobench_3di_best.fa
  common="--proteomes data/orthobench/Open_Orthobench/BENCHMARKS/Input --tree results/S1_speciestree/ob/ncbi.nwk --gene-order data/genomes/orthobench --preset accurate"
fi
for arm in "$@"; do
  case $arm in
    R1) opt="--family-rescue 3di --structure-3di $ss";;
    R2) opt="--family-rescue profile";;
    R3) opt="--family-rescue both --structure-3di $ss";;
    R4) opt="--family-rescue both --structure-3di $ss --family-rescue-guard";;
  esac
  d=work/rescue/$B/${arm}n
  mkdir -p $d/stages
  for s in input families scores; do [ -e $d/stages/$s ] || cp -al $base/$s $d/stages/$s; done
  [ -e $d/stages/rescue ] || cp -al work/rescue/$B/$arm/stages/rescue $d/stages/rescue
  [ -e $d/stages/rescue_hits ] || ln -s ../../hits_cache $d/stages/rescue_hits
  echo "[$(date +%T)] ${arm}n start" >> work/rescue/$B/chain.log
  /usr/bin/time -v python -m ancestra run $common --out $d --threads $T $opt > $d/run.log 2>&1
  rm -f $d/pairs.tsv $d/stages/output/pairs.tsv
  python scripts/rescue_eval_$B.py $d --label ${arm}n >> work/rescue/$B/chain.log 2>&1
  echo "[$(date +%T)] ${arm}n done" >> work/rescue/$B/chain.log
done
