#!/bin/bash
# S0 rescue grid on OrthoBench (DEV eval only), accurate preset + gene order, cached run_det S0/S2.
# usage: bash scripts/rescue_ob_arms.sh THREADS ARM...   (ARM in R0 R1 R2 R3 R4)
set -e
cd "$(dirname "$0")/.."
source /opt/conda/etc/profile.d/conda.sh && conda activate ./conda_envs/orth
T=$1; shift
for arm in "$@"; do
  case $arm in
    R0) opt="--family-rescue none";;
    R1) opt="--family-rescue 3di --structure-3di data/structure/orthobench_3di_best.fa";;
    R2) opt="--family-rescue profile";;
    R3) opt="--family-rescue both --structure-3di data/structure/orthobench_3di_best.fa";;
    R4) opt="--family-rescue both --structure-3di data/structure/orthobench_3di_best.fa --family-rescue-guard";;
  esac
  d=work/rescue/ob/$arm
  mkdir -p $d/stages
  for s in input families scores; do [ -e $d/stages/$s ] || cp -al work/core_v0_ob/run_det/stages/$s $d/stages/$s; done
  [ -e $d/stages/rescue_hits ] || ln -s ../../hits_cache $d/stages/rescue_hits
  echo "[$(date +%T)] $arm start" >> work/rescue/ob/chain.log
  /usr/bin/time -v python -m ancestra run --proteomes data/orthobench/Open_Orthobench/BENCHMARKS/Input \
     --tree results/S1_speciestree/ob/ncbi.nwk --gene-order data/genomes/orthobench --out $d \
     --threads $T --preset accurate $opt > $d/run.log 2>&1
  python scripts/rescue_eval_ob.py $d --label $arm >> work/rescue/ob/chain.log 2>&1
  echo "[$(date +%T)] $arm done" >> work/rescue/ob/chain.log
done
