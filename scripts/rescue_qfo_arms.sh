#!/bin/bash
# S0 rescue grid on QfO (TreeFam-A DEV only, scripts/qfo_treefam_dev_fast.py; never --allow-test),
# fast preset + gene order, on hard-linked cached det S0/S2 stages (work/rescue/qfo/base).
# usage: bash scripts/rescue_qfo_arms.sh THREADS ARM...
set -e
cd "$(dirname "$0")/.."
source /opt/conda/etc/profile.d/conda.sh && conda activate ./conda_envs/orth
T=$1; shift
mkdir -p work/rescue/qfo/hits_cache
for arm in "$@"; do
  case $arm in
    R0) opt="--family-rescue none";;
    R1) opt="--family-rescue 3di --structure-3di data/structure/qfo_3di_best.fa";;
    R2) opt="--family-rescue profile";;
    R3) opt="--family-rescue both --structure-3di data/structure/qfo_3di_best.fa";;
    R4) opt="--family-rescue both --structure-3di data/structure/qfo_3di_best.fa --family-rescue-guard";;
  esac
  d=work/rescue/qfo/$arm
  mkdir -p $d/stages
  for s in input families scores; do [ -e $d/stages/$s ] || cp -al work/rescue/qfo/base/stages/$s $d/stages/$s; done
  [ -e $d/stages/rescue_hits ] || ln -s ../../hits_cache $d/stages/rescue_hits
  echo "[$(date +%T)] $arm start" >> work/rescue/qfo/chain.log
  /usr/bin/time -v python -m ancestra run --proteomes data/inputs/qfo2022/proteomes \
     --tree data/inputs/qfo2022/guide_tree_ncbi.nwk --gene-order data/inputs/qfo2022/gene_order --out $d \
     --threads $T --preset fast $opt > $d/run.log 2>&1
  rm -f $d/pairs.tsv $d/stages/output/pairs.tsv
  python scripts/rescue_eval_qfo.py $d --label $arm >> work/rescue/qfo/chain.log 2>&1
  echo "[$(date +%T)] $arm done" >> work/rescue/qfo/chain.log
done
