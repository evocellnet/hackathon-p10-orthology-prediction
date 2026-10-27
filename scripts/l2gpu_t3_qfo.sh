#!/bin/bash
# L2-gpu task 3 on QfO: R4 rescue (fast preset + gene order) on the S0_rescue base stages (work/rescue/qfo/base).
# usage: l2gpu_t3_qfo.sh NAME 3DI.fa THREADS
set -u
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
source /opt/conda/etc/profile.d/conda.sh; conda activate $ROOT/conda_envs/orth
NAME=$1; TDI=$2; TH=$3
OUT=$ROOT/work/l2gpu/t3/qfo/$NAME
HC=$ROOT/work/l2gpu/t3/qfo_hits_cache
mkdir -p $OUT/stages $HC
[ -e $HC/profile_a13860c310bf11e6.npz ] || cp $ROOT/work/rescue/qfo/hits_cache/*.npz $HC/
for s in input families scores; do [ -e $OUT/stages/$s ] || cp -al $ROOT/work/rescue/qfo/base/stages/$s $OUT/stages/$s; done
[ -e $OUT/stages/rescue_hits ] || ln -s $HC $OUT/stages/rescue_hits
cd $ROOT/work/l2_gpu_wt
/usr/bin/time -v -o $OUT/time.txt nice -n 5 python -u -m ancestra run --proteomes $ROOT/data/inputs/qfo2022/proteomes \
  --tree $ROOT/data/inputs/qfo2022/guide_tree_ncbi.nwk --gene-order $ROOT/data/inputs/qfo2022/gene_order --out $OUT \
  --threads $TH --preset fast --family-rescue both --structure-3di $TDI --family-rescue-guard > $OUT/run.log 2>&1
echo "exit $?" >> $OUT/run.log
