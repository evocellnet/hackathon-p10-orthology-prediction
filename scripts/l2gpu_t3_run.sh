#!/bin/bash
# L2-gpu task 3: R4 rescue (3Di U profiles + RBH merge guard) on the frozen OrthoBench S0/S2 stages.
# usage: l2gpu_t3_run.sh NAME 3DI.fa THREADS
set -u
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
WT=$ROOT/work/l2_gpu_wt
source /opt/conda/etc/profile.d/conda.sh; conda activate $ROOT/conda_envs/orth
NAME=$1; TDI=$2; TH=$3
SRC=$ROOT/work/ob_final/main_accurate_go
OUT=$ROOT/work/l2gpu/t3/ob/$NAME
mkdir -p $OUT/stages $ROOT/work/l2gpu/t3/hits_cache
for s in input families scores; do [ -e $OUT/stages/$s ] || cp -al $SRC/stages/$s $OUT/stages/$s; done
[ -e $OUT/stages/rescue_hits ] || ln -s $ROOT/work/l2gpu/t3/hits_cache $OUT/stages/rescue_hits
cd $WT
/usr/bin/time -v -o $OUT/time.txt nice -n 5 python -u -m ancestra run --proteomes $ROOT/data/orthobench/Open_Orthobench/BENCHMARKS/Input \
  --tree $ROOT/results/S1_speciestree/ob/ncbi.nwk --gene-order $ROOT/data/genomes/orthobench --out $OUT --threads $TH \
  --preset accurate --family-rescue both --structure-3di $TDI --family-rescue-guard > $OUT/run.log 2>&1
echo "exit $?" >> $OUT/run.log
