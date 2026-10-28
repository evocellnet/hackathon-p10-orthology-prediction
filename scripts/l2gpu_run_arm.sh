#!/bin/bash
# L2-gpu: merge-only arm on an injected score table.
# usage: l2gpu_run_arm.sh ob|qfo NAME SCORES.npz CONFIG.json THREADS PRESET SRC_RUN
#   SRC_RUN: run dir providing the input + S0 families stages (hard-linked / copied)
set -u
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
WT=$ROOT/work/l2_gpu_wt
source /opt/conda/etc/profile.d/conda.sh; conda activate $ROOT/conda_envs/orth
DS=$1; NAME=$2; SC=$3; CFG=$4; TH=$5; PR=$6; SRC=$7
OUT=$ROOT/work/l2gpu/$DS/$NAME
mkdir -p $OUT
cd $WT
python scripts/gpu_inject_stages.py --src $SRC --out $OUT --scores $SC \
  --families-final $SRC/stages/families/families.npy --backend gpusw --config $CFG > $OUT/inject.log 2>&1 || exit 1
if [ $DS = ob ]; then
  ARGS="--proteomes $ROOT/data/orthobench/Open_Orthobench/BENCHMARKS/Input --tree $ROOT/results/S1_speciestree/ob/ncbi.nwk --gene-order $ROOT/data/genomes/orthobench"
else
  ARGS="--proteomes $ROOT/data/inputs/qfo2022/proteomes --tree $ROOT/data/inputs/qfo2022/guide_tree_ncbi.nwk --gene-order $ROOT/data/inputs/qfo2022/gene_order"
fi
/usr/bin/time -v -o $OUT/time.txt nice -n 5 python -u -m ancestra run $ARGS --out $OUT --threads $TH --preset $PR \
  --backend gpusw --config $CFG > $OUT/run.log 2>&1
echo "exit $?" >> $OUT/run.log
