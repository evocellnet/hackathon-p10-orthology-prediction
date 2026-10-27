#!/bin/bash
# L2-gpu task 2: accurate merge on the frozen OrthoBench S0/S2 stages (hard-linked) with a merge config.
# usage: l2gpu_t2_run.sh NAME CONFIG.json|none THREADS [ob|qfo]
set -u
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
WT=$ROOT/work/l2_gpu_wt
source /opt/conda/etc/profile.d/conda.sh; conda activate $ROOT/conda_envs/orth
export PYTHONPATH=$ROOT/work/l2gpu/pylib
NAME=$1; CFG=$2; TH=$3; DS=${4:-ob}
if [ $DS = ob ]; then
  SRC=$ROOT/work/ob_final/main_accurate_go
  ARGS="--proteomes $ROOT/data/orthobench/Open_Orthobench/BENCHMARKS/Input --tree $ROOT/results/S1_speciestree/ob/ncbi.nwk --gene-order $ROOT/data/genomes/orthobench"
else
  SRC=$ROOT/work/gpuv2/qfo_cpu_fast
  ARGS="--proteomes $ROOT/data/inputs/qfo2022/proteomes --tree $ROOT/data/inputs/qfo2022/guide_tree_ncbi.nwk --gene-order $ROOT/data/inputs/qfo2022/gene_order"
fi
OUT=$ROOT/work/l2gpu/t2/$DS/$NAME
mkdir -p $OUT/stages
for s in input families scores; do [ -e $OUT/stages/$s ] || cp -al $SRC/stages/$s $OUT/stages/$s; done
C=""; [ "$CFG" != none ] && C="--config $CFG"
cd $WT
/usr/bin/time -v -o $OUT/time.txt nice -n 5 python -u -m ancestra run $ARGS --out $OUT --threads $TH --preset accurate $C > $OUT/run.log 2>&1
echo "exit $?" >> $OUT/run.log
