#!/bin/bash
# L2-domains: run merge arms sequentially. scripts/L2d_queue.sh BENCH THREADS "arm:shuf arm:shuf ..." DOMTBL...
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
B=$1; T=$2; LIST=$3; shift 3
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=$T PATH=$ROOT/conda_envs/orth/bin:$PATH
cd $ROOT
for x in $LIST; do
  arm=${x%%:*}; sh=${x##*:}
  out=results/l2_domains/arms/${B}_${arm}_${sh}.json
  [ -s $out ] && continue
  nice -n 5 python work/l2_domains_wt/scripts/L2d_merge_arm.py --bench $B --arm $arm --shuffle $sh \
     --domtbl "$@" --out $out --threads $T > work/l2_domains/logs/${B}_${arm}_${sh}.log 2>&1
done
echo queue-done
