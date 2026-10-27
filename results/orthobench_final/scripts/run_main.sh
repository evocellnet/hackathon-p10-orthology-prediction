#!/bin/bash
# OrthoBench FINAL evaluation: fresh end-to-end Ancestra run, preset accurate + gene order (frozen freeze-v1).
# usage: bash results/orthobench_final/scripts/run_main.sh   (from repo root)
set -e
cd /scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
source /opt/conda/etc/profile.d/conda.sh && conda activate ./conda_envs/orth
T=10
d=work/ob_final/main_accurate_go
[ -e $d ] && { echo "$d exists; fresh run requires a new dir"; exit 1; }
mkdir -p $d
S=results/orthobench_final/scripts
( /usr/bin/time -v -o $d/time.txt python -u -m ancestra run \
    --proteomes data/orthobench/Open_Orthobench/BENCHMARKS/Input \
    --tree results/S1_speciestree/ob/ncbi.nwk --gene-order data/genomes/orthobench \
    --out $d --threads $T --preset accurate 2>&1 \
  | while IFS= read -r line; do printf '%s\t%s\n' "$(date +%s.%N)" "$line"; done > $d/run.log ) &
sleep 1
P=$(pgrep -f "time -v -o $d/time.txt" | head -1)
python $S/proc_sampler.py $P $d/proc_samples.tsv 2 &
wait
echo "exit done $(date -Iseconds)" >> $d/run.log
