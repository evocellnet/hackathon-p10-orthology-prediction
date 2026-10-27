#!/bin/bash
# OrthoBench FINAL: merge-only variants on the cache of the fresh main run (work/ob_final/main_accurate_go).
# Same-tree variants: stages input/families/scores hard-linked, unmodified CLI (as scripts/rescue_ob_arms.sh).
# Guide-tree arms: run_cached.py (same CLI; tree-independent S0/S2 checkpoints adopted under the new key).
# usage: bash results/orthobench_final/scripts/run_variants.sh [LABEL...]
set -e
cd /scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
source /opt/conda/etc/profile.d/conda.sh && conda activate ./conda_envs/orth
T=10
B=work/ob_final/main_accurate_go
S=results/orthobench_final/scripts
TR=results/S1_speciestree/ob
IN=data/orthobench/Open_Orthobench/BENCHMARKS/Input
GO="--gene-order data/genomes/orthobench"
LABELS=${@:-fast_go accurate_nogo fast_nogo accurate_go_rescue tree_astralpro tree_fastme tree_nni1_s11 tree_nni1_s23 tree_nni1_s47 tree_random_s11}
for lab in $LABELS; do
  tree=$TR/ncbi.nwk; preset=accurate; go="$GO"; extra=""; cached=0
  case $lab in
    fast_go) preset=fast;;
    accurate_nogo) go="";;
    fast_nogo) preset=fast; go="";;
    accurate_go_rescue) extra="--family-rescue both --structure-3di data/structure/orthobench_3di_best.fa --family-rescue-guard";;
    tree_astralpro) tree=$TR/armB_astralpro_ncbiroot.nwk; cached=1;;
    tree_fastme) tree=$TR/armA_fastme_ncbiroot.nwk; cached=1;;
    tree_nni1_s11|tree_nni1_s23|tree_nni1_s47) tree=$TR/perturbed_ncbi_nni1_${lab##*_}.nwk; cached=1;;
    tree_random_s11) tree=$TR/random_s11.nwk; cached=1;;
    *) echo "unknown $lab"; exit 1;;
  esac
  d=work/ob_final/$lab
  if [ -e $d/orthobench_ogs.txt ]; then echo "[$(date +%T)] $lab exists, skip" >> work/ob_final/variants.log; continue; fi
  mkdir -p $d/stages
  args="--proteomes $IN --tree $tree $go --out $d --threads $T --preset $preset $extra"
  if [ $cached = 1 ]; then
    cmd="python -u $S/run_cached.py $B -- $args"
  else
    for s in input families scores; do [ -e $d/stages/$s ] || cp -al $B/stages/$s $d/stages/$s; done
    cmd="python -u -m ancestra run $args"
  fi
  echo "$cmd" > $d/cmd.txt
  echo "[$(date +%T)] $lab start: $cmd" >> work/ob_final/variants.log
  ( /usr/bin/time -v -o $d/time.txt $cmd 2>&1 \
      | while IFS= read -r line; do printf '%s\t%s\n' "$(date +%s.%N)" "$line"; done > $d/run.log ) &
  sleep 1
  P=$(pgrep -f "time -v -o $d/time.txt" | head -1)
  python $S/proc_sampler.py $P $d/proc_samples.tsv 2 &
  wait
  echo "[$(date +%T)] $lab done ($(grep 'Elapsed' $d/time.txt | awk '{print $NF}'))" >> work/ob_final/variants.log
done
