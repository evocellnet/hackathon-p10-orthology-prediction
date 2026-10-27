#!/usr/bin/env bash
# Usage: launch.sh <arm> <omamer_db>  ; FastOMA 0.5.1 (local src checkout, tag v0.5.1 = 553f860) + gene-loss patch (patched_pkg via PYTHONPATH)
set -u
R=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
D=$R/data/orthobench/runs/fastoma_v2; A=$1; DB=$2
WORK=$3
source /opt/conda/etc/profile.d/conda.sh && conda activate $R/conda_envs/orth
export PYTHONPATH=$D/patched_pkg FASTOMA_RESCUE_LOG=$D/$A/rescued_subhogs.tsv NXF_HOME=$R/data/orthobench/runs/fastoma/nxf_home
ulimit -n 131072 2>/dev/null || ulimit -n 65536
cd $D/$A/run
EXTRA=""; [ -d $D/$A/in/hogmap_in ] && EXTRA="--hogmap_in $D/$A/in/hogmap_in"
exec bash $R/bench_repo/scripts/run_competitor_timed.sh $A $D/$A/logs -- nextflow run $R/data/orthobench/runs/fastoma/src/FastOMA.nf -c $D/$A/local.config \
  --input $D/$A/in --proteome_folder $D/$A/in/proteome --species_tree $D/$A/in/species_tree.nwk $EXTRA \
  --omamer_db $DB --output_folder $D/$A/run/Output -w $WORK -with-trace -with-report ${RESUME:-}
