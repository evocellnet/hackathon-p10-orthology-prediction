#!/usr/bin/env bash
# Usage: scale_large_run.sh N [PRESET=fast] [TAG=n{N}]
# Detached-safe wrapper: symlink input dir, /usr/bin/time -v, 60 s uptime+free log, RSS watchdog (kills run above LIMIT_GB).
set -u
N=$1; PRESET=${2:-fast}; TAG=${3:-n$N}
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
LIMIT_GB=${LIMIT_GB:-100}
cd $ROOT
source /opt/conda/etc/profile.d/conda.sh && conda activate $ROOT/conda_envs/orth
W=$ROOT/work/scale_large/$TAG
rm -rf $W; mkdir -p $W/input $W/logs
while read s; do ln -s $ROOT/data/scale1000/proteomes/$s.fasta.gz $W/input/$s.fasta.gz; done < $ROOT/results/scaling_large/subsets/n$N.txt
date -Iseconds > $W/logs/start.txt
(
  while true; do
    echo "== $(date -Iseconds)"; uptime; free -g | sed -n 1,2p
    # RSS (GB) of the python run tree = all processes whose cmdline mentions this work dir
    rss=$(ps -eo rss,args | grep -F "$W/" | grep -v grep | awk '{s+=$1} END{printf "%.1f", s/1048576}')
    echo "my_tree_rss_gb $rss"
    if awk "BEGIN{exit !($rss > $LIMIT_GB)}"; then echo "WATCHDOG: rss $rss > $LIMIT_GB, killing"; pkill -f "$W/" ; fi
    sleep 60
  done
) > $W/logs/monitor.log 2>&1 &
MON=$!
/usr/bin/time -v -o $W/logs/time_v.txt python -m ancestra run --preset $PRESET --proteomes $W/input \
   --tree $ROOT/results/scaling_large/trees/ncbi_n$N.nwk --out $W/out --threads 12 > $W/logs/run.log 2>&1
echo $? > $W/logs/exit_code.txt
date -Iseconds > $W/logs/end.txt
kill $MON
echo finished >> $W/logs/run.log
