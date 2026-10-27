#!/bin/bash
# L2-domains: Pfam-A (38.2) annotation with hmmsearch --cut_ga, NPROC single-thread processes
# over sequence chunks (hmmsearch --cpu scales poorly; one HMM parse per chunk).
#   scripts/L2d_hmmsearch.sh IN.fa OUTDIR [NPROC=8]
set -euo pipefail
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
IN=$1; OUT=$2; NP=${3:-8}
HMM=$ROOT/data/extra/l2_domains/Pfam-A.hmm
BIN=$ROOT/conda_envs/orth/bin
mkdir -p $OUT/chunks
# split by sequence into NP chunks, round-robin (balanced length)
awk -v np=$NP -v out=$OUT/chunks '/^>/{n++; f=sprintf("%s/c%02d.fa", out, n%np)} {print > f}' $IN
start=$(date +%s)
for f in $OUT/chunks/c*.fa; do
  b=$(basename $f .fa)
  ( /usr/bin/time -v nice -n 5 $BIN/hmmsearch --cut_ga --cpu 1 --noali -o /dev/null \
      --domtblout $OUT/chunks/$b.domtbl $HMM $f ) 2> $OUT/chunks/$b.time &
done
wait
cat $OUT/chunks/c*.domtbl | grep -v '^#' > $OUT/pfam.domtbl
echo "wall_s $(( $(date +%s) - start ))" > $OUT/runtime.txt
grep -h "User time\|System time" $OUT/chunks/*.time | awk '{s+=$NF} END {print "cpu_s", s}' >> $OUT/runtime.txt
md5sum $IN $HMM >> $OUT/runtime.txt
echo done >> $OUT/runtime.txt
