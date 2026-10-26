#!/usr/bin/env bash
# usage: calib_timing.sh proteome.fasta  (24 threads, all-vs-all self search timings)
source /opt/conda/etc/profile.d/conda.sh && conda activate /scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/conda_envs/orth
F=$1; B=$(basename $F .fasta); O=$(dirname $F)/out_$B; mkdir -p $O; : > $O/timings.tsv
t(){ n=$1; shift; s=$(date +%s.%N); "$@" > $O/$n.log 2>&1; rc=$?; e=$(date +%s.%N); python -c "print('$n\t%.1f\trc=$rc'%($e-$s))" | tee -a $O/timings.tsv; }
t diamond_makedb diamond makedb --in $F -d $O/db -p 24
for m in fast mid-sensitive sensitive very-sensitive ultra-sensitive; do t diamond_$m diamond blastp -q $F -d $O/db -o $O/d_$m.tsv --$m -p 24 -k 500 -e 1e-3 --outfmt 6 qseqid sseqid pident length evalue bitscore; done
t mmseqs_s5.7 mmseqs easy-search $F $F $O/mm57.tsv $O/tmp --threads 24 -s 5.7 --max-seqs 500; rm -rf $O/tmp
t mmseqs_s7.5 mmseqs easy-search $F $F $O/mm75.tsv $O/tmp --threads 24 -s 7.5 --max-seqs 500; rm -rf $O/tmp
wc -l $O/*.tsv
