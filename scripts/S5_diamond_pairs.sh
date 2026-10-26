#!/usr/bin/env bash
# Whole-proteome DIAMOND searches for the S5 species-pair benchmark (both directions).
# Output: work/S5/diamond/<A>__<B>.tsv  (qseqid sseqid bitscore evalue pident qcovhsp scovhsp)
set -euo pipefail
cd "$(dirname "$0")/.."
source /opt/conda/etc/profile.d/conda.sh; conda activate "$PWD/conda_envs/orth"
IN=data/orthobench/Open_Orthobench/BENCHMARKS/Input
D=work/S5/diamond; mkdir -p $D
T=${THREADS:-8}
PAIRS=${PAIRS:-"Homo_sapiens:Mus_musculus Homo_sapiens:Gallus_gallus Homo_sapiens:Danio_rerio Danio_rerio:Tetraodon_nigroviridis Homo_sapiens:Drosophila_melanogaster Homo_sapiens:Tetraodon_nigroviridis"}
for sp in $(echo $PAIRS | tr ' :' '\n\n' | sort -u); do
  [ -s $D/$sp.dmnd ] || diamond makedb --in $IN/$sp.*.pep.all.fa -d $D/$sp --quiet
done
for pr in $PAIRS; do
  A=${pr%%:*}; B=${pr##*:}
  for dir in "$A $B" "$B $A"; do
    set -- $dir
    out=$D/$1__$2.tsv
    [ -s $out ] && continue
    /usr/bin/time -v diamond blastp -q $IN/$1.*.pep.all.fa -d $D/$2.dmnd -o $out.tmp --sensitive \
      -k 100 --max-hsps 1 -e 1e-3 -p $T -b 2 -c 1 --quiet \
      --outfmt 6 qseqid sseqid bitscore evalue pident qcovhsp scovhsp 2> $D/$1__$2.time
    mv $out.tmp $out
  done
done
echo done
