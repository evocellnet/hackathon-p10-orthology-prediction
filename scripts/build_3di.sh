#!/usr/bin/env bash
# Reproduce the 3Di (Foldseek structural alphabet) tables for QfO 2022_02 and OrthoBench.
#   data/structure/qfo_3di_best.fa         (>UPID|acc      3Di)   + .source.tsv + .coverage.tsv
#   data/structure/orthobench_3di_best.fa  (>ensembl|acc   3Di)   + .source.tsv + .coverage.tsv
# Only proteins whose AFDB v6 sequence is IDENTICAL to ours are kept.  Source priority per protein:
#   1. fsexact  : true AFDB-v6 3Di from a Foldseek AFDB DB (Alphafold/UniProt if extracted, else Alphafold/Proteome),
#                 matched by exact sequence (own accession preferred, otherwise any AFDB entry with the same sequence)
#   2. foldcomp : 3Di from `foldseek createdb --input-format 5` on the Foldcomp AFDB-UniProt-v6 DB (lossy coordinate
#                 compression => ~97.5 % per-residue identity to the true 3Di, see docs/structure_feasibility.md)
# Steps are idempotent (skip if output exists).  Env vars: SKIP_DOWNLOAD=1, THREADS=16.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
S=$ROOT/data/structure; FC=$S/foldcomp; W=$FC/work; V6=$FC/v6
PY="python $ROOT/scripts/fc3di.py"
THREADS=${THREADS:-16}
source /opt/conda/etc/profile.d/conda.sh && conda activate "$ROOT/conda_envs/orth"
mkdir -p "$W"
OB_IN=$ROOT/data/orthobench/Open_Orthobench/BENCHMARKS/Input
V6_SIZE=1261454421092
log() { echo "[$(date +%H:%M:%S)] $*"; }

# 0. QfO canonical proteome FASTAs (UP*_taxid.fasta) and accession lists
if [ ! -d "$S/qfo_fasta" ] || [ "$(ls "$S/qfo_fasta" | wc -l)" -lt 78 ]; then
  log "extract QfO FASTAs"; mkdir -p "$S/qfo_fasta" "$W/qfo_tar"
  tar -xzf "$ROOT/data/qfo/proteomes_2022_02/QfO_release_2022_02.tar.gz" -C "$W/qfo_tar" --wildcards '*/UP*.fasta'
  find "$W/qfo_tar" -regextype posix-extended -regex '.*/UP[0-9]+_[0-9]+\.fasta' -exec cp {} "$S/qfo_fasta/" \;
fi
[ -f "$S/qfo_accessions.json" ] || python - <<EOF
import json,glob,os
d={}
for f in sorted(glob.glob("$S/qfo_fasta/*.fasta")):
    d[os.path.basename(f).split('_')[0]]=[l.split('|')[1] for l in open(f) if l.startswith('>')]
json.dump(d,open("$S/qfo_accessions.json",'w'))
EOF

# 1. Foldcomp AFDB-UniProt-v6 (1.26 TB; ~1.5-3 h at 150-350 MB/s; resumable)
if [ "${SKIP_DOWNLOAD:-0}" != 1 ]; then
  if [ ! -f "$V6/afdb_uniprot_v6" ] || [ -f "$V6/afdb_uniprot_v6.aria2" ] || [ "$(stat -c %s "$V6/afdb_uniprot_v6")" != $V6_SIZE ]; then
    log "download Foldcomp afdb_uniprot_v6"; "$FC/download_v6.sh"
  fi
fi
[ "$(stat -c %s "$V6/afdb_uniprot_v6")" = $V6_SIZE ] && [ ! -f "$V6/afdb_uniprot_v6.aria2" ] || { echo "Foldcomp v6 data file incomplete; re-run $FC/download_v6.sh"; exit 2; }

# 2. UniProt side inputs: taxon FASTAs of the 12 OrthoBench species (candidate Ensembl->UniProt by exact sequence), sec_ac
mkdir -p "$S/uniprot_taxon_fasta"
for t in 6239 9615 7719 7955 7227 9031 9606 13616 10090 9598 10116 99883; do
  f=$S/uniprot_taxon_fasta/tax_$t.fasta.gz
  [ -s "$f" ] && gzip -t "$f" 2>/dev/null || curl -sfL "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:$t&format=fasta&compressed=true" -o "$f"
done
[ -s "$FC/sec_ac.txt" ] || curl -sfL https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt -o "$FC/sec_ac.txt"

# 3. wanted accessions -> Foldcomp keys
[ -s "$W/ob_candidates.tsv" ] || $PY obcand "$OB_IN" "$W/ob_candidates.tsv" "$W/ob_cand_acc.txt" "$S"/uniprot_taxon_fasta/tax_*.fasta.gz "$S"/qfo_fasta/*.fasta
[ -s "$W/wanted_acc.txt" ] || $PY wanted "$S/qfo_accessions.json" "$FC/sec_ac.txt" "$W/wanted_acc.txt" "$W/qfo_alias.tsv" "$W/ob_cand_acc.txt"
[ -s "$W/keys_v6.txt" ] || $PY select "$V6/afdb_uniprot_v6.lookup" "$W/wanted_acc.txt" "$W/keys_v6.txt" "$W/key2acc_v6.tsv" --ver v6

# 4. sub-DB (soft link into the 1.2 TB file) and 3Di by foldseek (Foldcomp input = --input-format 5)
if [ ! -s "$W/sub_fs_ss" ]; then
  log "createsubdb + createdb"
  foldseek createsubdb --id-mode 0 --subdb-mode 1 "$W/keys_v6.txt" "$V6/afdb_uniprot_v6" "$W/sub" -v 2
  foldseek createdb "$W/sub" "$W/sub_fs" --input-format 5 --threads "$THREADS" -v 2
fi

# 5. Foldcomp-derived tables (exact AA sequence only; QfO also via UniProt merged/secondary accessions)
$PY qfo "$W/sub_fs" "$S/qfo_fasta" "$W/qfo_3di_foldcomp.fa" "$W/qfo_foldcomp.stats.tsv" --alias "$W/qfo_alias.tsv"
$PY ob "$W/sub_fs" "$OB_IN" "$W/ob_candidates.tsv" "$W/orthobench_3di_foldcomp.fa" "$W/ob_foldcomp.stats.tsv"

# 6. True-3Di tables from the Foldseek AFDB DB by exact sequence (Alphafold/UniProt preferred, else Alphafold/Proteome)
[ "${SKIP_FSDB:-0}" = 1 ] || [ -f "$S/afdb_uniprot_fs/ex/.extracted" ] || "$ROOT/scripts/get_foldseek_afdb_uniprot.sh"   # 491 GB download + parallel extract
FSDB=$S/afdb_uniprot_fs/ex/afdb; FSLAB=fsuniprot
if [ ! -f "$S/afdb_uniprot_fs/ex/.extracted" ]; then FSDB=$S/afdb_proteome_db/afdbp; FSLAB=fsproteome; fi
log "exact-sequence route on $FSDB"
$PY fsexact "$FSDB" "$S/qfo_fasta" "$OB_IN" "$W/$FSLAB" --procs "$THREADS"

# 7. merge (priority: exact Foldseek 3Di > Foldcomp 3Di)
$PY merge qfo "$S/qfo_fasta" "$S/qfo_3di_best" "$FSLAB=$W/$FSLAB.qfo.fa" "foldcomp=$W/qfo_3di_foldcomp.fa"
$PY merge ob  "$OB_IN"       "$S/orthobench_3di_best" "$FSLAB=$W/$FSLAB.ob.fa" "foldcomp=$W/orthobench_3di_foldcomp.fa"
log done; cat "$S/qfo_3di_best.coverage.tsv" | tail -1; tail -1 "$S/orthobench_3di_best.coverage.tsv"
