#!/usr/bin/env bash
# Foldseek "Alphafold/UniProt" DB (AFDB v6, 241 M entries, true 3Di). Replaces `foldseek databases Alphafold/UniProt`
# whose single-threaded `tar xvfz` runs at ~30-90 MB/s; here aria2c + parallel rapidgzip (~340 MB/s), skipping afdb_ca/afdb_h.
# Result: data/structure/afdb_uniprot_fs/ex/afdb{,_ss}{,.index,.dbtype} afdb.lookup afdb_mapping afdb_taxonomy
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; D=$ROOT/data/structure/afdb_uniprot_fs
source /opt/conda/etc/profile.d/conda.sh && conda activate "$ROOT/conda_envs/orth"   # needs aria2c, rapidgzip (pip install rapidgzip)
mkdir -p "$D/ex" "$D/dl"; cd "$D"
if [ ! -f dl/afdb.tar.gz ] || [ -f dl/afdb.tar.gz.aria2 ] || [ "$(stat -c %s dl/afdb.tar.gz)" != 491470141753 ]; then
  aria2c -x5 -s5 -c --file-allocation=none --summary-interval=300 -d dl -o afdb.tar.gz https://foldseek.steineggerlab.workers.dev/afdb.tar.gz
fi
[ -s ex/afdb.version ] || curl -sfL https://foldseek.steineggerlab.workers.dev/afdb.version -o ex/afdb.version
[ -f ex/.extracted ] || { rapidgzip -d -c -P "${THREADS:-12}" dl/afdb.tar.gz | tar -x -C ex --exclude='afdb_ca*' --exclude='afdb_h*' && touch ex/.extracted; }
