#!/usr/bin/env bash
# Create the small python env used by the benchmark wrappers (QfO python steps, OpenEBench fetchers)
# plus mafft (needed by the QfO G_STD2 tests). Default location: <project>/conda_envs/data
set -euo pipefail
source "$(dirname "$0")/config.sh"
ENV="${1:-$(cd "$BENCH_ROOT/.." && pwd)/conda_envs/data}"
CONDA="${CONDA:-$(command -v mamba || command -v conda || echo /opt/conda/bin/mamba)}"
"$CONDA" create -y -p "$ENV" -c conda-forge -c bioconda python=3.10 lxml biopython jsonschema \
    matplotlib pandas requests numpy scipy pip mafft
"$ENV/bin/pip" install dendropy ete3 tqdm
echo "env ready: $ENV  (export PYTHON=$ENV/bin/python)"
