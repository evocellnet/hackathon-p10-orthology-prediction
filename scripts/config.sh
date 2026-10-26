# Shared configuration; override any variable from the environment.
BENCH_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Raw data lives outside the benchmark repo (never committed). Default: ../data next to bench_repo
DATA_DIR="${DATA_DIR:-$(cd "$BENCH_ROOT/.." && pwd)/data}"
# python env with lxml dendropy biopython jsonschema pandas (see scripts/setup_env.sh)
PYTHON="${PYTHON:-$(cd "$BENCH_ROOT/.." && pwd)/conda_envs/data/bin/python}"
[[ -x "$PYTHON" ]] || PYTHON="$(command -v python3)"
DATA_DIR="$(readlink -f -m "$DATA_DIR")"
export BENCH_ROOT DATA_DIR PYTHON
