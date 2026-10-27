"""L2-plm shared helpers (benchmark/exploration code; never imported by the method)."""
import os, sys, json, hashlib, subprocess, time
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# data/work live in the main checkout (gitignored); the worktree sits at <main>/work/l2_plm_wt
DATA = HERE if os.path.exists(os.path.join(HERE, "work/core_v0_ob")) else os.path.abspath(os.path.join(HERE, "../.."))
RUN = os.path.join(DATA, "work/core_v0_ob/run_v2")
WORK = os.path.join(DATA, "work/l2_plm")
RES = os.path.join(HERE, "results/l2_plm")
HF = os.path.join(DATA, "data/extra/l2_plm/hf")
SEEDS = (11, 23, 47)
DEEP = {"vert-teleost", "vert-invert", "invert-invert"}
for p in (HERE, os.path.join(HERE, "scripts"), DATA):
    if p not in sys.path:
        sys.path.insert(0, p)


def git_commit():
    try:
        return subprocess.check_output(["git", "-C", HERE, "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return None


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def manifest(config, inputs, t0, extra=None):
    m = dict(git_commit=git_commit(), config=config, inputs=inputs, runtime_s=round(time.time() - t0, 1),
             finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"), host=os.uname().nodename)
    if extra:
        m.update(extra)
    return m


def gene_table():
    from ancestra.io import read_gene_table
    return read_gene_table(os.path.join(RUN, "stages/input/genes.tsv"))


def read_seqs(names):
    seqs, cur = {}, None
    want = set(names)
    for line in open(os.path.join(DATA, "work/S0/ob_all.fa")):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            cur = cur if cur in want else None
            if cur:
                seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    return ["".join(seqs[n]).rstrip("*") for n in names]
