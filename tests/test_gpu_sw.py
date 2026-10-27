"""v2-gpu: tests of the GPU Smith-Waterman kernels (ancestra.gpu_sw) and the gpusw backend key.

The kernel tests run the GPU code in the GPU python (conda_envs/gpu or $ANCESTRA_GPU_PYTHON) and are
skipped when it or a CUDA device is unavailable. The key test runs everywhere.
"""
import os
import subprocess
import sys
import textwrap

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def _gpu_python():
    from ancestra.scoring import gpu_python
    try:
        p = gpu_python()
    except RuntimeError:
        return None
    r = subprocess.run([p, "-c", "import cupy; cupy.cuda.runtime.getDeviceCount()"], capture_output=True)
    return p if r.returncode == 0 else None


GPU_PY = _gpu_python()
needs_gpu = pytest.mark.skipif(GPU_PY is None, reason="no GPU python with cupy / no CUDA device")


def _run_gpu(code: str) -> str:
    env = dict(os.environ, PYTHONPATH=ROOT)
    r = subprocess.run([GPU_PY, "-c", textwrap.dedent(code)], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    return r.stdout


@needs_gpu
def test_kernels_match_reference_and_each_other():
    out = _run_gpu("""
        import random, numpy as np
        import ancestra.gpu_sw as G
        rng = random.Random(5)
        AA = "ACDEFGHIKLMNPQRSTVWY"
        def mut(s, p):
            o = []
            for c in s:
                r = rng.random()
                if r < p: o.append(rng.choice(AA))
                elif r < 1.1 * p: continue
                elif r < 1.2 * p: o.append(c + rng.choice(AA))
                else: o.append(c)
            return "".join(o)
        seqs = []
        for _ in range(30):
            b = "".join(rng.choice(AA) for _ in range(rng.randint(3, 400)))
            seqs += [b, mut(b, 0.35)]
        seqs.append("XUBZ*acdOJ")
        seq, off, lens = G.encode(seqs)
        q, t = np.triu_indices(len(seqs))
        keep = np.array([rng.random() < 0.15 or a == b for a, b in zip(q, t)])
        q, t = q[keep], t[keep]
        r1 = G.align_pairs(seq, off, lens, q, t, evalue=1e30)
        bad = sum(G.sw_reference(seqs[a], seqs[b]) != s for a, b, s in zip(q, t, r1["score"]))
        G.LONG_CELLS = 0                      # everything through the warp kernel
        r2 = G.align_pairs(seq, off, lens, q, t, evalue=1e30, sw=G.GpuSW(w=16))
        same = all((r1[k] == r2[k]).all() for k in ("score", "qstart", "qend", "tstart", "tend"))
        sub = 0
        for k in np.flatnonzero(r1["score"] > 25)[:40]:
            a = seqs[q[k]][r1["qstart"][k]:r1["qend"][k] + 1]
            b = seqs[t[k]][r1["tstart"][k]:r1["tend"][k] + 1]
            sub += G.sw_reference(a, b) != r1["score"][k]
        print(len(q), bad, int(same), sub, r1["rev_mismatch"])
    """)
    n, bad, same, sub, rev = map(int, out.split())
    assert n > 100 and bad == 0 and same == 1 and sub == 0 and rev == 0


@needs_gpu
def test_family_mode_deterministic():
    out = _run_gpu("""
        import hashlib, random, numpy as np
        import ancestra.gpu_sw as G
        rng = random.Random(9)
        AA = "ACDEFGHIKLMNPQRSTVWY"
        seqs = ["".join(rng.choice(AA) for _ in range(rng.randint(20, 300))) for _ in range(300)]
        fam = np.array([i % 17 for i in range(300)])
        seq, off, lens = G.encode(seqs)
        h = []
        for chunk in (10**9, 500):
            r = G.score_families_gpu(seq, off, lens, fam, db_residues=5e7, chunk_pairs=chunk)
            m = hashlib.md5()
            for k in ("q", "t", "bits", "qcov", "tcov", "evalue"):
                m.update(np.ascontiguousarray(r[k]).tobytes())
            h.append(m.hexdigest())
        print(int(h[0] == h[1]), r["stats"]["n_pairs"])
    """)
    same, n = map(int, out.split())
    assert same == 1 and n == sum(k * (k + 1) // 2 for k in [18] * 11 + [17] * 6)


@needs_gpu
def test_composition_bias_kernels():
    out = _run_gpu("""
        import random, numpy as np
        import ancestra.gpu_sw as G
        rng = random.Random(3)
        AA = "ACDEFGHIKLMNPQRSTVWY"
        seqs = []
        for _ in range(25):
            b = "".join(rng.choice(AA if rng.random() < 0.7 else "QQQQSSPG") for _ in range(rng.randint(10, 200)))
            seqs += [b, "".join(c if rng.random() > 0.3 else rng.choice(AA) for c in b)]
        seq, off, lens = G.encode(seqs)
        bias = G.comp_bias(seq, off, lens)
        same_bias = bool((G.comp_bias_gpu(seq, off, lens, chunk=1000) == bias).all())
        q, t = np.triu_indices(len(seqs))
        sel = np.array([rng.random() < 0.12 for _ in q])
        q, t = q[sel], t[sel]
        r = G.align_pairs(seq, off, lens, q, t, evalue=1e30, bias=bias)
        bad = 0
        for k in range(len(q)):
            bq = bias[off[q[k]]:off[q[k]] + lens[q[k]]]
            bt = bias[off[t[k]]:off[t[k]] + lens[t[k]]]
            ref = max(G.sw_reference(seqs[q[k]], seqs[t[k]], ba=bq), G.sw_reference(seqs[q[k]], seqs[t[k]], bb=bt))
            bad += ref != r["score"][k]
        G.LONG_CELLS = 0
        r2 = G.align_pairs(seq, off, lens, q, t, evalue=1e30, bias=bias, sw=G.GpuSW(bias=True, w=16))
        same = all((r[k] == r2[k]).all() for k in ("score", "qstart", "qend", "tstart", "tend"))
        print(int(same_bias), bad, int(same), r["rev_mismatch"], int((bias != 0).any()))
    """)
    same_bias, bad, same, rev, nz = map(int, out.split())
    assert same_bias == 1 and bad == 0 and same == 1 and rev == 0 and nz == 1


def test_cpu_backend_keys_unchanged_by_gpu_fields():
    from ancestra.scoring import ScoringConfig
    for b in ("hybrid", "diamond", "mmseqs"):
        k = ScoringConfig(backend=b).key()
        assert not any(x.startswith("gpu_") for x in k), b
    k = ScoringConfig(backend="gpusw").key()
    assert {"gpu_max_family", "gpu_mask", "gpu_db_residues", "gpu_comp_bias"} <= set(k) and "gpu_device" not in k
