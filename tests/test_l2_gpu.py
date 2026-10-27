"""L2-gpu options (docs/L2_gpu.md): GPU-graph sparsification and batched local trees; both opt-in."""
import shutil

import numpy as np
import pytest

from ancestra.merge import MergeConfig
from ancestra.scoring import ScoringConfig, sparsify_mask


def test_keys_unchanged_at_defaults():
    assert "gpu_sparsify" not in ScoringConfig().key()
    assert "gpu_sparsify" not in ScoringConfig(backend="gpusw").key()
    assert ScoringConfig(backend="gpusw", gpu_sparsify="topk:5").key()["gpu_sparsify"] == "topk:5"
    k = MergeConfig().key()
    assert "tree_ft_batch" not in k and "tree_ft_inproc" not in k


def _toy():
    # genes 0..5; species: 0,0,1,1,2,2
    sp = np.array([0, 0, 1, 1, 2, 2])
    q = np.array([0, 0, 0, 0, 1, 1, 2, 3, 4])
    t = np.array([1, 2, 3, 4, 2, 3, 4, 5, 5])
    b = np.array([50., 100., 40., 80., 30., 90., 70., 20., 60.])
    selfb = np.full(6, 200.)
    return q, t, b, selfb, sp


def test_topk_either_direction():
    q, t, b, selfb, sp = _toy()
    m = sparsify_mask(q, t, b, selfb, sp, "topk:1")
    # 0's best in sp1 is 2 (100) -> (0,2) kept; 1's best in sp1 is 3 (90) -> (1,3) kept;
    # (0,3): 0 prefers 2, 3 prefers 1 in sp0 -> dropped; (1,2): 1 prefers 3, 2 prefers 0 -> dropped
    want = {(0, 2), (1, 3), (0, 1), (0, 4), (2, 4), (4, 5), (3, 5)}
    got = {(int(a), int(c)) for a, c, k in zip(q, t, m) if k}
    assert got == want
    assert sparsify_mask(q, t, b, selfb, sp, "topk:9").all()


def test_ratio_norm_rand():
    q, t, b, selfb, sp = _toy()
    r = sparsify_mask(q, t, b, selfb, sp, "ratio:0.9")
    assert r[1] and r[5] and not r[4]
    n = sparsify_mask(q, t, b, selfb, sp, "norm:0.3")
    assert (n == (b / 200 >= 0.3)).all()
    a = sparsify_mask(q, t, b, selfb, sp, "rand:0.5:11")
    assert (a == sparsify_mask(q, t, b, selfb, sp, "rand:0.5:11")).all()
    with pytest.raises(ValueError):
        sparsify_mask(q, t, b, selfb, sp, "bogus:1")


def test_topk_order_independent():
    q, t, b, selfb, sp = _toy()
    p = np.random.default_rng(3).permutation(len(q))
    m1 = sparsify_mask(q, t, b, selfb, sp, "topk:1")
    m2 = sparsify_mask(q[p], t[p], b[p], selfb, sp, "topk:1")
    assert (m1[p] == m2).all()


@pytest.mark.skipif(shutil.which("FastTree") is None or shutil.which("famsa") is None,
                    reason="FastTree/famsa not installed")
def test_fasttree_batch_equals_single(tmp_path):
    from ancestra.localtree import run_fasttree_jobs
    rng = np.random.default_rng(5)
    aa = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    jobs = []
    for j in range(6):
        base = rng.choice(aa, 120)
        seqs = []
        for i in range(4 + j % 3):
            s = base.copy()
            mut = rng.random(120) < 0.1 * (i + 1)
            s[mut] = rng.choice(aa, mut.sum())
            seqs.append("".join(s))
        jobs.append(([f"L{i}" for i in range(len(seqs))], seqs, str(tmp_path)))
    single = run_fasttree_jobs(jobs, 2, 0)
    batched = run_fasttree_jobs(jobs, 2, 4)
    assert single == batched


def test_mutual_topk():
    q, t, b, selfb, sp = _toy()
    m = sparsify_mask(q, t, b, selfb, sp, "mtopk:1")
    e = sparsify_mask(q, t, b, selfb, sp, "topk:1")
    assert (m <= e).all()
    got = {(int(a), int(c)) for a, c, k in zip(q, t, m) if k}
    assert (0, 2) in got and (1, 3) in got and (0, 3) not in got
