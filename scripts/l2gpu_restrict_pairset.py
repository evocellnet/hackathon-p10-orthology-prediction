#!/usr/bin/env python
"""L2-gpu diagnostic X1: GPU score table restricted to the CPU (DIAMOND) pair set in families with
> 1000 genes (GPU scores, CPU pair set); smaller families and self-scores unchanged.
    python scripts/l2gpu_restrict_pairset.py GPU.npz CPU.npz FAM.npy OUT.npz"""
import sys
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
from ancestra.scoring import ScoreTable  # noqa: E402
g, c, f, out = sys.argv[1:5]
fam = np.load(f)
sz = np.bincount(fam)
st = ScoreTable.load(g)
zc = np.load(c)
N = len(fam)
big = sz[fam[st.q]] > 1000
kc = zc["q"].astype(np.int64) * N + zc["t"]
kg = st.q.astype(np.int64) * N + st.t
keep = ~big | np.isin(kg, kc)
ScoreTable(*(getattr(st, k)[keep] for k in ("q", "t", "bits", "qcov", "tcov", "pident", "evalue")),
           st.self_bits).save(out)
print("kept", int(keep.sum()), "of", len(keep), "big kept", int((keep & big).sum()))
