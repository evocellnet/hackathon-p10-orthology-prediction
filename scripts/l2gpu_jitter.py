#!/usr/bin/env python
"""L2-gpu noise-floor control: multiply every non-self bit score of a table by exp(sigma * N(0,1)) (seeded);
pair set and self-scores unchanged.   python scripts/l2gpu_jitter.py IN.npz SIGMA SEED OUT.npz"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ancestra.scoring import ScoreTable  # noqa: E402

i, sig, seed, out = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
st = ScoreTable.load(i)
rng = np.random.default_rng(seed)
st.bits = (st.bits.astype(np.float64) * np.exp(sig * rng.standard_normal(len(st.bits)))).astype(np.float32)
st.save(out)
print("jittered", len(st.bits), "sigma", sig, "seed", seed)
