#!/usr/bin/env python
"""L2-gpu Pareto splice: GPU table for families <= T genes, CPU (DIAMOND pass-1) table for larger ones
(rows and self-scores). No cross-family rescue.
    python scripts/l2gpu_splice.py GPU.npz CPU.npz FAM.npy T OUT.npz"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ancestra.scoring import ScoreTable  # noqa: E402

g, c, f, T, out = sys.argv[1:6]
T = int(T)
fam = np.load(f)
big_gene = np.bincount(fam)[fam] > T
G, C = ScoreTable.load(g), ScoreTable.load(c)
mg, mc = ~big_gene[G.q], big_gene[C.q]
cols = {k: np.concatenate([getattr(G, k)[mg], getattr(C, k)[mc]]) for k in ("q", "t", "bits", "qcov", "tcov", "pident", "evalue")}
o = np.lexsort((cols["t"], cols["q"]))
sb = np.where(big_gene, C.self_bits, G.self_bits).astype(np.float32)
ScoreTable(**{k: v[o] for k, v in cols.items()}, self_bits=sb).save(out)
print("T", T, "rows gpu", int(mg.sum()), "rows cpu", int(mc.sum()))
