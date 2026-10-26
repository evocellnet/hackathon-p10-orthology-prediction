#!/usr/bin/env python
"""Derive a reversible 3Di substitution model (PAML format, for IQ-TREE ``-m FILE``) from
Foldseek's 3Di log-odds matrix (``mat3di.out``, embedded in the foldseek binary; extracted
to work/S4/mat3di.out with ``strings``; scale bit/2 with lambda 0.351568 nats per unit).

Joint p_ij ∝ f_i f_j exp(lambda S_ij) (f = Foldseek background), symmetric; pi = marginal;
P = p_ij / pi_i is a transition matrix at some divergence t*; Q = logm(P) (rates of a
reversible chain, as in the BLOSUM -> rate-matrix route); exchangeabilities r_ij = Q_ij / pi_j
(symmetrised; negative values from the log clipped to 1e-4 of the mean). IQ-TREE rescales Q to one
substitution per site, so t* does not matter. This is an approximation of a 3Di
model (the published 3DI model of Puente-Lelievre et al. 2024 was not available offline).

    python scripts/S4_3di_model.py  -> work/S4/3DI_foldseek.paml (+ .json diagnostics)
"""
import json, os, sys
import numpy as np
from scipy.linalg import logm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FS_ORDER = "ACDEFGHIKLMNPQRSTVWY"
PAML_ORDER = "ARNDCQEGHILKMFPSTWYV"


def read_mat(path):
    lines = [l for l in open(path) if l.strip()]
    bg = lam = None
    rows = {}
    for l in lines:
        if l.startswith("# Background"):
            bg = np.array([float(x) for x in l.split(":")[1].split()][:20])
        elif l.startswith("# Lambda"):
            lam = float(l.split(":")[1])
        elif l[0] in FS_ORDER and not l.startswith("#"):
            parts = l.split()
            rows[parts[0]] = [float(x) for x in parts[1:21]]
    S = np.array([rows[a] for a in FS_ORDER])
    return S, bg / bg.sum(), lam


def main():
    S, f, lam = read_mat(os.path.join(ROOT, "work/S4/mat3di.out"))
    J = np.outer(f, f) * np.exp(lam * S)
    J = (J + J.T) / 2
    J /= J.sum()
    pi = J.sum(1)
    P = J / pi[:, None]
    Q = np.real(logm(P))
    R = Q / pi[None, :]
    R = (R + R.T) / 2
    np.fill_diagonal(R, 0)
    off = R[~np.eye(20, dtype=bool)]
    n_neg = int((off <= 0).sum())
    floor = 1e-4 * off[off > 0].mean()
    R = np.where(R <= 0, floor, R)
    np.fill_diagonal(R, 0)
    idx = [FS_ORDER.index(a) for a in PAML_ORDER]
    Rp, pip = R[np.ix_(idx, idx)], pi[idx]
    Rp = Rp / Rp[np.tril_indices(20, -1)].mean()
    out = os.path.join(ROOT, "work/S4/3DI_foldseek.paml")
    with open(out, "w") as fh:
        for i in range(1, 20):
            fh.write(" ".join(f"{Rp[i, j]:.6f}" for j in range(i)) + "\n")
        fh.write("\n" + " ".join(f"{x:.6f}" for x in pip) + "\n")
    # FastTree -trans file: rate matrix Q transposed (columns = from-state, each sums to 0) + "*" column = pi
    Qp = Rp * pip[None, :]
    np.fill_diagonal(Qp, 0)
    np.fill_diagonal(Qp, -Qp.sum(1))
    Qp /= -(pip * np.diag(Qp)).sum()
    with open(out.replace(".paml", ".fasttree.tsv"), "w") as fh:
        fh.write("\t".join(list(PAML_ORDER) + ["*"]) + "\n")
        for i, a in enumerate(PAML_ORDER):
            fh.write(a + "\t" + "\t".join(f"{x:.8f}" for x in Qp[:, i]) + f"\t{pip[i]:.8f}\n")
    diag = dict(lambda_=lam, n_nonpositive_exchangeabilities_clipped=n_neg,
                expected_identity_at_tstar=float(np.trace(J)), pi_fs_order=dict(zip(FS_ORDER, pi.round(5).tolist())))
    json.dump(diag, open(out.replace(".paml", ".json"), "w"), indent=1)
    print(diag)


if __name__ == "__main__":
    main()
