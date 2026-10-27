"""Protein domain architectures (Pfam) as an optional merge layer (L2-domains, docs/L2_domains.md).

Off by default: nothing in the frozen presets imports this module.

* :func:`parse_domtbl` reads ``hmmsearch --domtblout`` (Pfam-A, ``--cut_ga``), resolves overlaps (domains
  sorted by i-Evalue; a domain is dropped if it overlaps a kept one by > ``max_overlap`` of the shorter
  envelope) and returns per protein the ordered list of domain labels (clan if the family has one, else the
  Pfam accession), consecutive repeats collapsed.
* :class:`DomainArch` maps genes (``GeneTable`` order) to architecture codes, optionally shuffled (control),
  and answers per-node questions for :func:`ancestra.merge.decide_node`: unit architecture = most frequent
  architecture among the unit's annotated members; ``jaccard`` / ``same`` per candidate edge.
"""
from __future__ import annotations

import gzip
import os

import numpy as np

__all__ = ["parse_domtbl", "load_clans", "DomainArch", "jaccard_sets"]


def load_clans(path: str) -> dict:
    """Pfam accession (no version) -> clan accession, from Pfam-A.clans.tsv(.gz)."""
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1]:
                out[f[0]] = f[1]
    return out


def parse_domtbl(path: str, clans: dict | None = None, max_overlap: float = 0.5,
                 use_clans: bool = True) -> dict:
    """protein -> tuple of domain labels in sequence order (repeats collapsed)."""
    hits: dict = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            prot, acc = f[0], f[4].split(".")[0]
            ie = float(f[12])
            s, e = int(f[19]), int(f[20])      # envelope coordinates
            hits.setdefault(prot, []).append((ie, s, e, acc))
    out = {}
    for prot, hs in hits.items():
        hs.sort(key=lambda x: (x[0], x[1], x[3]))
        kept = []
        for ie, s, e, acc in hs:
            ok = True
            for _, s2, e2, _ in kept:
                ov = min(e, e2) - max(s, s2) + 1
                if ov > 0 and ov > max_overlap * min(e - s + 1, e2 - s2 + 1):
                    ok = False
                    break
            if ok:
                kept.append((ie, s, e, acc))
        kept.sort(key=lambda x: (x[1], x[2]))
        lab = [(clans.get(a, a) if (use_clans and clans) else a) for _, _, _, a in kept]
        col = [x for i, x in enumerate(lab) if i == 0 or x != lab[i - 1]]
        out[prot] = tuple(col)
    return out


def jaccard_sets(a: frozenset, b: frozenset) -> float:
    if not a or not b:
        return float("nan")
    return len(a & b) / len(a | b)


class DomainArch:
    """Gene -> architecture code (-1 = unannotated / not in the annotated scope).

    ``shuffle_seed``: control; the architectures of *annotated* genes are permuted among the annotated genes
    of the same species (``scope="species"``) or of the same species and family (``scope="family"``, needs
    ``family``). Annotation status is unchanged.
    """

    def __init__(self, names, species_of, arch_by_name: dict, shuffle_seed=None, scope: str = "species",
                 family=None):
        archs: dict = {}
        code = np.full(len(names), -1, dtype=np.int64)
        for i, n in enumerate(names):
            a = arch_by_name.get(n)
            if a:
                code[i] = archs.setdefault(a, len(archs))
        self.archs = [None] * len(archs)
        for a, k in archs.items():
            self.archs[k] = a
        self.sets = [frozenset(a) for a in self.archs]
        if shuffle_seed is not None:
            rng = np.random.default_rng(shuffle_seed)
            sp = np.asarray(species_of, dtype=np.int64)
            key = sp if scope == "species" else sp * (int(np.max(family)) + 2) + np.asarray(family, np.int64)
            ann = np.flatnonzero(code >= 0)
            o = ann[np.argsort(key[ann], kind="stable")]
            ks = key[o]
            cuts = np.flatnonzero(np.diff(ks)) + 1
            new = code.copy()
            for grp in np.split(o, cuts):
                new[grp] = code[grp][rng.permutation(len(grp))]
            code = new
        self.code = code
        self.shuffle_seed, self.scope = shuffle_seed, scope
        self._jcache: dict = {}
        self.n_annotated = int((code >= 0).sum())

    @classmethod
    def from_domtbl(cls, gt, domtbl: str, clans_path: str | None = None, **kw):
        clans = load_clans(clans_path) if clans_path else None
        return cls(gt.names, gt.species_of, parse_domtbl(domtbl, clans), **kw)

    # ---------------------------------------------------------------- node level
    def unit_codes(self, ctx) -> np.ndarray:
        """Most frequent annotated architecture code per local unit (-1 if none)."""
        n = len(ctx.child)
        d = getattr(ctx.members, "__defaults__", None)
        if d and len(d) == 2:
            g_sorted, b = d
            uid = np.repeat(np.arange(n, dtype=np.int64), np.diff(b))
        else:
            parts = [np.asarray(ctx.members(i), np.int64) for i in range(n)]
            g_sorted = np.concatenate(parts) if parts else np.zeros(0, np.int64)
            uid = np.repeat(np.arange(n, dtype=np.int64), [len(p) for p in parts])
        c = self.code[g_sorted]
        m = c >= 0
        out = np.full(n, -1, dtype=np.int64)
        if not m.any():
            return out
        u, c = uid[m], c[m]
        nk = len(self.archs) + 1
        key, cnt = np.unique(u * nk + c, return_counts=True)
        uu, cc = key // nk, key % nk
        o = np.lexsort((cc, -cnt, uu))          # per unit: highest count, then smallest code
        uu, cc = uu[o], cc[o]
        first = np.r_[True, uu[1:] != uu[:-1]]
        out[uu[first]] = cc[first]
        return out

    def jaccard_codes(self, ca: np.ndarray, cb: np.ndarray) -> np.ndarray:
        """Jaccard of the clan sets of architecture codes (NaN if either is -1)."""
        out = np.full(len(ca), np.nan)
        ok = (ca >= 0) & (cb >= 0)
        if not ok.any():
            return out
        lo, hi = np.minimum(ca[ok], cb[ok]), np.maximum(ca[ok], cb[ok])
        nk = len(self.archs) + 1
        key = lo * nk + hi
        uk, inv = np.unique(key, return_inverse=True)
        vals = np.empty(len(uk))
        for i, k in enumerate(uk.tolist()):
            v = self._jcache.get(k)
            if v is None:
                a, b = divmod(k, nk)
                v = 1.0 if a == b else jaccard_sets(self.sets[a], self.sets[b])
                self._jcache[k] = v
            vals[i] = v
        out[ok] = vals[inv]
        return out

    def edge_sim(self, ctx, sel: np.ndarray | None = None):
        """(jaccard, same) for the candidate edges ``sel`` (index array; default all)."""
        uc = self.unit_codes(ctx)
        ca = ctx.cand_a if sel is None else ctx.cand_a[sel]
        cb = ctx.cand_b if sel is None else ctx.cand_b[sel]
        xa, xb = uc[ca], uc[cb]
        return self.jaccard_codes(xa, xb), (xa >= 0) & (xa == xb)

    # ---------------------------------------------------------------- gene level (benchmarks)
    def gene_sim(self, ia: np.ndarray, ib: np.ndarray):
        ca, cb = self.code[ia], self.code[ib]
        return self.jaccard_codes(ca, cb), (ca >= 0) & (ca == cb)


def default_paths(root: str):
    d = os.path.join(root, "data/extra/l2_domains")
    return os.path.join(d, "Pfam-A.clans.tsv.gz")
