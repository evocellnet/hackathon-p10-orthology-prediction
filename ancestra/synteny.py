"""S5 -- gene-order (synteny) context evidence.

Two levels, one kernel.

Leaf level (species pair benchmark, ``scripts/S5_pair_benchmark.py``)
    For genes a (species A) and b (species B), the *context score* is the number
    of anchor pairs (x, y) with x within ±w genes of a and y within ±w genes of b
    (same contig, ranks from :mod:`ancestra.genomes`). Anchors are genome-wide
    unambiguous reciprocal best hits (:func:`rbh_anchors`); anchors whose genes
    lie in the family of a or b are ignored (no circularity, no tandem copies).

Ancestral level (merge integration, DESIGN §3 S5)
    At node X with children C1..Ck the merge "units" are the ancestral genes of
    the children. Two units of the same child are *adjacent* if, in at least one
    species, a member gene of one lies within ±w genes of a member gene of the
    other (:meth:`AncestralAdjacency.unit_adjacency`; sparse, built from
    per-species neighbour lists precomputed once). For a candidate cross-child
    pair (A, B) and a set of anchor pairs (P, Q) (pass 1: unambiguous 1:1
    mutual-best matches between the children, :func:`select_anchors`), the
    support is

        support(A, B) = #{(P, Q) in anchors : P adj A and Q adj B,
                           fam(P) != fam(A), fam(Q) != fam(B)}

    (:func:`context_support`). At a leaf pair this reduces exactly to the
    species-pair score above. Three passes: (1) anchors from sequence only,
    (2) support for every candidate, (3) the merge decides with
    ``score + bonus(support)`` (:class:`SyntenyContext`). Anchor selection never
    uses synteny, so the evidence is not circular.

All index arrays are int64 numpy; the counting kernel is numba-compiled and
O(sum over candidates of deg(A) * (anchors per unit) * log deg(B)).
"""
from __future__ import annotations

import os

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*a, **k):
        def deco(f):
            return f
        return deco if not (a and callable(a[0])) else a[0]

from .genomes import neighbour_csr, shuffle_order

__all__ = ["gene_positions", "rbh_anchors", "pair_context_scores", "AncestralAdjacency",
           "select_anchors", "context_support", "SyntenyContext", "csr_from_pairs"]


# ----------------------------------------------------------------------------
# positions
# ----------------------------------------------------------------------------

def gene_positions(names, species_of, species, table_dir: str):
    """Integer gene IDs -> (contig_code, order) from ``table_dir/<species>.parquet``.

    ``names[g]`` protein ID, ``species_of[g]`` species index into ``species``.
    Contig codes are globally unique (species-specific); genes without a
    table row get contig_code -1. Returns (contig_code, order, coverage dict).
    """
    import pandas as pd
    n = len(names)
    cc = np.full(n, -1, np.int64)
    oo = np.zeros(n, np.int64)
    idx = {nm: i for i, nm in enumerate(names)}
    cov = {}
    base = 0
    for si, sp in enumerate(species):
        p = os.path.join(table_dir, f"{sp}.parquet")
        if not os.path.exists(p):
            cov[sp] = 0.0
            continue
        df = pd.read_parquet(p, columns=["protein_id", "contig", "order"])
        gi = df["protein_id"].map(idx)
        ok = gi.notna().values
        gi = gi[ok].astype(np.int64).values
        codes, uniq = pd.factorize(df["contig"][ok])
        cc[gi] = codes + base
        oo[gi] = df["order"].values[ok]
        base += len(uniq)
        nsp = int((np.asarray(species_of) == si).sum())
        cov[sp] = round(len(gi) / max(nsp, 1), 6)
    return cc, oo, cov


# ----------------------------------------------------------------------------
# anchors
# ----------------------------------------------------------------------------

def rbh_anchors(q, t, bits, amb: float = 0.9):
    """Unambiguous reciprocal best hits from directed hit lists.

    ``q, t, bits``: hits in both directions (query -> target, any species pair;
    pass the union of A->B and B->A searches, IDs global). (a, b) is an anchor
    iff b is a's best hit, a is b's best hit, and on both sides the best hit to
    any *other* gene scores < ``amb`` * best (unambiguous; ``amb=1`` = plain
    RBH). Returns (a, b) arrays with a < b.
    """
    q = np.asarray(q, np.int64); t = np.asarray(t, np.int64); bits = np.asarray(bits, float)
    # symmetric best score per unordered pair
    lo, hi = np.minimum(q, t), np.maximum(q, t)
    key = lo * (int(max(hi.max(initial=0), lo.max(initial=0))) + 1) + hi
    o = np.lexsort((-bits, key))
    first = np.r_[True, key[o][1:] != key[o][:-1]]
    sel = o[first]
    lo, hi, s = lo[sel], hi[sel], bits[sel]
    # per gene: best and second best partner over both endpoints
    g = np.concatenate([lo, hi]); p = np.concatenate([hi, lo]); ss = np.concatenate([s, s])
    o = np.lexsort((-ss, g))
    g, p, ss = g[o], p[o], ss[o]
    start = np.r_[True, g[1:] != g[:-1]]
    gi = np.flatnonzero(start)
    best_g, best_p, best_s = g[gi], p[gi], ss[gi]
    second = np.zeros(len(gi))
    has2 = np.r_[gi[1:], len(g)] - gi >= 2
    second[has2] = ss[gi[has2] + 1]
    ok = second < amb * best_s
    bp = dict(zip(best_g.tolist(), zip(best_p.tolist(), ok.tolist())))
    out = [(a, b) for a, (b, oka) in bp.items()
           if a < b and oka and bp.get(b, (-1, False))[0] == a and bp[b][1]]
    if not out:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    arr = np.array(out, np.int64)
    return arr[:, 0], arr[:, 1]


# ----------------------------------------------------------------------------
# kernel
# ----------------------------------------------------------------------------

def csr_from_pairs(n: int, a, b, symmetric: bool = True):
    """CSR (indptr, indices, sorted, deduplicated) of the pairs a->b (and b->a)."""
    a = np.asarray(a, np.int64); b = np.asarray(b, np.int64)
    if symmetric:
        a, b = np.concatenate([a, b]), np.concatenate([b, a])
    keep = a != b
    a, b = a[keep], b[keep]
    if len(a):
        key = a * n + b
        key = np.unique(key)
        a, b = key // n, key % n
    indptr = np.zeros(n + 1, np.int64)
    np.add.at(indptr, a + 1, 1)
    return np.cumsum(indptr), b.astype(np.int64)


@njit(cache=True)
def _support_kernel(ca, cb, adj_ptr, adj_idx, anc_ptr, anc_idx, fam):
    n = len(ca)
    out = np.zeros(n, np.int64)
    for k in range(n):
        a = ca[k]; b = cb[k]
        fa = fam[a]; fb = fam[b]
        lo_b = adj_ptr[b]; hi_b = adj_ptr[b + 1]
        c = 0
        for ii in range(adj_ptr[a], adj_ptr[a + 1]):
            p = adj_idx[ii]
            if fam[p] == fa or fam[p] == fb:
                continue
            for jj in range(anc_ptr[p], anc_ptr[p + 1]):
                q = anc_idx[jj]
                if fam[q] == fb or fam[q] == fa or q == b:
                    continue
                # binary search q in adj row of b
                lo = lo_b; hi = hi_b
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if adj_idx[mid] < q:
                        lo = mid + 1
                    else:
                        hi = mid
                if lo < hi_b and adj_idx[lo] == q:
                    c += 1
        out[k] = c
    return out


def context_support(candidate_pairs, anchors, node, fam=None):
    """Per-candidate synteny support at a node.

    ``candidate_pairs``: (ca, cb) local unit indices; ``anchors``: (pa, pb) local
    unit indices (1:1 or not; each anchor counted once per orientation that
    matches); ``node``: an object with ``adj_ptr``/``adj_idx`` (unit adjacency
    CSR, e.g. :class:`NodeAdjacency` from :meth:`AncestralAdjacency.at_node`)
    and optionally ``fam`` (unit -> family id; anchors in the family of either
    candidate unit are ignored). Returns int64 counts.
    """
    ca = np.ascontiguousarray(candidate_pairs[0], np.int64)
    cb = np.ascontiguousarray(candidate_pairs[1], np.int64)
    n = len(node.adj_ptr) - 1
    if fam is None:
        fam = getattr(node, "fam", None)
    if fam is None:
        fam = np.arange(n, dtype=np.int64)   # every unit its own family: only self excluded
    anc_ptr, anc_idx = csr_from_pairs(n, anchors[0], anchors[1], symmetric=True)
    return _support_kernel(ca, cb, node.adj_ptr, node.adj_idx, anc_ptr, anc_idx,
                           np.ascontiguousarray(fam, np.int64))


# ----------------------------------------------------------------------------
# species-pair scores (leaf level)
# ----------------------------------------------------------------------------

class _Node:
    def __init__(self, adj_ptr, adj_idx, fam=None):
        self.adj_ptr, self.adj_idx, self.fam = adj_ptr, adj_idx, fam


def pair_context_scores(cand_a, cand_b, anchor_a, anchor_b, contig_code, order, w, fam=None):
    """Leaf-level context score for candidate gene pairs (global gene IDs).

    Equivalent to :func:`context_support` with units = genes and adjacency =
    ±w neighbourhoods. ``fam``: gene -> family (anchors inside the family of
    a or b excluded).
    """
    ptr, idx = neighbour_csr(contig_code, order, w)
    return context_support((cand_a, cand_b), (anchor_a, anchor_b), _Node(ptr, idx), fam)


# ----------------------------------------------------------------------------
# ancestral adjacency (merge integration)
# ----------------------------------------------------------------------------

class NodeAdjacency:
    """Unit adjacency at one node: CSR over local units + unit families."""

    def __init__(self, adj_ptr, adj_idx, fam, n_species_support=None):
        self.adj_ptr, self.adj_idx, self.fam = adj_ptr, adj_idx, fam
        self.n_species_support = n_species_support

    @property
    def n_units(self):
        return len(self.adj_ptr) - 1

    @property
    def n_edges(self):
        return len(self.adj_idx) // 2


class AncestralAdjacency:
    """Per-species gene neighbour lists (precomputed once) -> unit adjacency per node.

    ``contig_code``/``order``: per global gene ID (:func:`gene_positions`);
    ``w``: window in genes (gap-tolerant: ranks count only table genes);
    ``family``: gene -> family id (S0), used for the anchor exclusion.
    """

    def __init__(self, contig_code, order, w: int = 5, family=None, shuffle_seed=None):
        contig_code = np.asarray(contig_code, np.int64)
        order = np.asarray(order, np.int64)
        if shuffle_seed is not None:      # negative control: permuted gene order
            contig_code, order = shuffle_order(contig_code, order, shuffle_seed)
        self.w = w
        self.N = len(order)
        self.gptr, self.gidx = neighbour_csr(contig_code, order, w)
        rows = np.repeat(np.arange(self.N, dtype=np.int64), np.diff(self.gptr))
        m = rows < self.gidx                     # each undirected leaf pair once
        self.pg, self.ph = rows[m], self.gidx[m]
        self.family = (np.asarray(family, np.int64) if family is not None
                       else np.zeros(self.N, np.int64) - 1)

    def unit_adjacency(self, gene_unit, n_units: int) -> NodeAdjacency:
        """Adjacency among units; ``gene_unit[g]`` = local unit of gene g or -1.

        Units u != v are adjacent iff some leaf neighbour pair (g, h) has
        gene_unit[g] = u, gene_unit[h] = v. Unit family = family of its genes
        (merges never cross families; if they do, the first gene's family).
        """
        gene_unit = np.asarray(gene_unit, np.int64)
        u, v = gene_unit[self.pg], gene_unit[self.ph]
        k = (u >= 0) & (v >= 0) & (u != v)
        ptr, idx = csr_from_pairs(n_units, u[k], v[k], symmetric=True)
        fam = np.full(n_units, -1, np.int64)
        g = np.flatnonzero(gene_unit >= 0)
        fam[gene_unit[g]] = self.family[g]
        # units without family info get unique negative ids (never excluded together)
        nf = fam < 0
        fam[nf] = -2 - np.flatnonzero(nf)
        return NodeAdjacency(ptr, idx, fam)

    def at_node(self, units, gene_anc, genes) -> NodeAdjacency:
        """Adjacency for the merge state of :func:`ancestra.merge.run_merge`.

        ``units``: global ancestral IDs of the children's units (local index =
        position), ``gene_anc``: gene -> current ancestral ID (before the node),
        ``genes``: the clade's gene IDs.
        """
        units = np.asarray(units, np.int64)
        pos = {int(x): i for i, x in enumerate(units)}
        gu = np.full(self.N, -1, np.int64)
        ga = np.asarray(gene_anc, np.int64)[genes]
        # vectorised lookup via sorting
        o = np.argsort(units)
        j = np.searchsorted(units[o], ga)
        j = np.clip(j, 0, len(units) - 1)
        ok = units[o][j] == ga
        gu[np.asarray(genes)[ok]] = o[j[ok]]
        del pos
        return self.unit_adjacency(gu, len(units))


def select_anchors(ca, cb, cs, child, amb: float = 0.9):
    """Pass 1: unambiguous 1:1 mutual-best cross-child pairs (sequence only).

    (A, B) with child(A)=i != child(B)=j is an anchor iff B is A's best partner
    in child j, A is B's best partner in child i, and the second-best partner
    score in that child is < ``amb`` * best on both sides.
    """
    ca = np.asarray(ca, np.int64); cb = np.asarray(cb, np.int64); cs = np.asarray(cs, float)
    child = np.asarray(child, np.int64)
    if len(ca) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    k = int(child.max()) + 1 if len(child) else 1
    # directed: (x, partner child) groups
    x = np.concatenate([ca, cb]); y = np.concatenate([cb, ca]); s = np.concatenate([cs, cs])
    grp = x * k + child[y]
    o = np.lexsort((-s, grp))
    grp, x, y, s = grp[o], x[o], y[o], s[o]
    st = np.flatnonzero(np.r_[True, grp[1:] != grp[:-1]])
    en = np.r_[st[1:], len(grp)]
    second = np.where(en - st >= 2, s[np.minimum(st + 1, len(s) - 1)], 0.0)
    ok = second < amb * s[st]
    best = {}
    for gg, yy, okk in zip(grp[st].tolist(), y[st].tolist(), ok.tolist()):
        best[gg] = (yy, okk)
    pa, pb = [], []
    for gg, (yy, okk) in best.items():
        xx = gg // k
        if not okk or xx > yy:
            continue
        back = best.get(yy * k + int(child[xx]))
        if back is not None and back[0] == xx and back[1]:
            pa.append(xx); pb.append(yy)
    return np.array(pa, np.int64), np.array(pb, np.int64)


def _unit_families(ctx, family):
    """Family id per local unit of a NodeContext (merges never cross families)."""
    n = len(ctx.child)
    fam = np.full(n, -1, np.int64)
    d = getattr(ctx.members, "__defaults__", None)
    if d and len(d) == 2:                      # fast path: members = lambda i, _g, _b
        g_sorted, b = d
        nz = np.diff(b) > 0
        fam[nz] = family[g_sorted[b[:-1][nz]]]
    elif ctx.members is not None:
        for i in range(n):
            m = ctx.members(i)
            if len(m):
                fam[i] = family[m[0]]
    neg = fam < 0
    fam[neg] = -2 - np.flatnonzero(neg)
    return fam


class SyntenyContext:
    """Context-evidence function for ``merge.run_merge(context=...)`` (DESIGN S5).

    Called once per node with a ``merge.NodeContext`` (``ctx``); returns one additive
    bonus per cross-child candidate. Three passes: anchors (pass 1, sequence only) ->
    support per candidate (:func:`context_support`) -> the merge rules run on
    ``score + bonus``.

    Adjacency source: if ``contig_code``/``order`` are given, unit adjacency is built
    here from per-gene neighbour lists (window ``w``; ``shuffle_seed`` = shuffled-order
    control); otherwise ``ctx.adjacency`` from the core (``run_merge(...,
    gene_order=GeneOrder, adjacency_w=w)``) is used, optionally requiring
    ``min_species`` supporting species per adjacency.
    Anchors: ``anchor="core"`` = ``ctx.anchors`` (strict mutual best, 1:1);
    ``anchor="amb"`` = :func:`select_anchors` (1:1 mutual best whose runner-up in the
    other child scores < ``amb`` x best; more conservative).
    Anchors in the family of either candidate unit are ignored (``family``: gene ->
    family id; None = only the candidate pair itself is excluded).

    ``mode``: ``off`` (zeros) | ``tiebreak`` (eps * support) | ``additive``
    (lam * min(support, cap)/cap) | ``multiplicative`` (score * lam * min(support,
    cap)/cap; lam=0.1: a supported candidate overtakes an unsupported one only within
    ~10% of its sequence score -- the rule that helped on OrthoBench dev) |
    ``anchors_only`` (lam if support >= min_support).
    """

    def __init__(self, contig_code=None, order=None, family=None, w: int = 5,
                 mode: str = "multiplicative", lam: float = 0.1, cap: int = 3, eps: float = 1e-6,
                 anchor: str = "amb", amb: float = 0.9, min_support: int = 1,
                 min_species: int = 1, shuffle_seed=None, near: float = 0.9, margin: int = 1,
                 delta: float = 0.01):
        self.near, self.margin, self.delta = near, margin, delta
        self.adj = (AncestralAdjacency(contig_code, order, w, family, shuffle_seed)
                    if contig_code is not None else None)
        self.family = np.asarray(family, np.int64) if family is not None else None
        self.mode, self.lam, self.cap, self.eps = mode, lam, cap, eps
        self.anchor, self.amb, self.min_support, self.min_species = anchor, amb, min_support, min_species
        self.stats = []

    def node_adjacency(self, ctx) -> NodeAdjacency:
        n = len(ctx.child)
        fam = (_unit_families(ctx, self.family) if self.family is not None
               else np.arange(n, dtype=np.int64))
        if self.adj is not None:
            gu = np.full(self.adj.N, -1, np.int64)
            d = getattr(ctx.members, "__defaults__", None)
            if d and len(d) == 2:
                g_sorted, b = d
                gu[g_sorted] = np.repeat(np.arange(n, dtype=np.int64), np.diff(b))
            else:
                for i in range(n):
                    gu[ctx.members(i)] = i
            nadj = self.adj.unit_adjacency(gu, n)
            nadj.fam = fam
            return nadj
        if ctx.adjacency is None:
            raise ValueError("no gene order: pass contig_code/order or run_merge(gene_order=...)")
        ua, ub, nsp, _ = ctx.adjacency
        k = nsp >= self.min_species
        ptr, idx = csr_from_pairs(n, ua[k], ub[k], symmetric=True)
        return NodeAdjacency(ptr, idx, fam)

    def anchors_of(self, ctx):
        if self.anchor == "core":
            return ctx.cand_a[ctx.anchors], ctx.cand_b[ctx.anchors]
        return select_anchors(ctx.cand_a, ctx.cand_b, ctx.cand_score, ctx.child, self.amb)

    def support(self, ctx):
        nadj = self.node_adjacency(ctx)
        pa, pb = self.anchors_of(ctx)
        sup = context_support((ctx.cand_a, ctx.cand_b), (pa, pb), nadj)
        self.stats.append(dict(node=ctx.node, n_cand=int(len(ctx.cand_a)), n_anchors=int(len(pa)),
                               n_adj_edges=int(nadj.n_edges), n_supported=int((sup > 0).sum())))
        return sup

    def __call__(self, ctx):
        cs = np.asarray(ctx.cand_score, float)
        if self.mode == "off":
            return np.zeros(len(cs))
        sup = self.support(ctx)
        if self.mode == "tiebreak":
            return self.eps * sup
        if self.mode == "multiplicative":
            return cs * self.lam * np.minimum(sup, self.cap) / self.cap
        if self.mode == "anchors_only":
            return np.where(sup >= self.min_support, self.lam, 0.0)
        if self.mode == "override":
            return override_bonus(ctx.cand_a, ctx.cand_b, cs, ctx.child, sup,
                                  self.near, self.margin, self.delta)
        return self.lam * np.minimum(sup, self.cap) / self.cap


def override_bonus(ca, cb, cs, child, sup, near=0.9, margin=1, delta=0.01):
    """Minimal bonus that lets synteny change only the *partner choice*.

    For endpoint x of candidate k=(a,b) and the other child c, let t = x's
    sequence-best candidate towards c. k *overrides* if, for a or for b,
    ``cs_k >= near * cs_t`` and ``sup_k - sup_t >= margin``. An overriding
    candidate gets the bonus that lifts it to (1 + delta) x the larger of the two
    endpoints' best scores (so it becomes mutual best); all others get 0, so the
    merge's in-paralog / repair thresholds are untouched elsewhere.
    """
    ca = np.asarray(ca, np.int64); cb = np.asarray(cb, np.int64)
    cs = np.asarray(cs, float); sup = np.asarray(sup, np.int64); child = np.asarray(child, np.int64)
    E = len(ca)
    if E == 0:
        return np.zeros(0)
    k = int(child.max()) + 1
    x = np.concatenate([ca, cb]); oc = np.concatenate([child[cb], child[ca]])
    e = np.concatenate([np.arange(E), np.arange(E)])
    grp = x * k + oc
    o = np.lexsort((-cs[e], grp))
    g_s = grp[o]
    st = np.r_[True, g_s[1:] != g_s[:-1]]
    first = np.flatnonzero(st)
    gid = np.cumsum(st) - 1
    top_e = e[o][first]                    # sequence-best edge per (unit, other child)
    best = cs[top_e][gid]; tsup = sup[top_e][gid]
    ee = e[o]
    ok = (cs[ee] >= near * best) & (sup[ee] - tsup >= margin)
    hit = np.zeros(E, bool)
    hit[ee[ok]] = True
    bestpos = np.zeros(2 * E)
    bestpos[o] = best
    bmax = np.maximum(bestpos[:E], bestpos[E:])
    return np.where(hit, np.maximum(0.0, (1 + delta) * bmax - cs), 0.0)
