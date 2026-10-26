"""Gene-tree utilities for the S4 tree-method benchmark (docs/S4_trees.md).

Generic, label-free code (no RefOG / benchmark data is read here):

* :class:`UTree` -- unrooted tree (adjacency + leaf names) from Newick or from NJ edges
  (``ancestra.localtree.nj`` output); pruning to a leaf subset, joining two trees, splits.
* :func:`nrf`, :func:`quartet_distance` -- topology distances on the common leaf set.
* :class:`SpeciesTree` + :func:`root_reconcile` -- species-tree-aware rooting: the root edge
  minimising LCA-reconciliation duplications, then losses, then |left - right| tip depth
  (midpoint-like). The same rule is used for every arm and for the reference trees.
* :func:`root_overlap` -- the S5 truth rooting (min species-overlap duplications, then balance).
* :func:`pair_labels` -- species-overlap labelling: a cross-species pair is an ortholog iff the
  children of its LCA have disjoint species sets.
* :func:`paired_bootstrap` -- percentile CI for a difference of pooled ratios, resampling units.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations

import numpy as np

from .newick import parse_newick

__all__ = ["UTree", "nrf", "quartet_distance", "SpeciesTree", "root_reconcile", "root_overlap",
           "pair_labels", "paired_bootstrap"]


class UTree:
    """Unrooted tree: ``adj[u][v] = length``; leaves are degree-1 nodes with ``name[u]``."""

    def __init__(self):
        self.adj: dict = defaultdict(dict)
        self.name: dict = {}
        self._next = 0

    # -- construction ------------------------------------------------------
    def _new(self, name=None):
        u = self._next
        self._next += 1
        self.adj[u]  # create
        if name is not None:
            self.name[u] = name
        return u

    def _link(self, a, b, L):
        L = 0.0 if L is None else float(L)
        self.adj[a][b] = L
        self.adj[b][a] = L

    @classmethod
    def from_newick(cls, s: str, rename=None) -> "UTree":
        root = parse_newick(s.strip())
        t = cls()
        ids = {}
        for n in root.preorder():
            nm = None
            if not n.children:
                nm = rename(n.name) if rename else n.name
            ids[id(n)] = t._new(nm)
        for n in root.preorder():
            for c in n.children:
                t._link(ids[id(n)], ids[id(c)], c.length)
        t._suppress_deg2()
        return t

    @classmethod
    def from_edges(cls, edges, n_nodes, leaf_names) -> "UTree":
        """From ``localtree.nj`` output; leaves 0..m-1 get ``leaf_names``."""
        t = cls()
        for i in range(n_nodes):
            t._new(leaf_names[i] if i < len(leaf_names) else None)
        for a, b, L in edges:
            t._link(a, b, L)
        t._suppress_deg2()
        return t

    def _suppress_deg2(self):
        for u in [u for u in list(self.adj) if len(self.adj[u]) == 2 and u not in self.name]:
            (a, la), (b, lb) = self.adj[u].items()
            del self.adj[a][u], self.adj[b][u], self.adj[u]
            self._link(a, b, la + lb)

    # -- queries -----------------------------------------------------------
    @property
    def leaves(self):
        return [u for u in self.adj if len(self.adj[u]) <= 1 and u in self.name]

    def leaf_names(self):
        return {self.name[u] for u in self.leaves}

    def edges(self):
        return [(a, b) for a in self.adj for b in self.adj[a] if a < b]

    def copy(self) -> "UTree":
        t = UTree()
        t._next = self._next
        for u, nb in self.adj.items():
            t.adj[u] = dict(nb)
        t.name = dict(self.name)
        return t

    def prune(self, keep) -> "UTree":
        """Keep only leaves whose name is in ``keep``; degree-2 nodes are suppressed."""
        keep = set(keep)
        t = self.copy()
        changed = True
        while changed:
            changed = False
            for u in list(t.adj):
                if u not in t.adj:
                    continue
                if len(t.adj[u]) <= 1 and (u not in t.name or t.name[u] not in keep):
                    for v in list(t.adj[u]):
                        del t.adj[v][u]
                    del t.adj[u]
                    t.name.pop(u, None)
                    changed = True
        t._suppress_deg2()
        return t

    def splits(self):
        """Non-trivial bipartitions as frozensets of leaf names (side without the min name)."""
        names = sorted(self.leaf_names())
        if len(names) < 4:
            return set()
        ref = names[0]
        out = set()
        for a, b in self.edges():
            side = self._side(b, a)
            if 1 < len(side) < len(names) - 1:
                s = frozenset(side)
                if ref in s:
                    s = frozenset(names) - s
                out.add(s)
        return out

    def _side(self, u, par):
        res, stack = [], [(u, par)]
        while stack:
            x, p = stack.pop()
            if x in self.name and len(self.adj[x]) <= 1:
                res.append(self.name[x])
            stack.extend((y, x) for y in self.adj[x] if y != p)
        return res

    def rooted(self, a, b):
        """Children dict of the tree rooted on edge (a, b); the virtual root is -1."""
        ch = {-1: [a, b]}
        stack = [(a, b), (b, a)]
        while stack:
            u, p = stack.pop()
            ch[u] = [v for v in self.adj[u] if v != p]
            stack.extend((v, u) for v in ch[u])
        return ch

    def tip_depth(self, ch, u):
        best, stack = 0.0, [(u, 0.0)]
        while stack:
            x, d = stack.pop()
            kids = ch.get(x, [])
            if not kids:
                best = max(best, d)
            stack.extend((k, d + self.adj[x][k]) for k in kids)
        return best

    @staticmethod
    def join(t1: "UTree", r1: tuple, t2: "UTree", r2: tuple) -> "UTree":
        """Join two trees by an edge between new nodes inserted on edge r1 of t1 and r2 of t2."""
        t = UTree()
        maps = []
        for tt in (t1, t2):
            m = {u: t._new(tt.name.get(u)) for u in tt.adj}
            for a, b in tt.edges():
                t._link(m[a], m[b], tt.adj[a][b])
            maps.append(m)
        mids = []
        for (a, b), m, tt in ((r1, maps[0], t1), (r2, maps[1], t2)):
            L = tt.adj[a][b]
            x = t._new()
            del t.adj[m[a]][m[b]], t.adj[m[b]][m[a]]
            t._link(m[a], x, L / 2)
            t._link(x, m[b], L / 2)
            mids.append(x)
        t._link(mids[0], mids[1], 1.0)
        return t


# ----------------------------------------------------------------------------
# topology distances
# ----------------------------------------------------------------------------

def nrf(t1: UTree, t2: UTree):
    """Normalised unrooted RF on the common leaf set: (|A-B| + |B-A|) / (2(n-3)).
    Returns (nrf, rf, max_rf, n)."""
    common = t1.leaf_names() & t2.leaf_names()
    n = len(common)
    if n < 4:
        return float("nan"), 0, 0, n
    a, b = t1.prune(common).splits(), t2.prune(common).splits()
    rf = len(a - b) + len(b - a)
    mx = 2 * (n - 3)
    return rf / mx, rf, mx, n


def _quartet_topos(splits, quartets, names_idx):
    """For each quartet (i, j, k, l) -> 0 (ij|kl), 1 (ik|jl), 2 (il|jk), -1 unresolved."""
    out = np.full(len(quartets), -1, dtype=np.int8)
    masks = []
    for s in splits:
        m = 0
        for x in s:
            m |= 1 << names_idx[x]
        masks.append(m)
    for qi, (i, j, k, l) in enumerate(quartets):
        for m in masks:
            bi, bj, bk, bl = (m >> i) & 1, (m >> j) & 1, (m >> k) & 1, (m >> l) & 1
            if bi == bj and bk == bl and bi != bk:
                out[qi] = 0
                break
            if bi == bk and bj == bl and bi != bj:
                out[qi] = 1
                break
            if bi == bl and bj == bk and bi != bj:
                out[qi] = 2
                break
    return out


def quartet_distance(t1: UTree, t2: UTree, max_quartets: int = 3000, seed: int = 11):
    """Fraction of quartets (common leaves) whose resolved topologies differ (an unresolved
    quartet in either tree counts as different). All quartets if <= max_quartets, else a
    random sample (``np.random.default_rng(seed)``). Returns (qd, n_quartets)."""
    common = sorted(t1.leaf_names() & t2.leaf_names())
    n = len(common)
    if n < 4:
        return float("nan"), 0
    idx = {x: i for i, x in enumerate(common)}
    from math import comb
    if comb(n, 4) <= max_quartets:
        qs = list(combinations(range(n), 4))
    else:
        rng = np.random.default_rng(seed)
        qs = [tuple(sorted(rng.choice(n, 4, replace=False).tolist())) for _ in range(max_quartets)]
    a = _quartet_topos(t1.prune(common).splits(), qs, idx)
    b = _quartet_topos(t2.prune(common).splits(), qs, idx)
    diff = (a != b) | (a < 0) | (b < 0)
    return float(diff.mean()), len(qs)


# ----------------------------------------------------------------------------
# species tree, rooting, labelling
# ----------------------------------------------------------------------------

class SpeciesTree:
    """Rooted species tree; clades as bitmasks over the leaf species."""

    def __init__(self, newick: str, rename=None):
        root = parse_newick(newick.strip())
        leaves = root.leaves()
        self.species = [rename(l.name) if rename else l.name for l in leaves]
        self.bit = {s: 1 << i for i, s in enumerate(self.species)}
        self.clades, self.depth = [], []
        mask = {}
        for n in root.postorder():
            if not n.children:
                mask[id(n)] = self.bit[rename(n.name) if rename else n.name]
            else:
                m = 0
                for c in n.children:
                    m |= mask[id(c)]
                mask[id(n)] = m
        dep = {id(root): 0}
        for n in root.preorder():
            for c in n.children:
                dep[id(c)] = dep[id(n)] + 1
        for n in root.preorder():
            self.clades.append(mask[id(n)])
            self.depth.append(dep[id(n)])
        # smallest clade containing a mask: scan clades by decreasing depth
        self._order = sorted(range(len(self.clades)), key=lambda i: -self.depth[i])
        self._cache = {}

    def lca(self, m: int) -> int:
        r = self._cache.get(m)
        if r is None:
            for i in self._order:
                if m & self.clades[i] == m:
                    r = i
                    break
            self._cache[m] = r
        return r


def _postorder(ch, root):
    order, stack = [], [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack.extend(ch.get(u, []))
    return order[::-1]


def reconcile_cost(ch, root, leaf_species: dict, st: SpeciesTree):
    """(duplications, losses) of the rooted gene tree under LCA reconciliation."""
    mask, mp = {}, {}
    dups = losses = 0
    for u in _postorder(ch, root):
        kids = ch.get(u, [])
        if not kids:
            mask[u] = st.bit[leaf_species[u]]
        else:
            m = 0
            for k in kids:
                m |= mask[k]
            mask[u] = m
        mp[u] = st.lca(mask[u])
        if kids:
            dup = any(mp[k] == mp[u] for k in kids)
            dups += dup
            du = st.depth[mp[u]]
            for k in kids:
                losses += st.depth[mp[k]] - du - (0 if dup else 1)
    return dups, losses


def _species_sets(ch, root, leaf_species):
    ss = {}
    for u in _postorder(ch, root):
        kids = ch.get(u, [])
        ss[u] = {leaf_species[u]} if not kids else set().union(*(ss[k] for k in kids))
    return ss


def _overlap_dups(ch, root, leaf_species):
    ss = _species_sets(ch, root, leaf_species)
    d = 0
    for u, kids in ch.items():
        if len(kids) >= 2:
            seen = set()
            for k in kids:
                if seen & ss[k]:
                    d += 1
                    break
                seen |= ss[k]
    return d


def _leaf_species(t: UTree, species_of_name):
    return {u: species_of_name(t.name[u]) for u in t.leaves}


def _balance(t, ch, a, b):
    """|max tip depth below a - below b| for the root on edge (a, b) (S5 tie-break)."""
    return abs(t.tip_depth(ch, a) - t.tip_depth(ch, b))


def root_reconcile(t: UTree, species_of_name, st: SpeciesTree):
    """Root on the edge minimising (LCA duplications, losses, tip-depth imbalance).
    Returns (children dict, root=-1, info)."""
    ls = _leaf_species(t, species_of_name)
    best = None
    for a, b in t.edges():
        ch = t.rooted(a, b)
        d, l = reconcile_cost(ch, -1, ls, st)
        key = (d, l, _balance(t, ch, a, b), a, b)
        if best is None or key < best[0]:
            best = (key, ch)
    (d, l, bal, a, b), ch = best
    return ch, -1, dict(dups=d, losses=l, root_edge=(a, b))


def root_overlap(t: UTree, species_of_name):
    """S5 truth rooting: min species-overlap duplications, ties -> most balanced."""
    ls = _leaf_species(t, species_of_name)
    best = None
    for a, b in t.edges():
        ch = t.rooted(a, b)
        key = (_overlap_dups(ch, -1, ls), _balance(t, ch, a, b), a, b)
        if best is None or key < best[0]:
            best = (key, ch)
    (d, bal, a, b), ch = best
    return ch, -1, dict(dups=d, root_edge=(a, b))


def pair_labels(t: UTree, ch, root, species_of_name, keep=None) -> dict:
    """{(name_a, name_b) sorted: 1 ortholog / 0 paralog} for all cross-species leaf pairs
    (restricted to ``keep`` names if given), by species overlap at their LCA."""
    ls = _leaf_species(t, species_of_name)
    ss = _species_sets(ch, root, ls)
    out = {}
    under = {}
    for u in _postorder(ch, root):
        kids = ch.get(u, [])
        if not kids:
            nm = t.name.get(u)
            under[u] = [u] if (nm is not None and (keep is None or nm in keep)) else []
            continue
        seen, dup = set(), False
        for k in kids:
            if seen & ss[k]:
                dup = True
            seen |= ss[k]
        lab = 0 if dup else 1
        for i in range(len(kids)):
            for j in range(i + 1, len(kids)):
                for x in under[kids[i]]:
                    for y in under[kids[j]]:
                        if ls[x] != ls[y]:
                            p = tuple(sorted((t.name[x], t.name[y])))
                            out[p] = lab
        under[u] = [x for k in kids for x in under[k]]
        for k in kids:
            del under[k]
    return out


# ----------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------

def paired_bootstrap(num_a, den_a, num_b, den_b, B: int = 2000, seed: int = 11):
    """Δ = sum(num_b)/sum(den_b) - sum(num_a)/sum(den_a) with a percentile 95% CI from
    resampling units (rows) with replacement. Returns (delta, lo, hi, p_le0)."""
    num_a, den_a, num_b, den_b = (np.asarray(x, dtype=np.float64) for x in (num_a, den_a, num_b, den_b))
    n = len(num_a)
    d0 = num_b.sum() / den_b.sum() - num_a.sum() / den_a.sum()
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(B, n))
    with np.errstate(invalid="ignore", divide="ignore"):
        db = num_b[idx].sum(1) / den_b[idx].sum(1) - num_a[idx].sum(1) / den_a[idx].sum(1)
    db = db[np.isfinite(db)]
    lo, hi = np.percentile(db, [2.5, 97.5])
    return float(d0), float(lo), float(hi), float((db <= 0).mean())
