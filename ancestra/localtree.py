"""Hybrid local-tree mode for ambiguous clusters (DESIGN §1(b); r1_external_pi must #1).

At a node X, after the score rules R1-R3, *ambiguous clusters* are re-decided with a
small gene tree:

1. **Detection.** Every unit u has best_all(u) = its best edge score (cross- or
   same-child). An edge (u, v) between two different rule groups is an *ambiguity edge*
   if S(u,v) >= tau * min(best_all(u), best_all(v)), i.e. it is (nearly) the best
   relation of one endpoint, yet the rules put u and v apart. Rule groups connected by
   ambiguity edges form a cluster; a cluster with >= 2 groups and <= ``max_units`` units
   goes to tree mode (larger ones keep the rule result and are counted).
2. **Tree.** Neighbour joining on unit distances d = -ln S, with S the linked (avg/max)
   member score; missing pairs get max(d) + 1. One *outgroup leaf* is added when
   possible: the sibling-clade gene o with the highest summed score to the cluster's
   units (``out_edges``), with d(u, o) = -ln S(u, o).
3. **Rooting.** On the outgroup's branch; otherwise (root node, or no outgroup scores)
   at the edge that minimises the number of duplications under species overlap
   (ties: longest edge, i.e. a midpoint-like choice).
4. **Labelling (species overlap on X's children).** A node whose two subtrees share a
   child of X is a duplication at or before X -> split into its subtrees (recursively).
   Any other subtree -- a speciation at X, or a duplication confined to one child
   (in-paralogs, duplication after X) -- is one ancestral gene of X.
NJ is used because it is rate-robust (unequal branch lengths are allowed), which is the
failure of the score rules (fast young duplicates, fast outgroup lineages such as
Ciona); it needs no MSA, so cost is O(m^3) per cluster with m <= max_units.
"""
from __future__ import annotations

import numpy as np

__all__ = ["nj", "resolve_ambiguous"]


def nj(D: np.ndarray):
    """Neighbour joining. Returns (edges list of (a, b, length), n_nodes); leaves 0..m-1."""
    m = D.shape[0]
    if m == 1:
        return [], 1
    if m == 2:
        return [(0, 1, max(D[0, 1], 0.0))], 2
    D = D.astype(np.float64).copy()
    active = list(range(m))
    node_id = list(range(m))
    nxt = m
    edges = []
    while len(active) > 2:
        k = len(active)
        sub = D[np.ix_(active, active)]
        r = sub.sum(axis=1)
        Q = (k - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        if i > j:
            i, j = j, i
        di = 0.5 * sub[i, j] + (r[i] - r[j]) / (2 * (k - 2))
        dj = sub[i, j] - di
        a, b = active[i], active[j]
        edges.append((node_id[a], nxt, max(di, 0.0)))
        edges.append((node_id[b], nxt, max(dj, 0.0)))
        # new row replaces a, b removed
        newd = 0.5 * (D[a, :] + D[b, :] - D[a, b])
        D[a, :] = newd
        D[:, a] = newd
        D[a, a] = 0.0
        node_id[a] = nxt
        nxt += 1
        active.pop(j)
    a, b = active
    edges.append((node_id[a], node_id[b], max(D[a, b], 0.0)))
    return edges, nxt


def _adjacency(edges, n_nodes):
    adj = [[] for _ in range(n_nodes)]
    for a, b, l in edges:
        adj[a].append((b, l))
        adj[b].append((a, l))
    return adj


def _rooted_children(adj, root_a, root_b):
    """Root on edge (root_a, root_b): returns children dict for a virtual root -1."""
    ch = {-1: [root_a, root_b]}
    stack = [(root_a, root_b), (root_b, root_a)]
    while stack:
        v, parent = stack.pop()
        kids = [w for w, _ in adj[v] if w != parent]
        ch[v] = kids
        for w in kids:
            stack.append((w, v))
    return ch


def _child_sets(ch, leaf_child, root=-1):
    """Post-order: set of X-children (frozenset) per node; leaves without a child label
    (outgroup, value -1) contribute nothing."""
    sets = {}
    order = []
    stack = [root]
    while stack:
        v = stack.pop()
        order.append(v)
        stack.extend(ch.get(v, []))
    for v in reversed(order):
        kids = ch.get(v, [])
        if not kids:
            c = leaf_child.get(v, -1)
            sets[v] = frozenset() if c < 0 else frozenset([c])
        else:
            s = frozenset()
            for w in kids:
                s = s | sets[w]
            sets[v] = s
    return sets


def _is_xdup(ch, sets, v):
    kids = ch.get(v, [])
    if len(kids) < 2:
        return False
    if len(sets[v]) < 2:          # duplication confined to one child of X: after X
        return False
    seen = set()
    for w in kids:
        if sets[w] & seen:
            return True
        seen |= sets[w]
    return False


def _count_dups(ch, sets):
    return sum(1 for v in ch if _is_xdup(ch, sets, v))


def _partition(ch, sets, root, leaves_of_interest):
    """Split at X-duplications; return list of leaf lists."""
    out = []
    stack = [root]
    while stack:
        v = stack.pop()
        kids = ch.get(v, [])
        if not kids:
            if v in leaves_of_interest:
                out.append([v])
            continue
        if _is_xdup(ch, sets, v) or _has_xdup_below(ch, sets, v):
            stack.extend(kids)
        else:
            lv = _leaves(ch, v)
            lv = [x for x in lv if x in leaves_of_interest]
            if lv:
                out.append(lv)
    return out


def _leaves(ch, v):
    res, stack = [], [v]
    while stack:
        x = stack.pop()
        kids = ch.get(x, [])
        if not kids:
            res.append(x)
        else:
            stack.extend(kids)
    return res


def _has_xdup_below(ch, sets, v):
    stack = list(ch.get(v, []))
    while stack:
        x = stack.pop()
        if _is_xdup(ch, sets, x):
            return True
        stack.extend(ch.get(x, []))
    return False


def tree_groups(D: np.ndarray, leaf_child: list, outgroup: int | None):
    """Groups (lists of leaf indices, outgroup excluded) of a cluster's distance matrix."""
    m = D.shape[0]
    edges, n_nodes = nj(D)
    if not edges:
        return [[0]]
    adj = _adjacency(edges, n_nodes)
    lc = {i: leaf_child[i] for i in range(m)}
    interest = set(range(m)) - ({outgroup} if outgroup is not None else set())
    if outgroup is not None:
        nb = adj[outgroup][0][0]
        ch = _rooted_children(adj, outgroup, nb)
        sets = _child_sets(ch, lc)
        return _partition(ch, sets, nb, interest)
    best = None
    for a, b, l in edges:
        ch = _rooted_children(adj, a, b)
        sets = _child_sets(ch, lc)
        key = (_count_dups(ch, sets), -l)
        if best is None or key < best[0]:
            best = (key, ch, sets)
    _, ch, sets = best
    return _partition(ch, sets, -1, interest)


def resolve_ambiguous(n, child, lu, lv, score, labels, out_edges=None, tau: float = 0.9,
                      max_units: int = 60, stats=None, require_outgroup: bool = False,
                      conservative: bool = False):
    """Re-decide ambiguous clusters with local trees; returns new labels (see module doc)."""
    stats = stats if stats is not None else {}
    labels = labels.copy()
    s = np.asarray(score, dtype=np.float64)
    best_all = np.full(n, -np.inf)
    np.maximum.at(best_all, lu, s)
    np.maximum.at(best_all, lv, s)
    amb = (labels[lu] != labels[lv]) & (s >= tau * np.minimum(best_all[lu], best_all[lv]))
    stats.update(tree_clusters=0, tree_units=0, tree_too_big=0, tree_changed=0)
    if not amb.any():
        return labels
    # union-find over labels via ambiguity edges
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(labels[lu[amb]].tolist(), labels[lv[amb]].tolist()):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    clusters: dict = {}
    for lab in list(parent):
        clusters.setdefault(find(lab), set()).add(lab)
    unit_by_label: dict = {}
    for i, lab in enumerate(labels.tolist()):
        if lab in parent:
            unit_by_label.setdefault(lab, []).append(i)
    # edge lookup for distance matrices
    key = np.minimum(lu, lv).astype(np.int64) * n + np.maximum(lu, lv)
    order = np.argsort(key)
    key_s, s_s = key[order], s[order]
    out_adj = None
    if out_edges is not None and len(out_edges[0]):
        out_adj = {}
        for u, o, sc in zip(out_edges[0].tolist(), out_edges[1].tolist(), out_edges[2].tolist()):
            out_adj.setdefault(u, {})[o] = sc
    next_label = int(labels.max()) + 1
    for root, labs in clusters.items():
        if len(labs) < 2:
            continue
        units = sorted(u for lab in labs for u in unit_by_label[lab])
        if len(units) > max_units:
            stats["tree_too_big"] += 1
            continue
        stats["tree_clusters"] += 1
        stats["tree_units"] += len(units)
        m = len(units)
        ua = np.array(units, dtype=np.int64)
        I, J = np.triu_indices(m, 1)
        kk = np.minimum(ua[I], ua[J]) * n + np.maximum(ua[I], ua[J])
        pos = np.searchsorted(key_s, kk)
        pos = np.minimum(pos, len(key_s) - 1)
        found = key_s[pos] == kk
        sim = np.where(found, s_s[pos], np.nan)
        d = -np.log(np.clip(sim, 1e-3, 1.0))
        dmax = np.nanmax(d) if np.isfinite(d).any() else 1.0
        d = np.where(np.isnan(d), dmax + 1.0, d)
        # outgroup leaf
        og = None
        if out_adj is not None:
            tot: dict = {}
            for u in units:
                for o, sc in out_adj.get(u, {}).items():
                    tot[o] = tot.get(o, 0.0) + sc
            if tot:
                og = max(tot, key=tot.get)
        if require_outgroup and og is None:
            stats["tree_no_outgroup"] = stats.get("tree_no_outgroup", 0) + 1
            continue
        size = m + (1 if og is not None else 0)
        D = np.zeros((size, size))
        D[I, J] = d
        D[J, I] = d
        leaf_child = [int(child[u]) for u in units]
        if og is not None:
            do = np.array([-np.log(min(max(out_adj.get(u, {}).get(og, np.nan), 1e-3), 1.0))
                           if og in out_adj.get(u, {}) else np.nan for u in units])
            do = np.where(np.isnan(do), max(dmax, np.nanmax(do) if np.isfinite(do).any() else dmax) + 1.0, do)
            D[m, :m] = do
            D[:m, m] = do
            leaf_child.append(-1)
        groups = tree_groups(D, leaf_child, m if og is not None else None)
        old = {u: labels[u] for u in units}
        new = {}
        for grp in groups:
            lab = next_label
            next_label += 1
            for li in grp:
                new[units[li]] = lab
        if conservative:
            # never fuse two multi-unit rule groups: a tree group may contain units of at
            # most one multi-unit rule group (plus leftover singletons); otherwise it is
            # split back along the rule groups
            sizes_old = {}
            for u in units:
                sizes_old[old[u]] = sizes_old.get(old[u], 0) + 1
            fixed = {}
            for grp in groups:
                us = [units[li] for li in grp]
                big = {old[u] for u in us if sizes_old[old[u]] > 1}
                if len(big) <= 1:
                    lab = next_label
                    next_label += 1
                    for u in us:
                        fixed[u] = lab
                else:
                    for b in big:
                        lab = next_label
                        next_label += 1
                        for u in us:
                            if old[u] == b:
                                fixed[u] = lab
                    for u in us:
                        if u not in fixed:
                            fixed[u] = next_label
                            next_label += 1
            new = fixed
        # count partition changes
        part_old = {frozenset(u for u in units if old[u] == l) for l in set(old.values())}
        part_new = {frozenset(u for u in units if new[u] == l) for l in set(new.values())}
        if part_old != part_new:
            stats["tree_changed"] += 1
            for u, lab in new.items():
                labels[u] = lab
    return labels
