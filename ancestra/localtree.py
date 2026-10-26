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
    """Groups (lists of leaf indices, outgroup excluded) of a cluster's distance matrix (NJ)."""
    m = D.shape[0]
    edges, n_nodes = nj(D)
    if not edges:
        return [[0]]
    return groups_from_edges(edges, n_nodes, m, leaf_child, outgroup)


def groups_from_edges(edges, n_nodes, m, leaf_child, outgroup):
    """Species-overlap labelling of an unrooted tree given as edges (leaves 0..m-1)."""
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


def newick_to_edges(nwk: str, leaf_index: dict):
    """FastTree Newick -> (edges, n_nodes) with leaves numbered by ``leaf_index[name]``."""
    from .newick import parse_newick
    root = parse_newick(nwk)
    m = len(leaf_index)
    ids = {}
    nxt = m
    for nd in root.preorder():
        if nd.children:
            ids[id(nd)] = nxt
            nxt += 1
        else:
            ids[id(nd)] = leaf_index[nd.name]
    edges = []
    for nd in root.preorder():
        for c in nd.children:
            edges.append((ids[id(nd)], ids[id(c)], max(c.length or 0.0, 0.0)))
    # a bifurcating root with a unary-like structure is fine for the unrooted labelling
    return edges, nxt


def _fasttree_job(job):
    """job = (leaf names, sequences, workdir) -> Newick string (FAMSA + FastTree -lg)."""
    import os
    import subprocess
    import tempfile
    names, seqs, wd = job
    with tempfile.TemporaryDirectory(dir=wd) as td:
        fa, aln = os.path.join(td, "in.fa"), os.path.join(td, "aln.fa")
        with open(fa, "w") as fh:
            for nm, sq in zip(names, seqs):
                fh.write(f">{nm}\n{sq}\n")
        subprocess.run(["famsa", "-t", "1", fa, aln], check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        r = subprocess.run(["FastTree", "-lg", "-quiet", "-nosupport", aln], check=True,
                           capture_output=True, text=True)
        return r.stdout.strip()


def resolve_ambiguous(n, child, lu, lv, score, labels, out_edges=None, tau: float = 0.9,
                      max_units: int = 60, stats=None, require_outgroup: bool = False,
                      conservative: bool = False, tree_method: str = "nj", members=None,
                      seqs=None, ft_mode: str = "medoid", ft_max_genes: int = 60, workers: int = 12,
                      workdir: str | None = None):
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
    state = {"next": int(labels.max()) + 1}
    ft_jobs = []

    def apply_groups(units, groups):
        """Relabel ``units`` by ``groups`` (lists of positions in ``units``) if the partition
        changed (``conservative``: never fuse two multi-unit rule groups)."""
        old = {u: labels[u] for u in units}
        new = {}
        for grp in groups:
            lab = state["next"]
            state["next"] += 1
            for li in grp:
                new[units[li]] = lab
        if conservative:
            sizes_old = {}
            for u in units:
                sizes_old[old[u]] = sizes_old.get(old[u], 0) + 1
            fixed = {}
            for grp in groups:
                us = [units[li] for li in grp]
                big = {old[u] for u in us if sizes_old[old[u]] > 1}
                if len(big) <= 1:
                    lab = state["next"]
                    state["next"] += 1
                    for u in us:
                        fixed[u] = lab
                else:
                    for b in big:
                        lab = state["next"]
                        state["next"] += 1
                        for u in us:
                            if old[u] == b:
                                fixed[u] = lab
                    for u in us:
                        if u not in fixed:
                            fixed[u] = state["next"]
                            state["next"] += 1
            new = fixed
        part_old = {frozenset(u for u in units if old[u] == l) for l in set(old.values())}
        part_new = {frozenset(u for u in units if new[u] == l) for l in set(new.values())}
        if part_old != part_new:
            stats["tree_changed"] += 1
            for u, lab in new.items():
                labels[u] = lab

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
        if tree_method == "fasttree" and members is not None and seqs is not None:
            leaves, leaf_unit = [], []
            mem = [np.asarray(members(u)) for u in units]
            allg = sum(len(x) for x in mem)
            use_all = ft_mode == "all" and allg <= ft_max_genes
            for ui, ms in enumerate(mem):
                if len(ms) == 0:
                    continue
                if use_all:
                    for g in ms.tolist():
                        leaves.append(int(g))
                        leaf_unit.append(ui)
                else:   # length-medoid member
                    o = sorted(ms.tolist(), key=lambda g: (len(seqs[g]), g))
                    leaves.append(int(o[len(o) // 2]))
                    leaf_unit.append(ui)
            if og is not None:
                leaves.append(int(og))
                leaf_unit.append(-1)
            if len(leaves) >= 4:
                ft_jobs.append((units, leaves, leaf_unit, og is not None, use_all))
                continue
        groups = tree_groups(D, leaf_child, m if og is not None else None)
        apply_groups(units, groups)
    if ft_jobs:
        import os
        import time as _t
        from concurrent.futures import ThreadPoolExecutor
        wd = workdir or os.environ.get("TMPDIR", "/tmp")
        t0 = _t.time()
        payload = [([f"L{i}" for i in range(len(lv_))], [seqs[g] for g in lv_], wd)
                   for (_, lv_, _, _, _) in ft_jobs]
        with ThreadPoolExecutor(workers) as ex:
            nwks = list(ex.map(_fasttree_job, payload))
        stats["tree_fasttree_jobs"] = len(ft_jobs)
        stats["tree_fasttree_s"] = round(_t.time() - t0, 2)
        for (units, lv_, leaf_unit, has_og, use_all), nwk in zip(ft_jobs, nwks):
            nl = len(lv_)
            edges, n_nodes = newick_to_edges(nwk, {f"L{i}": i for i in range(nl)})
            lc = [int(child[units[u]]) if u >= 0 else -1 for u in leaf_unit]
            ogl = nl - 1 if has_og else None
            lgroups = groups_from_edges(edges, n_nodes, nl, lc, ogl)
            # collapse leaves -> units (majority group per unit)
            votes: dict = {}
            for gi, grp in enumerate(lgroups):
                for li in grp:
                    u = leaf_unit[li]
                    if u >= 0:
                        votes.setdefault(u, {}).setdefault(gi, 0)
                        votes[u][gi] += 1
            ug: dict = {}
            for u, v in votes.items():
                gi = max(sorted(v), key=lambda k: v[k])
                ug.setdefault(gi, []).append(u)
            for u in range(len(units)):
                if u not in votes:
                    ug.setdefault(("solo", u), []).append(u)
            apply_groups(units, list(ug.values()))
    return labels
