"""S3 -- level-synchronous bottom-up merge over the guide tree (core v0, sequence only).

See docs/core_v0.md for the full description; summary:

State. After node ``X`` is processed, the genes of clade(X) are partitioned
into *ancestral genes of X* (global integer IDs). Each carries a nested
structure (leaf gene ID, or ``("O", level, children)`` / ``("P", level,
children)``) that is written as OrthoXML orthologGroup / paralogGroup. Leaf
genes are the ancestral genes of their species.

Similarity graph. Leaf-pair scores ``s(a, b)`` (S2, within-family only; bits or
normalised bits; pairs failing the bidirectional coverage guard are dropped)
are introduced exactly once, at the node where the pair *becomes* a
cross-child pair (the LCA of the two species; same-species pairs at the leaf).
After each node the graph is *contracted* along the merge map, keeping per
edge (sum, count, max) of the member leaf-pair scores, so

* max-link  S(A, B) = max_{a in A, b in B} s(a, b)
* avg-link  S(A, B) = mean over scored member pairs

are exact without re-searching. Cost per node is O(E_X log E_X) where E_X =
edges among the children's ancestral genes (all families at once).

Merge rule at X with children C1..Ck (k-way, multifurcations allowed):

R1 *seeds* -- for every unit A (ancestral gene of child i) and child j != i,
   best(A, j) = max_B in Cj S(A, B). An edge (A, B) is *mutual-best* if
   S >= tol * best(A, j) and S >= tol * best(B, i). Mutual-best edges are
   unioned greedily in decreasing score; a union that would put two units of
   the same child into one group is allowed only if they are in-paralogs:
   S_same(P_c, Q_c) >= inpar_factor * S(edge).
R2 *in-paralog absorption* -- an unassigned unit u (child c) whose best
   cross-child partner v lies in group G joins G if for some w in G of child c
   ``inpar_rule="ortholog"`` (InParanoid): S(u, w) >= inpar_factor * d(w), d(w) =
   w's best score to G's members in the other children; ``"self"``: S(u, w) >=
   inpar_factor * S(u, v) (u is closer to w than to its own best ortholog;
   robust to a fast-evolving duplicate u, since on an additive tree the path
   u-w is shorter than u-v whenever the duplication is after X); ``"either"``.
   ``"quartet"``: rate-robust four-point test with the closest outgroup gene o
   (a gene of a sibling clade of X scored against u, w and v): absorb iff
   S(u,w)*S(v,o) >= inpar_factor * S(w,v)*S(u,o), i.e. d(u,w)+d(v,o) <=
   d(w,v)+d(u,o) for d = -ln S, topology ((u,w),(v,o)); falls back to
   ``"ortholog"`` when no outgroup gene has scores to all three.
   => duplication after X on the branch to c.
R3 *fission repair* -- iterated mutual-best between remaining unassigned units
   and groups that have no member of the unit's child, accepted only if the
   score is >= repair_ratio * (the unit's best cross-child score), so out-
   paralogs of an existing group are not attached (hidden-paralogy guard).

Every group with >= 2 units becomes an ancestral gene of X:
``("O", X, [unit or ("P", child, units) per child])``. Unmatched units pass up
unchanged (their HOG at X contains one child's genes only).

Extension hooks (not used in v0): ``context`` -- a callable that returns an
additive bonus per candidate cross edge (synteny); ``resolver`` -- a callable
that may override the decision for an ambiguous connected cluster (local-tree
mode). Their signatures are documented on :func:`run_merge`.
"""
from __future__ import annotations

import time
from dataclasses import asdict, dataclass, replace

import numpy as np

from .io import GeneTable
from .newick import Node
from .orthoxml import GeneEntry, Group, OrthoXML, normalize_groups
from .scoring import ScoreTable

__all__ = ["MergeConfig", "MergeResult", "run_merge", "to_orthoxml", "decide_node", "ARMS"]


@dataclass
class MergeConfig:
    link: str = "max"            # max | avg
    score: str = "norm_max"      # bits | norm_min | norm_max | of (see leaf_scores)
    tol: float = 0.95            # mutual-best tolerance (1.0 = strict RBH)
    inparalogs: bool = True      # R2
    inpar_factor: float = 1.0
    inpar_rule: str = "ortholog"  # ortholog | self | either | quartet (see module docstring)
    repair: bool = True          # R3
    repair_ratio: float = 0.5
    min_cov: float = 0.5         # bidirectional coverage guard on leaf pairs
    min_score: float = 0.0       # minimum aggregated score for a candidate edge

    def key(self) -> dict:
        return asdict(self)


_D = MergeConfig()
ARMS = {   # named arms of the S3/S4 benchmark (docs/core_v0.md)
    "rbh": replace(_D, score="bits", tol=1.0, inparalogs=False, repair=False),
    "max": replace(_D, link="max", repair=False),
    "max+repair": replace(_D, link="max", repair=True),
    "avg": replace(_D, link="avg", repair=False),
    "avg+repair": replace(_D, link="avg", repair=True),
    "max+repair/bits": replace(_D, link="max", repair=True, score="bits"),
    "avg+repair/bits": replace(_D, link="avg", repair=True, score="bits"),
    "max+repair/quartet": replace(_D, link="max", repair=True, inpar_rule="quartet"),
    "avg+repair/quartet": replace(_D, link="avg", repair=True, inpar_rule="quartet"),
    "avg+repair/either": replace(_D, link="avg", repair=True, inpar_rule="either"),
}


@dataclass
class MergeResult:
    roots: list          # structures of the root's ancestral genes (int = singleton gene)
    node_stats: list     # per internal node dicts
    n_genes: int


# ----------------------------------------------------------------------------
# graph helpers
# ----------------------------------------------------------------------------

def _aggregate(u, v, ssum, scnt, smax, cap):
    """Merge duplicate undirected edges (u, v) -> unique (min, max) with sum/count/max."""
    if len(u) == 0:
        z = np.zeros(0)
        return np.zeros(0, np.int64), np.zeros(0, np.int64), z, z, z
    a = np.minimum(u, v)
    b = np.maximum(u, v)
    key = a.astype(np.int64) * cap + b
    order = np.argsort(key, kind="stable")
    key = key[order]
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    S = np.add.reduceat(ssum[order], starts)
    C = np.add.reduceat(scnt[order], starts)
    M = np.maximum.reduceat(smax[order], starts)
    k = key[starts]
    return k // cap, k % cap, S, C, M


class _UF:
    __slots__ = ("p", "members")

    def __init__(self, n, child):
        self.p = list(range(n))
        # members[root] = {child: [units]}
        self.members = [{child[i]: [i]} for i in range(n)]

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        ma, mb = self.members[ra], self.members[rb]
        if sum(map(len, ma.values())) < sum(map(len, mb.values())):
            ra, rb, ma, mb = rb, ra, mb, ma
        for c, lst in mb.items():
            ma.setdefault(c, []).extend(lst)
        self.members[rb] = None
        self.p[rb] = ra
        return ra

    def size(self, x):
        return sum(map(len, self.members[self.find(x)].values()))


# ----------------------------------------------------------------------------
# per-node decision (pure function of the local graph; unit-tested)
# ----------------------------------------------------------------------------

def decide_node(n: int, child: np.ndarray, lu: np.ndarray, lv: np.ndarray, score: np.ndarray,
                cfg: MergeConfig, context=None, resolver=None, node_name=None, out_edges=None):
    """Partition ``n`` local units into ancestral genes of the node.

    ``child[i]``: child index of unit i; edges ``(lu, lv, score)`` undirected,
    lu != lv. ``out_edges`` = (unit, outgroup_gene, score) arrays for the quartet
    rule. Returns (labels array length n, stats dict).
    """
    stats = dict(n_units=n, n_mb=0, r1_unions=0, r1_blocked=0, r2_absorbed=0, r3_repaired=0,
                 r2_quartet_tests=0)
    out_adj: list[dict] = [dict() for _ in range(n)] if out_edges is not None else []
    if out_edges is not None:
        for x, o, sc_ in zip(out_edges[0].tolist(), out_edges[1].tolist(), out_edges[2].tolist()):
            out_adj[x][o] = sc_
    child_l = child.tolist()
    cross = child[lu] != child[lv]
    ca, cb, cs = lu[cross], lv[cross], score[cross].astype(np.float64)
    if context is not None:          # hook: synteny / context evidence (additive bonus)
        cs = cs + context(node_name, ca, cb, cs)
    sa, sb, ss = lu[~cross], lv[~cross], score[~cross]
    same = {}
    for x, y, s in zip(sa.tolist(), sb.tolist(), ss.tolist()):
        same[(x, y) if x < y else (y, x)] = s
    k = int(child.max()) + 1 if n else 1
    E = len(ca)
    uf = _UF(n, child_l)
    if E == 0:
        return np.array([uf.find(i) for i in range(n)]), stats
    # best(A, j) for both endpoints
    dkey = np.concatenate([ca.astype(np.int64) * k + child[cb], cb.astype(np.int64) * k + child[ca]])
    dval = np.concatenate([cs, cs])
    uniq, inv = np.unique(dkey, return_inverse=True)
    best = np.full(len(uniq), -np.inf)
    np.maximum.at(best, inv, dval)
    best_a, best_b = best[inv[:E]], best[inv[E:]]
    if cfg.tol >= 1.0:
        mb = (cs >= best_a) & (cs >= best_b)
    else:
        mb = (cs >= cfg.tol * best_a) & (cs >= cfg.tol * best_b)
    mb &= cs > cfg.min_score
    stats["n_mb"] = int(mb.sum())
    # best cross score per unit (any child) and cross adjacency
    best_any = np.full(n, -np.inf)
    np.maximum.at(best_any, ca, cs)
    np.maximum.at(best_any, cb, cs)
    adj: list[dict] = [dict() for _ in range(n)]
    for x, y, s in zip(ca.tolist(), cb.tolist(), cs.tolist()):
        if s > cfg.min_score:
            adj[x][y] = s
            adj[y][x] = s

    def same_sim(x, y):
        return same.get((x, y) if x < y else (y, x), -np.inf)

    # ---- R1: greedy union of mutual-best edges with in-paralog constraint
    idx = np.flatnonzero(mb)
    idx = idx[np.argsort(-cs[idx], kind="stable")]
    for e in idx.tolist():
        x, y, s = int(ca[e]), int(cb[e]), float(cs[e])
        rx, ry = uf.find(x), uf.find(y)
        if rx == ry:
            continue
        mx, my = uf.members[rx], uf.members[ry]
        ok = True
        for c in mx.keys() & my.keys():
            sim = max(same_sim(p, q) for p in mx[c] for q in my[c])
            if sim < cfg.inpar_factor * s:
                ok = False
                break
        if ok:
            uf.union(x, y)
            stats["r1_unions"] += 1
        else:
            stats["r1_blocked"] += 1

    # ---- R2: in-paralog absorption
    def d_in_group(w, groot):
        mem = uf.members[groot]
        cw = child_l[w]
        best_d = -np.inf
        for c, lst in mem.items():
            if c == cw:
                continue
            for z in lst:
                s = adj[w].get(z)
                if s is not None and s > best_d:
                    best_d = s
        return best_d

    if cfg.inparalogs:
        for _ in range(5):
            changed = 0
            cand = []
            for u in range(n):
                if adj[u] and uf.size(u) == 1:
                    v, s = max(adj[u].items(), key=lambda kv: kv[1])
                    cand.append((s, u, v))
            cand.sort(reverse=True)
            for s, u, v in cand:
                if uf.size(u) != 1:
                    continue
                g = uf.find(v)
                mem = uf.members[g]
                cu = child_l[u]
                if cu not in mem or len(mem) < 2:
                    continue
                ok = False
                for w in mem[cu]:
                    sw = same_sim(u, w)
                    if sw == -np.inf:
                        continue
                    if cfg.inpar_rule == "quartet" and out_adj:
                        swv = adj[w].get(v)
                        common = out_adj[u].keys() & out_adj[w].keys() & out_adj[v].keys()
                        if swv is not None and common:
                            o = max(common, key=lambda z: out_adj[w][z] + out_adj[v][z])
                            stats["r2_quartet_tests"] += 1
                            if sw * out_adj[v][o] >= cfg.inpar_factor * swv * out_adj[u][o]:
                                ok = True
                                break
                            continue
                    if cfg.inpar_rule in ("self", "either") and sw >= cfg.inpar_factor * s:
                        ok = True       # u is closer to w than to its best ortholog in G
                        break
                    if cfg.inpar_rule in ("ortholog", "either", "quartet") and \
                            sw >= cfg.inpar_factor * d_in_group(w, g):
                        ok = True       # w is closer to u than to its orthologs in G
                        break
                if ok:
                    uf.union(u, v)
                    stats["r2_absorbed"] += 1
                    changed += 1
            if not changed:
                break

    # ---- R3: fission repair (mutual best among leftovers, hidden-paralogy guard)
    if cfg.repair:
        for _ in range(20):
            singles = [u for u in range(n) if adj[u] and uf.size(u) == 1]
            if not singles:
                break
            u_best = {}
            q_best: dict = {}
            for u in singles:
                cu = child_l[u]
                per_g: dict = {}
                for z, s in adj[u].items():
                    g = uf.find(z)
                    if cu in uf.members[g]:
                        continue
                    if s > per_g.get(g, -np.inf):
                        per_g[g] = s
                if not per_g:
                    continue
                g, s = max(per_g.items(), key=lambda kv: kv[1])
                if s < cfg.repair_ratio * best_any[u]:
                    continue
                u_best[u] = (g, s)
                for g2, s2 in per_g.items():
                    key = (g2, cu)
                    if s2 > q_best.get(key, (None, -np.inf))[1]:
                        q_best[key] = (u, s2)
            merged = 0
            for u, (g, s) in sorted(u_best.items(), key=lambda kv: -kv[1][1]):
                if uf.size(u) != 1 or uf.find(g) != g:
                    continue
                if q_best.get((g, child_l[u]), (None,))[0] != u:
                    continue
                mem = uf.members[g]
                if child_l[u] in mem:
                    continue
                if sum(map(len, mem.values())) == 1:      # singleton partner: guard both sides
                    q = next(iter(mem.values()))[0]
                    if q in u_best and u_best[q][0] != uf.find(u):
                        continue
                    if s < cfg.repair_ratio * best_any[q]:
                        continue
                uf.union(u, g)
                stats["r3_repaired"] += 1
                merged += 1
            if not merged:
                break
    labels = np.array([uf.find(i) for i in range(n)], dtype=np.int64)
    if resolver is not None:         # hook: local-tree mode for ambiguous clusters
        labels = resolver(node_name, child, lu, lv, score, labels)
    return labels, stats


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------

def _of_normalise(bits, x, sp_pair, min_hits=100, nbins=50, top=0.05):
    """OrthoFinder-style: per species pair fit log10(bits) ~ a + b*log10(len_q*len_t)
    on the top 5% hits of each length bin; score = bits / 10^(a + b x)."""
    y = np.log10(np.maximum(bits, 1e-3))

    def fit(xx, yy):
        if len(xx) < 10:
            return None
        edges = np.quantile(xx, np.linspace(0, 1, nbins + 1))
        bi = np.clip(np.searchsorted(edges, xx, side="right") - 1, 0, nbins - 1)
        xs, ys = [], []
        for b in range(nbins):
            m = bi == b
            if not m.any():
                continue
            yb, xb = yy[m], xx[m]
            k = max(1, int(np.ceil(top * m.sum())))
            o = np.argsort(-yb)[:k]
            xs.append(xb[o])
            ys.append(yb[o])
        xs, ys = np.concatenate(xs), np.concatenate(ys)
        if np.ptp(xs) == 0:
            return None
        b1, b0 = np.polyfit(xs, ys, 1)
        return b0, b1
    glob = fit(x, y) or (0.0, 0.0)
    out = np.empty_like(bits, dtype=np.float64)
    order = np.argsort(sp_pair, kind="stable")
    sp_sorted = sp_pair[order]
    cuts = np.flatnonzero(np.diff(sp_sorted)) + 1
    for grp in np.split(order, cuts):
        if len(grp) == 0:
            continue
        f = fit(x[grp], y[grp]) if len(grp) >= min_hits else None
        b0, b1 = f if f is not None else glob
        out[grp] = bits[grp] / 10 ** (b0 + b1 * x[grp])
    return out


def leaf_scores(gt: GeneTable, st: ScoreTable, mode: str) -> np.ndarray:
    """Leaf-pair similarity used by the merge.

    ``bits``: raw bitscore; ``norm_min`` / ``norm_max``: bits / min|max(self_a, self_b);
    ``of``: OrthoFinder-style length normalisation per species pair.
    """
    q, t = st.q.astype(np.int64), st.t.astype(np.int64)
    bits = st.bits.astype(np.float64)
    if mode == "bits":
        return bits
    if mode in ("norm_min", "norm_max"):
        sb = st.self_bits_filled(gt.lengths)
        f = np.minimum if mode == "norm_min" else np.maximum
        return bits / f(sb[q], sb[t])
    if mode == "of":
        L = gt.lengths.astype(np.float64)
        x = np.log10(L[q] * L[t])
        a, b = gt.species_of[q].astype(np.int64), gt.species_of[t].astype(np.int64)
        S = len(gt.species)
        sp_pair = np.minimum(a, b) * S + np.maximum(a, b)
        return _of_normalise(bits, x, sp_pair)
    raise ValueError(f"unknown score mode {mode!r}")


def _lca_matrix(tree: Node, species: list[str]):
    nodes = list(tree.postorder())
    nid = {nd.name: i for i, nd in enumerate(nodes)}
    sidx = {s: i for i, s in enumerate(species)}
    S = len(species)
    lca = np.full((S, S), -1, dtype=np.int32)
    clade: dict = {}
    for nd in nodes:
        if not nd.children:
            clade[nd.name] = np.array([sidx[nd.name]], dtype=np.int64)
            lca[sidx[nd.name], sidx[nd.name]] = nid[nd.name]
            continue
        cl = [clade[c.name] for c in nd.children]
        for i in range(len(cl)):
            for j in range(i + 1, len(cl)):
                lca[np.ix_(cl[i], cl[j])] = nid[nd.name]
                lca[np.ix_(cl[j], cl[i])] = nid[nd.name]
        clade[nd.name] = np.concatenate(cl)
    return nodes, nid, lca, clade


def run_merge(gt: GeneTable, tree: Node, st: ScoreTable, cfg: MergeConfig,
              context=None, resolver=None, log=None, trace: dict | None = None) -> MergeResult:
    """Run the bottom-up merge over all families at once.

    Hooks: ``context(node, a, b, s) -> bonus array`` (added to cross-child
    candidate scores); ``resolver(node, child, lu, lv, score, labels) -> labels``.
    ``trace`` (dict, diagnostics): filled with node -> (clade genes, their ancestral
    gene ID at the children, at the node).
    """
    N = gt.n
    cap = 2 * N + 2
    nodes, nid, lca, clade_sp = _lca_matrix(tree, gt.species)
    sp_genes = gt.genes_of_species()
    # leaf-pair scores
    q, t = st.q.astype(np.int64), st.t.astype(np.int64)
    keep = np.minimum(st.qcov, st.tcov) >= cfg.min_cov
    s = leaf_scores(gt, st, cfg.score)
    q, t, s = q[keep], t[keep], s[keep]
    plca = lca[gt.species_of[q], gt.species_of[t]]
    order = np.argsort(plca, kind="stable")
    q, t, s, plca = q[order], t[order], s[order], plca[order]
    bounds = np.searchsorted(plca, np.arange(len(nodes) + 1))
    struct: list = list(range(N))          # anc id -> structure
    size = np.zeros(cap, dtype=np.int64)
    size[:N] = 1
    next_id = N
    gene_anc = np.arange(N, dtype=np.int64)
    node_anc: dict = {}
    node_edges: dict = {}
    pos = np.full(cap, -1, dtype=np.int64)
    node_stats = []
    clade_genes = {}
    for ni, nd in enumerate(nodes):
        t0 = time.time()
        lo, hi = bounds[ni], bounds[ni + 1]
        nq, nt, ns = q[lo:hi], t[lo:hi], s[lo:hi]
        if not nd.children:
            g = sp_genes[gt.species.index(nd.name)]
            clade_genes[nd.name] = g
            node_anc[nd.name] = g.astype(np.int64)
            node_edges[nd.name] = _aggregate(nq, nt, ns, np.ones_like(ns), ns, cap)
            continue
        kids = [c.name for c in nd.children]
        units = np.concatenate([node_anc[c] for c in kids])
        child = np.concatenate([np.full(len(node_anc[c]), i, dtype=np.int64)
                                for i, c in enumerate(kids)])
        genes = np.concatenate([clade_genes.pop(c) for c in kids])
        clade_genes[nd.name] = genes
        newe = _aggregate(gene_anc[nq], gene_anc[nt], ns, np.ones_like(ns), ns, cap)
        parts = [node_edges.pop(c) for c in kids] + [newe]
        eu, ev, es, ec, em = (np.concatenate([p[i] for p in parts]) for i in range(5))
        n = len(units)
        pos[units] = np.arange(n)
        lu, lv = pos[eu], pos[ev]
        assert (lu >= 0).all() and (lv >= 0).all()
        score = em if cfg.link == "max" else es / np.maximum(ec, 1)
        out_edges = None
        if cfg.inpar_rule == "quartet" and nd.parent is not None:
            pi = nid[nd.parent.name]
            pq, pt, ps = q[bounds[pi]:bounds[pi + 1]], t[bounds[pi]:bounds[pi + 1]], s[bounds[pi]:bounds[pi + 1]]
            in_x = np.zeros(len(gt.species), bool)
            in_x[clade_sp[nd.name]] = True
            a_in, b_in = in_x[gt.species_of[pq]], in_x[gt.species_of[pt]]
            m1, m2 = a_in & ~b_in, b_in & ~a_in
            ug = np.concatenate([pq[m1], pt[m2]])
            og = np.concatenate([pt[m1], pq[m2]])
            os_ = np.concatenate([ps[m1], ps[m2]])
            if len(ug):
                ou, oo, osum, ocnt, omax = _aggregate(pos[gene_anc[ug]], og + cap, os_,
                                                      np.ones_like(os_), os_, 2 * cap + 2)
                # _aggregate orders each pair as (min, max): unit index < og + cap always
                osc = omax if cfg.link == "max" else osum / np.maximum(ocnt, 1)
                out_edges = (ou, oo - cap, osc)
        labels, stt = decide_node(n, child, lu, lv, score, cfg, context, resolver, nd.name,
                                  out_edges)
        pos[units] = -1
        # build ancestral genes
        newid = {}
        order_l = np.argsort(labels, kind="stable")
        sl = labels[order_l]
        cuts = np.flatnonzero(np.diff(sl)) + 1
        n_new = 0
        remap_units = np.empty(n, dtype=np.int64)
        for grp in np.split(order_l, cuts):
            if len(grp) == 1:
                remap_units[grp[0]] = units[grp[0]]
                continue
            by_child: dict = {}
            for i in grp.tolist():
                by_child.setdefault(int(child[i]), []).append(int(units[i]))
            parts_s = []
            for ci in sorted(by_child):
                us = by_child[ci]
                if len(us) == 1:
                    parts_s.append(struct[us[0]])
                else:
                    parts_s.append(("P", kids[ci], [struct[x] for x in us]))
            aid = next_id
            next_id += 1
            if len(parts_s) == 1:
                struct.append(parts_s[0])
            else:
                struct.append(("O", nd.name, parts_s))
            size[aid] = size[units[grp]].sum()
            remap_units[grp] = aid
            n_new += 1
            for i in grp.tolist():   # free memory of merged children's structures
                struct[int(units[i])] = None
        remap = np.arange(cap, dtype=np.int64)
        remap[units] = remap_units
        out_units = np.unique(remap_units)
        node_anc[nd.name] = out_units
        for c in kids:
            del node_anc[c]
        before = gene_anc[genes].copy() if trace is not None else None
        gene_anc[genes] = remap[gene_anc[genes]]
        if trace is not None:
            trace[nd.name] = (genes.copy(), before, gene_anc[genes].copy())
        cu_, cv_ = remap[eu], remap[ev]
        k = cu_ != cv_
        node_edges[nd.name] = _aggregate(cu_[k], cv_[k], es[k], ec[k], em[k], cap)
        stt.update(node=nd.name, n_children=len(kids), n_clade_genes=int(len(genes)),
                   n_edges=int(len(eu)), n_new_pairs=int(hi - lo), n_anc_out=int(len(out_units)),
                   n_merged_groups=n_new, wall_s=round(time.time() - t0, 4))
        node_stats.append(stt)
        if log:
            log(f"node {nd.name}: k={len(kids)} units={n} -> anc={len(out_units)} "
                f"(mb={stt['n_mb']} r1={stt['r1_unions']} blocked={stt['r1_blocked']} "
                f"r2={stt['r2_absorbed']} r3={stt['r3_repaired']}) {stt['wall_s']:.2f}s")
    root = nodes[-1].name
    roots = [struct[int(a)] for a in node_anc[root]]
    return MergeResult(roots=roots, node_stats=node_stats, n_genes=N)


def _to_group(x, names):
    if isinstance(x, (int, np.integer)):
        return names[int(x)]
    kind, level, kids = x
    return Group("ortholog" if kind == "O" else "paralog",
                 [_to_group(c, names) for c in kids], level=level)


def to_orthoxml(res: MergeResult, gt: GeneTable, tree: Node) -> OrthoXML:
    genes = {nm: GeneEntry(nm, gt.species[gt.species_of[g]]) for g, nm in enumerate(gt.names)}
    groups = [_to_group(x, gt.names) for x in res.roots if not isinstance(x, (int, np.integer))]
    return OrthoXML(genes=genes, groups=normalize_groups(groups), species_tree=tree)
