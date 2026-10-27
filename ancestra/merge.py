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

__all__ = ["MergeConfig", "MergeResult", "run_merge", "run_merge_constrained", "ConstrainedContext",
           "root_labels", "to_orthoxml", "decide_node", "ARMS", "NodeContext"]


@dataclass
class MergeConfig:
    link: str = "avg"            # avg | max   (default chosen by the pre-stated toy rule)
    score: str = "norm_max"      # bits | norm_min | norm_max | of (see leaf_scores)
    tol: float = 0.95            # mutual-best tolerance (1.0 = strict RBH)
    inparalogs: bool = True      # R2
    inpar_factor: float = 1.0
    inpar_rule: str = "ortholog"  # ortholog | self | either | quartet (see module docstring)
    repair: bool = True          # R3 (toy sweep: +0.0003 F1; results/core_v0/sweep.tsv)
    repair_ratio: float = 0.9
    min_cov: float = 0.0         # bidirectional coverage guard on leaf pairs (v0: 0.5; OrthoBench
                                 # dev: the guard removes fragments/odd gene models, -9 F)
    min_score: float = 0.0       # minimum aggregated score for a candidate edge
    # relative-divergence (RD) tests with the closest outgroup gene o (needs a norm_* score):
    # r(a,b) = d(a,b) / mean(d(a,o), d(b,o)), d = -ln S; calibrated per node on the R1 seeds
    rd_hidden: bool = False      # veto R1 seeds with r > median + rd_z * 1.4826 * MAD (hidden paralogy)
    rd_z: float = 3.0
    rd_inpar: bool = False       # R2: absorb u if r(u, w) <= rd_inpar_factor * median seed r
    rd_inpar_factor: float = 1.0
    engine: str = "numba"        # numba | python (reference implementation)
    tree_mode: bool = True       # hybrid local-tree mode for ambiguous clusters (ancestra.localtree)
    tree_tau: float = 0.9        # ambiguity edge: S(u,v) >= tau * min(best(u), best(v))
    tree_max_units: int = 60     # larger clusters keep the rule result
    tree_require_outgroup: bool = True   # only clusters rooted by an outgroup gene
    tree_conservative: bool = False      # tree may absorb singletons / split, never fuse groups
    tree_nodes: str = ""                 # comma-separated node names to restrict tree mode ("" = all)
    tree_method: str = "nj"              # nj (default) | fasttree (FAMSA + FastTree -lg on sequences)
    tree_ft_mode: str = "medoid"         # fasttree leaves: medoid (1 per unit) | all (genes, <=60)
    tree_ft_workers: int = 12
    tree_min_d: float = 0.5              # tree mode only at nodes whose median seed distance
                                         # -ln S >= this (divergent nodes; rules suffice shallow)
    ctx_mode: str = "r1"         # context evidence used only for R1 partner choice | "all" (additive)
    ctx_window: float = 0.75     # R1-only: bonus only for candidates >= window x best raw score
    ctx_r3: bool = False         # also use score+bonus for R3 (fission repair) partner choice
    ctx_veto: bool = False       # R2/R3 join of u blocked if its support is 0 while another
    ctx_veto_min: float = 2.0    #   candidate of u has support >= ctx_veto_min (needs .support)
    locus_repair: bool = False   # locus-aware fission repair (S8 proposal): a still-alone unit u
    locus_single_species: bool = True  # only units whose genes are all from one species
    locus_min_support: float = 1.0  # joins the group lacking child(u) with the highest context
                                 # support (>= locus_min_support) among in-family-scored partners
    so_struct: bool = False      # species-overlap labelling of same-child units inside a group:
                                 # species-disjoint units -> speciation (orthologGroup at their MRCA)
                                 # instead of a paralogGroup (QfO recall, missed speciations)
    so_absorb: bool = False      # a still-alone unit whose best partner's group holds only
                                 # species-disjoint units of its own child joins it (missed speciation)
    deep_min_d: float = 0.0      # "deep" node = median mutual-best seed distance -ln S >= this
                                 # (0 = off). QfO: Eukaryota 1.81, Opisthokonta 1.68 vs Bilateria 1.36
    deep_tol: float = 0.0        # at deep nodes: mutual-best tolerance (0 = keep tol)
    deep_so_absorb: bool = False # at deep nodes: species-overlap absorption (so_absorb)
    orphan_min: float = 0.0      # R4: minimum (normalised, linked) score of an orphan-merge edge
    orphan_merge: str = "off"    # R4: off | root | all -- merge partner-less same-child units
                                 # (no cross-child homolog at the node) connected by same-child
                                 # edges into one paralogGroup (parsimony: 1 dup vs dup + losses)

    def key(self) -> dict:
        d = asdict(self)
        # options added after the QfO freeze: omitted at their defaults (stable stage hashes)
        for k, v in (("tree_method", "nj"), ("tree_ft_mode", "medoid"), ("so_struct", False),
                     ("so_absorb", False), ("deep_min_d", 0.0), ("deep_tol", 0.0),
                     ("deep_so_absorb", False)):
            if d.get(k) == v:
                d.pop(k)
        d.pop("tree_ft_workers", None)   # results do not depend on the pool size (tested)
        return d


# v0 (sequence-only rules, coverage guard 0.5) is the base of the S3 arms of the toy grid;
# "default" is the current default (v0.2: min_cov 0 + local-tree mode at divergent nodes).
_D = MergeConfig(min_cov=0.5, tree_mode=False, tree_require_outgroup=False, tree_min_d=0.0)
ARMS = {   # named arms of the S3/S4 benchmark (docs/core_v0.md)
    "default": MergeConfig(),
    "v0_default": _D,
    "rbh": replace(_D, score="bits", tol=1.0, inparalogs=False, repair=False),
    "max": replace(_D, link="max", repair=False),
    "max+repair": replace(_D, link="max", repair=True, repair_ratio=0.5),
    "avg": replace(_D, link="avg", repair=False),
    "avg+repair": replace(_D, link="avg", repair=True, repair_ratio=0.5),
    "max+repair/bits": replace(_D, link="max", repair=True, score="bits", repair_ratio=0.5),
    "avg+repair/bits": replace(_D, link="avg", repair=True, score="bits", repair_ratio=0.5),
    "max+repair/quartet": replace(_D, link="max", repair=True, inpar_rule="quartet", repair_ratio=0.5),
    "avg+repair/quartet": replace(_D, link="avg", repair=True, inpar_rule="quartet", repair_ratio=0.5),
    "avg+repair/either": replace(_D, link="avg", repair=True, inpar_rule="either", repair_ratio=0.5),
    "rd_hidden": replace(_D, rd_hidden=True),
    "rd_inpar": replace(_D, rd_inpar=True),
    "rd_both": replace(_D, rd_hidden=True, rd_inpar=True),
    "tree": replace(_D, tree_mode=True, min_cov=0.0),
    "seq_ob": replace(_D, min_cov=0.0),
    "tree_ob": replace(_D, min_cov=0.0, tree_mode=True, tree_require_outgroup=True, tree_min_d=0.5),
    "tree+rd": replace(_D, tree_mode=True, rd_hidden=True, rd_inpar=True),
}


@dataclass
class NodeContext:
    """What a context-evidence function (synteny, S5) receives at one node X.

    All unit indices are *local* (0..n-1) to the node; a unit is an ancestral gene of one
    of X's children. The function returns a float array of length ``len(cand_a)``: an
    additive bonus (same scale as ``cand_score``; NaN = 0) applied before R1-R3.
    """
    node: str                        # species-tree node name X
    children: list | None            # child node names; child index = position in this list
    units: np.ndarray | None         # local unit -> global ancestral-gene ID
    child: np.ndarray                # local unit -> child index
    cand_a: np.ndarray               # cross-child candidate edges (local units, a in one child,
    cand_b: np.ndarray               #   b in another), every scored pair of units
    cand_score: np.ndarray           # sequence score of each candidate (avg/max link)
    anchors: np.ndarray              # bool per candidate: pass-1 anchor (strict mutual best, 1:1)
    adjacency: tuple | None          # (ua, ub, n_species, min_dist): within-child ancestral
                                     #   adjacencies between local units (genome.unit_adjacency,
                                     #   window w), None without gene order
    members: object = None           # callable: local unit -> array of leaf gene IDs


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
                cfg: MergeConfig, context=None, resolver=None, node_name=None, out_edges=None,
                node_info=None):
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
    ctx_support = None
    locus_sup = None
    child_l = child.tolist()
    cross = child[lu] != child[lv]
    ca, cb, cs = lu[cross], lv[cross], score[cross].astype(np.float64)
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
    uniq, inv = np.unique(dkey, return_inverse=True)

    def best_ab(sc):
        best = np.full(len(uniq), -np.inf)
        np.maximum.at(best, inv, np.concatenate([sc, sc]))
        return best[inv[:E]], best[inv[E:]]

    def mutual_best(sc, tol):
        best_a, best_b = best_ab(sc)
        if tol >= 1.0:
            m = (sc >= best_a) & (sc >= best_b)
        else:
            m = (sc >= tol * best_a) & (sc >= tol * best_b)
        return m & (sc > cfg.min_score)

    if context is not None:
        # three-pass context evidence (DESIGN S5): pass 1 = sequence-only anchors
        # (strict mutual best, 1:1), pass 2 = context function scores every candidate
        # against the anchors, pass 3 = the usual rules on score + bonus.
        strict = mutual_best(cs, 1.0)
        deg = np.bincount(np.concatenate([ca[strict], cb[strict]]), minlength=n)
        anchors = strict & (deg[ca] == 1) & (deg[cb] == 1)
        ctx = NodeContext(node=node_name, children=(node_info or {}).get("children"),
                          units=(node_info or {}).get("units"), child=child, cand_a=ca,
                          cand_b=cb, cand_score=cs.copy(), anchors=anchors,
                          adjacency=(node_info or {}).get("adjacency"),
                          members=(node_info or {}).get("members"))
        bonus = np.nan_to_num(np.asarray(context(ctx), dtype=np.float64))
        if bonus.shape != cs.shape:
            raise ValueError("context function must return one value per candidate edge")
        stats["n_anchors"] = int(anchors.sum())
        if (cfg.ctx_r3 or cfg.ctx_veto or cfg.locus_repair) and cfg.engine != "numba":
            raise NotImplementedError("ctx_r3 / ctx_veto are implemented in the numba engine only")
        if cfg.ctx_veto and hasattr(context, "support"):
            ctx_support = np.asarray(context.support(ctx), dtype=np.float64)
        if cfg.locus_repair and hasattr(context, "support"):
            fn = getattr(context, "locus_support", context.support)   # unmasked if constrained
            locus_sup = np.asarray(fn(ctx), dtype=np.float64)
        if cfg.ctx_mode == "all":
            cs = cs + bonus
            cs_sel = cs
        else:
            # R1-only: context may change the *partner choice* among candidates whose raw
            # sequence score is within ctx_window of the best for both endpoints; R1's
            # in-paralog constraint, R2, R3 and the RD tests keep the raw sequence scores.
            ba, bb = best_ab(cs)
            elig = (cs >= cfg.ctx_window * ba) & (cs >= cfg.ctx_window * bb)
            cs_sel = cs + np.where(elig, bonus, 0.0)
            stats["n_ctx_bonus"] = int(((bonus != 0) & elig).sum())
        cs_r3 = cs_sel if cfg.ctx_r3 else cs
    else:
        cs_sel = cs
        cs_r3 = cs
    mb = mutual_best(cs_sel, cfg.tol)
    stats["n_mb"] = int(mb.sum())
    if mb.any():   # node divergence: median distance -ln S of the mutual-best seeds
        stats["seed_median_d"] = float(np.median(-np.log(np.clip(cs[mb], 1e-6, 1.0))))
        if cfg.deep_min_d > 0 and stats["seed_median_d"] >= cfg.deep_min_d:
            stats["deep"] = 1
            if cfg.deep_tol > 0:     # relaxed mutual best at deep nodes only
                mb = mutual_best(cs_sel, cfg.deep_tol)
                stats["n_mb_deep"] = int(mb.sum())
    if cfg.engine == "numba":
        return _decide_numba(n, k, child, ca, cb, cs, sa, sb, ss, mb, cfg, out_edges, stats,
                             resolver, node_name, lu, lv, score, cs_sel, cs_r3, ctx_support, locus_sup,
                             (node_info or {}).get("unit_nsp"))
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

    # ---- relative divergence with the closest outgroup gene (RD tests)
    def _d(v):
        return -np.log(min(max(v, 1e-6), 1.0))

    def rel_div(a, b, s_ab):
        if not out_adj:
            return None
        oa, ob_ = out_adj[a], out_adj[b]
        if len(oa) > len(ob_):
            oa, ob_ = ob_, oa
        best_o, best_v = None, -1.0
        for z, va in oa.items():
            vb = ob_.get(z)
            if vb is not None and va + vb > best_v:
                best_o, best_v = z, va + vb
        if best_o is None:
            return None
        den = 0.5 * (_d(oa[best_o]) + _d(ob_[best_o]))
        return _d(s_ab) / den if den > 0 else None

    idx = np.flatnonzero(mb)
    idx = idx[np.argsort(-cs_sel[idx], kind="stable")]
    med_r = None
    veto = set()
    if (cfg.rd_hidden or cfg.rd_inpar) and out_adj:
        rs = {}
        for e in idx.tolist():
            r = rel_div(int(ca[e]), int(cb[e]), float(cs[e]))
            if r is not None:
                rs[e] = r
        if len(rs) >= 10:
            arr = np.fromiter(rs.values(), dtype=np.float64)
            med_r = float(np.median(arr))
            mad = float(np.median(np.abs(arr - med_r))) * 1.4826
            stats["rd_median"] = round(med_r, 4)
            stats["rd_mad"] = round(mad, 4)
            stats["rd_n"] = len(rs)
            if cfg.rd_hidden and mad > 0:
                thr = med_r + cfg.rd_z * mad
                veto = {e for e, r in rs.items() if r > thr}
    stats["r1_rd_veto"] = len(veto)

    # ---- R1: greedy union of mutual-best edges with in-paralog constraint
    for e in idx.tolist():
        if e in veto:
            continue
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
                    if cfg.rd_inpar and med_r is not None:
                        r = rel_div(u, w, sw)   # u-w split younger than X (rate-robust)
                        if r is not None and r <= cfg.rd_inpar_factor * med_r:
                            ok = True
                            stats["r2_rd_absorbed"] = stats.get("r2_rd_absorbed", 0) + 1
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


def _orphan_merge(n, child, lu, lv, labels, stats, score=None, min_score=0.0):
    """R4: units that are still alone and have no cross-child edge at this node are
    agglomerated with other such units of the same child by iterated *mutual best*
    (component similarity = best member edge), only through edges with score >=
    ``min_score``. Parsimony: one duplication inside the child instead of a duplication
    before X plus losses in every other child. Adds paralogGroups only (no ortholog pairs)."""
    cross = child[lu] != child[lv]
    has_cross = np.zeros(n, bool)
    has_cross[lu[cross]] = True
    has_cross[lv[cross]] = True
    alone = np.bincount(labels, minlength=n)[labels] == 1
    cand = alone & ~has_cross
    sc = np.ones(len(lu)) if score is None else score
    m = ~cross & cand[lu] & cand[lv] & (sc >= min_score)
    stats["r4_orphan_merged"] = 0
    if not m.any():
        return labels
    a, b, w = lu[m], lv[m], sc[m]
    comp = np.arange(n)
    merged = 0
    for _ in range(50):
        ca, cb = comp[a], comp[b]
        k = ca != cb
        if not k.any():
            break
        x, y, ww = ca[k], cb[k], w[k]
        best = np.full(n, -np.inf)
        np.maximum.at(best, x, ww)
        np.maximum.at(best, y, ww)
        # partner = component achieving the best score (ties: smallest id)
        part = np.full(n, n, dtype=np.int64)
        hx = ww >= best[x]
        np.minimum.at(part, x[hx], y[hx])
        hy = ww >= best[y]
        np.minimum.at(part, y[hy], x[hy])
        mutual = [(i, int(part[i])) for i in np.unique(np.r_[x, y])
                  if part[i] < n and part[part[i]] == i and i < part[i]]
        if not mutual:
            break
        for i, j in mutual:
            comp[comp == j] = i
            merged += 1
    stats["r4_orphan_merged"] = merged
    new = labels.copy()
    sel = cand & (np.bincount(comp, minlength=n)[comp] > 1)
    new[sel] = labels[comp[sel]]
    return new


def _so_absorb(n, child, lu, lv, score, labels, masks, stats):
    """Species-overlap absorption: a still-alone unit u joins the group of its best cross-child
    partner if every member of that group from u's child has a species set disjoint from u's
    (no duplication evidence -> missed speciation below X)."""
    labels = labels.copy()
    cross = child[lu] != child[lv]
    a, b, s = lu[cross], lv[cross], np.asarray(score, dtype=np.float64)[cross]
    gsize: dict = {}
    for l in labels.tolist():
        gsize[l] = gsize.get(l, 0) + 1
    best = {}
    for x, y, sc in zip(np.r_[a, b].tolist(), np.r_[b, a].tolist(), np.r_[s, s].tolist()):
        if gsize[int(labels[x])] == 1 and sc > best.get(x, (None, -1.0))[1]:
            best[x] = (y, sc)
    gmem: dict = {}
    for i, l in enumerate(labels.tolist()):
        gmem.setdefault(l, []).append(i)
    k = 0
    for u, (v, sc) in sorted(best.items(), key=lambda kv: -kv[1][1]):
        if gsize[int(labels[u])] != 1:
            continue
        g = int(labels[v])
        same = [w for w in gmem[g] if child[w] == child[u]]
        if not same or len({child[w] for w in gmem[g]}) < 2:
            continue
        if any(masks[w] & masks[u] for w in same):
            continue
        gsize[int(labels[u])] -= 1
        labels[u] = g
        gmem[g].append(u)
        gsize[g] += 1
        k += 1
    stats["so_absorbed"] = k
    return labels


def _locus_repair(n, child, ca, cb, cs, sup, labels, min_sup, stats, eligible=None):
    """Locus-aware fission repair: a unit that is still alone joins the group (lacking its
    child) to which it has the highest context support >= min_sup (ties: higher score),
    using any in-family-scored partner. Each group gains at most one unit per child."""
    labels = labels.copy()
    size = np.bincount(labels, minlength=int(labels.max()) + 1)
    alone = size[labels] == 1
    if eligible is not None:
        alone = alone & eligible
    m = (sup >= min_sup) & (alone[ca] | alone[cb])
    stats["locus_repaired"] = 0
    if not m.any():
        return labels
    a, b, sc, sp = ca[m], cb[m], cs[m], sup[m]
    u = np.where(alone[a], a, b)
    v = np.where(alone[a], b, a)
    order = np.lexsort((-sc, -sp))
    group_children = {}
    gsize = {}
    for i in range(n):
        group_children.setdefault(int(labels[i]), set()).add(int(child[i]))
        gsize[int(labels[i])] = gsize.get(int(labels[i]), 0) + 1
    done = set()
    for k in order.tolist():
        uu, vv = int(u[k]), int(v[k])
        if uu in done or labels[uu] == labels[vv] or gsize[int(labels[uu])] > 1:
            continue
        g = int(labels[vv])
        cu = int(child[uu])
        if cu in group_children[g]:
            continue
        gsize[int(labels[uu])] -= 1
        labels[uu] = g
        gsize[g] += 1
        group_children[g].add(cu)
        done.add(uu)
    stats["locus_repaired"] = len(done)
    return labels


_RULE_CODE = {"ortholog": 0, "self": 1, "either": 2, "quartet": 3}


def _decide_numba(n, k, child, ca, cb, cs, sa, sb, ss, mb, cfg, out_edges, stats,
                  resolver, node_name, lu, lv, score, cs_sel=None, cs_r3=None, ctx_support=None,
                  locus_sup=None, unit_nsp=None):
    from . import merge_numba as mn
    child = np.ascontiguousarray(child, dtype=np.int64)
    ca = np.ascontiguousarray(ca, dtype=np.int64)
    cb = np.ascontiguousarray(cb, dtype=np.int64)
    cs = np.ascontiguousarray(cs, dtype=np.float64)
    sa = np.ascontiguousarray(sa, dtype=np.int64)
    sb = np.ascontiguousarray(sb, dtype=np.int64)
    ss = np.ascontiguousarray(ss, dtype=np.float64)
    idx = np.flatnonzero(mb)
    order_sc = cs if cs_sel is None else cs_sel
    idx = idx[np.argsort(-order_sc[idx], kind="stable")].astype(np.int64)
    has_out = out_edges is not None and len(out_edges[0]) > 0
    if has_out:
        optr, onbr, osc = mn.csr_directed_sorted(n, out_edges[0].astype(np.int64),
                                                 out_edges[1].astype(np.int64),
                                                 out_edges[2].astype(np.float64))
    else:
        optr, onbr, osc = np.zeros(n + 1, np.int64), np.zeros(0, np.int64), np.zeros(0)
    med_r = np.nan
    veto = np.zeros(len(idx), np.bool_)
    if (cfg.rd_hidden or cfg.rd_inpar) and has_out:
        r = mn.rd_values(optr, onbr, osc, ca, cb, cs, idx)
        ok = ~np.isnan(r)
        if ok.sum() >= 10:
            med_r = float(np.median(r[ok]))
            mad = float(np.median(np.abs(r[ok] - med_r))) * 1.4826
            stats.update(rd_median=round(med_r, 4), rd_mad=round(mad, 4), rd_n=int(ok.sum()))
            if cfg.rd_hidden and mad > 0:
                veto = ok & (r > med_r + cfg.rd_z * mad)
    stats["r1_rd_veto"] = int(veto.sum())
    labels, cnt = mn.rules(n, k, child, ca, cb, cs, sa, sb, ss, idx, veto, float(cfg.min_score),
                           bool(cfg.inparalogs), float(cfg.inpar_factor), _RULE_CODE[cfg.inpar_rule],
                           bool(cfg.rd_inpar), float(cfg.rd_inpar_factor), med_r,
                           bool(cfg.repair), float(cfg.repair_ratio), optr, onbr, osc, has_out,
                           np.ascontiguousarray(cs if cs_r3 is None else cs_r3, dtype=np.float64),
                           np.ascontiguousarray(np.zeros(len(cs)) if ctx_support is None else ctx_support,
                                                dtype=np.float64),
                           bool(cfg.ctx_veto and ctx_support is not None), float(cfg.ctx_veto_min))
    if cfg.locus_repair and locus_sup is not None:
        eligible = (unit_nsp == 1) if (cfg.locus_single_species and unit_nsp is not None) else None
        labels = _locus_repair(n, child, ca, cb, cs, locus_sup, labels, cfg.locus_min_support, stats,
                               eligible)
    stats.update(r1_unions=int(cnt[0]), r1_blocked=int(cnt[1]), r2_absorbed=int(cnt[2]),
                 r2_rd_absorbed=int(cnt[3]), r2_quartet_tests=int(cnt[4]), r3_repaired=int(cnt[5]))
    if resolver is not None:
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
              context=None, resolver=None, log=None, trace: dict | None = None,
              gene_order=None, adjacency_w: int = 5) -> MergeResult:
    """Run the bottom-up merge over all families at once.

    Hooks: ``context(ctx: NodeContext) -> bonus array`` (added to the cross-child
    candidate scores, three-pass: anchors -> context -> rules; see NodeContext);
    ``resolver(node, child, lu, lv, score, labels) -> labels``. ``gene_order``
    (:class:`ancestra.genome.GeneOrder`) enables ``NodeContext.adjacency`` (window
    ``adjacency_w``).
    ``trace`` (dict, diagnostics): filled with node -> (clade genes, their ancestral
    gene ID at the children, at the node).
    """
    if (cfg.rd_hidden or cfg.rd_inpar) and cfg.score not in ("norm_max", "norm_min"):
        raise ValueError("relative-divergence tests need a normalised score (norm_max/norm_min)")
    N = gt.n
    cap = 2 * N + 2
    nodes, nid, lca, clade_sp = _lca_matrix(tree, gt.species)
    sp_genes = gt.genes_of_species()
    npairs = None
    if gene_order is not None and context is not None:
        from .genome import neighbor_pairs
        npairs = neighbor_pairs(gene_order, adjacency_w)
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
    amask: list = ([1 << int(x) for x in gt.species_of]
                   if (cfg.so_struct or cfg.so_absorb or cfg.deep_so_absorb) else None)
    if amask is not None:
        from .newick import TreeIndex
        _tix = TreeIndex(tree)
        _spn = list(gt.species)

        def mask_mrca(m: int) -> str:
            return _tix.mrca(_spn[i] for i in range(len(_spn)) if (m >> i) & 1)

        _sidx = {nm: i for i, nm in enumerate(_spn)}
        _cmask: dict = {}                     # node -> species mask of its clade
        for _nd in tree.postorder():
            _cmask[_nd.name] = (1 << _sidx[_nd.name]) if not _nd.children else \
                sum(_cmask[c.name] for c in _nd.children)
        _kids_of = {nd_.name: [c.name for c in nd_.children] for nd_ in tree.preorder()}

        def sp_mask(x) -> int:
            m, stack = 0, [x]
            while stack:
                y = stack.pop()
                if isinstance(y, (int, np.integer)):
                    m |= 1 << int(gt.species_of[int(y)])
                else:
                    stack.extend(y[2])
            return m

        def so_combine(items):
            """Species-disjoint structures -> one speciation structure, correctly nested:
            orthologGroup at the MRCA L; same-level orthologGroups are flattened, the
            pieces are bucketed by the child of L they lie in and combined recursively.
            Returns None if an item cannot lie strictly below L (caller keeps paralogs)."""
            ms = [sp_mask(x) for x in items]
            tot = 0
            for m_ in ms:
                tot |= m_
            L = mask_mrca(tot)
            flat = []
            for x, xm in zip(items, ms):
                if not isinstance(x, (int, np.integer)) and x[1] == L:
                    if x[0] == "O":
                        flat.extend((c, sp_mask(c)) for c in x[2])
                    else:              # paralogGroup at L: place at its members' MRCA if below L
                        lvl = mask_mrca(xm)
                        if lvl == L:
                            return None
                        flat.append((("P", lvl, x[2]), xm))
                else:
                    flat.append((x, xm))
            buckets: dict = {}
            for x, xm in flat:
                c = next((c for c in _kids_of[L] if not (xm & ~_cmask[c])), None)
                if c is None:
                    return None
                buckets.setdefault(c, []).append(x)
            out = []
            for c in _kids_of[L]:
                xs = buckets.get(c)
                if not xs:
                    continue
                if len(xs) == 1:
                    out.append(xs[0])
                else:
                    r = so_combine(xs)
                    if r is None:
                        return None
                    out.append(r)
            if len(out) == 1:
                return out[0]
            fixed = []
            for x in out:   # a paralogGroup re-bucketed below L must lie strictly below L:
                if not isinstance(x, (int, np.integer)) and x[0] == "P" and \
                        not _tix.is_strict_descendant(x[1], L):
                    x = ("P", mask_mrca(sp_mask(x)), x[2])    # latest-possible placement
                fixed.append(x)
            return ("O", L, fixed)
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
        if (cfg.inpar_rule == "quartet" or cfg.rd_hidden or cfg.rd_inpar or cfg.tree_mode) \
                and nd.parent is not None:
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
        node_info = None
        if context is not None or (cfg.tree_mode and cfg.tree_method == "fasttree"):
            node_info = {"children": kids, "units": units.copy()}
            # members: local unit -> leaf genes (gene_anc is still at the children level)
            ga_local = pos[gene_anc[genes]]
            o_m = np.argsort(ga_local, kind="stable")
            b_m = np.searchsorted(ga_local[o_m], np.arange(n + 1))
            g_sorted = genes[o_m]
            node_info["members"] = lambda i, _g=g_sorted, _b=b_m: _g[_b[i]:_b[i + 1]]
            if cfg.locus_repair:
                sp_sorted = gt.species_of[g_sorted].astype(np.int64)
                nsp = np.zeros(n, dtype=np.int64)
                for i in range(n):
                    seg = sp_sorted[b_m[i]:b_m[i + 1]]
                    nsp[i] = len(np.unique(seg)) if len(seg) else 0
                node_info["unit_nsp"] = nsp
            if npairs is not None:
                from .genome import unit_adjacency
                unit_of_gene = np.full(N, -1, dtype=np.int64)
                unit_of_gene[genes] = ga_local
                node_info["adjacency"] = unit_adjacency(npairs[0], npairs[1], npairs[2],
                                                        unit_of_gene, gt.species_of, n)
        labels, stt = decide_node(n, child, lu, lv, score, cfg, context, resolver, nd.name,
                                  out_edges, node_info)
        if cfg.tree_mode and (not cfg.tree_nodes or nd.name in cfg.tree_nodes.split(",")) \
                and stt.get("seed_median_d", 0.0) >= cfg.tree_min_d \
                and (not cfg.tree_require_outgroup or out_edges is not None):
            from .localtree import resolve_ambiguous
            t_tree = time.time()
            ft = cfg.tree_method == "fasttree"
            if ft and gt.seqs is None:
                raise ValueError("tree_method='fasttree' needs sequences (GeneTable.seqs)")
            labels = resolve_ambiguous(n, child, lu, lv, score, labels, out_edges, cfg.tree_tau,
                                       cfg.tree_max_units, stt, cfg.tree_require_outgroup,
                                       cfg.tree_conservative, cfg.tree_method,
                                       node_info["members"] if ft else None,
                                       gt.seqs if ft else None, cfg.tree_ft_mode, 60,
                                       cfg.tree_ft_workers)
            stt["tree_s"] = round(time.time() - t_tree, 3)
        if cfg.so_absorb or (cfg.deep_so_absorb and stt.get("deep")):
            labels = _so_absorb(n, child, lu, lv, score, labels, [amask[int(u)] for u in units], stt)
        if cfg.orphan_merge == "all" or (cfg.orphan_merge == "root" and nd.parent is None):
            labels = _orphan_merge(n, child, lu, lv, labels, stt, score, cfg.orphan_min)
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
                elif cfg.so_struct:
                    # species overlap: species-disjoint units are a missed speciation
                    clus = []          # [mask, [units]]
                    for x in sorted(us, key=lambda y: (-bin(amask[y]).count("1"), y)):
                        for cl in clus:
                            if not (cl[0] & amask[x]):
                                cl[0] |= amask[x]
                                cl[1].append(x)
                                break
                        else:
                            clus.append([amask[x], [x]])
                    sub = []
                    for cl in clus:
                        if len(cl[1]) == 1:
                            sub.append(struct[cl[1][0]])
                            continue
                        r = so_combine([struct[x] for x in cl[1]])
                        if r is None:          # cannot nest correctly: keep them as paralogs
                            sub.extend(struct[x] for x in cl[1])
                            stt["so_struct_fallback"] = stt.get("so_struct_fallback", 0) + 1
                        else:
                            sub.append(r)
                            stt["so_struct_merged"] = stt.get("so_struct_merged", 0) + len(cl[1]) - 1
                    parts_s.append(sub[0] if len(sub) == 1 else ("P", kids[ci], sub))
                else:
                    parts_s.append(("P", kids[ci], [struct[x] for x in us]))
            aid = next_id
            next_id += 1
            if len(parts_s) == 1:
                struct.append(parts_s[0])
            else:
                struct.append(("O", nd.name, parts_s))
            size[aid] = size[units[grp]].sum()
            if amask is not None:
                mm = 0
                for i in grp.tolist():
                    mm |= amask[int(units[i])]
                amask.append(mm)
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


def root_labels(res: "MergeResult", n_genes: int) -> np.ndarray:
    """gene -> index of its root ancestral gene (root HOG) in ``res.roots``."""
    lab = np.full(n_genes, -1, dtype=np.int64)
    for k, x in enumerate(res.roots):
        stack = [x]
        while stack:
            y = stack.pop()
            if isinstance(y, (int, np.integer)):
                lab[int(y)] = k
            else:
                stack.extend(y[2])
    return lab


class ConstrainedContext:
    """Two-pass constrained context evidence: the wrapped context function (e.g.
    ``ancestra.synteny.SyntenyContext``) acts only between units that lie in the same root
    HOG of a context-off run (``root_lab``: gene -> root HOG); bonus and support of other
    candidates are 0. OrthoBench dev: keeps synteny's pair-level gain without its OG-level
    fusions (docs/core_v0.md §7)."""

    def __init__(self, base, root_lab: np.ndarray):
        self.base, self.root_lab = base, root_lab

    def _same(self, ctx):
        lab = np.array([self.root_lab[m[0]] if len(m) else -1
                        for m in (ctx.members(i) for i in range(len(ctx.child)))])
        return lab[ctx.cand_a] == lab[ctx.cand_b]

    def __call__(self, ctx):
        return np.asarray(self.base(ctx), dtype=np.float64) * self._same(ctx)

    def support(self, ctx):
        return np.asarray(self.base.support(ctx), dtype=np.float64) * self._same(ctx)

    def locus_support(self, ctx):
        """Unmasked support (locus-aware fission repair may join across off-mode root HOGs)."""
        return np.asarray(self.base.support(ctx), dtype=np.float64)


def run_merge_constrained(gt, tree, st, cfg, context, off_cfg=None, **kw) -> "MergeResult":
    """Context-off merge (``off_cfg``, default ``cfg``) -> root HOGs -> merge with
    ``ConstrainedContext(context)`` under ``cfg``."""
    off = run_merge(gt, tree, st, off_cfg or cfg, **kw)
    return run_merge(gt, tree, st, cfg, context=ConstrainedContext(context, root_labels(off, gt.n)), **kw)


def _to_group(x, names):
    if isinstance(x, (int, np.integer)):
        return names[int(x)]
    kind, level, kids = x
    return Group("ortholog" if kind == "O" else "paralog",
                 [_to_group(c, names) for c in kids], level=level)


def _fix_paralog_levels(g: Group, sp: dict, tix) -> int:
    """Level consistency of paralogGroups (pairs unchanged): a paralogGroup directly inside a
    paralogGroup is flattened into it; a paralogGroup whose level does not contain its members
    or is not strictly below the enclosing orthologGroup gets its members' MRCA (latest-possible
    placement). Returns the number of changes."""
    n = 0
    stack = [(g, None)]               # (group, level of the enclosing orthologGroup)
    while stack:
        x, enc = stack.pop()
        if x.kind == "paralog":
            # flatten nested paralogGroups (repeat until none left)
            changed = True
            while changed:
                changed = False
                kids = []
                for c in x.children:
                    if isinstance(c, Group) and c.kind == "paralog":
                        kids.extend(c.children)
                        changed = True
                        n += 1
                    else:
                        kids.append(c)
                x.children = kids
            mem = {sp[y] for y in x.iter_genes()}
            ok = x.level in tix and all(tix.is_ancestor(x.level, m) for m in mem) and \
                (enc is None or tix.is_strict_descendant(x.level, enc))
            if not ok:
                x.level = tix.mrca(mem)
                n += 1
            sub_enc = enc
        else:
            sub_enc = x.level
        stack.extend((c, sub_enc) for c in x.children if isinstance(c, Group))
    return n


def to_orthoxml(res: MergeResult, gt: GeneTable, tree: Node) -> OrthoXML:
    from .newick import TreeIndex
    genes = {nm: GeneEntry(nm, gt.species[gt.species_of[g]]) for g, nm in enumerate(gt.names)}
    groups = [_to_group(x, gt.names) for x in res.roots if not isinstance(x, (int, np.integer))]
    groups = normalize_groups(groups)
    tix = TreeIndex(tree)
    sp = {nm: e.species for nm, e in genes.items()}
    for g in groups:
        _fix_paralog_levels(g, sp, tix)
    return OrthoXML(genes=genes, groups=normalize_groups(groups), species_tree=tree)
