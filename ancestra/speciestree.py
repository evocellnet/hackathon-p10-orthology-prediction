"""S1 -- guide species trees (DESIGN.md v2, section 3 row S1).

Arms (all produce a *rooted* Newick tree whose leaves are the proteome species
names, i.e. the FASTA file stems used by :mod:`ancestra.io`, and whose internal
nodes are named):

* **ncbi**    -- induced subtree of the NCBI taxdump (multifurcating; unary
  nodes collapsed; internal nodes named by NCBI scientific names, deduplicated).
* **A / besthit-NJ** -- a fixed random sample of proteins per species (seed 11)
  is searched with DIAMOND against each proteome separately; species distance =
  -ln(mean best-hit identity), symmetrised; BIONJ (+ optional FastME SPR);
  rooted on the edge that best agrees with the NCBI tree's clades (or midpoint).
* **B / ASTRAL-Pro** -- near-single-copy families (from S0 families or built
  from the sampled best hits) -> FAMSA -> FastTree gene trees -> ASTRAL-Pro3,
  plus an optional concatenation FastTree tree.
* **perturbed** -- k random NNIs (seeded); **random** -- uniform random rooted
  binary topology (seeded, negative control).

Data-derived trees get internal names from the NCBI tree: a node whose leaf set
equals an NCBI clade gets that clade's name, any other node gets
``<LCA name>~<k>`` (synthetic, dissolved on output by the pipeline, DESIGN §8).

Comparison utilities (``compare_trees``): unrooted nRF, recall of reference
splits, fraction of query splits *incompatible* with the reference (so that the
refinement of a reference polytomy is not counted as an error), and quartet
statistics.

Leakage (DESIGN §4): this module only reads proteomes, S0 families and the NCBI
taxdump.  Reference species trees from benchmark data are loaded by
``scripts/S1_evaluate.py`` for evaluation only.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .io import read_fasta, species_name
from .newick import Node, parse_newick, read_newick

__all__ = [
    "load_taxdump", "lookup_taxids_by_name", "ncbi_tree", "clades", "name_internal_nodes",
    "splits", "compare_trees", "quartet_compare", "nj", "midpoint_root", "root_by_reference",
    "reroot", "refine_polytomies", "nni_perturb", "random_topology", "prune_to", "collapse_unary",
    "besthit_distances", "sample_queries", "write_tree",
]


# ---------------------------------------------------------------------------
# generic tree helpers
# ---------------------------------------------------------------------------

def collapse_unary(root: Node) -> Node:
    """Remove nodes with exactly one child (lengths are added). Returns the root."""
    while len(root.children) == 1:
        ch = root.children[0]
        ch.parent = None
        root = ch
    for nd in list(root.postorder()):
        if nd.parent is not None and len(nd.children) == 1:
            ch, p = nd.children[0], nd.parent
            i = p.children.index(nd)
            p.children[i] = ch
            ch.parent = p
            if ch.length is not None and nd.length is not None:
                ch.length += nd.length
            elif nd.length is not None:
                ch.length = nd.length
    return root


def prune_to(root: Node, keep) -> Node:
    """Copy of the tree restricted to leaves in ``keep`` (unary nodes collapsed)."""
    keep = set(keep)
    new = copy_tree(root)
    for nd in list(new.postorder()):
        if nd.parent is None:
            continue
        if not nd.children and nd.name not in keep:
            nd.parent.children.remove(nd)
            nd.parent = None
    # drop internal nodes without leaves
    changed = True
    while changed:
        changed = False
        for nd in list(new.postorder()):
            if nd.parent is not None and not nd.children and nd.name not in keep:
                nd.parent.children.remove(nd)
                nd.parent = None
                changed = True
    return collapse_unary(new)


def copy_tree(root: Node) -> Node:
    mp = {}
    for nd in root.preorder():
        c = Node(nd.name, nd.length)
        c.data = dict(nd.data)
        mp[id(nd)] = c
        if nd.parent is not None:
            mp[id(nd.parent)].add_child(c)
    return mp[id(root)]


def write_tree(root: Node, path: str, lengths: bool = True) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(root.to_newick(lengths=lengths) + "\n")
    os.replace(tmp, path)


def clades(root: Node) -> dict[frozenset, str]:
    """Map leaf-set -> node name for every internal node."""
    out = {}
    below: dict[int, frozenset] = {}
    for nd in root.postorder():
        if not nd.children:
            below[id(nd)] = frozenset([nd.name])
        else:
            s = frozenset().union(*(below[id(c)] for c in nd.children))
            below[id(nd)] = s
            out.setdefault(s, nd.name)
    return out


def _leafsets(root: Node) -> dict[int, frozenset]:
    below: dict[int, frozenset] = {}
    for nd in root.postorder():
        below[id(nd)] = (frozenset([nd.name]) if not nd.children else
                         frozenset().union(*(below[id(c)] for c in nd.children)))
    return below


def name_internal_nodes(root: Node, ref: Node | None, prefix: str = "n") -> Node:
    """Name internal nodes after the reference (NCBI) clades.

    Exact leaf-set match -> reference name; otherwise ``<LCA name>~<k>`` where
    LCA is the smallest reference clade containing the node's leaves. Sets
    ``node.data['ref_match']``. Without a reference: ``<prefix><k>``.
    """
    below = _leafsets(root)
    refc = clades(ref) if ref is not None else {}
    # reference clades sorted by size for LCA lookup
    ref_items = sorted(refc.items(), key=lambda kv: len(kv[0]))
    used: set[str] = {n.name for n in root.leaves()}
    k = 0
    for nd in root.preorder():
        if not nd.children:
            continue
        s = below[id(nd)]
        name = refc.get(s)
        if name is not None and name not in used:
            nd.name = name
            nd.data["ref_match"] = True
        else:
            lca = None
            for cs, cn in ref_items:
                if s <= cs:
                    lca = cn
                    break
            k += 1
            base = f"{lca}~{k}" if lca else f"{prefix}{k}"
            while base in used:
                k += 1
                base = f"{lca}~{k}" if lca else f"{prefix}{k}"
            nd.name = base
            nd.data["ref_match"] = False
        used.add(nd.name)
    return root


# ---------------------------------------------------------------------------
# NCBI taxonomy arm
# ---------------------------------------------------------------------------

def load_taxdump(taxdump_dir: str, need_names_for: set[int] | None = None):
    """Parse nodes.dmp (+ merged.dmp) and scientific names from names.dmp.

    Returns (parent: dict[int,int], rank: dict[int,str], name: dict[int,str],
    merged: dict[int,int]). If ``need_names_for`` is given, only those names are
    kept (memory)."""
    parent, rank = {}, {}
    with open(os.path.join(taxdump_dir, "nodes.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t", 3)
            t = int(f[0])
            parent[t] = int(f[1])
            rank[t] = f[2]
    merged = {}
    mp = os.path.join(taxdump_dir, "merged.dmp")
    if os.path.exists(mp):
        with open(mp) as fh:
            for line in fh:
                f = line.split("\t|\t")
                merged[int(f[0])] = int(f[1].split("\t")[0])
    name = {}
    with open(os.path.join(taxdump_dir, "names.dmp")) as fh:
        for line in fh:
            if "scientific name" not in line:
                continue
            f = line.split("\t|\t")
            if f[3].startswith("scientific name"):
                t = int(f[0])
                if need_names_for is None or t in need_names_for:
                    name[t] = f[1]
    return parent, rank, name, merged


def lookup_taxids_by_name(taxdump_dir: str, names: list[str]) -> dict[str, int]:
    """Map organism names (e.g. 'Canis familiaris') to taxids (scientific
    names preferred over synonyms; case-insensitive)."""
    want = {n.lower(): n for n in names}
    found: dict[str, tuple[int, int]] = {}
    with open(os.path.join(taxdump_dir, "names.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t")
            key = f[1].lower()
            if key in want:
                pri = 0 if f[3].startswith("scientific name") else 1
                cur = found.get(want[key])
                if cur is None or pri < cur[1]:
                    found[want[key]] = (int(f[0]), pri)
    return {k: v[0] for k, v in found.items()}


def ncbi_tree(leaf_taxid: dict[str, int], taxdump_dir: str | None = None, tax=None) -> Node:
    """Induced NCBI taxonomy tree for the given leaves (name -> taxid).

    Keeps only nodes with >= 2 children (unary collapsed); internal nodes named
    by NCBI scientific name; a duplicated name gets ``_<taxid>`` appended.
    ``node.data['taxid']`` and ``['rank']`` are set."""
    if tax is None:
        tax = load_taxdump(taxdump_dir)
    parent, rank, name, merged = tax
    lineage = {}
    for leaf, t in leaf_taxid.items():
        t = merged.get(t, t)
        if t not in parent:
            raise ValueError(f"taxid {t} of {leaf} not in taxdump")
        path = [t]
        while path[-1] != 1 and parent[path[-1]] != path[-1]:
            path.append(parent[path[-1]])
        lineage[leaf] = path[::-1]          # root -> taxon
    nodes: dict[int, Node] = {}
    root = Node("root")
    root.data["taxid"] = 1
    nodes[1] = root
    for leaf, path in lineage.items():
        cur = root
        for t in path[1:]:
            nd = nodes.get(t)
            if nd is None:
                nd = cur.add_child(Node(name.get(t, str(t))))
                nd.data["taxid"] = t
                nd.data["rank"] = rank.get(t, "")
                nodes[t] = nd
            cur = nd
        lf = cur.add_child(Node(leaf))
        lf.data["taxid"] = path[-1]
        lf.data["rank"] = "leaf"
    root = collapse_unary(root)
    # the leaf's own taxon node may now be an internal node with >=2 children; fine.
    seen: Counter = Counter()
    for nd in root.preorder():
        if nd.children:
            seen[nd.name] += 1
    used = set(root.leaf_names())
    for nd in root.preorder():
        if nd.children:
            if seen[nd.name] > 1 or nd.name in used:
                nd.name = f"{nd.name}_{nd.data.get('taxid')}"
            used.add(nd.name)
    return root


def qfo_leaf_taxids(fasta_files: list[str]) -> dict[str, int]:
    """QfO file stems are UPxxxx_<taxid>."""
    out = {}
    for f in fasta_files:
        sp = species_name(f)
        out[sp] = int(sp.rsplit("_", 1)[1])
    return out


def ensembl_leaf_names(fasta_files: list[str]) -> dict[str, str]:
    """Ensembl pep file stem 'Homo_sapiens.GRCh38.pep.all' -> 'Homo sapiens'."""
    return {species_name(f): species_name(f).split(".")[0].replace("_", " ") for f in fasta_files}


# ---------------------------------------------------------------------------
# splits, nRF, quartets
# ---------------------------------------------------------------------------

def splits(root: Node, leaves: list[str]) -> set[int]:
    """Non-trivial unrooted bipartitions as bitmasks normalised to exclude leaf 0."""
    idx = {l: i for i, l in enumerate(leaves)}
    n = len(leaves)
    full = (1 << n) - 1
    out = set()
    mask: dict[int, int] = {}
    for nd in root.postorder():
        if not nd.children:
            m = 1 << idx[nd.name] if nd.name in idx else 0
        else:
            m = 0
            for c in nd.children:
                m |= mask[id(c)]
        mask[id(nd)] = m
        if m & 1:
            m = full ^ m
        c = bin(m).count("1")
        if 2 <= c <= n - 2:
            out.add(m)
    return out


def _compatible(a: int, b: int, full: int) -> bool:
    return (a & b) == 0 or (a & ~b & full) == 0 or (b & ~a & full) == 0 or ((full ^ a) & (full ^ b)) == 0


def _leaf_dist(root: Node, leaves: list[str]) -> np.ndarray:
    """Topological (unit-edge) leaf-to-leaf distances of the tree restricted to ``leaves``."""
    t = prune_to(root, leaves)
    idx = {l: i for i, l in enumerate(leaves)}
    n = len(leaves)
    depth = {}
    for nd in t.preorder():
        depth[id(nd)] = 0 if nd.parent is None else depth[id(nd.parent)] + 1
    D = np.zeros((n, n), dtype=np.int32)
    # for each internal node, leaves below via postorder; LCA depth via set merging
    below: dict[int, list[int]] = {}
    dep = np.zeros(n, dtype=np.int32)
    for nd in t.postorder():
        if not nd.children:
            below[id(nd)] = [idx[nd.name]]
            dep[idx[nd.name]] = depth[id(nd)]
        else:
            groups = [below.pop(id(c)) for c in nd.children]
            for a, b in itertools.combinations(range(len(groups)), 2):
                ga, gb = np.array(groups[a]), np.array(groups[b])
                D[np.ix_(ga, gb)] = -2 * depth[id(nd)]
            below[id(nd)] = [x for g in groups for x in g]
    D = D + D.T
    D = D + dep[:, None] + dep[None, :]
    np.fill_diagonal(D, 0)
    return D


def _quartet_topologies(D: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """0: ab|cd, 1: ac|bd, 2: ad|bc, -1: unresolved (four-point condition)."""
    a, b, c, d = quads.T
    s = np.stack([D[a, b] + D[c, d], D[a, c] + D[b, d], D[a, d] + D[b, c]], axis=1)
    mn = s.min(axis=1)
    arg = s.argmin(axis=1)
    nmin = (s == mn[:, None]).sum(axis=1)
    arg[nmin > 1] = -1
    return arg


def quartet_compare(t1: Node, t2: Node, leaves: list[str], max_quartets: int = 3_000_000,
                    seed: int = 11) -> dict:
    """Quartet agreement between two trees on ``leaves`` (all quartets if
    C(n,4) <= max_quartets, else a seeded random sample)."""
    n = len(leaves)
    if n < 4:
        return {"n_quartets": 0}
    D1, D2 = _leaf_dist(t1, leaves), _leaf_dist(t2, leaves)
    total = math.comb(n, 4)
    if total <= max_quartets:
        quads = np.array(list(itertools.combinations(range(n), 4)), dtype=np.int32)
        sampled = False
    else:
        rng = np.random.default_rng(seed)
        quads = np.sort(np.array([rng.choice(n, 4, replace=False) for _ in range(max_quartets)],
                                 dtype=np.int32), axis=1)
        sampled = True
    q1, q2 = _quartet_topologies(D1, quads), _quartet_topologies(D2, quads)
    both = (q1 >= 0) & (q2 >= 0)
    differ = both & (q1 != q2)
    return {
        "n_quartets": int(len(quads)), "sampled": sampled,
        "resolved_t1": float((q1 >= 0).mean()), "resolved_t2": float((q2 >= 0).mean()),
        "n_both_resolved": int(both.sum()),
        # quartet distance restricted to quartets resolved in both trees
        "qd_resolved": float(differ.sum() / max(1, both.sum())),
        # classic normalised quartet distance (unresolved vs resolved counts as different)
        "qd_all": float((q1 != q2).mean()),
    }


def compare_trees(query: Node, ref: Node, quartets: bool = True, leaves=None) -> dict:
    """Unrooted comparison of ``query`` vs ``ref`` on common leaves.

    nRF = (|Q\\R| + |R\\Q|) / (|Q| + |R|) (non-trivial splits; polytomies give
    fewer splits). ``ref_split_recall`` = |Q∩R|/|R| (resolved reference splits
    recovered). ``conflict_frac`` = share of query splits incompatible with at
    least one reference split (refining a reference polytomy is not a conflict).
    """
    common = sorted(set(query.leaf_names()) & set(ref.leaf_names()))
    if leaves is not None:
        common = sorted(set(common) & set(leaves))
    n = len(common)
    out = {"n_leaves": n}
    if n < 4:
        return out
    q, r = prune_to(query, common), prune_to(ref, common)
    sq, sr = splits(q, common), splits(r, common)
    full = (1 << n) - 1
    inter = len(sq & sr)
    rf = len(sq - sr) + len(sr - sq)
    conf = sum(1 for a in sq if any(not _compatible(a, b, full) for b in sr))
    rconf = sum(1 for b in sr if any(not _compatible(a, b, full) for a in sq))
    out.update({
        "n_splits_query": len(sq), "n_splits_ref": len(sr), "max_splits": n - 3,
        "rf": rf, "nrf": rf / max(1, len(sq) + len(sr)),
        "nrf_binary_norm": rf / (2 * (n - 3)) if n > 3 else 0.0,
        "ref_split_recall": inter / max(1, len(sr)),
        "query_split_precision": inter / max(1, len(sq)),
        "conflict_frac": conf / max(1, len(sq)),
        "n_conflicting_query_splits": conf,
        "ref_splits_contradicted": rconf / max(1, len(sr)),
    })
    # rooted clades (both trees taken as rooted; only meaningful for rooted inputs)
    cq = {c for c in clades(q) if 1 < len(c) < n}
    cr = {c for c in clades(r) if 1 < len(c) < n}
    out["rooted_clade_recall"] = len(cq & cr) / max(1, len(cr))
    out["root_split"] = sorted(len(c.leaves()) for c in q.children)
    if quartets:
        out["quartet"] = quartet_compare(q, r, common)
    return out


def conflicting_splits(query: Node, ref: Node) -> list[list[str]]:
    """Query splits incompatible with the reference, each as the sorted leaf
    names of its smaller side."""
    common = sorted(set(query.leaf_names()) & set(ref.leaf_names()))
    n = len(common)
    q, r = prune_to(query, common), prune_to(ref, common)
    sq, sr = splits(q, common), splits(r, common)
    full = (1 << n) - 1
    out = []
    for a in sq:
        if any(not _compatible(a, b, full) for b in sr):
            side = a if bin(a).count("1") <= n / 2 else full ^ a
            out.append(sorted(common[i] for i in range(n) if side >> i & 1))
    return sorted(out, key=lambda x: (len(x), x))


# ---------------------------------------------------------------------------
# distance trees, rooting
# ---------------------------------------------------------------------------

def nj(D: np.ndarray, names: list[str]) -> Node:
    """Plain neighbour joining (Saitou & Nei 1987); returns an unrooted tree as a
    root with 3 children. O(n^3); fine for n <= a few hundred."""
    D = np.array(D, dtype=float)
    nodes = [Node(n) for n in names]
    active = list(range(len(names)))
    while len(active) > 3:
        m = len(active)
        sub = D[np.ix_(active, active)]
        r = sub.sum(axis=1)
        Q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        a, b = active[i], active[j]
        li = 0.5 * sub[i, j] + (r[i] - r[j]) / (2 * (m - 2))
        lj = sub[i, j] - li
        u = Node()
        na, nb = nodes[a], nodes[b]
        na.length, nb.length = max(li, 0.0), max(lj, 0.0)
        u.add_child(na)
        u.add_child(nb)
        newrow = 0.5 * (D[a] + D[b] - D[a, b])
        D = np.vstack([D, newrow[None, :]])
        D = np.hstack([D, np.append(newrow, 0.0)[:, None]])
        nodes.append(u)
        active = [x for x in active if x not in (a, b)] + [len(nodes) - 1]
    root = Node()
    if len(active) == 3:
        a, b, c = active
        la = 0.5 * (D[a, b] + D[a, c] - D[b, c])
        lb = D[a, b] - la
        lc = D[a, c] - la
        for x, l in ((a, la), (b, lb), (c, lc)):
            nodes[x].length = max(l, 0.0)
            root.add_child(nodes[x])
    else:
        for x in active:
            nodes[x].length = D[active[0], active[-1]] / 2
            root.add_child(nodes[x])
    return root


def write_phylip_dist(D: np.ndarray, names: list[str], path: str) -> dict[str, str]:
    """Write a PHYLIP square distance matrix with short safe names; returns short->name."""
    short = {f"S{i:05d}": n for i, n in enumerate(names)}
    with open(path, "w") as fh:
        fh.write(f"{len(names)}\n")
        for i, s in enumerate(short):
            fh.write(s + " " + " ".join(f"{x:.6f}" for x in D[i]) + "\n")
    return short


def fastme_tree(D: np.ndarray, names: list[str], workdir: str, method: str = "B",
                spr: bool = True, threads: int = 1) -> tuple[Node, str]:
    """BIONJ (method 'B'/'I'... as in fastme -m) [+ SPR/NNI balanced ME] via fastme."""
    os.makedirs(workdir, exist_ok=True)
    dpath = os.path.join(workdir, "dist.phy")
    short = write_phylip_dist(D, names, dpath)
    out = os.path.join(workdir, "fastme.nwk")
    cmd = ["fastme", "-i", dpath, "-o", out, "-m", method, "-T", str(threads),
           "-I", os.path.join(workdir, "fastme.info")]
    if spr:
        cmd += ["-s", "-n"]   # SPR + NNI on the balanced ME criterion
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t = read_newick(out)
    for lf in t.leaves():
        lf.name = short[lf.name]
    return t, " ".join(cmd)


def _undirected(root: Node):
    adj: dict[int, list[tuple[int, float | None]]] = defaultdict(list)
    obj: dict[int, Node] = {}
    for nd in root.preorder():
        obj[id(nd)] = nd
        if nd.parent is not None:
            adj[id(nd)].append((id(nd.parent), nd.length))
            adj[id(nd.parent)].append((id(nd), nd.length))
    return adj, obj


def reroot(root: Node, target: Node, frac: float = 0.5) -> Node:
    """New rooted tree with the root placed on the edge (target, target.parent),
    at fraction ``frac`` of its length from ``target``. Degree-2 old root is
    removed. Node names/data are kept (on copies)."""
    if target.parent is None:
        return copy_tree(root)
    if len(root.children) == 1:          # unary root would become a fake leaf
        if target.parent is root:
            return copy_tree(root.children[0])
        root2 = copy_tree(root)
        # locate the copy of target by preorder position
        pos = [id(n) for n in root.preorder()].index(id(target))
        tgt2 = list(root2.preorder())[pos]
        root2 = collapse_unary(root2)
        return reroot(root2, tgt2, frac)
    adj, obj = _undirected(root)
    tid, pid = id(target), id(target.parent)
    L = target.length
    newroot = Node("")
    made: dict[int, Node] = {}

    def build(start: int, came_from: int, length):
        # iterative DFS copy away from came_from
        nd0 = obj[start]
        c0 = Node(nd0.name, length)
        c0.data = dict(nd0.data)
        stack = [(start, came_from, c0)]
        while stack:
            u, frm, cu = stack.pop()
            for v, l in adj[u]:
                if v == frm:
                    continue
                nv = obj[v]
                cv = Node(nv.name, l)
                cv.data = dict(nv.data)
                cu.add_child(cv)
                stack.append((v, u, cv))
        return c0

    la = None if L is None else L * frac
    lb = None if L is None else L * (1 - frac)
    newroot.add_child(build(tid, pid, la))
    newroot.add_child(build(pid, tid, lb))
    return collapse_unary_internal(newroot)


def collapse_unary_internal(root: Node) -> Node:
    """Collapse degree-2 internal (non-root) nodes created by rerooting."""
    for nd in list(root.postorder()):
        if nd.parent is not None and len(nd.children) == 1:
            ch, p = nd.children[0], nd.parent
            i = p.children.index(nd)
            p.children[i] = ch
            ch.parent = p
            if ch.length is not None and nd.length is not None:
                ch.length += nd.length
    return root


def midpoint_root(root: Node) -> Node:
    adj, obj = _undirected(root)
    leaves = [id(l) for l in root.leaves()]

    def far(src):
        dist = {src: 0.0}
        prev = {src: None}
        stack = [src]
        while stack:
            u = stack.pop()
            for v, l in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + (l or 0.0)
                    prev[v] = u
                    stack.append(v)
        best = max(leaves, key=lambda x: dist[x])
        return best, dist, prev

    a, _, _ = far(leaves[0])
    b, dist, prev = far(a)
    half = dist[b] / 2
    # walk from b towards a
    u = b
    while prev[u] is not None and dist[prev[u]] >= half:
        u = prev[u]
    # edge (u, prev[u]) contains the midpoint; dist decreasing towards a
    p = prev[u]
    if p is None:
        return copy_tree(root)
    # orient: find node object for edge in original rooted tree
    nu, np_ = obj[u], obj[p]
    if nu.parent is np_:
        tgt, frac_from_tgt = nu, (dist[u] - half) / max(1e-12, dist[u] - dist[p])
    else:
        tgt, frac_from_tgt = np_, (half - dist[p]) / max(1e-12, dist[u] - dist[p])
    return reroot(root, tgt, frac=min(max(frac_from_tgt, 0.0), 1.0))


def root_by_reference(root: Node, ref: Node) -> tuple[Node, dict]:
    """Root on the edge that maximises the number of reference (NCBI) clades that
    appear as clades of the rerooted tree (ties: more balanced root, then
    first in preorder). Uses only the reference topology (allowed input)."""
    refc = {s for s in clades(ref) if 1 < len(s) < len(ref.leaf_names())}
    best, best_score, best_info = None, None, {}
    nall = len(root.leaf_names())
    for nd in root.preorder():
        if nd.parent is None:
            continue
        t = reroot(root, nd)
        cs = set(clades(t))
        hit = len(cs & refc)
        sizes = sorted(len(_leafsets(t)[id(c)]) for c in t.children)
        bal = min(sizes) / nall
        score = (hit, bal)
        if best_score is None or score > best_score:
            best, best_score = t, score
            best_info = {"ref_clades_matched": hit, "ref_clades_total": len(refc),
                         "root_split_sizes": sizes}
    return best, best_info


def refine_polytomies(ref: Node, data: Node, refine_root: bool = False) -> tuple[Node, int]:
    """Hybrid arm: resolve polytomies of ``ref`` (NCBI) with the rooted clades of
    ``data`` that are compatible with every ``ref`` clade (nested or disjoint).
    Clades are inserted largest first; each new node groups a subset (>= 2, not
    all) of the children of the smallest node containing it. New nodes are named
    ``<parent>~r<k>``. A root polytomy is left unresolved unless ``refine_root``
    (the data trees' root position is itself chosen by NCBI agreement, so any
    resolution there is a rooting artefact). Returns (tree, number of clades added)."""
    t = copy_tree(ref)
    leaves = set(t.leaf_names())
    data = prune_to(data, leaves)
    cand = sorted((c for c in clades(data) if 1 < len(c) < len(leaves)), key=len, reverse=True)
    added = 0
    for c in cand:
        below = _leafsets(t)
        # smallest node containing c
        best = None
        for nd in t.preorder():
            if nd.children and c <= below[id(nd)]:
                if best is None or len(below[id(nd)]) < len(below[id(best)]):
                    best = nd
        if best is None or below[id(best)] == c or (best is t and not refine_root):
            continue
        inside = [ch for ch in best.children if below[id(ch)] <= c]
        covered = frozenset().union(*(below[id(ch)] for ch in inside)) if inside else frozenset()
        if covered != c:
            continue          # c cuts through a child clade -> incompatible with ref
        if len(inside) < 2 or len(inside) == len(best.children):
            continue
        added += 1
        u = Node(f"{best.name}~r{added}")
        u.data["ref_match"] = False
        for ch in inside:
            best.children.remove(ch)
            u.add_child(ch)
        best.add_child(u)
    return t, added


# ---------------------------------------------------------------------------
# perturbation / negative control
# ---------------------------------------------------------------------------

def nni_perturb(root: Node, k: int, seed: int) -> tuple[Node, int]:
    """Apply k random NNIs (seeded). An NNI on internal edge (p, u): swap a
    random child of u with a random other child of p. Internal names are
    cleared (rename with :func:`name_internal_nodes`). Returns (tree, n_done)."""
    rng = random.Random(seed)
    t = copy_tree(root)
    done = 0
    tries = 0
    while done < k and tries < 100 * max(1, k):
        tries += 1
        edges = [nd for nd in t.preorder() if nd.parent is not None and nd.children
                 and len(nd.parent.children) >= 2]
        if not edges:
            break
        u = rng.choice(edges)
        p = u.parent
        sibs = [s for s in p.children if s is not u]
        c = rng.choice(u.children)
        s = rng.choice(sibs)
        # swap c (child of u) with s (child of p)
        iu = u.children.index(c)
        ip = p.children.index(s)
        u.children[iu] = s
        s.parent = u
        p.children[ip] = c
        c.parent = p
        done += 1
    for nd in t.preorder():
        if nd.children:
            nd.name = ""
            nd.data.pop("ref_match", None)
    return t, done


def random_topology(leaves: list[str], seed: int) -> Node:
    """Random rooted binary tree: random sequential joining of a shuffled list
    (coalescent-like; seeded). Branch lengths 1."""
    rng = random.Random(seed)
    pool = [Node(l, 1.0) for l in leaves]
    while len(pool) > 1:
        i, j = sorted(rng.sample(range(len(pool)), 2))
        b, a = pool.pop(j), pool.pop(i)
        u = Node("", 1.0)
        u.add_child(a)
        u.add_child(b)
        pool.append(u)
    pool[0].length = None
    return pool[0]


# ---------------------------------------------------------------------------
# arm A: sampled best-hit distances
# ---------------------------------------------------------------------------

def list_proteomes(paths: list[str]) -> list[str]:
    out = []
    for p in paths:
        if os.path.isdir(p):
            for root_, _, files in os.walk(p):
                for f in files:
                    if f.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")) and not any(
                            x in f for x in ("_DNA", "_additional")):
                        out.append(os.path.join(root_, f))
        else:
            out.append(p)
    out.sort(key=species_name)
    return out


def sample_queries(files: list[str], n_per_species, seed: int, out_fa: str,
                   min_len: int = 50) -> dict:
    """Seeded random sample of proteins (length >= min_len) per species.
    ``n_per_species``: int or dict species -> int. IDs become ``s<idx>__<pid>``."""
    rng = np.random.default_rng(seed)
    info = {}
    with open(out_fa, "w") as fh:
        for si, f in enumerate(files):
            sp = species_name(f)
            recs = [(pid, s) for pid, s in read_fasta(f) if len(s) >= min_len]
            n = n_per_species[sp] if isinstance(n_per_species, dict) else n_per_species
            n = min(n, len(recs))
            pick = np.sort(rng.choice(len(recs), n, replace=False))
            for k in pick:
                pid, s = recs[k]
                fh.write(f">s{si}__{pid}\n{s}\n")
            info[sp] = {"n_proteins": len(recs), "n_sampled": int(n)}
    return info


def build_species_dbs(files: list[str], dbdir: str, threads: int = 4) -> list[str]:
    os.makedirs(dbdir, exist_ok=True)
    out = []
    for si, f in enumerate(files):
        db = os.path.join(dbdir, f"s{si}.dmnd")
        if not os.path.exists(db):
            fa = os.path.join(dbdir, f"s{si}.fa")
            with open(fa, "w") as fh:
                for pid, s in read_fasta(f):
                    fh.write(f">s{si}__{pid}\n{s}\n")
            subprocess.run(["diamond", "makedb", "--in", fa, "-d", db[:-5], "-p", str(threads),
                            "--quiet"], check=True)
            os.remove(fa)
        out.append(db)
    return out


def search_per_species(queries: str, dbs: list[str], outdir: str, threads: int,
                       sens: str = "--mid-sensitive", k: int = 5, evalue: float = 1e-5) -> dict:
    """DIAMOND blastp of the sample against each species DB separately (top-k
    per species, so best and second-best hit per species are known)."""
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    cpu0 = os.times()
    for db in dbs:
        out = os.path.join(outdir, os.path.basename(db)[:-5] + ".tsv")
        if os.path.exists(out):
            continue
        cmd = ["diamond", "blastp", "-q", queries, "-d", db, "-o", out + ".tmp", sens,
               "-k", str(k), "--max-hsps", "1", "-e", str(evalue), "-p", str(threads), "--quiet",
               "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen",
               "evalue", "bitscore"]
        subprocess.run(cmd, check=True)
        os.replace(out + ".tmp", out)
    c1 = os.times()
    return {"wall_s": time.time() - t0,
            "cpu_children_s": (c1.children_user - cpu0.children_user) + (c1.children_system - cpu0.children_system),
            "cmd_template": f"diamond blastp {sens} -k {k} --max-hsps 1 -e {evalue} -p {threads}"}


def load_hits(outdir: str, n_species: int, min_qcov: float = 0.0):
    """Concatenate per-species hit tables -> list of numpy record arrays.
    Columns: q_sp, query, t_sp, target, pident, length, qlen, slen, bits."""
    import pandas as pd
    frames = []
    for si in range(n_species):
        p = os.path.join(outdir, f"s{si}.tsv")
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            continue
        df = pd.read_csv(p, sep="\t", header=None,
                         names=["query", "target", "pident", "length", "qlen", "slen", "evalue", "bits"])
        df["t_sp"] = si
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["q_sp"] = df["query"].str.split("__", n=1).str[0].str[1:].astype(int)
    if min_qcov > 0:
        df = df[df["length"] / df["qlen"] >= min_qcov]
    return df


def besthit_distances(df, n_species: int, self_exclude: bool = True, min_hits: int = 5,
                      stat: str = "mean_identity"):
    """Species distance matrix from best hits.

    For queries of species i hitting species j: best hit (max bits; the query
    itself is excluded when i == j). d_ij = -ln(mean fractional identity over
    the best hits) (Poisson correction), symmetrised as the mean of d_ij and
    d_ji. Pairs with < ``min_hits`` hits get the maximum observed distance x 1.5.
    Also returns the hit fraction matrix (share of queries with a hit)."""
    d = df[df["q_sp"] != df["t_sp"]] if self_exclude else df
    best = d.sort_values("bits", ascending=False).drop_duplicates(["query", "t_sp"])
    g = best.groupby(["q_sp", "t_sp"])
    ident = g["pident"].mean() / 100.0
    cnt = g.size()
    M = np.full((n_species, n_species), np.nan)
    C = np.zeros((n_species, n_species))
    for (i, j), v in ident.items():
        if cnt[(i, j)] >= min_hits:
            M[i, j] = -math.log(max(v, 1e-3))
        C[i, j] = cnt[(i, j)]
    S = np.where(np.isnan(M), M.T, np.where(np.isnan(M.T), M, (M + M.T) / 2))
    mx = np.nanmax(S) if np.isfinite(np.nanmax(S)) else 1.0
    S = np.where(np.isnan(S), 1.5 * mx, S)
    np.fill_diagonal(S, 0.0)
    return S, C


# ---------------------------------------------------------------------------
# arm B: near-single-copy families, gene trees, ASTRAL-Pro
# ---------------------------------------------------------------------------

def families_from_s0(fam_tsv: str, gene_species: dict[str, str], n_species: int,
                     min_occ: float = 0.5, max_copy: int = 2, max_mean_copy: float = 1.3):
    """Near-single-copy families from an S0 families.tsv (gene_id, family_id)."""
    fams: dict[str, list[str]] = defaultdict(list)
    with open(fam_tsv) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            g, f = line.rstrip("\n").split("\t")[:2]
            fams[f].append(g)
    out = {}
    for f, genes in fams.items():
        c = Counter(gene_species[g] for g in genes)
        occ = len(c) / n_species
        if occ < min_occ:
            continue
        if max(c.values()) > max_copy or len(genes) / len(c) > max_mean_copy:
            continue
        out[f] = genes
    return out


def families_from_hits(df, species: list[str], min_occ: float = 0.5, ratio: float = 0.9,
                       min_qcov: float = 0.5, max_overlap: float = 0.2, max_fams: int = 500,
                       seed: int = 11):
    """Seed families from the sampled best hits (QfO, where S0 is not run).

    For a seed query q and target species t: include t's best hit h1 if
    query coverage >= min_qcov and (no second hit, or bits(h2) < ratio*bits(h1))
    (clear best hit = near single copy in t). The seed is included for its own
    species (self-hit excluded). Keep seeds with occupancy >= min_occ, then
    greedily (by occupancy, ties random-seeded) accept families sharing < max_overlap
    of members with any accepted family. Returns dict fam -> list of (species, protein)."""
    n = len(species)
    d = df[(df["length"] / df["qlen"] >= min_qcov)]
    d = d[d["query"] != d["target"]]
    d = d.sort_values(["query", "t_sp", "bits"], ascending=[True, True, False])
    top2 = d.groupby(["query", "t_sp"]).head(2)
    rank = top2.groupby(["query", "t_sp"]).cumcount()
    top2 = top2.assign(rank=rank.values)
    h1 = top2[top2["rank"] == 0].set_index(["query", "t_sp"])
    h2 = top2[top2["rank"] == 1].set_index(["query", "t_sp"])["bits"]
    h1 = h1.join(h2.rename("bits2"), how="left")
    ok = h1["bits2"].isna() | (h1["bits2"] < ratio * h1["bits"])
    h1 = h1[ok].reset_index()
    # own species: seed itself replaces the best (non-self) own-species hit
    h1 = h1[h1["t_sp"] != h1["q_sp"]]
    cand = {}
    for q, grp in h1.groupby("query"):
        qsp = int(q.split("__", 1)[0][1:])
        mem = {qsp: q}
        for t, tg in zip(grp["t_sp"].values, grp["target"].values):
            mem[int(t)] = tg
        if len(mem) / n >= min_occ:
            cand[q] = mem
    rng = random.Random(seed)
    keys = list(cand)
    rng.shuffle(keys)
    keys.sort(key=lambda q: -len(cand[q]))
    owner: dict[str, int] = {}
    fams = {}
    for q in keys:
        mem = cand[q]
        ids = set(mem.values())
        ov = Counter(owner[x] for x in ids if x in owner)
        if ov and max(ov.values()) / len(ids) >= max_overlap:
            continue
        if any(x in owner for x in ids):
            # drop members already used elsewhere, keep the rest
            mem = {s: g for s, g in mem.items() if g not in owner}
            if len(mem) / n < min_occ:
                continue
        fid = len(fams)
        fams[f"H{fid:05d}"] = [(species[s], g) for s, g in sorted(mem.items())]
        for g in mem.values():
            owner[g] = fid
        if len(fams) >= max_fams:
            break
    return fams, len(cand)


def _famsa_fasttree(args):
    fam, fa, outdir, ft_threads = args
    aln = os.path.join(outdir, fam + ".aln")
    tre = os.path.join(outdir, fam + ".nwk")
    if os.path.exists(tre) and os.path.getsize(tre) > 0:
        return fam, tre
    subprocess.run(["famsa", "-t", "1", fa, aln], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(tre + ".tmp", "w") as fh:
        subprocess.run(["FastTree", "-lg", "-gamma", "-quiet", "-nopr", aln], check=True,
                       stdout=fh, stderr=subprocess.DEVNULL)
    os.replace(tre + ".tmp", tre)
    return fam, tre


def gene_trees(fams: dict, seqs: dict, outdir: str, threads: int):
    """Write per-family FASTA (gene ids prefixed for uniqueness), FAMSA + FastTree
    (LG+Gamma). ``fams``: fam -> list of (species, protein). Returns (gene->species
    map, list of tree paths, cpu seconds of children)."""
    os.makedirs(outdir, exist_ok=True)
    gmap = {}
    jobs = []
    for fam, mem in fams.items():
        fa = os.path.join(outdir, fam + ".fa")
        with open(fa, "w") as fh:
            for k, (sp, pid) in enumerate(mem):
                gid = f"{fam}_g{k}"
                gmap[gid] = sp
                fh.write(f">{gid}\n{seqs[pid]}\n")
        jobs.append((fam, fa, outdir, 1))
    c0 = os.times()
    with ThreadPoolExecutor(threads) as ex:
        res = list(ex.map(_famsa_fasttree, jobs))
    c1 = os.times()
    cpu = (c1.children_user - c0.children_user) + (c1.children_system - c0.children_system)
    return gmap, [r[1] for r in res], cpu


def run_astral_pro(tree_paths: list[str], gmap: dict, outdir: str, threads: int,
                   binary: str = "astral-pro3") -> tuple[Node, dict]:
    gt = os.path.join(outdir, "genetrees.nwk")
    with open(gt, "w") as fh:
        for p in tree_paths:
            s = open(p).read().strip()
            if s:
                fh.write(s + "\n")
    mp = os.path.join(outdir, "gene2species.tsv")
    with open(mp, "w") as fh:
        for g, sp in gmap.items():
            fh.write(f"{g}\t{sp}\n")
    out = os.path.join(outdir, "astral_pro.nwk")
    t0 = time.time()
    c0 = os.times()
    cmd = [binary, "-t", str(threads), "-a", mp, "-o", out, "-i", gt, "--seed", "11"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=open(os.path.join(outdir, "astral.log"), "w"))
    c1 = os.times()
    t = read_newick(out)
    for nd in t.preorder():
        if nd.children:
            nd.data["support"] = nd.name
            nd.name = ""
    return t, {"wall_s": time.time() - t0,
               "cpu_s": (c1.children_user - c0.children_user) + (c1.children_system - c0.children_system),
               "cmd": " ".join(cmd)}


def concat_tree(fams: dict, outdir: str, gmap: dict, threads: int, species: list[str]):
    """Concatenate single-copy alignments (one sequence per species; families with
    copies keep the first) and run FastTree -lg -gamma on the supermatrix."""
    blocks = []
    for fam in fams:
        aln = os.path.join(outdir, fam + ".aln")
        if not os.path.exists(aln):
            continue
        rec = dict(read_fasta(aln))
        per = {}
        for gid, s in rec.items():
            per.setdefault(gmap[gid], s)
        L = len(next(iter(rec.values())))
        # drop columns with > 50% gaps
        arr = np.array([list(s) for s in per.values()])
        keep = (arr != "-").mean(axis=0) >= 0.5
        per = {sp: "".join(np.array(list(s))[keep]) for sp, s in per.items()}
        L = int(keep.sum())
        blocks.append((per, L))
    short = {sp: f"S{i:05d}" for i, sp in enumerate(species)}
    sm = os.path.join(outdir, "supermatrix.fa")
    tot = 0
    with open(sm, "w") as fh:
        for sp in species:
            parts = [per.get(sp, "-" * L) for per, L in blocks]
            seq = "".join(parts)
            tot = len(seq)
            fh.write(f">{short[sp]}\n{seq}\n")
    out = os.path.join(outdir, "concat_fasttree.nwk")
    t0 = time.time()
    c0 = os.times()
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(out, "w") as fh:
        subprocess.run(["FastTreeMP", "-lg", "-gamma", "-quiet", "-nopr", sm], check=True,
                       stdout=fh, stderr=subprocess.DEVNULL, env=env)
    c1 = os.times()
    t = read_newick(out)
    inv = {v: k for k, v in short.items()}
    for lf in t.leaves():
        lf.name = inv[lf.name]
    for nd in t.preorder():
        if nd.children:
            nd.name = ""
    return t, {"wall_s": time.time() - t0,
               "cpu_s": (c1.children_user - c0.children_user) + (c1.children_system - c0.children_system),
               "n_sites": tot, "n_blocks": len(blocks)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def finalize(tree: Node, ncbi: Node | None, rooting: str = "ncbi") -> tuple[Node, dict]:
    info = {"rooting": rooting}
    if rooting == "ncbi" and ncbi is not None:
        tree, ri = root_by_reference(tree, ncbi)
        info.update(ri)
    elif rooting == "midpoint":
        tree = midpoint_root(tree)
    name_internal_nodes(tree, ncbi)
    info["n_internal_named_as_ncbi"] = sum(1 for n in tree.preorder() if n.children and n.data.get("ref_match"))
    info["n_internal"] = sum(1 for n in tree.preorder() if n.children)
    return tree, info


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m ancestra.speciestree", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("ncbi", help="NCBI taxdump tree for proteome files")
    p.add_argument("--proteomes", nargs="+", required=True)
    p.add_argument("--taxdump", required=True)
    p.add_argument("--leaf-naming", choices=["qfo", "ensembl"], default="qfo")
    p.add_argument("--out", required=True)

    p = sub.add_parser("perturb", help="k random NNIs (seeded) of a tree")
    p.add_argument("--tree", required=True)
    p.add_argument("--k", type=int, required=True)
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--name-ref", help="tree used to name internal nodes (e.g. NCBI)")
    p.add_argument("--out", required=True)

    p = sub.add_parser("random", help="random binary topology (negative control)")
    p.add_argument("--tree", required=True, help="tree providing the leaf set")
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--name-ref")
    p.add_argument("--out", required=True)

    p = sub.add_parser("compare", help="nRF / quartet comparison")
    p.add_argument("query")
    p.add_argument("ref")

    a = ap.parse_args(argv)
    if a.cmd == "ncbi":
        files = list_proteomes(a.proteomes)
        if a.leaf_naming == "qfo":
            lt = qfo_leaf_taxids(files)
        else:
            nm = ensembl_leaf_names(files)
            tx = lookup_taxids_by_name(a.taxdump, list(nm.values()))
            lt = {leaf: tx[n] for leaf, n in nm.items()}
        t = ncbi_tree(lt, a.taxdump)
        write_tree(t, a.out, lengths=False)
        print(json.dumps({"n_leaves": len(t.leaves()),
                          "n_internal": sum(1 for n in t.preorder() if n.children)}))
    elif a.cmd in ("perturb", "random"):
        t = read_newick(a.tree)
        ref = read_newick(a.name_ref) if a.name_ref else None
        if a.cmd == "perturb":
            t2, done = nni_perturb(t, a.k, a.seed)
        else:
            t2 = random_topology(sorted(t.leaf_names()), a.seed)
        name_internal_nodes(t2, ref)
        write_tree(t2, a.out, lengths=False)
    elif a.cmd == "compare":
        print(json.dumps(compare_trees(read_newick(a.query), read_newick(a.ref)), indent=1))


if __name__ == "__main__":
    main()
