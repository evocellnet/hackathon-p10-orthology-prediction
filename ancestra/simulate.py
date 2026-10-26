"""Toy genome-evolution simulator with exact orthology truth.

Usage::

    python -m ancestra.simulate --out DIR --seed 11 [--species 8 --genes 300 ...]

Model (see ``docs/simulator.md`` for the full description and the list of
simplifications):

* **Species tree**: a random Yule (pure-birth) ultrametric tree with ``S``
  leaves, rescaled to a root-to-tip depth of ``tree_height`` expected
  substitutions per site, or a user Newick. Internal nodes ``N0`` (root),
  ``N1``... in pre-order; leaves ``sp0``... (user names are kept).
* **Root genome**: ``G`` genes ("families") placed in order on
  ``n_chromosomes`` chromosomes, each with a random root protein of length
  U[len_min, len_max] drawn from the WAG stationary amino-acid frequencies.
* **Gene-content events** along branches, as a continuous-time process over
  *all branches simultaneously* (global clock, needed for HGT): per gene and
  per unit branch length: tandem/dispersed duplication, loss, inversion and
  translocation of gene blocks, HGT (off by default). Optional whole-genome
  duplication on one branch, followed by fractionation losses.
* **Sequence evolution**: WAG exchangeabilities and frequencies (Whelan &
  Goldman 2001), discrete-gamma among-site rate heterogeneity (``K``
  categories, mean-of-category rates), per-family rate multipliers
  (Gamma(``family_rate_shape``) with mean 1), optional temporary rate boost
  of new duplicates (neofunctionalisation), and small insertions/deletions.
  Sequences are evolved lazily: the expected number of substitutions per
  site accumulated along a gene lineage is applied exactly (``P(t)=exp(Qt)``)
  whenever the lineage hits a node (speciation, duplication, transfer, leaf).

Truth (exact, from the simulated history):

* gene trees with speciation (S), duplication (D) and transfer (T) nodes;
* true HOGs (OrthoXML, nested per species-tree node);
* true ortholog pairs: pairs whose gene-tree LCA is a speciation node and
  whose paths to the LCA do not cross a horizontal transfer (xenologs are
  not orthologs). A transferred subtree becomes its own HOG family.
"""
from __future__ import annotations

import argparse
import heapq
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass, field

import numpy as np

from .newick import Node, parse_newick
from .orthoxml import GeneEntry, Group, OrthoXML, normalize_groups, write_orthoxml

__all__ = ["SimParams", "SimResult", "simulate", "write_simulation", "AAModel", "main"]

AA = "ARNDCQEGHILKMFPSTWYV"
_AA_BYTES = np.frombuffer(AA.encode(), dtype=np.uint8)

# WAG (Whelan & Goldman 2001, MBE 18:691). Lower-triangular exchangeabilities
# s_ij (rows 2..20, order ARNDCQEGHILKMFPSTWYV) followed by frequencies pi.
# Copied from the PAML-format wag.dat distributed with HMMER/easel
# (share/easel/formats/wag.dat in the project conda env).
_WAG = """
0.551571
0.509848 0.635346
0.738998 0.147304 5.429420
1.027040 0.528191 0.265256 0.0302949
0.908598 3.035500 1.543640 0.616783 0.0988179
1.582850 0.439157 0.947198 6.174160 0.021352 5.469470
1.416720 0.584665 1.125560 0.865584 0.306674 0.330052 0.567717
0.316954 2.137150 3.956290 0.930676 0.248972 4.294110 0.570025 0.249410
0.193335 0.186979 0.554236 0.039437 0.170135 0.113917 0.127395 0.0304501 0.138190
0.397915 0.497671 0.131528 0.0848047 0.384287 0.869489 0.154263 0.0613037 0.499462 3.170970
0.906265 5.351420 3.012010 0.479855 0.0740339 3.894900 2.584430 0.373558 0.890432 0.323832 0.257555
0.893496 0.683162 0.198221 0.103754 0.390482 1.545260 0.315124 0.174100 0.404141 4.257460 4.854020 0.934276
0.210494 0.102711 0.0961621 0.0467304 0.398020 0.0999208 0.0811339 0.049931 0.679371 1.059470 2.115170 0.088836 1.190630
1.438550 0.679489 0.195081 0.423984 0.109404 0.933372 0.682355 0.243570 0.696198 0.0999288 0.415844 0.556896 0.171329 0.161444
3.370790 1.224190 3.974230 1.071760 1.407660 1.028870 0.704939 1.341820 0.740169 0.319440 0.344739 0.967130 0.493905 0.545931 1.613280
2.121110 0.554413 2.030060 0.374866 0.512984 0.857928 0.822765 0.225833 0.473307 1.458160 0.326622 1.386980 1.516120 0.171903 0.795384 4.378020
0.113133 1.163920 0.0719167 0.129767 0.717070 0.215737 0.156557 0.336983 0.262569 0.212483 0.665309 0.137505 0.515706 1.529640 0.139405 0.523742 0.110864
0.240735 0.381533 1.086000 0.325711 0.543833 0.227710 0.196303 0.103604 3.873440 0.420170 0.398618 0.133264 0.428437 6.454280 0.216046 0.786993 0.291148 2.485390
2.006010 0.251849 0.196246 0.152335 1.002140 0.301281 0.588731 0.187247 0.118358 7.821300 1.800340 0.305434 2.058450 0.649892 0.314887 0.232739 1.388230 0.365369 0.314730
0.0866279 0.043972 0.0390894 0.0570451 0.0193078 0.0367281 0.0580589 0.0832518 0.0244313 0.048466 0.086209 0.0620286 0.0195027 0.0384319 0.0457631 0.0695179 0.0610127 0.0143859 0.0352742 0.0708956
"""


def wag_matrices() -> tuple[np.ndarray, np.ndarray]:
    """Return (exchangeabilities S[20,20], frequencies pi[20]) of WAG."""
    rows = [r.split() for r in _WAG.strip().splitlines()]
    S = np.zeros((20, 20))
    for i in range(1, 20):
        vals = rows[i - 1]
        assert len(vals) == i
        for j in range(i):
            S[i, j] = S[j, i] = float(vals[j])
    pi = np.array([float(x) for x in rows[19]])
    return S, pi / pi.sum()


def discrete_gamma_rates(shape: float, k: int) -> np.ndarray:
    """Mean-of-category discrete gamma rates (Yang 1994), mean 1."""
    if k <= 1 or not np.isfinite(shape):
        return np.ones(1)
    from scipy.special import gammainc
    from scipy.stats import gamma
    q = gamma.ppf(np.arange(1, k) / k, a=shape, scale=1.0 / shape)
    cdf = np.concatenate([[0.0], gammainc(shape + 1.0, q * shape), [1.0]])
    r = k * np.diff(cdf)
    return r / r.mean()


class AAModel:
    """Reversible amino-acid substitution model (WAG) with discrete gamma.

    ``Q_ij = s_ij * pi_j`` normalised to one expected substitution per unit
    time at stationarity. ``P(t)`` from the eigendecomposition of the
    symmetrised matrix.
    """

    def __init__(self, gamma_shape: float = 1.0, n_categories: int = 4):
        S, pi = wag_matrices()
        Q = S * pi[None, :]
        np.fill_diagonal(Q, 0.0)
        np.fill_diagonal(Q, -Q.sum(1))
        Q /= -(pi * np.diag(Q)).sum()
        self.Q, self.pi = Q, pi
        sq = np.sqrt(pi)
        A = sq[:, None] * Q / sq[None, :]
        A = (A + A.T) / 2
        lam, U = np.linalg.eigh(A)
        self.lam = lam
        self.left = U / sq[:, None]
        self.right = U.T * sq[None, :]
        self.rates = discrete_gamma_rates(gamma_shape, n_categories)
        self.ncat = len(self.rates)
        self.pi_cdf = np.cumsum(pi)

    def P(self, t: float) -> np.ndarray:
        """Transition matrices for each rate category, shape (ncat, 20, 20)."""
        E = np.exp(self.lam[None, :] * (self.rates[:, None] * t))
        P = (self.left[None, :, :] * E[:, None, :]) @ self.right[None, :, :]
        np.clip(P, 0.0, None, out=P)
        P /= P.sum(2, keepdims=True)
        return P

    def random_seq(self, n: int, rng: np.random.Generator) -> np.ndarray:
        return np.minimum(np.searchsorted(self.pi_cdf, rng.random(n)), 19).astype(np.uint8)

    def substitute(self, seq: np.ndarray, cats: np.ndarray, t: float,
                   rng: np.random.Generator) -> np.ndarray:
        if t <= 0 or len(seq) == 0:
            return seq
        C = np.cumsum(self.P(t), axis=2)
        rows = C[cats, seq]
        u = rng.random(len(seq))
        new = (u[:, None] > rows).sum(1)
        return np.minimum(new, 19).astype(np.uint8)


@dataclass
class SimParams:
    """All simulator parameters. Rates are per gene per unit branch length
    (branch lengths = expected substitutions/site for an average family)."""
    n_species: int = 8
    n_genes: int = 300
    n_chromosomes: int = 2
    newick: str | None = None          # Newick string of a fixed species tree
    tree_height: float | None = None   # root-to-tip depth (default 1.0 for random trees)
    dup_rate: float = 0.2
    loss_rate: float = 0.2
    inversion_rate: float = 0.05
    translocation_rate: float = 0.02
    hgt_rate: float = 0.0
    dispersed_dup_frac: float = 0.2
    block_mean: float = 4.0            # mean rearranged block length (genes, geometric)
    wgd_branch: str | None = None      # species-node name below the WGD branch, or "random";
                                       # a list "N3@0.3,N3@0.7" gives several WGDs (name@position)
    wgd_position: float = 0.5          # relative position of the WGD on that branch
    wgd_fractionation: float = 0.5     # P(one copy of a WGD pair is later lost)
    gamma_shape: float = 1.0           # among-site rate heterogeneity
    gamma_categories: int = 4
    family_rate_shape: float = 2.0     # per-family rate multiplier ~ Gamma(shape, 1/shape)
    indel_rate: float = 0.02           # indel events per site per substitution unit
    indel_mean_len: float = 3.0
    indel_max_len: int = 20
    neo_factor: float = 1.0            # rate multiplier of a new duplicate copy (1 = off)
    lineage_rates: str | None = None   # "sp9:3,N5:2": substitution-rate multiplier for every gene
                                       # entering that branch (inherited by the whole lineage below)
    neo_duration: float = 0.05         # duration of the boost (branch-length units)
    len_min: int = 150
    len_max: int = 600
    min_protein_len: int = 50


@dataclass
class SimResult:
    params: SimParams
    seed: int
    species_tree: Node
    proteomes: dict            # species -> list[(prot_id, sequence)]
    layout: dict               # species -> list[(chrom, list[(prot_id, strand, prot_len)])]
    gene_trees: dict           # family -> Node (pruned, unary-collapsed, with HGT marks)
    orthoxml: OrthoXML
    ortholog_pairs: list       # (gene_a, gene_b, level, family)
    events: list               # dicts
    gene_family: dict          # prot_id -> (family, hog_family)
    runtime_s: float = 0.0
    stats: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------
# species tree
# ----------------------------------------------------------------------------

def random_species_tree(n: int, rng: np.random.Generator, height: float = 1.0) -> Node:
    """Yule tree with ``n`` leaves, ultrametric, root-to-tip depth ``height``."""
    if n < 2:
        raise ValueError("need at least 2 species")
    root = Node()
    root.data["t"] = 0.0
    lineages = [root.add_child(Node()), root.add_child(Node())]
    t = 0.0
    for nd in lineages:
        nd.data["t0"] = 0.0
    while True:
        t += rng.exponential(1.0 / len(lineages))
        if len(lineages) == n:
            break
        k = int(rng.integers(len(lineages)))
        nd = lineages.pop(k)
        nd.data["t"] = t
        a, b = nd.add_child(Node()), nd.add_child(Node())
        a.data["t0"] = b.data["t0"] = t
        lineages.extend([a, b])
    for nd in lineages:
        nd.data["t"] = t
    for nd in root.preorder():
        if nd.parent is not None:
            nd.length = (nd.data["t"] - nd.parent.data["t"]) * height / t
    for nd in root.preorder():
        nd.data.clear()
    _name_tree(root)
    return root


def _name_tree(root: Node) -> None:
    ni = li = 0
    used = {n.name for n in root.preorder() if n.name}
    for nd in root.preorder():
        if nd.name:
            continue
        if nd.children:
            while f"N{ni}" in used:
                ni += 1
            nd.name = f"N{ni}"
            ni += 1
        else:
            while f"sp{li}" in used:
                li += 1
            nd.name = f"sp{li}"
            li += 1
        used.add(nd.name)


def _depths(root: Node) -> dict:
    d = {}
    for nd in root.preorder():
        d[nd.name] = 0.0 if nd.parent is None else d[nd.parent.name] + (nd.length or 0.0)
    return d


# ----------------------------------------------------------------------------
# simulation core
# ----------------------------------------------------------------------------

class _Gene:
    __slots__ = ("seq", "cats", "pending", "t_last", "node", "fam", "rate",
                 "boost_until", "strand", "alive")

    def clone(self, node: int) -> "_Gene":
        g = _Gene()
        g.seq, g.cats, g.pending, g.t_last = self.seq, self.cats, 0.0, self.t_last
        g.node, g.fam, g.rate = node, self.fam, self.rate
        g.boost_until, g.strand, g.alive = self.boost_until, self.strand, True
        return g


class _Branch:
    __slots__ = ("name", "node", "t_start", "t_end", "chroms", "chrom_names", "n", "version")


class _Simulator:
    KINDS = ("dup", "loss", "inversion", "translocation", "hgt")

    def __init__(self, params: SimParams, seed: int):
        self.p = params
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.model = AAModel(params.gamma_shape, params.gamma_categories)
        # gene-tree node table (parent id < child id always)
        self.gt_parent: list[int] = []
        self.gt_kind: list[str] = []
        self.gt_sp: list[str] = []
        self.gt_time: list[float] = []
        self.gt_len: list[float] = []
        self.gt_fam: list[int] = []
        self.gt_extra: list[str] = []
        self.leaf_name: dict[int, str] = {}
        self.events: list[dict] = []
        self.heap: list = []
        self.counter = 0
        self.active: dict[str, _Branch] = {}
        self.leaf_genomes: dict[str, tuple[list, list]] = {}
        rates = np.array([params.dup_rate, params.loss_rate, params.inversion_rate,
                          params.translocation_rate, params.hgt_rate], dtype=float)
        if (rates < 0).any():
            raise ValueError("rates must be >= 0")
        self.rate_tot = float(rates.sum())
        self.kind_cdf = np.cumsum(rates) / self.rate_tot if self.rate_tot > 0 else None

    # -- helpers ------------------------------------------------------------
    def new_node(self, kind, parent, sp, t, length, fam, extra=""):
        self.gt_parent.append(parent)
        self.gt_kind.append(kind)
        self.gt_sp.append(sp)
        self.gt_time.append(t)
        self.gt_len.append(length)
        self.gt_fam.append(fam)
        self.gt_extra.append(extra)
        return len(self.gt_parent) - 1

    def log(self, t, etype, branch, fam, node, details=""):
        self.events.append(dict(event_id=len(self.events), time=round(t, 6), type=etype,
                                branch=branch, family=fam, gene_tree_node=node, details=details))

    def push(self, t, kind, payload):
        self.counter += 1
        heapq.heappush(self.heap, (t, self.counter, kind, payload))

    def accumulate(self, g: _Gene, t: float) -> None:
        dt = t - g.t_last
        if dt > 0:
            extra = 0.0
            if self.p.neo_factor != 1.0 and g.boost_until > g.t_last:
                extra = (self.p.neo_factor - 1.0) * (min(t, g.boost_until) - g.t_last)
            g.pending += g.rate * (dt + extra)
            g.t_last = t

    def materialise(self, g: _Gene, t: float) -> float:
        """Apply the accumulated evolution to g's sequence; return edge length."""
        self.accumulate(g, t)
        d = g.pending
        if d > 0:
            g.seq = self.model.substitute(g.seq, g.cats, d, self.rng)
            if self.p.indel_rate > 0:
                self._indels(g, d)
        g.pending = 0.0
        return d

    def _indels(self, g: _Gene, d: float) -> None:
        p, rng = self.p, self.rng
        n_ev = rng.poisson(p.indel_rate * d * len(g.seq))
        for _ in range(int(n_ev)):
            L = len(g.seq)
            ln = min(int(rng.geometric(1.0 / p.indel_mean_len)), p.indel_max_len)
            pos = int(rng.integers(L + 1))
            if rng.random() < 0.5:
                ins = self.model.random_seq(ln, rng)
                cats = rng.integers(self.model.ncat, size=ln).astype(np.uint8)
                g.seq = np.concatenate([g.seq[:pos], ins, g.seq[pos:]])
                g.cats = np.concatenate([g.cats[:pos], cats, g.cats[pos:]])
            elif L - ln >= p.min_protein_len - 1:
                pos = min(pos, L - ln)
                g.seq = np.concatenate([g.seq[:pos], g.seq[pos + ln:]])
                g.cats = np.concatenate([g.cats[:pos], g.cats[pos + ln:]])

    def pick(self, br: _Branch) -> tuple[int, int]:
        r = int(self.rng.integers(br.n))
        for ci, ch in enumerate(br.chroms):
            if r < len(ch):
                return ci, r
            r -= len(ch)
        raise RuntimeError("gene count out of sync")

    def schedule_next(self, br: _Branch, t: float) -> None:
        br.version += 1
        if self.rate_tot <= 0 or br.n == 0:
            return
        tn = t + self.rng.exponential(1.0 / (self.rate_tot * br.n))
        if tn < br.t_end:
            self.push(tn, "ev", (br.name, br.version))

    # -- main ---------------------------------------------------------------
    def run(self, tree: Node) -> None:
        p, rng = self.p, self.rng
        self.tree = tree
        self.depth = _depths(tree)
        self.nodes = {n.name: n for n in tree.preorder()}
        if p.wgd_branch == "random":
            cands = [n.name for n in tree.preorder() if n.parent is not None and (n.length or 0) > 0]
            self.wgd_branch = cands[int(rng.integers(len(cands)))]
        else:
            self.wgd_branch = p.wgd_branch
        self.wgd_list = []          # [(branch, position)]
        if self.wgd_branch is not None:
            for item in str(self.wgd_branch).split(","):
                nm, _, pos = item.partition("@")
                self.wgd_list.append((nm, float(pos) if pos else p.wgd_position))
        for nm, _ in self.wgd_list:
            if nm not in self.nodes:
                raise ValueError(f"wgd_branch {nm!r} not in species tree")
        self.lineage_rates = {}
        if p.lineage_rates:
            for item in p.lineage_rates.split(","):
                nm, _, f = item.partition(":")
                if nm not in self.nodes:
                    raise ValueError(f"lineage_rates node {nm!r} not in species tree")
                self.lineage_rates[nm] = float(f)
        # root genome
        genes = []
        for f in range(p.n_genes):
            g = _Gene()
            L = int(rng.integers(p.len_min, p.len_max + 1)) - 1   # 'M' is prepended at output
            g.seq = self.model.random_seq(L, rng)
            g.cats = rng.integers(self.model.ncat, size=L).astype(np.uint8)
            g.pending, g.t_last, g.node, g.fam = 0.0, 0.0, -1, f
            g.rate = float(rng.gamma(p.family_rate_shape, 1.0 / p.family_rate_shape)) \
                if np.isfinite(p.family_rate_shape) else 1.0
            g.boost_until, g.strand, g.alive = -1.0, (1 if rng.random() < 0.5 else -1), True
            genes.append(g)
        nchr = max(1, min(p.n_chromosomes, p.n_genes))
        cuts = np.sort(rng.choice(np.arange(1, p.n_genes), size=nchr - 1, replace=False)) \
            if nchr > 1 else np.array([], dtype=int)
        bounds = [0, *cuts.tolist(), p.n_genes]
        chroms = [genes[bounds[i]:bounds[i + 1]] for i in range(nchr)]
        names = [f"chr{i + 1}" for i in range(nchr)]
        self.speciate(tree, chroms, names, 0.0)
        while self.heap:
            t, _, kind, payload = heapq.heappop(self.heap)
            if kind == "end":
                self.end_branch(self.active.pop(payload), t)
            elif kind == "ev":
                name, ver = payload
                br = self.active.get(name)
                if br is None or br.version != ver:
                    continue
                self.random_event(br, t)
                self.schedule_next(br, t)
            elif kind == "wgd":
                br = self.active[payload]
                self.wgd(br, t)
                self.schedule_next(br, t)
            elif kind == "sched_loss":
                name, g = payload
                br = self.active.get(name)
                if br is None or not g.alive:
                    continue
                for ci, ch in enumerate(br.chroms):
                    if g in ch:   # identity check via __eq__ default
                        self.loss(br, ci, ch.index(g), t, reason="wgd_fractionation")
                        break
                self.schedule_next(br, t)

    def speciate(self, node: Node, chroms, names, t: float) -> None:
        new_ids = []
        for ch in chroms:
            ids = []
            for g in ch:
                d = self.materialise(g, t)
                s = self.new_node("S", g.node, node.name, t, d, g.fam)
                g.alive = False
                ids.append(s)
            new_ids.append(ids)
        for c in node.children:
            br = _Branch()
            br.name, br.node, br.t_start = c.name, c, t
            br.t_end = t + max(0.0, c.length or 0.0)
            br.chroms = [[g.clone(s) for g, s in zip(ch, ids)] for ch, ids in zip(chroms, new_ids)]
            br.chrom_names = list(names)
            br.n = sum(len(ch) for ch in br.chroms)
            br.version = 0
            self.active[c.name] = br
            self.push(br.t_end, "end", c.name)
            lr = self.lineage_rates.get(c.name)
            if lr is not None:
                for ch in br.chroms:
                    for g in ch:
                        g.rate *= lr
            for nm, pos in self.wgd_list:
                if nm == c.name:
                    self.push(t + pos * (br.t_end - t), "wgd", c.name)
            self.schedule_next(br, t)

    def end_branch(self, br: _Branch, t: float) -> None:
        br.version += 1
        if br.node.children:
            self.speciate(br.node, br.chroms, br.chrom_names, t)
        else:
            out = []
            for ch, cname in zip(br.chroms, br.chrom_names):
                lst = []
                for g in ch:
                    d = self.materialise(g, t)
                    leaf = self.new_node("leaf", g.node, br.name, t, d, g.fam)
                    g.alive = False
                    lst.append((leaf, g))
                out.append((cname, lst))
            self.leaf_genomes[br.name] = out

    # -- events -------------------------------------------------------------
    def random_event(self, br: _Branch, t: float) -> None:
        k = self.KINDS[int(np.searchsorted(self.kind_cdf, self.rng.random(), side="right"))]
        ci, idx = self.pick(br)
        getattr(self, k)(br, ci, idx, t)

    def dup(self, br: _Branch, ci: int, idx: int, t: float, extra: str = "") -> _Gene:
        rng, p = self.rng, self.p
        g = br.chroms[ci][idx]
        d = self.materialise(g, t)
        D = self.new_node("D", g.node, br.name, t, d, g.fam, extra)
        g.node = D
        g2 = g.clone(D)
        if self.p.neo_factor != 1.0:
            g2.boost_until = t + p.neo_duration
        if extra == "wgd":
            return g2
        if rng.random() < p.dispersed_dup_frac:
            cj = self.pick_chrom(br)
            pos = int(rng.integers(len(br.chroms[cj]) + 1))
            g2.strand = 1 if rng.random() < 0.5 else -1
            br.chroms[cj].insert(pos, g2)
            mode = f"dispersed:{br.chrom_names[cj]}"
        else:
            br.chroms[ci].insert(idx + 1, g2)
            mode = "tandem"
        br.n += 1
        self.log(t, "duplication", br.name, g.fam, D, mode)
        return g2

    def pick_chrom(self, br: _Branch) -> int:
        sizes = np.array([len(c) for c in br.chroms], dtype=float) + 1.0
        return int(self.rng.choice(len(sizes), p=sizes / sizes.sum()))

    def loss(self, br: _Branch, ci: int, idx: int, t: float, reason: str = "") -> None:
        g = br.chroms[ci].pop(idx)
        self.accumulate(g, t)
        L = self.new_node("L", g.node, br.name, t, g.pending, g.fam)
        g.alive = False
        br.n -= 1
        self.log(t, "loss", br.name, g.fam, L, reason)

    def _block(self, br: _Branch, ci: int, idx: int) -> int:
        k = int(self.rng.geometric(1.0 / max(1.0, self.p.block_mean)))
        return max(1, min(k, len(br.chroms[ci]) - idx))

    def inversion(self, br: _Branch, ci: int, idx: int, t: float) -> None:
        k = self._block(br, ci, idx)
        ch = br.chroms[ci]
        blk = ch[idx:idx + k][::-1]
        for g in blk:
            g.strand = -g.strand
        ch[idx:idx + k] = blk
        self.log(t, "inversion", br.name, -1, -1, f"{br.chrom_names[ci]}:{idx}:{k}")

    def translocation(self, br: _Branch, ci: int, idx: int, t: float) -> None:
        k = self._block(br, ci, idx)
        ch = br.chroms[ci]
        blk = ch[idx:idx + k]
        del ch[idx:idx + k]
        inv = self.rng.random() < 0.5
        if inv:
            blk = blk[::-1]
            for g in blk:
                g.strand = -g.strand
        cj = int(self.rng.integers(len(br.chroms)))
        pos = int(self.rng.integers(len(br.chroms[cj]) + 1))
        br.chroms[cj][pos:pos] = blk
        self.log(t, "translocation", br.name, -1, -1,
                 f"{br.chrom_names[ci]}:{idx}:{k}->{br.chrom_names[cj]}:{pos}:{'inv' if inv else 'fwd'}")

    def hgt(self, br: _Branch, ci: int, idx: int, t: float) -> None:
        others = [b for b in self.active.values() if b is not br]
        if not others:
            return
        rec = others[int(self.rng.integers(len(others)))]
        g = br.chroms[ci][idx]
        d = self.materialise(g, t)
        T = self.new_node("T", g.node, br.name, t, d, g.fam, rec.name)
        g.node = T
        R = self.new_node("R", T, rec.name, t, 0.0, g.fam, br.name)
        g2 = g.clone(R)
        g2.strand = 1 if self.rng.random() < 0.5 else -1
        cj = self.pick_chrom(rec)
        rec.chroms[cj].insert(int(self.rng.integers(len(rec.chroms[cj]) + 1)), g2)
        rec.n += 1
        self.log(t, "hgt", br.name, g.fam, T, f"to:{rec.name}")
        self.schedule_next(rec, t)

    def wgd(self, br: _Branch, t: float) -> None:
        self.log(t, "wgd", br.name, -1, -1, f"n_genes={br.n}")
        new_chroms = []
        pairs = []
        for ci, ch in enumerate(br.chroms):
            copy = []
            for idx in range(len(ch)):
                g2 = self.dup(br, ci, idx, t, extra="wgd")
                copy.append(g2)
                pairs.append((ch[idx], g2))
                self.log(t, "duplication", br.name, g2.fam, g2.node, "wgd")
            new_chroms.append(copy)
        br.chroms.extend(new_chroms)
        br.chrom_names.extend(f"{n}b" for n in list(br.chrom_names))
        br.n *= 2
        f = self.p.wgd_fractionation
        for a, b in pairs:
            if self.rng.random() < f:
                victim = a if self.rng.random() < 0.5 else b
                self.push(t + self.rng.uniform(0, br.t_end - t), "sched_loss", (br.name, victim))


# ----------------------------------------------------------------------------
# truth extraction
# ----------------------------------------------------------------------------

def _build_truth(sim: _Simulator, fam_names: list[str]):
    N = len(sim.gt_parent)
    parent, kind = sim.gt_parent, sim.gt_kind
    alive = [False] * N
    children: list[list[int]] = [[] for _ in range(N)]
    for n in range(N - 1, -1, -1):
        if kind[n] == "leaf":
            alive[n] = True
        if alive[n] and parent[n] >= 0:
            alive[parent[n]] = True
    for n in range(N):
        if alive[n] and parent[n] >= 0:
            children[parent[n]].append(n)
    eff: dict[int, Node] = {}
    roots: dict[int, Node] = {}
    for n in range(N - 1, -1, -1):
        if not alive[n]:
            continue
        k = kind[n]
        if k == "leaf":
            obj = Node(sim.leaf_name[n], sim.gt_len[n])
            obj.data.update(kind="leaf", species=sim.gt_sp[n])
        else:
            ch = [eff.pop(c) for c in children[n]]
            if len(ch) == 1:
                obj = ch[0]
                obj.length = (obj.length or 0.0) + sim.gt_len[n]
                if k == "R":
                    obj.data["xfer"] = f"{sim.gt_extra[n]}>{sim.gt_sp[n]}"
            else:
                label = f"{k}_{sim.gt_sp[n]}" if k != "T" else f"T_{sim.gt_sp[n]}_to_{sim.gt_extra[n]}"
                obj = Node(label, sim.gt_len[n])
                obj.data.update(kind=k, species=sim.gt_sp[n], wgd=(sim.gt_extra[n] == "wgd"))
                for c in ch:
                    obj.add_child(c)
        if parent[n] < 0:
            obj.length = None
            roots[sim.gt_fam[n]] = obj
        else:
            eff[n] = obj
    gene_trees = {fam_names[f]: roots[f] for f in sorted(roots)}

    # cut at transfers -> HOG families; assign hog family to each node
    hog_groups: list[Group] = []
    gene_family: dict[str, tuple[str, str]] = {}
    pairs: list[tuple[str, str, str, str]] = []
    for fam, root in gene_trees.items():
        xcount = 0
        # assign hog-family ids top-down
        hogfam: dict[int, str] = {}
        for nd in root.preorder():
            if nd.parent is None:
                hogfam[id(nd)] = fam
            elif "xfer" in nd.data:
                xcount += 1
                hogfam[id(nd)] = f"{fam}.x{xcount}"
            else:
                hogfam[id(nd)] = hogfam[id(nd.parent)]
        # leaves per node restricted to same hog family (post-order)
        same: dict[int, list[str]] = {}
        for nd in root.postorder():
            hf = hogfam[id(nd)]
            if nd.is_leaf():
                gene_family[nd.name] = (fam, hf)
                same[id(nd)] = [nd.name]
                continue
            child_lists = [same[id(c)] if hogfam[id(c)] == hf else [] for c in nd.children]
            if nd.data["kind"] == "S":
                for i in range(len(child_lists)):
                    for j in range(i + 1, len(child_lists)):
                        for a in child_lists[i]:
                            for b in child_lists[j]:
                                pa, pb = (a, b) if a < b else (b, a)
                                pairs.append((pa, pb, nd.data["species"], hf))
            same[id(nd)] = [x for cl in child_lists for x in cl]
            for c in nd.children:
                same.pop(id(c), None)
        # HOG groups: convert each hog-family subtree into nested Groups
        tops: dict[str, object] = {}
        conv: dict[int, object] = {}
        for nd in root.postorder():
            hf = hogfam[id(nd)]
            if nd.is_leaf():
                obj = nd.name
            else:
                kids = [conv.pop(id(c)) for c in nd.children if hogfam[id(c)] == hf]
                k = nd.data["kind"]
                if k == "S":
                    obj = Group("ortholog", kids, level=nd.data["species"])
                elif k == "D":
                    obj = Group("paralog", kids, level=nd.data["species"])
                else:  # T node: the transferred child was cut; keep donor side
                    obj = Group("paralog", kids, level=nd.data["species"])  # unary -> collapsed
            if nd.parent is None or hogfam[id(nd.parent)] != hf:
                tops[hf] = obj
            else:
                conv[id(nd)] = obj
        for hf in sorted(tops):
            # collapse unary nodes (losses, cut transfers); singleton families vanish
            for g in normalize_groups([tops[hf]]):
                g.id = hf
                hog_groups.append(g)
    return gene_trees, hog_groups, pairs, gene_family


def _nhx(root: Node) -> str:
    """Newick with NHX annotations (Ev=S/D/T, S=species node, XFER=donor>recipient)."""
    parts = []
    stack: list = [(root, 0)]
    while stack:
        nd, st = stack.pop()
        if st == 2:
            parts.append(",")
            continue
        if st == 0 and nd.children:
            parts.append("(")
            stack.append((nd, 1))
            for i, c in enumerate(reversed(nd.children)):
                stack.append((c, 0))
                if i < len(nd.children) - 1:
                    stack.append((None, 2))
            continue
        if st == 1:
            parts.append(")")
        parts.append(nd.name)
        if nd.length is not None:
            parts.append(f":{nd.length:.6g}")
        ann = []
        if nd.children:
            ann.append(f"Ev={nd.data['kind']}")
            ann.append(f"S={nd.data['species']}")
            if nd.data.get("wgd"):
                ann.append("WGD=1")
        else:
            ann.append(f"S={nd.data['species']}")
        if "xfer" in nd.data:
            ann.append(f"XFER={nd.data['xfer']}")
        parts.append("[&&NHX:" + ":".join(ann) + "]")
    return "".join(parts) + ";"


# ----------------------------------------------------------------------------
# public API
# ----------------------------------------------------------------------------

def simulate(params: SimParams | None = None, seed: int = 11) -> SimResult:
    """Run the simulator and return all outputs in memory."""
    p = params or SimParams()
    t0 = time.time()
    rng_tree = np.random.default_rng([seed, 1])
    if p.newick:
        tree = parse_newick(p.newick)
        _name_tree(tree)
        if any(n.length is None for n in tree.preorder() if n.parent is not None):
            raise ValueError("user species tree must have branch lengths")
        if p.tree_height is not None:
            h = max(_depths(tree).values())
            for n in tree.preorder():
                if n.length is not None:
                    n.length *= p.tree_height / h
        tree.length = None
    else:
        tree = random_species_tree(p.n_species, rng_tree, p.tree_height or 1.0)
    sim = _Simulator(p, seed)
    sim.run(tree)
    # name proteins by genomic order
    width = max(4, len(str(max([sum(len(l) for _, l in v) for v in sim.leaf_genomes.values()] + [1]))))
    proteomes, layout = {}, {}
    leaf_order = tree.leaf_names()
    for sp in leaf_order:
        prots, lay = [], []
        k = 0
        for cname, lst in sim.leaf_genomes[sp]:
            entries = []
            for leaf, g in lst:
                k += 1
                pid = f"{sp}_g{k:0{width}d}"
                sim.leaf_name[leaf] = pid
                seq = "M" + _AA_BYTES[g.seq].tobytes().decode()
                prots.append((pid, seq))
                entries.append((pid, g.strand, len(seq)))
            if entries:
                lay.append((cname, entries))
        proteomes[sp], layout[sp] = prots, lay
    fam_names = [f"fam{f:0{max(4, len(str(p.n_genes - 1)))}d}" for f in range(p.n_genes)]
    gene_trees, hog_groups, pairs, gene_family = _build_truth(sim, fam_names)
    genes = {}
    for sp in leaf_order:
        for pid, _ in proteomes[sp]:
            genes[pid] = GeneEntry(pid, sp, gene_id=pid)
    ox = OrthoXML(genes=genes, groups=hog_groups, species_tree=tree, origin="ancestra.simulate",
                  origin_version=_version(), database="ancestra_sim", database_version=str(seed))
    for e in sim.events:
        if e["family"] >= 0:
            e["family"] = fam_names[e["family"]]
        else:
            e["family"] = ""
    stats = dict(n_proteins=len(genes), n_families_surviving=len(gene_trees),
                 n_true_ortholog_pairs=len(pairs),
                 n_events={k: sum(1 for e in sim.events if e["type"] == k)
                           for k in ("duplication", "loss", "inversion", "translocation", "hgt", "wgd")},
                 wgd_branch=sim.wgd_branch)
    res = SimResult(p, seed, tree, proteomes, layout, gene_trees, ox, pairs, sim.events,
                    gene_family, runtime_s=time.time() - t0, stats=stats)
    return res


def _version() -> str:
    try:
        from . import __version__
        return __version__
    except Exception:  # pragma: no cover
        return "0"


def write_simulation(res: SimResult, out: str) -> None:
    """Write all simulator outputs into directory ``out``."""
    for d in ("proteomes", "gff", "truth"):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    tree = res.species_tree
    with open(os.path.join(out, "species_tree.nwk"), "w") as fh:
        fh.write(tree.to_newick() + "\n")
    for sp, prots in res.proteomes.items():
        with open(os.path.join(out, "proteomes", f"{sp}.fa"), "w") as fh:
            for pid, seq in prots:
                fh.write(f">{pid}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
    rng = np.random.default_rng([res.seed, 2])
    for sp, lay in res.layout.items():
        lines = ["##gff-version 3"]
        body = []
        for cname, entries in lay:
            pos = int(rng.integers(1000, 10000))
            for pid, strand, plen in entries:
                start = pos + 1
                end = start + 3 * (plen + 1) - 1          # CDS incl. stop codon, single exon
                s = "+" if strand > 0 else "-"
                body.append(f"{cname}\tancestra_sim\tgene\t{start}\t{end}\t.\t{s}\t.\t"
                            f"ID=gene:{pid};Name={pid};biotype=protein_coding")
                body.append(f"{cname}\tancestra_sim\tmRNA\t{start}\t{end}\t.\t{s}\t.\t"
                            f"ID=transcript:{pid}_t1;Parent=gene:{pid};biotype=protein_coding")
                body.append(f"{cname}\tancestra_sim\tCDS\t{start}\t{end}\t.\t{s}\t0\t"
                            f"ID=CDS:{pid};Parent=transcript:{pid}_t1;protein_id={pid}")
                pos = end + int(rng.integers(200, 5000))
            lines.append(f"##sequence-region {cname} 1 {pos + 1000}")
        with open(os.path.join(out, "gff", f"{sp}.gff3"), "w") as fh:
            fh.write("\n".join(lines + body) + "\n")
    td = os.path.join(out, "truth")
    for fam, root in res.gene_trees.items():
        with open(os.path.join(td, f"gene_tree_{fam}.nwk"), "w") as fh:
            fh.write(_nhx(root) + "\n")
    write_orthoxml(res.orthoxml, os.path.join(td, "hogs.orthoxml"))
    with open(os.path.join(td, "ortholog_pairs.tsv"), "w") as fh:
        fh.write("gene_a\tgene_b\tlevel\thog_family\n")
        for a, b, lvl, hf in sorted(res.ortholog_pairs):
            fh.write(f"{a}\t{b}\t{lvl}\t{hf}\n")
    with open(os.path.join(td, "events.tsv"), "w") as fh:
        cols = ["event_id", "time", "type", "branch", "family", "gene_tree_node", "details"]
        fh.write("\t".join(cols) + "\n")
        for e in res.events:
            fh.write("\t".join(str(e[c]) for c in cols) + "\n")
    with open(os.path.join(td, "gene_families.tsv"), "w") as fh:
        fh.write("protein\tspecies\tfamily\thog_family\n")
        for sp, prots in res.proteomes.items():
            for pid, _ in prots:
                fam, hf = res.gene_family[pid]
                fh.write(f"{pid}\t{sp}\t{fam}\t{hf}\n")
    meta = dict(seed=res.seed, params=asdict(res.params), version=_version(),
                runtime_s=round(res.runtime_s, 3), stats=res.stats)
    with open(os.path.join(out, "params.json"), "w") as fh:
        json.dump(meta, fh, indent=1)


def _parser() -> argparse.ArgumentParser:
    d = SimParams()
    ap = argparse.ArgumentParser(prog="python -m ancestra.simulate",
                                 description="Toy genome-evolution simulator with exact orthology truth.")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--species", "-S", type=int, default=d.n_species, help="number of species (random tree)")
    ap.add_argument("--genes", "-G", type=int, default=d.n_genes, help="number of root genes / families")
    ap.add_argument("--chromosomes", type=int, default=d.n_chromosomes)
    ap.add_argument("--tree", help="Newick file of a fixed rooted species tree (with branch lengths)")
    ap.add_argument("--tree-height", type=float, default=None,
                    help="root-to-tip depth in subst/site (default 1.0; rescales a --tree if given)")
    for name in ("dup_rate", "loss_rate", "inversion_rate", "translocation_rate", "hgt_rate",
                 "dispersed_dup_frac", "block_mean", "wgd_position", "wgd_fractionation",
                 "gamma_shape", "family_rate_shape", "indel_rate", "indel_mean_len",
                 "neo_factor", "neo_duration"):
        ap.add_argument("--" + name.replace("_", "-"), type=float, default=getattr(d, name))
    for name in ("gamma_categories", "indel_max_len", "len_min", "len_max", "min_protein_len"):
        ap.add_argument("--" + name.replace("_", "-"), type=int, default=getattr(d, name))
    ap.add_argument("--lineage-rates", default=None,
                    help="substitution-rate multipliers per lineage, e.g. 'sp9:3,N5:2'")
    ap.add_argument("--wgd-branch", default=None,
                    help="species-tree node below the branch carrying a WGD, or 'random'")
    return ap


def main(argv=None) -> int:
    a = _parser().parse_args(argv)
    newick = None
    if a.tree:
        with open(a.tree) as fh:
            newick = fh.read()
    p = SimParams(n_species=a.species, n_genes=a.genes, n_chromosomes=a.chromosomes, newick=newick,
                  tree_height=a.tree_height, wgd_branch=a.wgd_branch, lineage_rates=a.lineage_rates,
                  **{k: getattr(a, k) for k in (
                      "dup_rate", "loss_rate", "inversion_rate", "translocation_rate", "hgt_rate",
                      "dispersed_dup_frac", "block_mean", "wgd_position", "wgd_fractionation",
                      "gamma_shape", "family_rate_shape", "indel_rate", "indel_mean_len",
                      "neo_factor", "neo_duration", "gamma_categories", "indel_max_len",
                      "len_min", "len_max", "min_protein_len")})
    t0 = time.time()
    res = simulate(p, a.seed)
    write_simulation(res, a.out)
    s = res.stats
    print(f"[ancestra.simulate] seed={a.seed} species={len(res.proteomes)} proteins={s['n_proteins']} "
          f"families={s['n_families_surviving']} true_pairs={s['n_true_ortholog_pairs']} "
          f"events={s['n_events']} sim={res.runtime_s:.1f}s total={time.time() - t0:.1f}s -> {a.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
