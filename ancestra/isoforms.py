"""S7: isoform-aware pair scoring (canonical vs best-isoform-pair).

OrthoBench / QfO inputs carry one protein per gene. Ensembl ``pep.all`` files list
every translated isoform with ``gene:`` and ``transcript:`` tags. Given all-vs-all
hits between isoforms, a gene pair's score is

* ``canonical``: max bitscore of (canonical a, canonical b) over both directions;
* ``best_iso``:  max bitscore over all isoform pairs of the two genes (>= canonical
  whenever the canonical pair is reported);
* control: ``decoy`` isoform sets in which a gene's non-canonical isoforms are
  replaced by the same number of non-canonical isoforms drawn from *other* genes of
  the same family (:func:`decoy_isoform_sets`).
"""
from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd

_TAG = re.compile(r"(\S+):(\S+)")


def parse_ensembl_pep(path: str, keep_seq: bool = True) -> pd.DataFrame:
    """Ensembl pep.all FASTA -> protein_id, gene_id, transcript_id, length, seq.

    Ensembl IDs (ENS*) and gene IDs lose their version suffix; FlyBase / WormBase
    protein IDs (FBpp*, C. elegans ``Y110A7A.10.1``) are kept as is.
    """
    rows = []
    op = gzip.open if path.endswith(".gz") else open
    pid = gid = tid = None
    seq = []

    def flush():
        if pid is not None:
            s = "".join(seq)
            rows.append((pid, gid, tid, len(s), s if keep_seq else None))
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                parts = line[1:].split()
                pid = parts[0].split(".")[0] if parts[0].startswith("ENS") else parts[0]
                tags = dict(m.groups() for m in map(_TAG.fullmatch, parts[1:]) if m)
                gid = tags.get("gene", "").split(".")[0] or None
                tid = tags.get("transcript", "") or None
                if tid and tid.startswith("ENS"):
                    tid = tid.split(".")[0]
                seq = []
            else:
                seq.append(line.strip())
    flush()
    return pd.DataFrame(rows, columns=["protein_id", "gene_id", "transcript_id", "length", "seq"])


def gene_pair_scores(q, t, bits, prot_gene: dict, canonical: set, pairs) -> pd.DataFrame:
    """Per gene pair (ga, gb): canonical / best-isoform bits, and the best isoform pair.

    q, t, bits: protein-level hits (any direction); prot_gene: protein -> gene;
    canonical: set of canonical protein IDs; pairs: iterable of (gene_a, gene_b).
    Missing hits score 0.
    """
    h = pd.DataFrame({"q": q, "t": t, "bits": np.asarray(bits, float)})
    h["ga"] = h.q.map(prot_gene); h["gb"] = h.t.map(prot_gene)
    h = h[h.ga.notna() & h.gb.notna() & (h.ga != h.gb)]
    # orient every hit as (lower gene, higher gene) so both directions pool
    sw = h.ga > h.gb
    h["g1"] = np.where(sw, h.gb, h.ga); h["g2"] = np.where(sw, h.ga, h.gb)
    h["p1"] = np.where(sw, h.t, h.q); h["p2"] = np.where(sw, h.q, h.t)
    h["canon"] = h.p1.isin(canonical) & h.p2.isin(canonical)
    best = h.sort_values("bits", ascending=False).drop_duplicates(["g1", "g2"])
    can = h[h.canon].groupby(["g1", "g2"]).bits.max()
    best = best.set_index(["g1", "g2"])
    out = []
    for ga, gb in pairs:
        k = (ga, gb) if ga < gb else (gb, ga)
        c = float(can.get(k, 0.0))
        if k in best.index:
            r = best.loc[k]
            b = float(r.bits)
            pa, pb = (r.p1, r.p2) if ga < gb else (r.p2, r.p1)
        else:
            b, pa, pb = 0.0, None, None
        out.append((ga, gb, c, max(b, c), pa, pb))
    df = pd.DataFrame(out, columns=["gene_a", "gene_b", "canonical", "best_iso", "best_prot_a", "best_prot_b"])
    df["best_is_canonical"] = df.best_iso <= df.canonical
    return df


def decoy_isoform_sets(iso: dict, canon_of: dict, family_of: dict, species_of: dict, seed: int = 11):
    """Control isoform sets: canonical + k decoys from other genes of the same family.

    iso: gene -> list of protein IDs (incl. canonical); canon_of: gene -> canonical protein;
    family_of / species_of: gene -> family / species. Decoys are non-canonical isoforms of
    other same-family genes of the same species (fallback: any other same-family gene).
    Returns gene -> list of protein IDs (canonical first).
    """
    rng = np.random.default_rng(seed)
    pool_fs, pool_f = {}, {}
    for g, ps in iso.items():
        alts = [p for p in ps if p != canon_of.get(g)]
        for p in alts:
            pool_fs.setdefault((family_of.get(g), species_of.get(g)), []).append((g, p))
            pool_f.setdefault(family_of.get(g), []).append((g, p))
    out = {}
    for g in sorted(iso):
        k = len(iso[g]) - 1
        c = canon_of.get(g)
        res = [c] if c else []
        if k > 0:
            cand = [p for h, p in pool_fs.get((family_of.get(g), species_of.get(g)), []) if h != g]
            if not cand:
                cand = [p for h, p in pool_f.get(family_of.get(g), []) if h != g]
            if cand:
                res += list(rng.choice(cand, size=min(k, len(cand)), replace=False))
        out[g] = res
    return out


def set_pair_scores(hits: pd.DataFrame, sets: dict, pairs) -> np.ndarray:
    """Max bits over protein pairs (sets[ga] x sets[gb]); hits has columns q, t, bits."""
    d = {}
    for a, b, s in zip(hits.q.values, hits.t.values, hits.bits.values):
        k = (a, b) if a < b else (b, a)
        if s > d.get(k, 0.0):
            d[k] = s
    out = np.zeros(len(pairs))
    for i, (ga, gb) in enumerate(pairs):
        best = 0.0
        for pa in sets.get(ga, []):
            for pb in sets.get(gb, []):
                if pa == pb:
                    continue
                k = (pa, pb) if pa < pb else (pb, pa)
                v = d.get(k, 0.0)
                if v > best:
                    best = v
        out[i] = best
    return out
