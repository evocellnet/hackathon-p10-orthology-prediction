"""S9: expression-profile evidence from Bgee expression calls (exploratory layer).

Input: Bgee ``*_expr_simple.tsv.gz`` (gene x anatomical entity calls, present/absent,
with an expression score 0-100). Anatomical entities are Uberon / CL IDs (or
"CL:x ∩ UBERON:y" composites) that are shared across species by construction, so
two species' profiles are compared over the entities present in *both* species.

Bgee is method *input* (like a proteome), never a label source; the module reads
only a declared expression file and gene IDs.

Main objects
------------
* :class:`ExpressionProfiles` - dense gene x entity matrices (score, call) for one species.
* :func:`pair_similarity` - Spearman (expression score) / Jaccard (present sets) over
  shared entities in which both genes have a call; NaN when < ``min_shared`` entities.
* :func:`shuffled_profiles` - negative control: permute gene <-> profile within species.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass

import numpy as np
import pandas as pd

QUALITIES = ("gold quality", "silver quality")


@dataclass
class ExpressionProfiles:
    genes: np.ndarray          # gene IDs (row order)
    entities: np.ndarray       # anatomical entity IDs (column order)
    score: np.ndarray          # float32 [n_genes, n_entities], NaN = no call
    present: np.ndarray        # bool    [n_genes, n_entities], True = present call
    species: str = ""

    def __post_init__(self):
        self.gene_index = {g: i for i, g in enumerate(self.genes)}

    def rows(self, gene_ids):
        """Row index per gene ID (-1 when the gene has no Bgee calls)."""
        return np.array([self.gene_index.get(g, -1) for g in gene_ids], np.int64)

    def permuted(self, seed: int) -> "ExpressionProfiles":
        return shuffled_profiles(self, seed)


def read_bgee(path: str, qualities=QUALITIES) -> pd.DataFrame:
    """Long table (gene, entity, present, score) from a Bgee ``expr_simple`` file."""
    df = pd.read_csv(path, sep="\t", usecols=[0, 2, 4, 5, 7], header=0,
                     names=["gene", "entity", "call", "quality", "score"],
                     dtype={"gene": str, "entity": str, "call": "category", "quality": "category",
                            "score": np.float32},
                     compression="gzip" if path.endswith(".gz") else None)
    df = df[df["quality"].isin(qualities)]
    df["present"] = (df["call"] == "present").values
    return df[["gene", "entity", "present", "score"]].reset_index(drop=True)


def build_profiles(df: pd.DataFrame, species: str = "") -> ExpressionProfiles:
    """Dense matrices from the long table (one row per gene with >= 1 call)."""
    gcode, genes = pd.factorize(df["gene"], sort=True)
    ecode, ents = pd.factorize(df["entity"], sort=True)
    score = np.full((len(genes), len(ents)), np.nan, np.float32)
    present = np.zeros((len(genes), len(ents)), bool)
    score[gcode, ecode] = df["score"].values
    present[gcode, ecode] = df["present"].values
    return ExpressionProfiles(np.asarray(genes), np.asarray(ents), score, present, species)


def load_profiles(path: str, species: str = "", qualities=QUALITIES) -> ExpressionProfiles:
    return build_profiles(read_bgee(path, qualities), species)


def shared_columns(pa: ExpressionProfiles, pb: ExpressionProfiles):
    """Column indices (in a, in b) of entities present in both species, and their IDs."""
    common = np.intersect1d(pa.entities, pb.entities)
    ia = {e: i for i, e in enumerate(pa.entities)}
    ib = {e: i for i, e in enumerate(pb.entities)}
    return (np.array([ia[e] for e in common], np.int64), np.array([ib[e] for e in common], np.int64),
            common)


def _rank(x):
    """Average ranks (ties averaged), 1-D."""
    return pd.Series(x).rank(method="average").values


def spearman(x, y):
    if len(x) < 3:
        return np.nan
    rx, ry = _rank(x), _rank(y)
    sx, sy = rx.std(), ry.std()
    if sx == 0 or sy == 0:
        return np.nan
    return float(((rx - rx.mean()) * (ry - ry.mean())).mean() / (sx * sy))


def pair_similarity(pa: ExpressionProfiles, pb: ExpressionProfiles, ra, rb, min_shared: int = 5,
                    cols=None):
    """Expression similarity for row pairs (ra[i], rb[i]) (row -1 = no data -> NaN).

    Over the entities shared by both species *and* called in both genes:
    ``spearman`` of Bgee expression scores, ``jaccard`` of present-call sets
    (NaN if neither gene is present anywhere there), ``n_shared`` entities.
    """
    ca, cb, _ = cols if cols is not None else shared_columns(pa, pb)
    ra = np.asarray(ra, np.int64); rb = np.asarray(rb, np.int64)
    n = len(ra)
    rho = np.full(n, np.nan); jac = np.full(n, np.nan); nsh = np.zeros(n, np.int64)
    for i in range(n):
        if ra[i] < 0 or rb[i] < 0:
            continue
        xa = pa.score[ra[i], ca]; xb = pb.score[rb[i], cb]
        m = ~np.isnan(xa) & ~np.isnan(xb)
        k = int(m.sum()); nsh[i] = k
        if k < min_shared:
            continue
        rho[i] = spearman(xa[m], xb[m])
        qa = pa.present[ra[i], ca][m]; qb = pb.present[rb[i], cb][m]
        u = int((qa | qb).sum())
        jac[i] = (qa & qb).sum() / u if u else np.nan
    return pd.DataFrame({"spearman": rho, "jaccard": jac, "n_shared": nsh})


def shuffled_profiles(p: ExpressionProfiles, seed: int) -> ExpressionProfiles:
    """Negative control: gene IDs keep their slot, profiles are permuted among genes."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(p.genes))
    return ExpressionProfiles(p.genes, p.entities, p.score[perm], p.present[perm], p.species + f"_shuf{seed}")


def protein_to_gene(table_parquet: str) -> dict:
    """OrthoBench protein ID -> Ensembl gene ID from a gene-order table (data/genomes/orthobench)."""
    d = pd.read_parquet(table_parquet, columns=["protein_id", "gene_id"])
    return dict(zip(d.protein_id, d.gene_id.astype(str).str.split(".").str[0]))
