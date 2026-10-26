"""Input layer: proteomes -> integer gene IDs, species table, guide tree, families.

Conventions (DESIGN.md §8)
--------------------------
* One FASTA per species in a directory; the species name is the file stem
  (``Homo_sapiens.fa`` -> ``Homo_sapiens``; ``.gz`` allowed; extensions
  ``.fa .fasta .faa .fas .pep``).
* Species are sorted by name; genes keep their file order. Gene ``g`` gets the
  global integer ID ``g`` (0-based, contiguous), so ``species_of[g]`` and
  ``names[g]`` are plain array lookups and every downstream file (FASTA for the
  search tools, score tables, merge state) uses the integer only.
* The protein ID is the first whitespace-delimited token of the header; a
  UniProt ``sp|ACC|NAME`` / ``tr|ACC|NAME`` header is reduced to ``ACC``.
  IDs must be unique over all species.
* Guide tree: rooted Newick whose leaf names are species names and whose
  internal names are taxonomic level names (they become ``TaxRange`` in the
  OrthoXML). Multifurcations are kept (k-way merge). Unnamed internal nodes
  get ``node<k>``; unary nodes are collapsed; leaves without a proteome are
  pruned. Every proteome species must be a leaf (error otherwise).
"""
from __future__ import annotations

import gzip
import os
from dataclasses import dataclass

import numpy as np

from .newick import Node, TreeIndex, parse_newick

__all__ = ["GeneTable", "load_proteomes", "write_fasta", "write_gene_table", "read_gene_table",
           "load_guide_tree", "load_families", "families_from_mapping"]

FASTA_EXT = (".fa", ".fasta", ".faa", ".fas", ".pep")


@dataclass
class GeneTable:
    """All genes of all species with contiguous integer IDs."""
    names: list            # gene id -> protein ID string
    species: list          # species index -> species name
    species_of: np.ndarray  # int32, gene id -> species index
    lengths: np.ndarray     # int32, gene id -> protein length
    seqs: list | None = None  # gene id -> sequence (None when loaded from the table only)

    @property
    def n(self) -> int:
        return len(self.names)

    def index(self) -> dict:
        return {nm: i for i, nm in enumerate(self.names)}

    def genes_of_species(self) -> list[np.ndarray]:
        order = np.argsort(self.species_of, kind="stable")
        bounds = np.searchsorted(self.species_of[order], np.arange(len(self.species) + 1))
        return [order[bounds[i]:bounds[i + 1]] for i in range(len(self.species))]


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _clean_id(header: str) -> str:
    tok = header.split()[0]
    parts = tok.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr"):
        return parts[1]
    return tok


def read_fasta(path: str):
    """Yield (id, sequence) from a (gzipped) FASTA."""
    name, buf = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = _clean_id(line[1:].strip()), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def species_name(path: str) -> str:
    b = os.path.basename(path)
    if b.endswith(".gz"):
        b = b[:-3]
    return os.path.splitext(b)[0]


def load_proteomes(directory: str) -> GeneTable:
    files = sorted(f for f in os.listdir(directory)
                   if f.removesuffix(".gz").endswith(FASTA_EXT))
    if not files:
        raise FileNotFoundError(f"no FASTA files ({', '.join(FASTA_EXT)}) in {directory}")
    files.sort(key=species_name)
    names, seqs, sp_of, species = [], [], [], []
    seen: set = set()
    for si, f in enumerate(files):
        sp = species_name(f)
        if sp in species:
            raise ValueError(f"two proteome files map to species {sp!r}")
        species.append(sp)
        for pid, s in read_fasta(os.path.join(directory, f)):
            if pid in seen:
                raise ValueError(f"duplicate protein ID {pid!r} (species {sp})")
            seen.add(pid)
            names.append(pid)
            seqs.append(s.rstrip("*").replace("*", "X"))
            sp_of.append(si)
    lengths = np.fromiter((len(s) for s in seqs), dtype=np.int32, count=len(seqs))
    return GeneTable(names, species, np.asarray(sp_of, dtype=np.int32), lengths, seqs)


def write_fasta(gt: GeneTable, path: str, genes=None, names: bool = False) -> None:
    """FASTA with integer headers (gene IDs), or protein IDs with ``names=True``;
    optionally a subset of genes."""
    tmp = path + ".tmp"
    it = range(gt.n) if genes is None else genes
    with open(tmp, "w") as fh:
        for g in it:
            h = gt.names[g] if names else int(g)
            fh.write(f">{h}\n{gt.seqs[g]}\n")
    os.replace(tmp, path)


def write_gene_table(gt: GeneTable, path: str) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("#gid\tprotein\tspecies\tlength\n")
        for g in range(gt.n):
            fh.write(f"{g}\t{gt.names[g]}\t{gt.species[gt.species_of[g]]}\t{gt.lengths[g]}\n")
    os.replace(tmp, path)


def read_gene_table(path: str) -> GeneTable:
    names, sps, lens = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            g, nm, sp, ln = line.rstrip("\n").split("\t")
            assert int(g) == len(names), "gene table must be ordered by gid"
            names.append(nm)
            sps.append(sp)
            lens.append(int(ln))
    species = list(dict.fromkeys(sps))
    sidx = {s: i for i, s in enumerate(species)}
    return GeneTable(names, species, np.array([sidx[s] for s in sps], dtype=np.int32),
                     np.array(lens, dtype=np.int32))


# ----------------------------------------------------------------------------
# guide tree
# ----------------------------------------------------------------------------

def load_guide_tree(path_or_newick: str, species: list[str]) -> Node:
    """Read, name, prune and collapse the guide tree (see module docstring)."""
    if os.path.exists(path_or_newick):
        with open(path_or_newick) as fh:
            root = parse_newick(fh.read())
    else:
        root = parse_newick(path_or_newick)
    leaves = set(root.leaf_names())
    missing = [s for s in species if s not in leaves]
    if missing:
        raise ValueError(f"{len(missing)} proteome species missing from the guide tree, "
                         f"e.g. {missing[:5]}")
    keep = set(species)
    # prune leaves without proteome (post-order)
    for nd in list(root.postorder()):
        if nd.parent is None:
            continue
        if (not nd.children and nd.name not in keep):
            nd.parent.children.remove(nd)
            nd.parent = None
    # collapse unary nodes and drop internal nodes that lost all children
    changed = True
    while changed:
        changed = False
        for nd in list(root.postorder()):
            if nd.children or nd.name in keep:
                if len(nd.children) == 1:
                    ch = nd.children[0]
                    if nd.parent is None:
                        ch.parent = None
                        root = ch
                    else:
                        p = nd.parent
                        i = p.children.index(nd)
                        p.children[i] = ch
                        ch.parent = p
                        if ch.length is not None and nd.length is not None:
                            ch.length += nd.length
                    changed = True
            elif nd.parent is not None:      # internal node with no children left
                nd.parent.children.remove(nd)
                changed = True
    k = 0
    used = {n.name for n in root.preorder() if n.name}
    for nd in root.preorder():
        if not nd.name:
            while f"node{k}" in used:
                k += 1
            nd.name = f"node{k}"
            used.add(nd.name)
    TreeIndex(root)   # checks unique names
    return root


# ----------------------------------------------------------------------------
# families
# ----------------------------------------------------------------------------

def families_from_mapping(gt: GeneTable, gene2fam: dict) -> np.ndarray:
    """int32 family index per gene; genes absent from the mapping become singletons."""
    fam_idx: dict = {}
    out = np.empty(gt.n, dtype=np.int32)
    n_missing = 0
    for g, nm in enumerate(gt.names):
        f = gene2fam.get(nm)
        if f is None:
            f = gene2fam.get(str(g))
        if f is None:
            f = ("__singleton__", g)
            n_missing += 1
        out[g] = fam_idx.setdefault(f, len(fam_idx))
    return out


def load_families(path: str, gt: GeneTable, gene_col: int = 0, fam_col: int = 1) -> np.ndarray:
    """Read a gene->family TSV (lines starting with '#' and an unknown-gene header are skipped).

    Genes may be given by protein ID or by integer gene ID.
    """
    idx = gt.index()
    m: dict = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            g, f = parts[gene_col], parts[fam_col]
            if i == 0 and g not in idx and not g.isdigit():
                continue      # header
            m[g] = f
    unknown = [g for g in m if g not in idx and not (g.isdigit() and int(g) < gt.n)]
    if len(unknown) > 0.5 * max(1, len(m)):
        raise ValueError(f"{len(unknown)} of {len(m)} genes in {path} are unknown (e.g. {unknown[:3]})")
    return families_from_mapping(gt, m)
