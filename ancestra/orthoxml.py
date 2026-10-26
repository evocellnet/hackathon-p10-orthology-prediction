"""OrthoXML (v0.3/0.4/0.5) writer, reader, and flattening utilities.

In-memory model
---------------
* :class:`GeneEntry` -- one gene/protein (``prot_id`` is the unique key, and is
  what ends up in the ``protId`` attribute; QfO requires UniProt accessions
  there).
* :class:`Group` -- a nested group. ``kind`` is ``"ortholog"`` (speciation /
  HOG node) or ``"paralog"`` (duplication node). ``children`` holds nested
  :class:`Group` objects and/or gene ``prot_id`` strings.
  ``level`` is the taxonomic level:

  - for an orthologGroup: the species-tree node of the speciation (written as
    ``<property name="TaxRange" value=...>`` and, in v0.5 when a taxonomy is
    present, as ``taxonId``);
  - for a paralogGroup: the species-tree node **directly below the branch on
    which the duplication happened** (i.e. the first level at which the
    copies are separate HOGs). Written only as the v0.5 ``taxonId`` attribute
    (never as a TaxRange property: pyham attaches properties of a
    paralogGroup to the enclosing orthologGroup, which would corrupt it).
    If unknown, readers fall back to the MRCA of the species below the group
    (latest-possible placement, the same convention pyham uses).
* :class:`OrthoXML` -- genes + top-level groups + optional species tree.

Flattening semantics (QfO benchmark, ``map_relations.py``)
---------------------------------------------------------
The QfO benchmark service (qfo/benchmark-webservice, ``map_relations.py``,
``PairwiseOrthologRelationExtractor.extract_pairwise_relations``) converts
HOGs to pairs as follows: for every ``orthologGroup`` (at any depth, including
orthologGroups nested in a top-level paralogGroup), for every two *distinct*
children of that group, every gene of the first child is paired with every
gene of the second child **if the two genes belong to different species**.
Children of a ``paralogGroup`` are never paired with each other. Hence:

* ortholog pair  <=> the deepest group containing both genes is an
  orthologGroup and the genes are from different species;
* paralog pair   <=> the deepest group containing both genes is a
  paralogGroup (any species);
* a same-species pair whose deepest common group is an orthologGroup is
  dropped by QfO (it is an inconsistent prediction) -- we report it as
  neither;
* genes in different top-level groups are unrelated in this file.

:func:`hogs_to_pairs` implements exactly this (tested against the QfO code).

Groups with fewer than two children are invalid in the XSD (``minOccurs=2``);
the writer collapses single-child groups into their child and omits
top-level singletons (the genes are still listed in the species blocks).
"""
from __future__ import annotations

import io
import itertools
import os
from dataclasses import dataclass, field
from typing import IO, Iterable, Iterator, Union
from xml.sax.saxutils import escape, quoteattr

from .newick import Node, TreeIndex, parse_newick

__all__ = [
    "GeneEntry", "Group", "OrthoXML", "write_orthoxml", "read_orthoxml",
    "hogs_to_pairs", "hogs_at_level", "to_orthobench_ogs", "validate_xsd",
    "check_pyham", "normalize_groups", "ORTHOXML_NS",
]

ORTHOXML_NS = "http://orthoXML.org/2011/"
SCHEMA_DIR = os.path.join(os.path.dirname(__file__), "schemas")

Child = Union["Group", str]


@dataclass
class GeneEntry:
    prot_id: str
    species: str
    gene_id: str | None = None
    transcript_id: str | None = None


@dataclass
class Group:
    kind: str                       # "ortholog" | "paralog"
    children: list = field(default_factory=list)
    level: str | None = None
    id: str | None = None
    properties: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.kind not in ("ortholog", "paralog"):
            raise ValueError(f"Group.kind must be 'ortholog' or 'paralog', got {self.kind!r}")

    def iter_genes(self) -> Iterator[str]:
        stack: list = [self]
        while stack:
            x = stack.pop()
            if isinstance(x, Group):
                stack.extend(reversed(x.children))
            else:
                yield x

    def genes(self) -> list[str]:
        return list(self.iter_genes())

    def iter_groups(self) -> Iterator["Group"]:
        stack: list = [self]
        while stack:
            x = stack.pop()
            if isinstance(x, Group):
                yield x
                stack.extend(c for c in reversed(x.children) if isinstance(c, Group))


@dataclass
class OrthoXML:
    genes: dict                     # prot_id -> GeneEntry (insertion-ordered)
    groups: list                    # top-level Group or gene ids
    species_tree: Node | None = None
    origin: str = "ancestra"
    origin_version: str = "0.1"
    database: str = "ancestra"
    database_version: str = "1"
    species_taxid: dict = field(default_factory=dict)   # species -> NCBI taxid (optional)

    @property
    def species(self) -> list[str]:
        seen: dict[str, None] = {}
        for g in self.genes.values():
            seen.setdefault(g.species, None)
        if self.species_tree is not None:
            order = [n for n in self.species_tree.leaf_names() if n in seen]
            return order + [s for s in seen if s not in set(order)]
        return list(seen)

    def species_of(self, prot_id: str) -> str:
        return self.genes[prot_id].species

    def tree_index(self) -> TreeIndex | None:
        return TreeIndex(self.species_tree) if self.species_tree is not None else None


# ----------------------------------------------------------------------------
# normalisation
# ----------------------------------------------------------------------------

def _normalize(x: Child) -> Child | None:
    """Collapse groups with <2 children (post-order, iterative)."""
    if not isinstance(x, Group):
        return x
    # iterative post-order rebuild
    stack: list[tuple[Group, int]] = [(x, 0)]
    results: dict[int, Child | None] = {}
    while stack:
        g, st = stack.pop()
        if st == 0:
            stack.append((g, 1))
            for c in g.children:
                if isinstance(c, Group):
                    stack.append((c, 0))
        else:
            new_children = []
            for c in g.children:
                r = results.pop(id(c)) if isinstance(c, Group) else c
                if r is not None:
                    new_children.append(r)
            if len(new_children) == 0:
                results[id(g)] = None
            elif len(new_children) == 1:
                results[id(g)] = new_children[0]
            else:
                ng = Group(g.kind, new_children, g.level, g.id, dict(g.properties))
                results[id(g)] = ng
    return results[id(x)]


def normalize_groups(groups: Iterable[Child]) -> list[Group]:
    """Collapse single-child groups and drop top-level singletons."""
    out = []
    for g in groups:
        r = _normalize(g)
        if isinstance(r, Group):
            out.append(r)
    return out


# ----------------------------------------------------------------------------
# writer
# ----------------------------------------------------------------------------

def write_orthoxml(ox: OrthoXML, dest: Union[str, IO[str]], version: str = "0.5",
                   include_taxonomy: bool = True, group_id_prefix: str = "HOG",
                   pyham_compat: bool = False) -> None:
    """Write an :class:`OrthoXML` object.

    Parameters
    ----------
    version : "0.5" (default; adds ``<taxonomy>`` and ``taxonId`` when a species
        tree is attached), "0.4" or "0.3" (no taxonomy; v0.3 forbids top-level
        paralogGroups, which are then written as their children).
    pyham_compat : pyham (1.1.x) crashes on a top-level paralogGroup that
        directly contains a geneRef. With ``True``, every top-level
        paralogGroup is replaced by its child groups (recursively; direct
        genes become singletons). Ortholog pairs (QfO flattening) are
        unchanged; only paralog relations across the split are lost. Also
        implied by ``version="0.3"``.
    """
    if version not in ("0.3", "0.4", "0.5"):
        raise ValueError("version must be 0.3, 0.4 or 0.5")
    own = isinstance(dest, str)
    fh = open(dest, "w") if own else dest
    try:
        _write(ox, fh, version, include_taxonomy, group_id_prefix, pyham_compat)
    finally:
        if own:
            fh.close()


def _split_toplevel_paralogs(groups: list[Group]) -> list[Group]:
    out: list[Group] = []
    stack = list(reversed(groups))
    while stack:
        g = stack.pop()
        if g.kind == "paralog":
            stack.extend(reversed([c for c in g.children if isinstance(c, Group)]))
        else:
            out.append(g)
    return out


def _write(ox: OrthoXML, fh: IO[str], version: str, include_taxonomy: bool, prefix: str,
           pyham_compat: bool = False) -> None:
    groups = normalize_groups(ox.groups)
    if version == "0.3" or pyham_compat:
        groups = _split_toplevel_paralogs(groups)
    use_tax = version == "0.5" and include_taxonomy and ox.species_tree is not None
    taxid: dict[str, int] = {}
    if use_tax:
        for k, nd in enumerate(ox.species_tree.preorder(), start=1):
            if not nd.name:
                raise ValueError("species tree nodes must all be named to write a taxonomy")
            taxid[nd.name] = k
    gene_num: dict[str, int] = {}
    by_sp: dict[str, list[GeneEntry]] = {}
    for ge in ox.genes.values():
        by_sp.setdefault(ge.species, []).append(ge)
    w = fh.write
    w('<?xml version="1.0" encoding="UTF-8"?>\n')
    w(f'<orthoXML xmlns="{ORTHOXML_NS}" origin={quoteattr(ox.origin)} '
      f'originVersion={quoteattr(ox.origin_version)} version="{version}">\n')
    k = 0
    for sp in ox.species:
        attrs = f"name={quoteattr(sp)}"
        if sp in ox.species_taxid:
            attrs += f' NCBITaxId="{int(ox.species_taxid[sp])}"'
        elif version != "0.5":
            attrs += ' NCBITaxId="0"'   # required in 0.3/0.4
        if use_tax and sp in taxid:
            attrs += f' taxonId="{taxid[sp]}"'
        w(f"  <species {attrs}>\n")
        w(f"    <database name={quoteattr(ox.database)} version={quoteattr(ox.database_version)}>\n")
        w("      <genes>\n")
        for ge in by_sp[sp]:
            k += 1
            gene_num[ge.prot_id] = k
            a = f'id="{k}" protId={quoteattr(ge.prot_id)}'
            if ge.gene_id:
                a += f" geneId={quoteattr(ge.gene_id)}"
            if ge.transcript_id:
                a += f" transcriptId={quoteattr(ge.transcript_id)}"
            w(f"        <gene {a}/>\n")
        w("      </genes>\n    </database>\n  </species>\n")
    if use_tax:
        w("  <taxonomy>\n")
        _write_taxonomy(ox.species_tree, taxid, w)
        w("  </taxonomy>\n")
    w("  <groups>\n")
    used_ids: set[str] = set()
    for gi, g in enumerate(groups, start=1):
        top_id = g.id or f"{prefix}{gi:07d}"
        _write_group(g, top_id, gene_num, taxid if use_tax else None, w, used_ids)
    w("  </groups>\n</orthoXML>\n")


def _write_taxonomy(root: Node, taxid: dict[str, int], w) -> None:
    stack: list[tuple[Node, int, int]] = [(root, 0, 2)]
    while stack:
        nd, st, ind = stack.pop()
        pad = " " * (2 * ind)
        if st == 0:
            if nd.children:
                w(f'{pad}<taxon id="{taxid[nd.name]}" name={quoteattr(nd.name)}>\n')
                stack.append((nd, 1, ind))
                for c in reversed(nd.children):
                    stack.append((c, 0, ind + 1))
            else:
                w(f'{pad}<taxon id="{taxid[nd.name]}" name={quoteattr(nd.name)}/>\n')
        else:
            w(f"{pad}</taxon>\n")


def _write_group(top: Group, top_id: str, gene_num, taxid, w, used_ids: set[str]) -> None:
    counter = itertools.count(1)
    stack: list[tuple[Child, int, int, str | None]] = [(top, 0, 2, top_id)]
    while stack:
        x, st, ind, gid = stack.pop()
        pad = " " * (2 * ind)
        if not isinstance(x, Group):
            if x not in gene_num:
                raise KeyError(f"gene {x!r} referenced in a group but not declared in genes")
            w(f'{pad}<geneRef id="{gene_num[x]}"/>\n')
            continue
        tag = "orthologGroup" if x.kind == "ortholog" else "paralogGroup"
        if st == 1:
            w(f"{pad}</{tag}>\n")
            continue
        if gid is None:
            gid = x.id or f"{top_id}.{next(counter)}"
        while gid in used_ids:
            gid = f"{gid}_{next(counter)}"
        used_ids.add(gid)
        attrs = f"id={quoteattr(gid)}"
        if taxid is not None and x.level is not None and x.level in taxid:
            attrs += f' taxonId="{taxid[x.level]}"'
        w(f"{pad}<{tag} {attrs}>\n")
        props = dict(x.properties)
        if x.kind == "ortholog" and x.level is not None:
            props = {"TaxRange": x.level, **{k: v for k, v in props.items() if k != "TaxRange"}}
        elif x.kind == "paralog":
            props.pop("TaxRange", None)   # see module docstring (pyham)
        for kname, val in props.items():
            w(f"{pad}  <property name={quoteattr(str(kname))} value={quoteattr(str(val))}/>\n")
        stack.append((x, 1, ind, gid))
        for c in reversed(x.children):
            stack.append((c, 0, ind + 1, None))


# ----------------------------------------------------------------------------
# reader
# ----------------------------------------------------------------------------

def _local(tag) -> str:
    if not isinstance(tag, str):
        return ""
    return tag.rsplit("}", 1)[-1]


def read_orthoxml(src: Union[str, IO[bytes]]) -> OrthoXML:
    """Parse an OrthoXML file (any of v0.3/0.4/0.5, gz allowed) into :class:`OrthoXML`.

    Gene keys are ``protId`` (falling back to ``geneId``, then the internal
    ``id``). If the file has a v0.5 ``<taxonomy>`` it becomes ``species_tree``.
    Group levels come from the ``TaxRange`` property (orthologGroups) or the
    ``taxonId`` attribute (mapped through the taxonomy).
    """
    from lxml import etree

    if isinstance(src, str) and src.endswith(".gz"):
        import gzip
        fh = gzip.open(src, "rb")
        own = True
    elif isinstance(src, str):
        fh = open(src, "rb")
        own = True
    else:
        fh, own = src, False
    genes: dict[str, GeneEntry] = {}
    internal: dict[str, str] = {}
    species_taxon: dict[str, str] = {}
    species_taxid: dict[str, int] = {}
    taxa: dict[str, Node] = {}
    tax_stack: list[Node] = []
    tax_root: Node | None = None
    groups: list = []
    gstack: list[Group] = []
    group_taxonid: dict[int, str] = {}
    cur_species = None
    root_attrs: dict = {}
    db_attrs: dict = {}
    try:
        for event, el in etree.iterparse(fh, events=("start", "end")):
            tag = _local(el.tag)
            if event == "start":
                if tag == "orthoXML":
                    root_attrs = dict(el.attrib)
                elif tag == "species":
                    cur_species = el.get("name")
                    if el.get("taxonId") is not None:
                        species_taxon[cur_species] = el.get("taxonId")
                    if el.get("NCBITaxId") not in (None, "0"):
                        species_taxid[cur_species] = int(el.get("NCBITaxId"))
                elif tag == "database" and not db_attrs:
                    db_attrs = dict(el.attrib)
                elif tag == "taxon":
                    nd = Node(el.get("name") or el.get("id"))
                    nd.data["taxon_id"] = el.get("id")
                    taxa[el.get("id")] = nd
                    if tax_stack:
                        tax_stack[-1].add_child(nd)
                    else:
                        tax_root = nd
                    tax_stack.append(nd)
                elif tag in ("orthologGroup", "paralogGroup"):
                    g = Group("ortholog" if tag == "orthologGroup" else "paralog", id=el.get("id"))
                    if el.get("taxonId") is not None:
                        group_taxonid[id(g)] = el.get("taxonId")
                    if gstack:
                        gstack[-1].children.append(g)
                    gstack.append(g)
            else:  # end
                if tag == "gene":
                    pid = el.get("protId") or el.get("geneId") or el.get("id")
                    if pid in genes:
                        raise ValueError(f"duplicate gene identifier {pid!r} in OrthoXML")
                    genes[pid] = GeneEntry(pid, cur_species, el.get("geneId"), el.get("transcriptId"))
                    internal[el.get("id")] = pid
                    el.clear()
                elif tag == "species":
                    cur_species = None
                    el.clear()
                elif tag == "taxon":
                    tax_stack.pop()
                elif tag == "geneRef" and gstack:
                    gstack[-1].children.append(internal[el.get("id")])
                elif tag == "property" and gstack:
                    # properties belong to the innermost open group
                    gstack[-1].properties[el.get("name")] = el.get("value")
                elif tag in ("orthologGroup", "paralogGroup"):
                    g = gstack.pop()
                    if not gstack:
                        groups.append(g)
                        el.clear()
    finally:
        if own:
            fh.close()
    # resolve levels
    for top in groups:
        for g in top.iter_groups():
            tr = g.properties.get("TaxRange")
            if tr is not None:
                g.level = tr
            elif id(g) in group_taxonid and group_taxonid[id(g)] in taxa:
                g.level = taxa[group_taxonid[id(g)]].name
    species_tree = None
    if tax_root is not None:
        # rename leaf taxa to species names when linked through taxonId
        rev = {v: k for k, v in species_taxon.items()}
        for tid, nd in taxa.items():
            if not nd.children and tid in rev:
                nd.name = rev[tid]
        species_tree = tax_root
    return OrthoXML(genes=genes, groups=groups, species_tree=species_tree,
                    origin=root_attrs.get("origin", ""),
                    origin_version=root_attrs.get("originVersion", ""),
                    database=db_attrs.get("name", "unknown"),
                    database_version=db_attrs.get("version", ""),
                    species_taxid=species_taxid)


def _as_ox(x: Union[OrthoXML, str]) -> OrthoXML:
    return x if isinstance(x, OrthoXML) else read_orthoxml(x)


# ----------------------------------------------------------------------------
# flattening to pairs (QfO semantics)
# ----------------------------------------------------------------------------

def hogs_to_pairs(ox: Union[OrthoXML, str], with_paralogs: bool = False,
                  with_level: bool = False):
    """Flatten nested HOGs into pairwise relations exactly like the QfO benchmark.

    Returns a set of ``(gene_a, gene_b)`` tuples with ``gene_a < gene_b``
    (orthologs). With ``with_level=True`` a dict pair -> level (the level of
    the orthologGroup that induced the pair, may be None) is returned instead.
    With ``with_paralogs=True`` a tuple ``(orthologs, paralogs)`` is returned,
    where paralogs are pairs whose deepest common group is a paralogGroup.
    See module docstring for the exact rules.
    """
    ox = _as_ox(ox)
    sp = {k: v.species for k, v in ox.genes.items()}
    orth: dict | set = {} if with_level else set()
    para: set = set()
    for top in ox.groups:
        if not isinstance(top, Group):
            continue
        # iterative post-order; for each group compute gene list of each child
        stack: list[tuple[Group, int]] = [(top, 0)]
        memo: dict[int, list[str]] = {}
        while stack:
            g, st = stack.pop()
            if st == 0:
                stack.append((g, 1))
                for c in g.children:
                    if isinstance(c, Group):
                        stack.append((c, 0))
                continue
            child_sets = []
            for c in g.children:
                child_sets.append(memo.pop(id(c)) if isinstance(c, Group) else [c])
            if g.kind == "ortholog":
                for a_set, b_set in itertools.combinations(child_sets, 2):
                    for a in a_set:
                        sa = sp[a]
                        for b in b_set:
                            if sa != sp[b]:
                                p = (a, b) if a < b else (b, a)
                                if with_level:
                                    orth[p] = g.level
                                else:
                                    orth.add(p)
            elif with_paralogs:
                for a_set, b_set in itertools.combinations(child_sets, 2):
                    for a in a_set:
                        for b in b_set:
                            if a != b:
                                para.add((a, b) if a < b else (b, a))
            memo[id(g)] = [x for cs in child_sets for x in cs]
    if with_paralogs:
        return orth, para
    return orth


# ----------------------------------------------------------------------------
# HOGs at a given level
# ----------------------------------------------------------------------------

def hogs_at_level(ox: Union[OrthoXML, str], level: str | None = None,
                  species_tree: Node | None = None, include_singletons: bool = True) -> list[set]:
    """Return the HOGs (gene sets) at a taxonomic ``level``.

    ``level=None`` returns one set per top-level group (the gene "families").
    Otherwise a species tree is needed (argument, or the file's taxonomy).
    Only genes of species below ``level`` are returned. Rules (``Y`` = group
    level; if missing or unknown in the tree, ``Y`` = MRCA of the group's
    species):

    * orthologGroup with ``Y`` at or below ``level`` -> one HOG;
    * paralogGroup with ``Y`` strictly below ``level`` (duplication after the
      speciation at ``level``) -> one HOG;
    * otherwise descend into the children separately.

    With ``include_singletons``, genes of the clade not in any multi-gene HOG
    are returned as singleton sets, so the result partitions the clade's genes.
    """
    ox = _as_ox(ox)
    if level is None:
        out = [set(g.genes()) for g in ox.groups if isinstance(g, Group)]
        if include_singletons:
            seen = set().union(*out) if out else set()
            out.extend({gid} for gid in ox.genes if gid not in seen)
        return out
    tree = species_tree if species_tree is not None else ox.species_tree
    if tree is None:
        raise ValueError("hogs_at_level with a level requires a species tree")
    tix = TreeIndex(tree)
    if level not in tix:
        raise KeyError(f"level {level!r} not in species tree")
    clade = tix.clade(level)
    sp = {k: v.species for k, v in ox.genes.items()}
    out: list[set] = []

    def group_level(g: Group, genes: list[str]) -> str:
        if g.level is not None and g.level in tix:
            return g.level
        return tix.mrca(sp[x] for x in genes)

    for top in ox.groups:
        if not isinstance(top, Group):
            continue
        stack: list[Child] = [top]
        while stack:
            x = stack.pop()
            if not isinstance(x, Group):
                if sp[x] in clade:
                    out.append({x})
                continue
            genes = x.genes()
            in_clade = [gname for gname in genes if sp[gname] in clade]
            if not in_clade:
                continue
            y = group_level(x, genes)
            if x.kind == "ortholog":
                one = tix.is_ancestor(level, y)
            else:
                one = tix.is_strict_descendant(y, level)
            if one:
                out.append(set(in_clade))
            else:
                stack.extend(x.children)
    if include_singletons:
        seen = set().union(*out) if out else set()
        out.extend({gid} for gid in ox.genes if sp[gid] in clade and gid not in seen)
    else:
        out = [s for s in out if len(s) > 1]
    return out


def to_orthobench_ogs(ox: Union[OrthoXML, str], level: str | None = None,
                      out: Union[str, IO[str], None] = None, species_tree: Node | None = None,
                      include_singletons: bool = True) -> list[list[str]]:
    """Export HOGs at ``level`` (default: top-level groups) as OrthoBench OGs.

    Output format (Open_Orthobench ``benchmark.py``, both the default regex
    reader and the ``-b`` basic reader accept it): a header line
    ``Orthogroup<TAB>Genes``, then one OG per line,
    ``OG0000000:<TAB>gene1<TAB>gene2...``. OGs are sorted by decreasing size.
    Each gene appears in exactly one OG (benchmark.py rejects duplicates).
    """
    ogs = hogs_at_level(ox, level, species_tree=species_tree, include_singletons=include_singletons)
    ogs = [sorted(s) for s in ogs if include_singletons or len(s) > 1]
    ogs.sort(key=lambda s: (-len(s), s[0]))
    if out is not None:
        own = isinstance(out, str)
        fh = open(out, "w") if own else out
        try:
            fh.write("Orthogroup\tGenes\n")
            for i, og in enumerate(ogs):
                fh.write(f"OG{i:07d}:\t" + "\t".join(og) + "\n")
        finally:
            if own:
                fh.close()
    return ogs


# ----------------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------------

_SCHEMA_CACHE: dict = {}


def validate_xsd(path_or_bytes: Union[str, bytes], version: str | None = None) -> tuple[bool, list[str]]:
    """Validate an OrthoXML document against the official XSD.

    The schema is chosen from the document's ``version`` attribute unless
    ``version`` is given. Returns ``(ok, error_messages)``.
    """
    from lxml import etree

    if isinstance(path_or_bytes, bytes):
        doc = etree.parse(io.BytesIO(path_or_bytes))
    elif path_or_bytes.endswith(".gz"):
        import gzip
        with gzip.open(path_or_bytes, "rb") as fh:
            doc = etree.parse(fh)
    else:
        doc = etree.parse(path_or_bytes)
    v = version or doc.getroot().get("version", "0.5")
    v = {"0.3": "0.3", "0.4": "0.4", "0.5": "0.5"}.get(str(v), "0.5")
    if v not in _SCHEMA_CACHE:
        _SCHEMA_CACHE[v] = etree.XMLSchema(etree.parse(os.path.join(SCHEMA_DIR, f"orthoxml-{v}.xsd")))
    schema = _SCHEMA_CACHE[v]
    ok = schema.validate(doc)
    return ok, [f"line {e.line}: {e.message}" for e in schema.error_log]


def check_pyham(orthoxml_path: str, newick_path: str) -> int:
    """Load an OrthoXML with pyham (species-tree internal names = TaxRange).

    Returns the number of top-level HOGs pyham built. Raises on failure.
    """
    import logging

    import pyham
    logging.getLogger("pyham").setLevel(logging.WARNING)
    ham = pyham.Ham(newick_path, orthoxml_path, use_internal_name=True, tree_format="newick")
    return len(ham.get_list_top_level_hogs())
