#!/usr/bin/env python3
"""Verify a pipeline OrthoXML output (DESIGN.md §6).

Checks (exit status 1 if any fails; JSON report with ``--json``):

``xsd``                  XSD validity (``ancestra.orthoxml.validate_xsd``, schema from the
                         document version).
``generef_integrity``    gene ``id`` and ``protId`` unique; every ``geneRef`` points to a
                         declared gene; group ids unique.
``genes_once``           every gene is referenced by at most one ``geneRef`` (genes in no
                         group are singletons); with the input gene table (``--genes``,
                         default ``<run>/stages/input/genes.tsv``) the declared genes are
                         exactly the input proteins, each in its own species.
``taxrange_lca``         every orthologGroup's TaxRange == LCA of its member species in the
                         guide tree (``--tree``, default ``<run>/stages/input/guide_tree.nwk``)
                         and agrees with its ``taxonId``; nested orthologGroups lie strictly
                         below their enclosing orthologGroup; a paralogGroup's level contains
                         all its members and lies strictly below the enclosing orthologGroup.
``pairs_count``          number of ortholog pairs implied by the HOGs (QfO flattening,
                         counted, not stored) == lines of ``pairs.tsv``.
``flattening_vs_qfo``    on a seeded sample of top-level groups (+ the largest ones): pairs
                         from the QfO benchmark's own ``map_relations.py``
                         (``PairwiseOrthologRelationExtractor`` + ``parse_orthoxml`` run on
                         this file) == ``ancestra.orthoxml.hogs_to_pairs`` == ``pairs.tsv``.
                         Uses the QfO 2022 ``mapping.json.gz`` if it covers every protId
                         (QfO run), else a synthetic mapping (toy / OrthoBench).
``pyham``                pyham loads a ``pyham_compat`` copy (full file, or the sample for
                         large inputs, ``--pyham``) with the guide tree; #top-level HOGs ==
                         #top-level orthologGroups of that copy.
``orthobench_export``    ``orthobench_ogs.txt`` read by the official Open_Orthobench
                         ``benchmark.py`` reader (``read_orthogroups``, ``get_n_col_skip``):
                         no gene twice, every gene known, OGs == multi-gene top-level groups.
                         For OrthoBench inputs the official ``check_orthogroups`` is run too.

    python scripts/verify_orthoxml.py RUN_DIR [--json out.json] [--sample 3000] \
        [--pyham auto|full|sample|off] [--seed 11]
"""
from __future__ import annotations

import argparse
import contextlib
import gzip
import importlib.util
import io
import json
import logging
import os
import random
import sys
import tempfile
import time
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

QFO_DIR = os.path.join(ROOT, "data/qfo/benchmark-webservice")
QFO_MAPPING = os.path.join(ROOT, "data/qfo/reference_data/2022.1/mapping.json.gz")
OB_DIR = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS")
NS = "{http://orthoXML.org/2011/}"


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ----------------------------------------------------------------------------
# streaming structure checks
# ----------------------------------------------------------------------------

def check_structure(path: str, gene_table: str | None) -> dict:
    """generef_integrity + genes_once in one streaming pass."""
    from lxml import etree
    ids, prot, species_of = {}, set(), {}
    dup_id, dup_prot, dangling, gids = [], [], [], Counter()
    refcount: Counter = Counter()
    cur_sp = None
    for ev, el in etree.iterparse(path, events=("start", "end")):
        tag = el.tag.replace(NS, "")
        if ev == "start":
            if tag == "species":
                cur_sp = el.get("name")
            elif tag in ("orthologGroup", "paralogGroup"):
                gids[el.get("id")] += 1
            continue
        if tag == "gene":
            i, p = el.get("id"), el.get("protId")
            if i in ids:
                dup_id.append(i)
            if p in prot:
                dup_prot.append(p)
            ids[i] = p
            prot.add(p)
            species_of[p] = cur_sp
            el.clear()
        elif tag == "geneRef":
            i = el.get("id")
            if i not in ids:
                dangling.append(i)
            refcount[i] += 1
        elif tag in ("orthologGroup", "paralogGroup") and el.getparent() is not None \
                and el.getparent().tag == NS + "groups":
            el.clear()
    dup_gid = [g for g, c in gids.items() if c > 1]
    multi = [(ids.get(i), c) for i, c in refcount.items() if c > 1]
    out = {"generef_integrity": {"ok": not (dup_id or dup_prot or dangling or dup_gid),
                                 "n_genes_declared": len(ids), "n_generefs": sum(refcount.values()),
                                 "n_groups": sum(gids.values()), "dup_gene_ids": dup_id[:10],
                                 "dup_protids": dup_prot[:10], "dangling_generefs": dangling[:10],
                                 "dup_group_ids": dup_gid[:10]},
           "genes_once": {"n_genes_in_groups": len(refcount),
                          "n_singletons": len(ids) - len(refcount),
                          "n_genes_referenced_more_than_once": len(multi),
                          "examples": multi[:10]}}
    ok = not multi
    if gene_table and os.path.exists(gene_table):
        exp = {}
        with open(gene_table) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                _, p, sp, _ = line.rstrip("\n").split("\t")
                exp[p] = sp
        missing = [p for p in exp if p not in prot]
        extra = [p for p in prot if p not in exp]
        wrong_sp = [p for p in exp if p in species_of and species_of[p] != exp[p]]
        out["genes_once"].update(gene_table=gene_table, n_input_genes=len(exp),
                                 n_missing=len(missing), n_extra=len(extra),
                                 n_wrong_species=len(wrong_sp),
                                 examples_missing=missing[:5], examples_extra=extra[:5],
                                 examples_wrong_species=wrong_sp[:5])
        ok = ok and not (missing or extra or wrong_sp)
    out["genes_once"]["ok"] = ok
    return out


# ----------------------------------------------------------------------------
# levels
# ----------------------------------------------------------------------------

def check_levels(ox, tree_root, path: str) -> dict:
    from ancestra.newick import TreeIndex
    from ancestra.orthoxml import Group
    from lxml import etree
    tix = TreeIndex(tree_root)
    sp = {k: v.species for k, v in ox.genes.items()}
    # taxonId attributes of orthologGroups (streamed, by group id) vs taxonomy names
    taxname, og_taxid = {}, {}
    for ev, el in etree.iterparse(path, events=("end",)):
        tag = el.tag.replace(NS, "")
        if tag == "taxon":
            taxname[el.get("id")] = el.get("name")
        elif tag == "orthologGroup" and el.get("taxonId") is not None:
            og_taxid[el.get("id")] = el.get("taxonId")
        elif tag == "gene":
            el.clear()
    bad_lca, bad_taxid, bad_nest, bad_par, unknown = [], [], [], [], []
    n_og = n_pg = 0
    per_level = Counter()
    for top in ox.groups:
        if not isinstance(top, Group):
            continue
        # iterative post-order with species sets; parent orthologGroup level passed down
        stack = [(top, None, 0)]
        memo = {}
        while stack:
            g, parent_lvl, st = stack.pop()
            if st == 0:
                stack.append((g, parent_lvl, 1))
                lvl_down = g.level if g.kind == "ortholog" else parent_lvl
                for c in g.children:
                    if isinstance(c, Group):
                        stack.append((c, lvl_down, 0))
                continue
            s = set()
            for c in g.children:
                s |= memo.pop(id(c)) if isinstance(c, Group) else {sp[c]}
            memo[id(g)] = s
            lca = tix.mrca(s)
            if g.kind == "ortholog":
                n_og += 1
                per_level[g.level] += 1
                if g.level not in tix:
                    unknown.append((g.id, g.level))
                    continue
                if g.level != lca:
                    bad_lca.append((g.id, g.level, lca))
                t = og_taxid.get(g.id)
                if t is not None and taxname.get(t) != g.level:
                    bad_taxid.append((g.id, g.level, taxname.get(t)))
                if parent_lvl is not None and not tix.is_strict_descendant(g.level, parent_lvl):
                    bad_nest.append((g.id, g.level, parent_lvl))
            else:
                n_pg += 1
                if g.level is None or g.level not in tix:
                    continue        # unknown -> readers use the MRCA (always consistent)
                if not tix.is_ancestor(g.level, lca) or (
                        parent_lvl is not None and not tix.is_strict_descendant(g.level, parent_lvl)):
                    bad_par.append((g.id, g.level, lca, parent_lvl))
        memo.clear()
    return {"ok": not (bad_lca or bad_taxid or bad_nest or bad_par or unknown),
            "n_orthologGroups": n_og, "n_paralogGroups": n_pg,
            "n_taxrange_ne_lca": len(bad_lca), "examples_taxrange_ne_lca": bad_lca[:10],
            "n_taxonid_ne_taxrange": len(bad_taxid), "examples_taxonid": bad_taxid[:10],
            "n_nesting_violations": len(bad_nest), "examples_nesting": bad_nest[:10],
            "n_paralog_level_violations": len(bad_par), "examples_paralog": bad_par[:10],
            "n_unknown_levels": len(unknown), "examples_unknown": unknown[:10],
            "orthologGroups_per_level": dict(per_level.most_common())}


def count_pairs(ox) -> int:
    """#ortholog pairs of the QfO flattening, counted from species compositions."""
    from ancestra.orthoxml import Group
    sp = {k: v.species for k, v in ox.genes.items()}
    total = 0
    for top in ox.groups:
        if not isinstance(top, Group):
            continue
        stack = [(top, 0)]
        memo = {}
        while stack:
            g, st = stack.pop()
            if st == 0:
                stack.append((g, 1))
                stack.extend((c, 0) for c in g.children if isinstance(c, Group))
                continue
            kids = [memo.pop(id(c)) if isinstance(c, Group) else Counter({sp[c]: 1})
                    for c in g.children]
            if g.kind == "ortholog":
                acc = Counter()
                for k in kids:
                    n_k = sum(k.values())
                    n_acc = sum(acc.values())
                    total += n_acc * n_k - sum(acc[s] * c for s, c in k.items())
                    acc.update(k)
                memo[id(g)] = acc
            else:
                acc = Counter()
                for k in kids:
                    acc.update(k)
                memo[id(g)] = acc
    return total


# ----------------------------------------------------------------------------
# flattening vs the QfO benchmark code
# ----------------------------------------------------------------------------

def _synthetic_mapping(ox) -> dict:
    by = {}
    for p, ge in ox.genes.items():
        by.setdefault(ge.species, []).append(p)
    mapping, goff, species = {}, [0], []
    k = 0
    for s in sorted(by):
        for p in sorted(by[s]):
            k += 1
            mapping[p] = k
        goff.append(k)
        species.append(s[:5])
    return {"mapping": mapping, "Goff": goff, "species": species, "excluded_ids": []}


def qfo_pairs_for(path: str, mapping: dict, sample_ids: set, prot_ids) -> tuple[set, bool, dict]:
    """Run the QfO map_relations extractor on ``path``, keeping only top-level groups in
    ``sample_ids`` (by id). Returns (pairs as sorted protId tuples, mapping_ok, stats)."""
    if QFO_DIR not in sys.path:
        sys.path.insert(0, QFO_DIR)
    mr = _load_module("qfo_map_relations", os.path.join(QFO_DIR, "map_relations.py"))
    pairs_nr: set = set()

    class Collect:
        def add_orthologs(self, a, b):
            pairs_nr.add((a, b) if a < b else (b, a))

    class Sampled(mr.PairwiseOrthologRelationExtractor):
        def extract_pairwise_relations(self, node):
            if node.get("id") in sample_ids:
                super().extract_pairwise_relations(node)
            else:
                self.processed_stats["processed_toplevel"] += 1

    proc = Sampled(mapping, Collect())
    logging.getLogger("relations-processor").setLevel(logging.ERROR)
    ok = True
    try:
        with open(path, "rb") as fh:
            mr.parse_orthoxml(fh, proc)
    except SystemExit:
        ok = False
    # several xrefs map to one protein number: invert through the file's own protIds
    inv = {mapping["mapping"][p]: p for p in prot_ids if p in mapping["mapping"]}
    pairs = set()
    for a, b in pairs_nr:
        x, y = inv[a], inv[b]
        pairs.add((x, y) if x < y else (y, x))
    return pairs, ok, {"n_toplevel_seen": proc.processed_stats["processed_toplevel"]}


def check_flattening(path: str, ox, pairs_tsv: str, n_sample: int, n_largest: int, seed: int,
                     mapping_path: str) -> tuple[dict, list]:
    from ancestra.orthoxml import Group, OrthoXML, hogs_to_pairs
    from lxml import etree
    # top-level orthologGroups as the QfO parser sees them (orthologGroups at og depth 0,
    # incl. those nested in top-level paralogGroups) -> id, genes
    units = []
    for top in ox.groups:
        stack = [top]
        while stack:
            g = stack.pop()
            if not isinstance(g, Group):
                continue
            if g.kind == "ortholog":
                units.append(g)
            else:
                stack.extend(g.children)
    # ids in the file: our reader keeps the 'id' attribute
    rng = random.Random(seed)
    by_size = sorted(range(len(units)), key=lambda i: -len(units[i].genes()))
    largest = by_size[:n_largest]
    rest = [i for i in range(len(units)) if i not in set(largest)]
    pick = sorted(set(largest) | set(rng.sample(rest, min(n_sample, len(rest)))))
    sample = [units[i] for i in pick]
    sample_ids = {g.id for g in sample}
    if None in sample_ids:
        raise RuntimeError("orthologGroups without id attribute")
    genes = {x for g in sample for x in g.genes()}

    # ours
    sub = OrthoXML(genes=ox.genes, groups=sample, species_tree=ox.species_tree)
    ours = hogs_to_pairs(sub)
    # pairs.tsv restricted to the sampled genes (top-level groups partition the genes)
    from_tsv = set()
    with open(pairs_tsv) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            a, b = line.split("\t", 2)[:2]
            if a in genes:
                from_tsv.add((a, b) if a < b else (b, a))
    # QfO
    with gzip.open(mapping_path, "rt") if os.path.exists(mapping_path) else contextlib.nullcontext() as fh:
        mp = json.load(fh) if fh is not None else None
    used = "qfo2022 mapping.json.gz"
    if mp is None or any(p not in mp["mapping"] for p in ox.genes):
        mp, used = _synthetic_mapping(ox), "synthetic (protIds not in the QfO 2022 mapping)"
    qfo, map_ok, st = qfo_pairs_for(path, mp, sample_ids, ox.genes)
    res = {"ok": map_ok and ours == qfo == from_tsv, "mapping": used, "qfo_mapping_checks_ok": map_ok,
           "n_toplevel_units": len(units), "n_sampled_units": len(sample),
           "n_largest_included": len(largest),
           "largest_unit_genes": len(units[largest[0]].genes()) if largest else 0,
           "n_sample_genes": len(genes), "n_pairs_qfo": len(qfo), "n_pairs_ours": len(ours),
           "n_pairs_tsv": len(from_tsv), "ours_minus_qfo": len(ours - qfo),
           "qfo_minus_ours": len(qfo - ours), "tsv_sym_diff_ours": len(from_tsv ^ ours), **st}
    return res, sample


# ----------------------------------------------------------------------------
# pyham, OrthoBench
# ----------------------------------------------------------------------------

def check_pyham(ox, groups, tree_path: str, workdir: str) -> dict:
    from ancestra.orthoxml import OrthoXML, _split_toplevel_paralogs, normalize_groups, write_orthoxml
    import pyham
    logging.getLogger("pyham").setLevel(logging.ERROR)
    used = {x for g in groups for x in g.genes()}
    sub = OrthoXML(genes={k: v for k, v in ox.genes.items() if k in used}, groups=groups,
                   species_tree=ox.species_tree, species_taxid=ox.species_taxid)
    q = os.path.join(workdir, "pyham_compat.orthoxml")
    write_orthoxml(sub, q, pyham_compat=True)
    expected = len(_split_toplevel_paralogs(normalize_groups(groups)))
    t0 = time.time()
    ham = pyham.Ham(tree_path, q, use_internal_name=True, tree_format="newick")
    n = len(ham.get_list_top_level_hogs())
    return {"ok": n == expected, "n_top_level_hogs_pyham": n, "n_expected": expected,
            "n_genes": len(sub.genes), "load_s": round(time.time() - t0, 1)}


def check_orthobench(ogs_path: str, ox) -> dict:
    from ancestra.orthoxml import Group
    ob = _load_module("ob_benchmark", os.path.join(OB_DIR, "benchmark.py"))
    ob.__file__ = os.path.join(OB_DIR, "benchmark.py")
    exp_genes = set(ox.genes)
    ncs = ob.get_n_col_skip(ogs_path)
    ogs = ob.read_orthogroups(ogs_path, exp_genes, ncs)
    allg = [g for og in ogs for g in og]
    dups = [g for g, c in Counter(allg).items() if c > 1]
    unknown = [g for g in set(allg) if g not in exp_genes]
    tops = {frozenset(t.genes()) for t in ox.groups if isinstance(t, Group) and len(t.genes()) > 1}
    got = {frozenset(o) for o in ogs}
    res = {"n_ogs_read": len(ogs), "n_col_skip": ncs, "n_duplicate_genes": len(dups),
           "n_unknown_genes": len(unknown), "n_multi_gene_toplevel_groups": len(tops),
           "ogs_equal_toplevel_groups": got == tops}
    ok = not dups and not unknown and got == tops
    try:
        ob_genes = ob.get_expected_genes()
    except Exception:
        ob_genes = set()
    if ob_genes and ob_genes == exp_genes:
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                ob.check_orthogroups(ob.read_orthogroups_smart(ogs_path), ob_genes)
            res["official_check_orthogroups"] = buf.getvalue().strip().splitlines()[-1:]
            ok = ok and "ERROR" not in buf.getvalue()
        except (Exception, SystemExit) as e:
            res["official_check_orthogroups"] = f"failed: {e!r}"
            ok = False
    else:
        res["official_check_orthogroups"] = "skipped (input is not the OrthoBench gene set)"
    res["ok"] = ok
    return res


# ----------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", help="pipeline output dir (hogs.orthoxml, pairs.tsv, orthobench_ogs.txt)")
    ap.add_argument("--orthoxml")
    ap.add_argument("--tree")
    ap.add_argument("--genes")
    ap.add_argument("--mapping", default=QFO_MAPPING)
    ap.add_argument("--sample", type=int, default=3000, help="random top-level groups for the QfO check")
    ap.add_argument("--largest", type=int, default=5, help="also include the N largest groups")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--pyham", default="auto", choices=["auto", "full", "sample", "off"],
                    help="auto = full if <= 300k genes else sample")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    from ancestra.newick import read_newick
    from ancestra.orthoxml import read_orthoxml, validate_xsd

    rd = a.run_dir
    x = a.orthoxml or os.path.join(rd, "hogs.orthoxml")
    tree = a.tree or os.path.join(rd, "stages/input/guide_tree.nwk")
    genes = a.genes or os.path.join(rd, "stages/input/genes.tsv")
    rep = {"orthoxml": os.path.abspath(x), "tree": os.path.abspath(tree), "checks": {},
           "timing_s": {}}
    ck, tm = rep["checks"], rep["timing_s"]

    def timed(name, fn):
        t0 = time.time()
        try:
            r = fn()
        except Exception as e:      # a crash is a failed check, not a crash of the verifier
            r = {"ok": False, "error": repr(e)}
        tm[name] = round(time.time() - t0, 1)
        return r

    def xsd():
        ok, errs = validate_xsd(x)
        return {"ok": ok, "errors": errs[:10]}
    ck["xsd"] = timed("xsd", xsd)
    r = timed("structure", lambda: check_structure(x, genes))
    if "error" in r:
        ck["generef_integrity"] = ck["genes_once"] = r
    else:
        ck.update(r)
    tm["generef_integrity"] = tm["genes_once"] = tm.pop("structure")
    t0 = time.time()
    try:
        ox = read_orthoxml(x)
    except Exception as e:
        ox = None
        ck["read_orthoxml"] = {"ok": False, "error": repr(e)}
    tm["read_orthoxml"] = round(time.time() - t0, 1)
    troot = read_newick(tree)
    if ox is None:      # the remaining checks need the parsed file
        return _finish(rep, a.json)
    ck["taxrange_lca"] = timed("taxrange_lca", lambda: check_levels(ox, troot, x))
    pairs_tsv = os.path.join(rd, "pairs.tsv")

    def pc():
        n = count_pairs(ox)
        with open(pairs_tsv) as fh:
            m = sum(1 for line in fh if not line.startswith("#"))
        return {"ok": n == m, "n_pairs_hogs": n, "n_lines_pairs_tsv": m}
    ck["pairs_count"] = timed("pairs_count", pc)
    sample: list = []

    def fl():
        nonlocal sample
        r, sample = check_flattening(x, ox, pairs_tsv, a.sample, a.largest, a.seed, a.mapping)
        return r
    ck["flattening_vs_qfo"] = timed("flattening_vs_qfo", fl)
    mode = a.pyham
    if mode == "auto":
        mode = "full" if len(ox.genes) <= 300_000 else "sample"
    if mode != "off":
        with tempfile.TemporaryDirectory() as td:
            groups = ox.groups if mode == "full" else (sample or ox.groups[:a.sample])
            ck["pyham"] = timed("pyham", lambda: {**check_pyham(ox, groups, tree, td), "mode": mode})
    ck["orthobench_export"] = timed("orthobench_export",
                                    lambda: check_orthobench(os.path.join(rd, "orthobench_ogs.txt"), ox))
    return _finish(rep, a.json)


def _finish(rep: dict, json_path: str | None) -> int:
    ck, tm = rep["checks"], rep["timing_s"]
    rep["ok"] = all(c.get("ok", False) for c in ck.values())
    for name, c in ck.items():
        brief = {k: v for k, v in c.items() if k not in ("ok", "orthologGroups_per_level")
                 and not (isinstance(v, (list, dict)) and not v)}
        print(f"[verify_orthoxml] {'PASS' if c.get('ok') else 'FAIL'} {name} ({tm.get(name, '')}s): "
              f"{json.dumps(brief, default=str)[:500]}")
    print(f"[verify_orthoxml] overall: {'PASS' if rep['ok'] else 'FAIL'}")
    if json_path:
        os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)
        with open(json_path, "w") as fh:
            json.dump(rep, fh, indent=1, default=str)
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
