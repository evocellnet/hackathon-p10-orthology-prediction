"""Tests for the toy simulator: determinism, file consistency, exact truth, model sanity."""
import glob
import hashlib
import itertools
import os
import time
from collections import defaultdict

import numpy as np
import pytest

from ancestra.newick import read_newick
from ancestra.orthoxml import hogs_to_pairs, read_orthoxml
from ancestra.simulate import AAModel, SimParams, main, simulate, write_simulation


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = []
        elif line:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def dir_digest(d):
    h = hashlib.sha256()
    for f in sorted(glob.glob(os.path.join(d, "**", "*"), recursive=True)):
        if os.path.isfile(f) and not f.endswith("params.json"):
            h.update(os.path.relpath(f, d).encode())
            h.update(open(f, "rb").read())
    return h.hexdigest()


# ---------------------------------------------------------------- substitution model
def test_aa_model_properties():
    m = AAModel(gamma_shape=1.0, n_categories=4)
    assert np.allclose(m.rates.mean(), 1.0)
    assert np.allclose(m.P(0.0), np.eye(20)[None].repeat(4, 0), atol=1e-10)
    for t in (0.05, 0.5, 3.0):
        Pt = m.P(t)
        assert np.allclose(Pt.sum(2), 1.0)
        for k in range(4):
            assert np.allclose(m.pi @ Pt[k], m.pi, atol=1e-8)                     # stationarity
            F = m.pi[:, None] * Pt[k]
            assert np.allclose(F, F.T, atol=1e-8)                                 # reversibility
    # Q normalised to 1 expected substitution per unit time
    assert np.isclose(-(m.pi * np.diag(m.Q)).sum(), 1.0)


def test_substitution_rate_matches_expectation():
    m = AAModel(gamma_shape=1.0, n_categories=4)
    rng = np.random.default_rng(1)
    n, d = 200_000, 0.3
    seq = m.random_seq(n, rng)
    cats = rng.integers(4, size=n).astype(np.uint8)
    new = m.substitute(seq, cats, d, rng)
    Pd = m.P(d)
    expected = np.mean([(m.pi * (1 - np.diag(Pd[k]))).sum() for k in range(4)])
    assert abs((new != seq).mean() - expected) < 0.005


# ---------------------------------------------------------------- determinism + CLI
def test_deterministic_and_seed_dependent(tmp_path):
    p = SimParams(n_species=5, n_genes=30)
    for i, seed in enumerate((11, 11, 47)):
        write_simulation(simulate(p, seed), str(tmp_path / f"r{i}"))
    d = [dir_digest(str(tmp_path / f"r{i}")) for i in range(3)]
    assert d[0] == d[1] and d[0] != d[2]


def test_cli_and_given_tree(tmp_path):
    nwk = tmp_path / "t.nwk"
    nwk.write_text("((a:0.2,b:0.2):0.3,(c:0.1,d:0.1)cd:0.4);")
    out = tmp_path / "o"
    assert main(["--out", str(out), "--seed", "3", "--genes", "20", "--tree", str(nwk)]) == 0
    t = read_newick(str(out / "species_tree.nwk"))
    assert sorted(t.leaf_names()) == ["a", "b", "c", "d"]
    assert all(n.name for n in t.preorder())
    assert t.find("cd") is not None
    assert sorted(os.listdir(out / "proteomes")) == ["a.fa", "b.fa", "c.fa", "d.fa"]


# ---------------------------------------------------------------- file consistency
@pytest.mark.parametrize("which", ["basic", "complex"])
def test_files_consistent(which, request):
    res, out = request.getfixturevalue(f"sim_{which}")
    tree = read_newick(os.path.join(out, "species_tree.nwk"))
    species = tree.leaf_names()
    all_ids = []
    for sp in species:
        seqs = read_fasta(os.path.join(out, "proteomes", f"{sp}.fa"))
        assert all(k.startswith(sp + "_g") for k in seqs)
        assert all(set(s) <= set("ACDEFGHIKLMNPQRSTVWY") and s[0] == "M" for s in seqs.values())
        all_ids.extend(seqs)
        # GFF: same genes, same order, CDS length = 3*(len+1), no overlaps
        gff = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "gff", f"{sp}.gff3"))
               if not l.startswith("#")]
        cds = [r for r in gff if r[2] == "CDS"]
        pids = [dict(kv.split("=") for kv in r[8].split(";"))["protein_id"] for r in cds]
        assert pids == list(seqs)
        last_end = defaultdict(int)
        for r, pid in zip(cds, pids):
            s, e = int(r[3]), int(r[4])
            assert e - s + 1 == 3 * (len(seqs[pid]) + 1)
            assert s > last_end[r[0]]
            last_end[r[0]] = e
        assert {r[2] for r in gff} == {"gene", "mRNA", "CDS"}
    assert len(all_ids) == len(set(all_ids)) == res.stats["n_proteins"]
    # every protein is a leaf of exactly one gene tree
    leaves = []
    for f in glob.glob(os.path.join(out, "truth", "gene_tree_*.nwk")):
        leaves.extend(read_newick(f).leaf_names())
    assert sorted(leaves) == sorted(all_ids)
    # truth orthoxml covers all genes; pairs file == QfO flattening of the truth HOGs
    ox = read_orthoxml(os.path.join(out, "truth", "hogs.orthoxml"))
    assert set(ox.genes) == set(all_ids)
    tsv = {tuple(l.split("\t")[:2]) for l in open(os.path.join(out, "truth", "ortholog_pairs.tsv"))
           .read().splitlines()[1:]}
    assert tsv == hogs_to_pairs(ox)


def test_species_tree_shape(sim_basic):
    res, out = sim_basic
    t = read_newick(os.path.join(out, "species_tree.nwk"))
    assert len(t.leaf_names()) == 8
    depths = {}
    for n in t.preorder():
        depths[n.name] = 0.0 if n.parent is None else depths[n.parent.name] + n.length
    leaf_depths = [depths[x] for x in t.leaf_names()]
    assert np.allclose(leaf_depths, 1.0)       # ultrametric, default height
    assert all(n.name.startswith("N") for n in t.preorder() if n.children)


# ---------------------------------------------------------------- exact truth from gene trees
def _lca_pairs_from_gene_trees(out):
    """Independent truth: pairs whose gene-tree LCA is labelled S_* (no HGT case)."""
    pairs = set()
    for f in glob.glob(os.path.join(out, "truth", "gene_tree_*.nwk")):
        root = read_newick(f)
        below = {}
        for nd in root.postorder():
            if nd.is_leaf():
                below[id(nd)] = [nd.name]
                continue
            kids = [below.pop(id(c)) for c in nd.children]
            if nd.name.startswith("S_"):
                for a_set, b_set in itertools.combinations(kids, 2):
                    for a in a_set:
                        for b in b_set:
                            pairs.add(tuple(sorted((a, b))))
            else:
                assert nd.name.startswith(("D_", "T_"))
            below[id(nd)] = [x for k in kids for x in k]
    return pairs


def test_truth_pairs_equal_gene_tree_lca(sim_basic):
    res, out = sim_basic
    assert res.stats["n_events"]["hgt"] == 0
    tsv = {tuple(l.split("\t")[:2]) for l in open(os.path.join(out, "truth", "ortholog_pairs.tsv"))
           .read().splitlines()[1:]}
    assert tsv == _lca_pairs_from_gene_trees(out)
    # every S-labelled LCA is at a species-tree node that is the MRCA of the two species
    levels = {tuple(l.split("\t")[:2]): l.split("\t")[2]
              for l in open(os.path.join(out, "truth", "ortholog_pairs.tsv")).read().splitlines()[1:]}
    from ancestra.newick import TreeIndex
    tix = TreeIndex(read_newick(os.path.join(out, "species_tree.nwk")))
    for (a, b), lvl in levels.items():
        assert tix.mrca([a.split("_")[0], b.split("_")[0]]) == lvl


def test_hgt_wgd_truth(sim_complex):
    res, out = sim_complex
    ev = [l.split("\t") for l in open(os.path.join(out, "truth", "events.tsv")).read().splitlines()[1:]]
    types = defaultdict(int)
    for e in ev:
        types[e[2]] += 1
    assert types["hgt"] > 0 and types["wgd"] == 1
    wgd_dups = sum(1 for e in ev if e[2] == "duplication" and e[6] == "wgd")
    assert wgd_dups > 0
    lca = _lca_pairs_from_gene_trees(out)
    truth = {tuple(l.split("\t")[:2]) for l in open(os.path.join(out, "truth", "ortholog_pairs.tsv"))
             .read().splitlines()[1:]}
    # xenologs (paths crossing a transfer) are excluded from orthologs: truth is a subset
    assert truth <= lca
    hogfams = {l.split("\t")[3] for l in open(os.path.join(out, "truth", "gene_families.tsv"))
               .read().splitlines()[1:]}
    n_x = sum(1 for h in hogfams if ".x" in h)
    assert n_x > 0
    if truth != lca:
        # every excluded pair involves a gene of a transferred hog family
        fam = dict(l.split("\t")[::3] for l in open(os.path.join(out, "truth", "gene_families.tsv"))
                   .read().splitlines()[1:])
        for a, b in lca - truth:
            assert fam[a] != fam[b]


def test_sequence_divergence_tracks_tree(sim_basic):
    """3-mer sharing of true orthologs decreases with divergence and exceeds random pairs."""
    res, out = sim_basic
    seqs = {}
    for f in glob.glob(os.path.join(out, "proteomes", "*.fa")):
        seqs.update(read_fasta(f))
    km = {k: {s[i:i + 3] for i in range(len(s) - 2)} for k, s in seqs.items()}
    sim = lambda a, b: len(km[a] & km[b]) / min(len(km[a]), len(km[b]))  # noqa: E731
    from ancestra.newick import TreeIndex
    tix = TreeIndex(res.species_tree)
    by_depth = defaultdict(list)
    for a, b, lvl, _ in res.ortholog_pairs:
        by_depth[tix.depth[lvl]].append(sim(a, b))
    shallow, deep = max(by_depth), min(by_depth)
    rng = np.random.default_rng(0)
    keys = list(seqs)
    rand = [sim(*rng.choice(keys, 2, replace=False)) for _ in range(300)]
    assert np.mean(by_depth[shallow]) > np.mean(by_depth[deep]) > np.mean(rand)


# ---------------------------------------------------------------- speed
@pytest.mark.slow
def test_speed_default_size(tmp_path):
    t0 = time.time()
    res = simulate(SimParams(n_species=8, n_genes=300), seed=47)
    write_simulation(res, str(tmp_path))
    assert time.time() - t0 < 15.0
