#!/usr/bin/env python3
"""Prepare the QfO 2022 prediction inputs (DESIGN.md §4) -- idempotent.

Creates (all under data/, gitignored):
  data/inputs/qfo2022/proteomes/UPxxxxxxxxx_<taxid>.fasta
        symlinks to the 78 *canonical* FASTAs of
        data/qfo/proteomes_2022_02_updZEBRAFISH/{Archaea,Bacteria,Eukaryota}/
        (QfO 2022_02_with_updated_UP000000437). Only files matching
        ``UP[0-9]+_[0-9]+.fasta`` are linked: no ``_additional`` (isoforms), no
        ``_DNA``. The link name (= file stem) is the species name used by
        ancestra.io and is identical to the leaf names of the S1 QfO trees.
  data/inputs/qfo2022/guide_tree_ncbi.nwk
        copy of results/S1_speciestree/qfo/ncbi.nwk (NCBI taxdump tree, S1 default).
  data/inputs/qfo2022/subsets/sub<k>_s<seed>/
        symlinks for a seeded, domain-stratified species subset (cost extrapolation).
  data/inputs/qfo2022/inputs.json
        md5 of every source FASTA (checked against the release md5sums.txt), tree md5,
        leaf-name check, subset membership.

    python scripts/prepare_qfo_inputs.py [--subset 10 --seed 11]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import re
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

SRC = os.path.join(ROOT, "data/qfo/proteomes_2022_02_updZEBRAFISH")
OUT = os.path.join(ROOT, "data/inputs/qfo2022")
TREE_SRC = os.path.join(ROOT, "results/S1_speciestree/qfo/ncbi.nwk")
CANON = re.compile(r"^UP[0-9]+_[0-9]+\.fasta$")
DOMAINS = ("Archaea", "Bacteria", "Eukaryota")


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def canonical_files() -> dict[str, tuple[str, str]]:
    """species -> (domain, absolute path) of the canonical FASTA."""
    out = {}
    for dom in DOMAINS:
        for f in sorted(os.listdir(os.path.join(SRC, dom))):
            if CANON.match(f):
                out[f[:-len(".fasta")]] = (dom, os.path.join(SRC, dom, f))
    return out


def release_md5() -> dict[str, str]:
    m = {}
    with open(os.path.join(SRC, "md5sums.txt")) as fh:
        for line in fh:
            h, p = line.split(None, 1)
            m[os.path.basename(p.strip())] = h
    return m


def link_dir(d: str, files: dict[str, str]) -> None:
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):        # remove stale entries
        if f.removesuffix(".fasta") not in files:
            os.remove(os.path.join(d, f))
    for sp, path in files.items():
        dst = os.path.join(d, sp + ".fasta")
        if os.path.islink(dst) and os.readlink(dst) == path:
            continue
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(path, dst)


def stratified_subset(cf: dict, k: int, seed: int) -> list[str]:
    """k species, domain-stratified (proportional, >=1 per domain), seeded."""
    by = {d: sorted(s for s, (dd, _) in cf.items() if dd == d) for d in DOMAINS}
    n = len(cf)
    quota = {d: max(1, round(k * len(v) / n)) for d, v in by.items()}
    while sum(quota.values()) > k:
        quota[max(quota, key=quota.get)] -= 1
    while sum(quota.values()) < k:
        quota["Eukaryota"] += 1
    rng = random.Random(seed)
    out = []
    for d in DOMAINS:
        out += rng.sample(by[d], quota[d])
    return sorted(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--subset", type=int, default=10)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args(argv)

    cf = canonical_files()
    if len(cf) != 78:
        raise SystemExit(f"expected 78 canonical proteomes, found {len(cf)}")
    rel = release_md5()
    src_md5 = {}
    for sp, (_, p) in cf.items():
        h = md5(p)
        if rel.get(os.path.basename(p)) != h:
            raise SystemExit(f"md5 mismatch vs release md5sums.txt: {p}")
        src_md5[sp] = h
    link_dir(os.path.join(OUT, "proteomes"), {sp: p for sp, (_, p) in cf.items()})

    # guide tree: copy (frozen input) + leaf-name check through the pipeline loader
    tree_dst = os.path.join(OUT, "guide_tree_ncbi.nwk")
    shutil.copyfile(TREE_SRC, tree_dst + ".tmp")
    os.replace(tree_dst + ".tmp", tree_dst)
    from ancestra.io import load_guide_tree
    from ancestra.newick import parse_newick
    raw = parse_newick(open(tree_dst).read())
    leaves = set(raw.leaf_names())
    if leaves != set(cf):
        raise SystemExit(f"tree leaves != proteome stems: only-tree {sorted(leaves - set(cf))[:5]}, "
                         f"only-proteomes {sorted(set(cf) - leaves)[:5]}")
    t = load_guide_tree(tree_dst, sorted(cf))
    n_int = sum(1 for nd in t.preorder() if nd.children)

    sub = stratified_subset(cf, a.subset, a.seed)
    sub_name = f"sub{a.subset}_s{a.seed}"
    link_dir(os.path.join(OUT, "subsets", sub_name), {sp: cf[sp][1] for sp in sub})

    info = {"source_release": "qfo_release-2022_02_with_updated_UP000000437",
            "source_dir": os.path.relpath(SRC, ROOT),
            "n_species": len(cf), "proteome_md5": src_md5,
            "domain": {sp: d for sp, (d, _) in cf.items()},
            "guide_tree": {"source": os.path.relpath(TREE_SRC, ROOT), "md5": md5(tree_dst),
                           "n_leaves": len(leaves), "n_internal_after_load": n_int,
                           "leaf_names_match_proteome_stems": True},
            "subsets": {sub_name: {"seed": a.seed, "k": a.subset, "species": sub,
                                   "domains": {d: sum(cf[s][0] == d for s in sub) for d in DOMAINS}}}}
    old = os.path.join(OUT, "inputs.json")
    if os.path.exists(old):
        prev = json.load(open(old))
        info["subsets"] = {**prev.get("subsets", {}), **info["subsets"]}
    with open(old + ".tmp", "w") as fh:
        json.dump(info, fh, indent=1)
    os.replace(old + ".tmp", old)
    print(f"[prepare] 78 proteomes linked, md5 OK; tree leaves match ({n_int} internal nodes); "
          f"subset {sub_name}: {sub}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
