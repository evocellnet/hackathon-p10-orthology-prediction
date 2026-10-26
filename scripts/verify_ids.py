#!/usr/bin/env python3
"""Verify the QfO prediction input before every QfO run (DESIGN.md §4).

Checks (exit status 1 if any fails; report as JSON with ``--json``):

1. ``species_78``      exactly 78 proteome files; stems ``UP<digits>_<taxid>``.
2. ``canonical_files`` every file resolves to a canonical ``UPxxx_<taxid>.fasta`` of
                       qfo_release-2022_02_with_updated_UP000000437 (not ``_additional`` /
                       ``_DNA``), and its md5 equals the release ``md5sums.txt``.
3. ``bare_accessions`` every protein ID, as the pipeline reads it (``ancestra.io.read_fasta``:
                       ``sp|ACC|NAME`` -> ``ACC``), is a bare UniProt accession (no ``-n``
                       isoform suffix, no ``|``, no version), and IDs are unique.
4. ``no_isoforms``     no input ID occurs in any release ``*_additional.fasta`` (isoforms).
5. ``coverage``        >= 99.5% of the proteins of the QfO 2022 benchmark mapping
                       (``mapping.json.gz``, evaluation data, read only here) are present
                       in the input; <= 0.5% of input IDs are unknown to it; every input
                       protein maps to one benchmark species, species map 1:1, and no two
                       input IDs map to the same protein.
6. ``tree_leaves``     (with ``--tree``) guide-tree leaves == proteome species.
7. ``gene_order``      (with ``--gene-order DIR``) one ``<species>.parquet`` per species; every
                       ``protein_id`` is an input protein of that species, listed once, and
                       (contig, order) positions are unique. Coverage per species is reported
                       (S5 expectation: 0.931 overall), not thresholded.

    python scripts/verify_ids.py [--proteomes data/inputs/qfo2022/proteomes] \
        [--tree data/inputs/qfo2022/guide_tree_ncbi.nwk] [--json out.json]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
from bisect import bisect_right
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

RELEASE = os.path.join(ROOT, "data/qfo/proteomes_2022_02_updZEBRAFISH")
MAPPING = os.path.join(ROOT, "data/qfo/reference_data/2022.1/mapping.json.gz")
# UniProt accession (https://www.uniprot.org/help/accession_numbers), fully anchored
ACC = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
STEM = re.compile(r"^UP[0-9]+_[0-9]+$")
CANON = re.compile(r"^UP[0-9]+_[0-9]+\.fasta$")
MIN_COVERAGE = 0.995


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fasta_ids(path: str):
    from ancestra.io import read_fasta
    for pid, _ in read_fasta(path):
        yield pid


def raw_header_ids(path: str):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                tok = line[1:].split()[0]
                p = tok.split("|")
                yield p[1] if len(p) >= 3 else tok


def check_gene_order(go_dir: str, ids_by_sp: dict) -> dict:
    import pandas as pd
    missing, foreign, dup_ids, dup_pos, cov = [], {}, {}, {}, {}
    n_in = n_all = 0
    for sp, ids in ids_by_sp.items():
        p = os.path.join(go_dir, sp + ".parquet")
        if not os.path.exists(p):
            missing.append(sp)
            continue
        df = pd.read_parquet(p, columns=["protein_id", "contig", "order"])
        idset = set(ids)
        f = int((~df["protein_id"].isin(idset)).sum())
        d = int(df["protein_id"].duplicated().sum())
        dp = int(df.duplicated(["contig", "order"]).sum())
        if f:
            foreign[sp] = f
        if d:
            dup_ids[sp] = d
        if dp:
            dup_pos[sp] = dp
        n = int(df["protein_id"].isin(idset).sum()) - d
        cov[sp] = round(n / max(1, len(ids)), 4)
        n_in += n
        n_all += len(ids)
    return {"ok": not (missing or foreign or dup_ids or dup_pos), "dir": go_dir,
            "missing_species": missing, "foreign_protein_ids": foreign,
            "duplicate_protein_ids": dup_ids, "duplicate_positions": dup_pos,
            "coverage_overall": round(n_in / max(1, n_all), 4),
            "min_species_coverage": min(cov.values()) if cov else None,
            "n_species_coverage_lt_0.5": sum(v < 0.5 for v in cov.values()),
            "coverage_per_species": cov}


def run(proteomes: str, tree: str | None, mapping_path: str, release: str,
        gene_order: str | None = None) -> dict:
    from ancestra.io import FASTA_EXT, species_name
    res: dict = {"proteomes": os.path.abspath(proteomes), "release": os.path.abspath(release),
                 "checks": {}}
    ck = res["checks"]
    files = sorted(f for f in os.listdir(proteomes) if f.removesuffix(".gz").endswith(FASTA_EXT))
    species = [species_name(f) for f in files]

    # 1. species
    bad_stem = [s for s in species if not STEM.match(s)]
    ck["species_78"] = {"ok": len(species) == 78 and not bad_stem, "n_species": len(species),
                        "bad_stems": bad_stem}

    # 2. canonical release files, md5
    rel_md5 = {}
    with open(os.path.join(release, "md5sums.txt")) as fh:
        for line in fh:
            h, p = line.split(None, 1)
            rel_md5[os.path.basename(p.strip())] = h
    rel_real = os.path.realpath(release)
    probs = []
    for f in files:
        rp = os.path.realpath(os.path.join(proteomes, f))
        b = os.path.basename(rp)
        if not rp.startswith(rel_real + os.sep):
            probs.append((f, "outside release dir"))
        elif not CANON.match(b):
            probs.append((f, f"not a canonical FASTA: {b}"))
        elif b != f:
            probs.append((f, f"name differs from target {b}"))
        elif rel_md5.get(b) != md5(rp):
            probs.append((f, "md5 differs from release md5sums.txt"))
    ck["canonical_files"] = {"ok": not probs, "problems": probs}

    # 3. bare accessions, unique
    ids_by_sp: dict[str, list[str]] = {}
    bad_ids, n_raw_mismatch = [], 0
    for f, sp in zip(files, species):
        ids = list(fasta_ids(os.path.join(proteomes, f)))
        raw = list(raw_header_ids(os.path.join(proteomes, f)))
        n_raw_mismatch += sum(a != b for a, b in zip(ids, raw)) + abs(len(ids) - len(raw))
        bad_ids += [(sp, i) for i in ids if not ACC.match(i)]
        ids_by_sp[sp] = ids
    all_ids = [i for v in ids_by_sp.values() for i in v]
    dups = [i for i, c in Counter(all_ids).items() if c > 1]
    ck["bare_accessions"] = {"ok": not bad_ids and not dups and n_raw_mismatch == 0,
                             "n_proteins": len(all_ids), "n_not_accession": len(bad_ids),
                             "examples_not_accession": bad_ids[:10], "n_duplicate_ids": len(dups),
                             "examples_duplicates": dups[:10],
                             "n_pipeline_vs_header_mismatch": n_raw_mismatch,
                             "per_species": {s: len(v) for s, v in ids_by_sp.items()}}

    # 4. isoforms from _additional files
    id_set = set(all_ids)
    iso_hits, n_iso = [], 0
    for dom in ("Archaea", "Bacteria", "Eukaryota"):
        d = os.path.join(release, dom)
        for f in sorted(os.listdir(d)):
            if f.endswith("_additional.fasta"):
                for i in raw_header_ids(os.path.join(d, f)):
                    n_iso += 1
                    if i in id_set:
                        iso_hits.append((f, i))
    with_dash = [i for i in all_ids if "-" in i]
    ck["no_isoforms"] = {"ok": not iso_hits and not with_dash, "n_isoform_entries_scanned": n_iso,
                         "n_input_ids_in_additional": len(iso_hits), "examples": iso_hits[:10],
                         "n_ids_with_isoform_suffix": len(with_dash)}

    # 5. coverage vs the benchmark mapping
    with gzip.open(mapping_path, "rt") as fh:
        mp = json.load(fh)
    mapping, goff, msp = mp["mapping"], mp["Goff"], mp["species"]
    n_ref = goff[-1]
    covered = set()
    not_in_mapping = []
    sp2codes: dict[str, Counter] = defaultdict(Counter)
    nr2ids: dict[int, list] = defaultdict(list)
    for sp, ids in ids_by_sp.items():
        for i in ids:
            nr = mapping.get(i)
            if nr is None:
                not_in_mapping.append((sp, i))
                continue
            covered.add(nr)
            nr2ids[nr].append(i)
            sp2codes[sp][msp[bisect_right(goff, nr - 1) - 1]] += 1
    multi = {nr: v for nr, v in nr2ids.items() if len(v) > 1}
    sp_bad = {sp: dict(c) for sp, c in sp2codes.items() if len(c) != 1}
    code_of = {sp: next(iter(c)) for sp, c in sp2codes.items() if len(c) == 1}
    code_dup = [c for c, n in Counter(code_of.values()).items() if n > 1]
    per_code_ref = {msp[k]: goff[k + 1] - goff[k] for k in range(len(msp))}
    per_sp_cov = {sp: round(sum(1 for i in ids_by_sp[sp] if i in mapping) / per_code_ref[c], 5)
                  for sp, c in code_of.items()}
    cov = len(covered) / n_ref
    extra = len(not_in_mapping) / max(1, len(all_ids))
    ck["coverage"] = {"ok": cov >= MIN_COVERAGE and extra <= 1 - MIN_COVERAGE and not sp_bad
                      and not code_dup and not multi and len(code_of) == len(msp),
                      "frac_input_not_in_mapping": round(extra, 6),
                      "n_benchmark_proteins": n_ref, "n_covered": len(covered),
                      "coverage": round(cov, 6), "threshold": MIN_COVERAGE,
                      "n_input_not_in_mapping": len(not_in_mapping),
                      "examples_not_in_mapping": not_in_mapping[:10],
                      "n_benchmark_species": len(msp), "n_species_mapped_1to1": len(code_of),
                      "species_mapping_to_several_codes": sp_bad, "codes_hit_by_several_species": code_dup,
                      "n_benchmark_proteins_hit_twice": len(multi),
                      "min_species_coverage": min(per_sp_cov.values()) if per_sp_cov else None,
                      "species_code": code_of}

    # 6. tree leaves
    if tree:
        from ancestra.newick import parse_newick
        leaves = set(parse_newick(open(tree).read()).leaf_names())
        ck["tree_leaves"] = {"ok": leaves == set(species), "tree": tree,
                             "only_in_tree": sorted(leaves - set(species)),
                             "only_in_proteomes": sorted(set(species) - leaves)}
    if gene_order:
        ck["gene_order"] = check_gene_order(gene_order, ids_by_sp)
    res["ok"] = all(c["ok"] for c in ck.values())
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--proteomes", default=os.path.join(ROOT, "data/inputs/qfo2022/proteomes"))
    ap.add_argument("--tree", default=None)
    ap.add_argument("--mapping", default=MAPPING)
    ap.add_argument("--release", default=RELEASE)
    ap.add_argument("--gene-order", help="S5 gene-order table dir (<species>.parquet)")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    res = run(a.proteomes, a.tree, a.mapping, a.release, a.gene_order)
    for name, c in res["checks"].items():
        brief = {k: v for k, v in c.items() if k not in ("per_species", "species_code", "ok", "coverage_per_species")
                 and not (isinstance(v, (list, dict)) and not v)}
        print(f"[verify_ids] {'PASS' if c['ok'] else 'FAIL'} {name}: "
              f"{json.dumps(brief, default=str)[:400]}")
    print(f"[verify_ids] overall: {'PASS' if res['ok'] else 'FAIL'}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(res, fh, indent=1, default=str)
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
