#!/usr/bin/env python3
"""Create the pre-registered OrthoBench dev/test split (DESIGN.md §4), v2: homology-aware.

v1 (random 35/35 over RefOG ids) leaked through paralogous RefOGs (review r2_stats_ai_failure #1).
v2 procedure (fixed, do not change):
 1. Homology components. All 1,945 RefOG sequences
    (Supporting_Data/Data_for_RefOGs/refog_sequences/*.fa) are searched all-vs-all with
    `diamond blastp --sensitive -e 1e-5 -k 0`. Two RefOGs are linked if any hit with
    e-value < 1e-5 joins a gene of one to a gene of the other, or if they share a gene
    (FBpp0309618 is listed in both RefOG021 and RefOG068). Units = connected components.
 2. Size stratification. Units are sorted by (mean RefOG size of the unit, smallest RefOG id)
    and cut into consecutive strata of 2 units (the last stratum may hold 1).
 3. Assignment with random.Random(11). Strata are visited in a random order (rng.shuffle).
    Inside a stratum, the two units are put in random order (rng.shuffle); if dev and test
    currently hold the same number of RefOGs, the first unit goes to the subset chosen by
    rng.random() < 0.5 (dev) else test; otherwise the unit with more RefOGs goes to the subset
    that currently has fewer RefOGs (keeps ~35/35 even with the 7-RefOG component). The
    other unit goes to the other subset. A single-unit stratum goes to the smaller subset.
Writes splits/orthobench_split.json, splits/orthobench_refog_components.tsv and their
lines in splits/SHA256SUMS.
"""
import argparse, hashlib, json, os, random, re, subprocess, sys, tempfile
from collections import defaultdict

SEED = 11
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OB = os.path.join(ROOT, "data/orthobench/Open_Orthobench")


def refog_ids(refog_dir):
    ids = sorted(f[:-4] for f in os.listdir(refog_dir) if re.fullmatch(r"RefOG\d{3}\.txt", f))
    if len(ids) != 70:
        sys.exit(f"expected 70 RefOGs, found {len(ids)}")
    return ids


def homology_components(seq_dir, ids, threads=24, evalue=1e-5, hits_out=None):
    gene2ref = defaultdict(set)
    seen = set()
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "refogs.fa")
        with open(fa, "w") as fo:
            for rid in ids:
                for line in open(os.path.join(seq_dir, rid + ".fa")):
                    if line.startswith(">"):
                        g = line[1:].split()[0]
                        gene2ref[g].add(rid)
                        w = g not in seen
                        seen.add(g)
                        line = f">{g}\n"
                    if w:
                        fo.write(line)
        db = os.path.join(td, "refogs")
        hits = hits_out or os.path.join(td, "hits.tsv")
        subprocess.run(["diamond", "makedb", "--in", fa, "-d", db, "--quiet"], check=True)
        subprocess.run(["diamond", "blastp", "-q", fa, "-d", db, "-o", hits, "--sensitive", "-k", "0",
                        "-e", str(evalue), "--threads", str(threads), "--quiet",
                        "--outfmt", "6", "qseqid", "sseqid", "evalue", "bitscore"], check=True)
        parent = {r: r for r in ids}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        links = defaultdict(int)

        def link(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)

        for g, rs in gene2ref.items():
            rs = sorted(rs)
            for b in rs[1:]:
                links[(rs[0], b)] += 1
                link(rs[0], b)
        for line in open(hits):
            q, s, e, _ = line.rstrip("\n").split("\t")
            if float(e) >= evalue:
                continue
            for a in gene2ref[q]:
                for b in gene2ref[s]:
                    if a != b:
                        links[tuple(sorted((a, b)))] += 1
                        link(a, b)
    comps = defaultdict(list)
    for r in ids:
        comps[find(r)].append(r)
    return sorted((sorted(c) for c in comps.values()), key=lambda c: c[0]), dict(links)


def assign(units, sizes, seed=SEED):
    rng = random.Random(seed)
    key = lambda u: (sum(sizes[r] for r in u) / len(u), u[0])
    order = sorted(units, key=key)
    strata = [order[i:i + 2] for i in range(0, len(order), 2)]
    rng.shuffle(strata)
    dev, test = [], []
    for st in strata:
        st = list(st)
        rng.shuffle(st)
        if len(st) == 1:
            (dev if len(dev) <= len(test) else test).extend(st[0])
            continue
        u1, u2 = st
        if len(dev) == len(test):
            first_dev = rng.random() < 0.5
        else:
            big, small = (u1, u2) if len(u1) >= len(u2) else (u2, u1)
            u1, u2 = big, small
            first_dev = len(dev) < len(test)
        if first_dev:
            dev.extend(u1); test.extend(u2)
        else:
            test.extend(u1); dev.extend(u2)
    return sorted(dev), sorted(test)


def _update_sums(sums, path):
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    name = os.path.basename(path)
    lines = [l for l in (open(sums).read().splitlines() if os.path.exists(sums) else []) if not l.endswith("  " + name)]
    lines.append(f"{digest}  {name}")
    open(sums, "w").write("\n".join(sorted(lines, key=lambda l: l.split()[1])) + "\n")
    return digest


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--refogs", default=os.path.join(OB, "BENCHMARKS/RefOGs"))
    ap.add_argument("--seqs", default=os.path.join(OB, "Supporting_Data/Data_for_RefOGs/refog_sequences"))
    ap.add_argument("--outdir", default=os.path.join(ROOT, "splits"))
    ap.add_argument("--threads", type=int, default=24)
    a = ap.parse_args()
    ids = refog_ids(a.refogs)
    sizes = {r: sum(1 for l in open(os.path.join(a.refogs, r + ".txt")) if l.strip()) for r in ids}
    units, links = homology_components(a.seqs, ids, a.threads)
    dev, test = assign(units, sizes)
    os.makedirs(a.outdir, exist_ok=True)
    comp_path = os.path.join(a.outdir, "orthobench_refog_components.tsv")
    with open(comp_path, "w") as fh:
        fh.write("#component\tn_refogs\trefogs\tsizes\tsubset\n")
        for i, u in enumerate(units, 1):
            sub = "dev" if u[0] in dev else "test"
            assert all((r in dev) == (sub == "dev") for r in u)
            fh.write(f"C{i:02d}\t{len(u)}\t{','.join(u)}\t{','.join(str(sizes[r]) for r in u)}\t{sub}\n")
    split = {
        "version": 2,
        "description": "OrthoBench RefOG dev/test split (DESIGN.md §4), homology-aware v2 (whole homology components "
                       "assigned together, size-stratified). Tune on dev only; test is touched once per frozen config.",
        "rule": "see scripts/make_orthobench_split.py docstring (diamond --sensitive e<1e-5 RefOG-RefOG links -> "
                "components; strata of 2 by mean RefOG size; random.Random(11); balance #RefOGs)",
        "seed": SEED,
        "source": "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs + Supporting_Data/Data_for_RefOGs/refog_sequences",
        "n_components": len(units),
        "multi_refog_components": [u for u in units if len(u) > 1],
        "refog_links_n_hits": {f"{a_}-{b_}": n for (a_, b_), n in sorted(links.items())},
        "dev": dev,
        "test": test,
        "dev_genes": sum(sizes[r] for r in dev),
        "test_genes": sum(sizes[r] for r in test),
    }
    out = os.path.join(a.outdir, "orthobench_split.json")
    with open(out, "w") as fh:
        json.dump(split, fh, indent=1)
        fh.write("\n")
    sums = os.path.join(a.outdir, "SHA256SUMS")
    print(_update_sums(sums, out), out)
    print(_update_sums(sums, comp_path), comp_path)
    print(f"components={len(units)} multi={split['multi_refog_components']}")
    print(f"dev={len(dev)} RefOGs/{split['dev_genes']} genes; test={len(test)} RefOGs/{split['test_genes']} genes")
