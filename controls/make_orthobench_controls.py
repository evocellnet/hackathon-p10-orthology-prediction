#!/usr/bin/env python3
"""Generate positive/negative control orthogroup files for OrthoBench.

positive: every RefOG as one orthogroup (all other genes as singletons) -> expect F=100%
negative: random partition of all genes with the same group-size distribution as the
          positive control's non-singleton groups, and random partition with groups of size ~10
          -> expect near-zero F-score.
Usage: make_orthobench_controls.py <OrthoBench BENCHMARKS dir> <outdir> [--seed 42]
"""
import argparse, glob, os, random

p = argparse.ArgumentParser()
p.add_argument("bench_dir"); p.add_argument("outdir"); p.add_argument("--seed", type=int, default=42)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
genes = []
for fn in sorted(glob.glob(os.path.join(a.bench_dir, "Input", "*.fa"))):
    with open(fn) as f:
        genes += [l[1:].strip() for l in f if l.startswith(">")]
refogs = []
for i in range(1, 71):
    with open(os.path.join(a.bench_dir, "RefOGs", "RefOG%03d.txt" % i)) as f:
        refogs.append([g.strip() for g in f if g.strip()])
# FBpp0309618 occurs in two RefOGs upstream; a partition can only place it once, so keep
# the first occurrence (the positive control is therefore ~99.9%, not exactly 100%).
seen = set(); dedup = []
for r in refogs:
    dedup.append([g for g in r if g not in seen]); seen.update(r)
refogs_part = [r for r in dedup if r]
inref = seen

def write(name, groups):
    with open(os.path.join(a.outdir, name), "w") as f:
        f.write("Orthogroup\tGenes\n")
        for i, g in enumerate(groups):
            f.write("OG%07d: %s\n" % (i, ", ".join(g)))

write("positive_refogs.txt", refogs_part + [[g] for g in genes if g not in inref])
rng = random.Random(a.seed)
sh = genes[:]; rng.shuffle(sh)
# random partition into groups of size 10
write("negative_random10.txt", [sh[i:i + 10] for i in range(0, len(sh), 10)])
# random partition with same size distribution as positive control
sizes = [len(r) for r in refogs] + [1] * (len(genes) - len(inref))
rng.shuffle(sh); out, k = [], 0
for s in sizes:
    out.append(sh[k:k + s]); k += s
write("negative_shuffled_refog_sizes.txt", out)
# all genes in a single group (recall=100, precision tiny)
write("negative_one_group.txt", [genes])
print("wrote controls to", a.outdir, "n_genes", len(genes))
# intermediate control: each RefOG split into two halves (precision 100%, recall ~50%)
halves = []
for r in refogs_part:
    k = len(r) // 2
    halves += [r[:k], r[k:]] if k else [r]
write("intermediate_split_halves.txt", [h for h in halves if h] + [[g] for g in genes if g not in inref])
# intermediate control: each RefOG merged with its successor (recall 100%, precision ~50%)
merged = [refogs_part[i] + refogs_part[i + 1] for i in range(0, len(refogs_part) - 1, 2)]
if len(refogs_part) % 2: merged.append(refogs_part[-1])
write("intermediate_merged_pairs.txt", merged + [[g] for g in genes if g not in inref])
