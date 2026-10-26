#!/usr/bin/env python3
"""Create the pre-registered OrthoBench dev/test split (DESIGN.md §4).

Rule (fixed, do not change): take the sorted RefOG ids (RefOG001..RefOG070, from the
files in Open_Orthobench/BENCHMARKS/RefOGs), shuffle a copy with random.Random(11),
first 35 -> dev, remaining 35 -> test (each list re-sorted). Writes
splits/orthobench_split.json and (re)writes its line in splits/SHA256SUMS.
"""
import argparse, hashlib, json, os, random, re, sys

SEED = 11

def make_split(refog_dir):
    ids = sorted(f[:-4] for f in os.listdir(refog_dir) if re.fullmatch(r"RefOG\d{3}\.txt", f))
    if len(ids) != 70:
        sys.exit(f"expected 70 RefOGs, found {len(ids)}")
    shuffled = list(ids)
    random.Random(SEED).shuffle(shuffled)
    return {
        "description": "OrthoBench RefOG dev/test split (DESIGN.md §4). Tune on dev only; test is touched once per frozen config.",
        "rule": "ids=sorted(RefOG ids); s=list(ids); random.Random(11).shuffle(s); dev=sorted(s[:35]); test=sorted(s[35:])",
        "seed": SEED,
        "source": "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs",
        "dev": sorted(shuffled[:35]),
        "test": sorted(shuffled[35:]),
    }

if __name__ == "__main__":
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser()
    ap.add_argument("--refogs", default=os.path.join(root, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"))
    ap.add_argument("--out", default=os.path.join(root, "splits/orthobench_split.json"))
    a = ap.parse_args()
    split = make_split(a.refogs)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(split, fh, indent=1)
        fh.write("\n")
    digest = hashlib.sha256(open(a.out, "rb").read()).hexdigest()
    sums = os.path.join(os.path.dirname(a.out), "SHA256SUMS")
    name = os.path.basename(a.out)
    lines = [l for l in (open(sums).read().splitlines() if os.path.exists(sums) else []) if not l.endswith("  " + name)]
    lines.append(f"{digest}  {name}")
    open(sums, "w").write("\n".join(lines) + "\n")
    print(digest, a.out)
