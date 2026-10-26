#!/usr/bin/env python3
"""Summarise QfO assessment stubs (OpenEBench 'assessment' JSON written by the upstream
benchmark code) into a TSV: benchmark, challenge_id, metric, value, stderr.
By default only aggregate challenges are printed (e.g. 'SwissTrees', 'GO', 'G_STD2_Luca');
use --all to include per-family sub-challenges (e.g. 'SwissTrees-APP')."""
import argparse, glob, json, os
p = argparse.ArgumentParser(); p.add_argument("assess_dir"); p.add_argument("--all", action="store_true")
a = p.parse_args()
print("benchmark\tchallenge_id\tmetric\tvalue\tstderr")
for fn in sorted(glob.glob(os.path.join(a.assess_dir, "*.json"))):
    b = os.path.basename(fn)[:-5]
    try:
        data = json.load(open(fn))
    except Exception as e:
        print(f"{b}\t\tERROR\t{e}\t"); continue
    for x in (data if isinstance(data, list) else [data]):
        if x.get("type") != "assessment":
            continue
        cid = x.get("challenge_id", "")
        if not a.all and "-" in cid.replace(b, "", 1) and cid != b:
            continue
        m = x.get("metrics", {})
        print(f"{b}\t{cid}\t{m.get('metric_id')}\t{m.get('value')}\t{m.get('stderr', '')}")
