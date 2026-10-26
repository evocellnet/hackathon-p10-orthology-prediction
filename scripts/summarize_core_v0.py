"""Markdown tables from results/core_v0/toy_grid.jsonl (mean +- sd over seeds).

    python scripts/summarize_core_v0.py [results/core_v0/toy_grid.jsonl] > results/core_v0/toy_tables.md
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def ms(xs, f="{:.3f}"):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "-"
    if len(xs) == 1:
        return f.format(xs[0])
    return f.format(np.mean(xs)) + "±" + f.format(np.std(xs, ddof=1))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results/core_v0/toy_grid.jsonl")
    recs = [json.loads(l) for l in open(path)]
    # de-duplicate (keep last record per key)
    d = {}
    for r in recs:
        d[(r["sim"], r["families"], r["backend"], r["evalue"], r["arm"])] = r
    recs = list(d.values())
    settings = sorted({(r["setting"], r["n_species"]) for r in recs}, key=lambda x: (x[1], x[0]))
    arms = []
    for r in recs:
        if r["arm"] not in arms:
            arms.append(r["arm"])
    print("## Merge arms (DIAMOND blocks, e=10); mean±sd over seeds 11/23/47\n")
    for fams in ("truth", "s0"):
        print(f"### families = {fams}\n")
        print("| setting | arm | P | R | F1 | FP paralog | root exact | level exact | level Jacc | deep exact | fissions | fusions |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for st in settings:
            for arm in arms:
                rs = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == fams
                      and r["backend"] == "diamond" and r["evalue"] == 10.0 and r["arm"] == arm]
                if not rs:
                    continue
                print(f"| {st[0]}{st[1]} | {arm} | {ms([r['P'] for r in rs])} | {ms([r['R'] for r in rs])} | "
                      f"{ms([r['F1'] for r in rs])} | {ms([r['fp_paralogs'] for r in rs], '{:.0f}')} | "
                      f"{ms([r['root_exact'] for r in rs])} | {ms([r['lvl_exact'] for r in rs])} | "
                      f"{ms([r['lvl_jacc'] for r in rs])} | {ms([r['deep_exact'] for r in rs])} | "
                      f"{ms([r['lvl_fissions'] for r in rs], '{:.0f}')} | {ms([r['lvl_fusions'] for r in rs], '{:.0f}')} |")
        print()
    print("### Overall mean over the 12 truth-family conditions (selection rule)\n")
    print("| arm | P | R | F1 |")
    print("|---|---|---|---|")
    rows = []
    for arm in arms:
        rs = [r for r in recs if r["families"] == "truth" and r["backend"] == "diamond"
              and r["evalue"] == 10.0 and r["arm"] == arm]
        if rs:
            rows.append((np.mean([r["F1"] for r in rs]), np.mean([r["P"] for r in rs]),
                         np.mean([r["R"] for r in rs]), arm, len(rs)))
    for f1, p, rr, arm, n in sorted(rows, reverse=True):
        print(f"| {arm} (n={n}) | {p:.4f} | {rr:.4f} | {f1:.4f} |")
    print("\n## S2 arms (truth families)\n")
    print("| setting | merge arm | backend | e-value | P | R | F1 | scored pairs | score wall s |")
    print("|---|---|---|---|---|---|---|---|---|")
    for st in settings:
        for arm in ("max+repair", "avg+repair"):
            for be, ev in (("diamond", 10.0), ("diamond", 0.001), ("mmseqs", 10.0)):
                rs = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == "truth"
                      and r["backend"] == be and r["evalue"] == ev and r["arm"] == arm]
                if rs:
                    print(f"| {st[0]}{st[1]} | {arm} | {be} | {ev:g} | {ms([r['P'] for r in rs])} | "
                          f"{ms([r['R'] for r in rs])} | {ms([r['F1'] for r in rs])} | "
                          f"{ms([r['n_scored_pairs'] for r in rs], '{:.0f}')} | {ms([r['score_s'] for r in rs], '{:.1f}')} |")
    print("\n## S0 default clustering vs truth families\n")
    print("| setting | S0 families | truth families exactly recovered | fissions | fusions | S0 wall s |")
    print("|---|---|---|---|---|---|")
    for st in settings:
        rs = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == "s0"
              and r["arm"] == arms[0]]
        if rs:
            print(f"| {st[0]}{st[1]} | {ms([r['s0_n_families'] for r in rs], '{:.0f}')} | "
                  f"{ms([r['s0_exact'] for r in rs])} | {ms([r['s0_fissions'] for r in rs], '{:.0f}')} | "
                  f"{ms([r['s0_fusions'] for r in rs], '{:.0f}')} | {ms([r['s0_wall_s'] for r in rs], '{:.0f}')} |")
    print("\n## Merge wall time (s, all nodes, one process)\n")
    print("| setting | genes | merge s (max+repair) | merge s (avg+repair/quartet) |")
    print("|---|---|---|---|")
    for st in settings:
        a = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == "truth"
             and r["backend"] == "diamond" and r["evalue"] == 10.0]
        print(f"| {st[0]}{st[1]} | {ms([r['n_genes'] for r in a if r['arm'] == 'max+repair'], '{:.0f}')} | "
              f"{ms([r['merge_s'] for r in a if r['arm'] == 'max+repair'], '{:.2f}')} | "
              f"{ms([r['merge_s'] for r in a if r['arm'] == 'avg+repair/quartet'], '{:.2f}')} |")


if __name__ == "__main__":
    main()
