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
    print("## Headline: default configuration (arm `default`); mean±sd over seeds 11/23/47\n")
    print("Families: `truth` = simulator root families; `s0b` = current S0 default `2sA_dc30_cc_sh30_RA`; "
          "`s0` = earlier S0 default (dc30 reps-vs-reps CC, min-cov 50%).\n")
    print("Scoring `hybrid+rescue` = final default (DIAMOND blocks + MMseqs2 rescue of missed in-family "
          "pairs + cross-family S0 rescue); `hybrid` = without the family rescue.\n")
    print("| setting | families | scoring | P | R | F1 | FP paralog | root exact | level exact | level Jacc | deep exact (d<=1) | shallow exact (d>=2) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for st in settings:
        for fams, be in (("truth", "hybrid+rescue"), ("truth", "hybrid"), ("s0b", "hybrid+rescue"), ("s0b", "hybrid"),
                         ("s0", "hybrid+rescue"), ("s0", "hybrid")):
            rs = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == fams
                  and r["backend"] == be and r["arm"] == "default"]
            if rs:
                print(f"| {st[0]}{st[1]} | {fams} | {be} | {ms([r['P'] for r in rs])} | {ms([r['R'] for r in rs])} | "
                      f"{ms([r['F1'] for r in rs])} | {ms([r['fp_paralogs'] for r in rs], '{:.0f}')} | "
                      f"{ms([r['root_exact'] for r in rs])} | {ms([r['lvl_exact'] for r in rs])} | "
                      f"{ms([r['lvl_jacc'] for r in rs])} | {ms([r['deep_exact'] for r in rs])} | "
                      f"{ms([r['shallow_exact'] for r in rs])} |")
    print()
    # per-depth HOG metrics from the saved evaluation JSONs
    print("### Per-level HOG metrics by depth (default, hybrid+rescue, truth families; pooled over seeds)\n")
    print("| setting | depth | n levels | truth HOGs | exact | mean best Jaccard | pair F1 at level | fissions/level | fusions/level |")
    print("|---|---|---|---|---|---|---|---|---|")
    for st in settings:
        rows = defaultdict(list)
        for r in recs:
            if (r["setting"], r["n_species"]) == st and r["families"] == "truth" and r["backend"] == "hybrid+rescue" \
                    and r["arm"] == "default":
                ej = os.path.join(ROOT, "work/core_v0/runs", r["sim"], "truth", "eval_hybrid+rescue_e10_default.json")
                if not os.path.exists(ej):
                    continue
                ev = json.load(open(ej))
                for lvl, m in ev["hogs"].items():
                    if lvl == "families":
                        continue
                    dbin = str(m["depth"]) if m["depth"] <= 2 else ("3-5" if m["depth"] <= 5 else "6+")
                    rows[dbin].append(m)
        for dbin in sorted(rows, key=lambda x: int(x.split("-")[0].rstrip("+"))):
            ms_ = rows[dbin]
            print(f"| {st[0]}{st[1]} | {dbin} | {len(ms_)} | {sum(m['n_truth_hogs'] for m in ms_)} | "
                  f"{np.mean([m['exact'] for m in ms_]):.3f} | {np.mean([m['mean_best_jaccard'] for m in ms_]):.3f} | "
                  f"{np.mean([m['pair_f1'] for m in ms_]):.3f} | {np.mean([m['fissions'] for m in ms_]):.1f} | "
                  f"{np.mean([m['fusions'] for m in ms_]):.1f} |")
    print()
    for be in ("diamond", "hybrid", "hybrid+rescue"):
      print(f"## Merge arms ({be} scoring, e=10); mean±sd over seeds 11/23/47\n")
      for fams in ("truth", "s0"):
        print(f"### families = {fams}\n")
        print("| setting | arm | P | R | F1 | FP paralog | root exact | level exact | level Jacc | deep exact | fissions | fusions |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for st in settings:
            for arm in arms:
                rs = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == fams
                      and r["backend"] == be and r["evalue"] == 10.0 and r["arm"] == arm]
                if not rs:
                    continue
                print(f"| {st[0]}{st[1]} | {arm} | {ms([r['P'] for r in rs])} | {ms([r['R'] for r in rs])} | "
                      f"{ms([r['F1'] for r in rs])} | {ms([r['fp_paralogs'] for r in rs], '{:.0f}')} | "
                      f"{ms([r['root_exact'] for r in rs])} | {ms([r['lvl_exact'] for r in rs])} | "
                      f"{ms([r['lvl_jacc'] for r in rs])} | {ms([r['deep_exact'] for r in rs])} | "
                      f"{ms([r['lvl_fissions'] for r in rs], '{:.0f}')} | {ms([r['lvl_fusions'] for r in rs], '{:.0f}')} |")
        print()
    print("### Overall mean over the 12 truth-family conditions, DIAMOND e=10 (selection rule)\n")
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
            for be, ev in (("diamond", 10.0), ("diamond", 0.001), ("mmseqs", 10.0), ("hybrid", 10.0)):
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
      for fsrc in ("s0", "s0b"):
        rs = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == fsrc
              and r.get("s0_exact") is not None]
        best = {}
        for r in rs:   # keep the record that actually ran S0 (cached reruns report ~0 s)
            if r["sim"] not in best or (r.get("s0_wall_s") or 0) > (best[r["sim"]].get("s0_wall_s") or 0):
                best[r["sim"]] = r
        rs = list(best.values())
        if rs:
            print(f"| {st[0]}{st[1]} {fsrc} | {ms([r['s0_n_families'] for r in rs], '{:.0f}')} | "
                  f"{ms([r['s0_exact'] for r in rs])} | {ms([r['s0_fissions'] for r in rs], '{:.0f}')} | "
                  f"{ms([r['s0_fusions'] for r in rs], '{:.0f}')} | {ms([r['s0_wall_s'] for r in rs], '{:.0f}')} |")
    print("\n## Merge wall time (s, all nodes, one process)\n")
    print("| setting | genes | scoring s: diamond | mmseqs | hybrid | merge s: default | avg+repair/quartet |")
    print("|---|---|---|---|---|---|---|")
    for st in settings:
        a = [r for r in recs if (r["setting"], r["n_species"]) == st and r["families"] == "truth"]
        g = lambda be, arm, k: ms([r[k] for r in a if r["backend"] == be and r["evalue"] == 10.0  # noqa: E731
                                   and r["arm"] == arm], "{:.2f}")
        print(f"| {st[0]}{st[1]} | {ms([r['n_genes'] for r in a if r['arm'] == 'max+repair' and r['backend'] == 'diamond'], '{:.0f}')} | "
              f"{g('diamond', 'max+repair', 'score_s')} | {g('mmseqs', 'max+repair', 'score_s')} | "
              f"{g('hybrid', 'default', 'score_s')} | {g('hybrid', 'default', 'merge_s')} | "
              f"{g('diamond', 'avg+repair/quartet', 'merge_s')} |")
    print("\n(Wall times on a shared 24-CPU node, load average 35-80 during the runs; 3-4 threads per job.)")


if __name__ == "__main__":
    main()
