#!/usr/bin/env python3
"""Test/dev/all scores of the OMA comparators and paired bootstrap dF = F(Ancestra accurate) - F(OMA comparator)
(RefOGs resampled, B=2000, seed 11; non-inferiority margin 3 pp as DESIGN 5). Writes a markdown table to stdout.
Usage: oma_comparators_vs_ancestra.py <ancestra score.json> <split.json> <method_dir> [<method_dir> ...]"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orthobench_eval import bootstrap_delta_f, ni_verdict
anc = json.load(open(sys.argv[1])); split = json.load(open(sys.argv[2]))
rows = []
print("| comparator | test F | P | R | exact | dev F | all-70 F | ΔF test (Ancestra − comparator) [95% CI] | SE | NI verdict (Ancestra vs comparator) | ΔF all-70 [95% CI] | CPU-h | wall-h | peak RSS GB |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for d in sys.argv[3:]:
    s = json.load(open(f"{d}/score.json")); r = json.load(open(f"{d}/runtime.json"))
    bt = bootstrap_delta_f(anc["per_refog"], s["per_refog"], split["test"], 2000, 11)
    ba = bootstrap_delta_f(anc["per_refog"], s["per_refog"], None, 2000, 11)
    t = s["test"]
    print(f"| {os.path.basename(d.rstrip('/'))} | {t['F']:.1f} | {t['precision']:.1f} | {t['recall']:.1f} | {t['exact_refogs']}/36 | {s['dev']['F']:.1f} | {s['all']['F']:.1f} | "
          f"{bt['delta_point']:+.1f} [{bt['ci95'][0]:+.1f}, {bt['ci95'][1]:+.1f}] | {bt['se']:.1f} | {ni_verdict(bt['ci95'][0], 3.0)} | "
          f"{ba['delta_point']:+.1f} [{ba['ci95'][0]:+.1f}, {ba['ci95'][1]:+.1f}] | {r['cpu_h']:.2f} | {r['wall_h']:.2f} | {r['max_rss_gb']:.1f} |")
