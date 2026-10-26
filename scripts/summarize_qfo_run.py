#!/usr/bin/env python3
"""Summarise a pipeline run directory (QfO run engineer report; no benchmark scoring).

Reads only the run's own outputs: ``stages/*/stage.json``, ``stages/families/families.npy``,
``stages/scores/families_final.npy``, ``manifest.json``, ``hogs.orthoxml``,
``stage_accounting.json`` (scripts/qfo_pipeline_driver.py) and optional ``/usr/bin/time -v``
logs. Writes JSON (+ a short Markdown table with ``--md``).

Reported: S0 families (count, singletons, largest, Σg²), families after the S2 cross-family
rescue, scoring statistics, root HOGs (count, largest), orthologGroups per TaxRange, HOGs per
level (``ancestra.orthoxml.hogs_at_level``, multi-gene HOGs), per-stage wall / CPU / RSS /
node load, and per-node merge statistics.

    python scripts/summarize_qfo_run.py RUN_DIR --json out.json [--md out.md] [--timelog log ...]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def fam_stats(fam: np.ndarray, top: int = 10) -> dict:
    sizes = np.bincount(fam)
    sizes = sizes[sizes > 0]
    srt = np.sort(sizes)[::-1]
    return {"n_families": int(len(sizes)), "n_multi_gene": int((sizes > 1).sum()),
            "n_singletons": int((sizes == 1).sum()),
            "frac_genes_singleton": round(float((sizes == 1).sum() / sizes.sum()), 4),
            "largest": [int(x) for x in srt[:top]],
            "n_over_5000": int((sizes > 5000).sum()), "n_over_7500": int((sizes > 7500).sum()),
            "sum_g2": int((sizes.astype(np.int64) ** 2).sum()),
            "sum_g2_over_N2": round(float((sizes.astype(np.float64) ** 2).sum() / sizes.sum() ** 2), 6)}


def parse_timelog(path: str) -> dict:
    out = {}
    pats = {"wall": r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)",
            "user_s": r"User time \(seconds\): (\S+)", "sys_s": r"System time \(seconds\): (\S+)",
            "cpu_pct": r"Percent of CPU this job got: (\S+)%",
            "max_rss_kb": r"Maximum resident set size \(kbytes\): (\d+)",
            "exit": r"Exit status: (\d+)"}
    txt = open(path).read()
    for k, p in pats.items():
        m = re.findall(p, txt)
        if m:
            out[k] = m[-1]
    if "wall" in out:
        parts = [float(x) for x in out["wall"].split(":")]
        out["wall_s"] = round(sum(x * 60 ** i for i, x in enumerate(reversed(parts))), 1)
    if "max_rss_kb" in out:
        out["max_rss_gb"] = round(int(out["max_rss_kb"]) / 1024 ** 2, 2)
    if "user_s" in out:
        out["cpu_h"] = round((float(out["user_s"]) + float(out["sys_s"])) / 3600, 3)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--json", required=True)
    ap.add_argument("--md")
    ap.add_argument("--timelog", nargs="*", default=[])
    ap.add_argument("--no-levels", action="store_true", help="skip HOGs-per-level (slow on big runs)")
    a = ap.parse_args(argv)
    from ancestra.newick import TreeIndex, read_newick
    from ancestra.orthoxml import Group, hogs_at_level, read_orthoxml

    rd = a.run_dir
    sd = os.path.join(rd, "stages")
    S: dict = {"run_dir": os.path.abspath(rd)}
    st = {s: json.load(open(os.path.join(sd, s, "stage.json")))
          for s in ("input", "families", "scores", "merge", "output")
          if os.path.exists(os.path.join(sd, s, "stage.json"))}
    S["stage_stats"] = {k: v.get("stats") for k, v in st.items()}
    S["input"] = st.get("input", {}).get("stats")
    if os.path.exists(os.path.join(sd, "families", "families.npy")):
        S["families_S0"] = fam_stats(np.load(os.path.join(sd, "families", "families.npy")))
    if os.path.exists(os.path.join(sd, "scores", "families_final.npy")):
        S["families_after_rescue"] = fam_stats(np.load(os.path.join(sd, "scores", "families_final.npy")))
    acc = os.path.join(rd, "stage_accounting.json")
    if os.path.exists(acc):
        runtime = {}
        for run in json.load(open(acc)):
            for s, e in run["stages"].items():
                if not e.get("skipped"):
                    runtime[s] = {k: e.get(k) for k in ("wall_s", "cpu_total_s", "cpu_self_s",
                                                        "cpu_children_s", "peak_rss_mb",
                                                        "children_maxrss_mb_so_far", "load_start",
                                                        "load_end")}
                    runtime[s]["threads"] = run.get("threads")
        S["runtime_per_stage"] = runtime
    S["timelogs"] = {os.path.basename(p): parse_timelog(p) for p in a.timelog if os.path.exists(p)}

    x = os.path.join(rd, "hogs.orthoxml")
    if os.path.exists(x):
        ox = read_orthoxml(x)
        tree = read_newick(os.path.join(sd, "input", "guide_tree.nwk"))
        tix = TreeIndex(tree)
        sp = {k: v.species for k, v in ox.genes.items()}
        tops = [g for g in ox.groups if isinstance(g, Group)]
        sizes = sorted(((len(g.genes()), len({sp[y] for y in g.genes()}), g.level or "?", g.kind)
                        for g in tops), reverse=True)
        in_groups = sum(s[0] for s in sizes)
        S["root_hogs"] = {"n_toplevel_groups": len(tops),
                          "n_toplevel_paralogGroups": sum(g.kind == "paralog" for g in tops),
                          "n_genes": len(ox.genes), "n_genes_in_groups": in_groups,
                          "n_singleton_genes": len(ox.genes) - in_groups,
                          "largest": [{"genes": a_, "species": b, "level": c, "kind": d}
                                      for a_, b, c, d in sizes[:10]],
                          "toplevel_by_level": dict(Counter(g.level for g in tops).most_common())}
        og = Counter()
        for t in tops:
            for g in t.iter_groups():
                if g.kind == "ortholog":
                    og[g.level] += 1
        S["orthologGroups_by_taxrange"] = dict(og.most_common())
        if not a.no_levels:
            per = {}
            for lvl in tix.internal:
                h = hogs_at_level(ox, lvl, tree, include_singletons=False)
                per[lvl] = {"depth": tix.depth[lvl], "n_species": len(tix.clade(lvl)),
                            "n_hogs": len(h), "n_genes_in_hogs": sum(len(s) for s in h),
                            "max_hog": max((len(s) for s in h), default=0)}
            S["hogs_per_level"] = per
    man = os.path.join(rd, "manifest.json")
    if os.path.exists(man):
        m = json.load(open(man))
        S["merge_nodes"] = m.get("nodes")
        S["total_wall_s_pipeline"] = m.get("total_wall_s")
        S["config"] = m.get("config")
        S["git"] = m.get("git")
    os.makedirs(os.path.dirname(os.path.abspath(a.json)), exist_ok=True)
    with open(a.json, "w") as fh:
        json.dump(S, fh, indent=1, default=str)

    if a.md:
        L = [f"# Run summary: `{os.path.relpath(rd, ROOT)}`", ""]
        if "runtime_per_stage" in S:
            L += ["| stage | wall s | CPU s (self+children) | peak RSS MB (own / max child) | load start→end |",
                  "|---|---|---|---|---|"]
            for s, e in S["runtime_per_stage"].items():
                L.append(f"| {s} | {e['wall_s']} | {e['cpu_total_s']} | {e['peak_rss_mb']} / "
                         f"{e['children_maxrss_mb_so_far']} | {e['load_start'][0] if e['load_start'] else ''}"
                         f"→{e['load_end'][0] if e['load_end'] else ''} |")
            L.append("")
        for k in ("families_S0", "families_after_rescue"):
            if k in S:
                f = S[k]
                L.append(f"- **{k}**: {f['n_families']} families ({f['n_multi_gene']} multi-gene, "
                         f"{f['n_singletons']} singletons = {100 * f['frac_genes_singleton']:.1f}% of genes); "
                         f"largest {f['largest']}; Σg²/N² = {f['sum_g2_over_N2']}")
        if "root_hogs" in S:
            r = S["root_hogs"]
            L.append(f"- **root HOGs**: {r['n_toplevel_groups']} top-level groups "
                     f"({r['n_genes_in_groups']} genes; {r['n_singleton_genes']} singleton genes); largest "
                     + ", ".join(f"{d['genes']} genes/{d['species']} sp @{d['level']}" for d in r["largest"][:5]))
        if "hogs_per_level" in S:
            L += ["", "| level | depth | #species | #HOGs (≥2 genes) | genes in HOGs | max HOG | orthologGroups with this TaxRange |",
                  "|---|---|---|---|---|---|---|"]
            for lvl, e in sorted(S["hogs_per_level"].items(), key=lambda kv: (kv[1]["depth"], kv[0])):
                L.append(f"| {lvl} | {e['depth']} | {e['n_species']} | {e['n_hogs']} | {e['n_genes_in_hogs']} | "
                         f"{e['max_hog']} | {S['orthologGroups_by_taxrange'].get(lvl, 0)} |")
        with open(a.md, "w") as fh:
            fh.write("\n".join(L) + "\n")
    print(f"[summarize] wrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
