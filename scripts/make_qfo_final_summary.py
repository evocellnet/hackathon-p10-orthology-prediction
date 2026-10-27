#!/usr/bin/env python3
"""Final QfO summary (freeze-v1): main table, ablation table, Pareto plot.

Inputs
* ours (full suite): work/qfo_final/cand2_{accurate,fast}/assess/*.json (bench_repo/scripts/run_qfo_local.sh,
  + FAS via scripts/qfo_final_fas.sh)
* ablations (fast tests): work/qfo_final/ablation/<arm>/assess/*.json (scripts/qfo_final_ablations.sh);
  baseline = cand2_fast (same fast tests, from its full-suite run)
* published QfO 2022 (OpenEBench OEBE0020000003): bench_repo/reference_results/openebench_QfO_2022.tsv
  (all 20 participants; the local pipeline reproduces their SwissTree/TreeFam-A/VGNC TPR/PPV exactly,
  results/qfo_peek1/competitor_reproduction.tsv)

Values are as reported by the QfO code: "±" = the reported stderr field, which for SwissTree/TreeFam-A/VGNC
is a 95% CI half-width (1.96*sqrt(var)) and for GO/EC/FAS/STD the standard error of the mean.
STD/G_STD2: FRAC_INCORRECT_TREES and RF_DISTANCE with n = NR_COMPLETED_TREE_SAMPLINGS (random tree
samplings, no seed exposed: stochastic). Lower is better for STD metrics.

Outputs: results/qfo_final/summary.{tsv,md}, ablation.tsv, pareto.png (+ pareto_data.tsv = the plotted values).
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PUB = os.path.join(ROOT, "bench_repo/reference_results/openebench_QfO_2022.tsv")
FIN = os.path.join(ROOT, "work/qfo_final")
OUT = os.path.join(ROOT, "results/qfo_final")
# published challenge -> our benchmark name
CH = {"SwissTree": "SwissTrees", "TreeFam-A": "TreeFam-A", "VGNC": "VGNC", "GOtest": "GO", "ECtest": "EC",
      "FAS": "FAS"}
STD = ["STD_Eukaryota", "STD_Fungi", "STD_Bacteria", "G_STD2_Luca", "G_STD2_Eukaryota", "G_STD2_Fungi",
       "G_STD2_Vertebrata"]
PRIMARY = ["OrthoFinder", "FastOMA_v0.1.4_nov", "Hieranoid_2", "InParanoid5", "SonicParanoid2-sens-g"]
# OMA resources are always shown explicitly (user request): FastOMA is a primary comparator; the others here
OMA = ["OMA_HOGs", "OMA_Pairs", "OMA_Groups"]
OURS = [("Ancestra accurate (cand2)", "cand2_accurate"), ("Ancestra fast (cand2)", "cand2_fast")]
ABL = [("a_nogo", "no gene order"), ("b_rescue", "+ family rescue (3Di+profile) + guard"),
       ("c_astral", "guide tree ASTRAL-Pro"), ("d_fastme", "guide tree FastME"),
       ("e_nni5", "NCBI + 5 random NNIs (seed 11)"), ("f_random", "random tree (negative control)"),
       ("g_accurate_nogo", "accurate preset, no gene order")]


def published() -> dict:
    out: dict = {}
    for r in csv.DictReader(open(PUB), delimiter="\t"):
        ch = r["challenge"]
        if "_agg_" not in r["dataset"]:
            continue
        b = CH.get(ch, ch)
        d = out.setdefault(r["tool"], {})
        if b in STD and r["x_metric"] == "NR_ORTHOLOGS":
            d[(b, "NR_ORTHOLOGS")] = (float(r["x"]), 0.0)
            continue
        d[(b, r["x_metric"])] = (float(r["x"]), float(r["stderr_x"] or 0))
        d[(b, r["y_metric"])] = (float(r["y"]), float(r["stderr_y"] or 0))
    return out


def ours(run: str) -> dict:
    d: dict = {}
    for fn in glob.glob(os.path.join(run, "assess", "*.json")):
        b = os.path.basename(fn)[:-5]
        try:
            data = json.load(open(fn))
        except Exception:
            continue
        for x in data if isinstance(data, list) else [data]:
            if x.get("type") != "assessment":
                continue
            cid = x.get("challenge_id", "")
            if cid != b and not (b == "SP" and cid == "SwissProtIDs"):
                continue
            m = x["metrics"]
            d[(b, m["metric_id"])] = (float(m["value"]), float(m.get("stderr") or 0))
    return d


def fmt(v, nd=3, n=False):
    if v is None:
        return ""
    val, se = v
    if n:
        return f"{int(round(val)):,}"
    return f"{val:.{nd}f}±{se:.{nd}f}" if se else f"{val:.{nd}f}"


def main_row(name: str, kind: str, d: dict) -> dict:
    g = d.get
    r = {"method": name, "kind": kind}
    for b, lab in (("SwissTrees", "SwissTree"), ("TreeFam-A", "TreeFam-A"), ("VGNC", "VGNC")):
        r[f"{lab} TPR"] = fmt(g((b, "TPR")))
        r[f"{lab} PPV"] = fmt(g((b, "PPV")))
    for b in ("GO", "EC"):
        r[f"{b} avg Schlicker"] = fmt(g((b, "avg Schlicker")), 4)
        r[f"{b} n"] = fmt(g((b, "NR_ORTHOLOGS")), n=True)
    r["FAS"] = fmt(g(("FAS", "FAS")), 4)
    r["FAS n"] = fmt(g(("FAS", "NR_ORTHOLOGS")), n=True)
    for b in STD:
        r[f"{b} frac incorrect"] = fmt(g((b, "FRAC_INCORRECT_TREES")), 4)
        r[f"{b} RF"] = fmt(g((b, "RF_DISTANCE")), 4)
        r[f"{b} n trees"] = fmt(g((b, "NR_COMPLETED_TREE_SAMPLINGS")), n=True)
    return r


def pareto(points: dict) -> set:
    """Names not dominated in (TPR, PPV) (>= both, > one)."""
    nd = set()
    for a, (x, y) in points.items():
        if not any((x2 >= x and y2 >= y and (x2 > x or y2 > y)) for b, (x2, y2) in points.items() if b != a):
            nd.add(a)
    return nd


def place_labels(ax, labels, ink, ink2, muted, gap=0.052) -> None:
    """Unambiguous callouts: labels sit in a column just right of the axes, ordered by the point's y and
    spread to a minimum gap; each is joined to its point by a hairline leader (dense clusters stay legible)."""
    y0, y1 = ax.get_ylim()
    items = sorted(((y - y0) / (y1 - y0), x, y, t, b) for x, y, t, b in labels)
    pos = [it[0] for it in items]
    for k in range(1, len(pos)):
        pos[k] = max(pos[k], pos[k - 1] + gap)
    top = 0.98
    for k in range(len(pos) - 1, -1, -1):
        lim = top if k == len(pos) - 1 else pos[k + 1] - gap
        pos[k] = min(pos[k], lim)
    for (yf, x, y, txt, bold), p in zip(items, pos):
        ax.annotate(txt, xy=(x, y), xycoords="data", xytext=(1.03, p), textcoords="axes fraction",
                    ha="left", va="center", fontsize=8.5 if bold else 8, color=ink if bold else ink2,
                    fontweight="bold" if bold else "normal", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=muted, lw=0.6, shrinkA=1, shrinkB=4, alpha=0.8))


def plot(pub: dict, our: dict, path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
    # validated 3-slot set (dataviz validate_palette.js --pairs all, light): blue / orange / aqua; aqua is
    # < 3:1 on the surface -> relief rule: every aqua point is direct-labelled (+ square marker, + table)
    C_OURS, C_PRIM, C_OMA = "#2a78d6", "#eb6834", "#1baf7a"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
                         "ytick.color": MUTED, "text.color": INK, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), facecolor=SURF)
    rows = []
    for ax, (b, title) in zip(axes, (("TreeFam-A", "TreeFam-A"), ("SwissTrees", "SwissTree"))):
        ax.set_facecolor(SURF)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        labels = []
        pts = {}
        for tool, d in pub.items():
            if (b, "TPR") in d:
                pts[tool] = (d[(b, "TPR")][0], d[(b, "PPV")][0])
        for name, d in our.items():
            if (b, "TPR") in d:
                pts[name] = (d[(b, "TPR")][0], d[(b, "PPV")][0])
        front = sorted((p for n, p in pts.items() if n in pareto(pts)), key=lambda p: p[0])
        if len(front) > 1:
            ax.plot([p[0] for p in front], [p[1] for p in front], color=MUTED, linewidth=1, zorder=1)
        for tool, d in pub.items():
            if (b, "TPR") not in d:
                continue
            (x, sx), (y, sy) = d[(b, "TPR")], d[(b, "PPV")]
            prim, oma = tool in PRIMARY, tool in OMA
            col = C_PRIM if prim else C_OMA if oma else "#b5b3ab"
            ax.errorbar(x, y, xerr=sx, yerr=sy, fmt="s" if oma else "o", ms=6 if (prim or oma) else 4.5,
                        color=col, ecolor=col if (prim or oma) else GRID,
                        elinewidth=0.8, mec=SURF, mew=1.2, zorder=3 if (prim or oma) else 2)
            if prim or oma:
                labels.append((x, y, tool.replace("_v0.1.4_nov", "").replace("_", " "), False))
            rows.append((title, tool, "primary comparator" if prim else "OMA resource" if oma else "other published",
                         x, sx, y, sy))
        for name, d in our.items():
            if (b, "TPR") not in d:
                continue
            (x, sx), (y, sy) = d[(b, "TPR")], d[(b, "PPV")]
            ax.errorbar(x, y, xerr=sx, yerr=sy, fmt="D" if "accurate" in name else "^", ms=6.5,
                        color=C_OURS, ecolor=C_OURS, elinewidth=1, mec=SURF, mew=1.0, zorder=4, alpha=0.95)
            labels.insert(0, (x, y, name.replace(" (cand2)", ""), True))
            rows.append((title, name, "ours", x, sx, y, sy))
        place_labels(ax, labels, INK, INK2, MUTED)
        ax.set_xlabel("TPR (recall)")
        ax.set_ylabel("PPV (precision)")
        ax.set_title(f"{title}: TPR vs PPV (95% CI); line = Pareto front", fontsize=10, color=INK, loc="left")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="D", ls="", color=C_OURS, label="Ancestra accurate"),
               Line2D([], [], marker="^", ls="", color=C_OURS, label="Ancestra fast"),
               Line2D([], [], marker="o", ls="", color=C_PRIM, label="primary comparators (incl. FastOMA)"),
               Line2D([], [], marker="s", ls="", color=C_OMA, label="OMA HOGs / Pairs / Groups"),
               Line2D([], [], marker="o", ls="", color="#b5b3ab", label="other QfO 2022 participants")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.06, right=0.86, bottom=0.17, top=0.92, wspace=0.55)
    fig.savefig(path, dpi=160, facecolor=SURF)
    with open(path.replace(".png", "_data.tsv"), "w") as fh:
        fh.write("benchmark\tmethod\tgroup\tTPR\tTPR_ci\tPPV\tPPV_ci\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    pub = published()
    our = {name: ours(os.path.join(FIN, d)) for name, d in OURS}
    rows = [main_row(n, "ours (final, freeze-v1)", d) for n, d in our.items()]
    order = PRIMARY + OMA + sorted(t for t in pub if t not in PRIMARY and t not in OMA)
    kind = lambda t: ("primary comparator" + (" (OMA resource)" if "OMA" in t else "")) if t in PRIMARY \
        else "OMA resource (published)" if t in OMA else "published QfO 2022"
    missing = [t for t in PRIMARY + OMA if t not in pub]
    if missing:
        raise SystemExit(f"published rows missing for {missing}")
    rows += [main_row(t, kind(t), pub[t]) for t in order]
    # Pareto flags on TreeFam-A and SwissTree over all methods
    for b, lab in (("TreeFam-A", "TreeFam-A"), ("SwissTrees", "SwissTree")):
        pts = {n: (d[(b, "TPR")][0], d[(b, "PPV")][0]) for n, d in {**pub, **our}.items() if (b, "TPR") in d}
        front = pareto(pts)
        for r in rows:
            r[f"{lab} Pareto"] = "yes" if r["method"] in front else ""
    cols = list(rows[0])
    with open(os.path.join(OUT, "summary.tsv"), "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ablations vs cand2_fast
    base = our["Ancestra fast (cand2)"]
    arows = []
    for arm, desc in [("cand2_fast", "baseline: fast preset, gene order, NCBI tree")] + ABL:
        d = base if arm == "cand2_fast" else ours(os.path.join(FIN, "ablation", arm))
        if not d:
            arows.append({"arm": arm, "description": desc, "status": "not available"})
            continue
        r = {"arm": arm, "description": desc, "status": "ok"}
        for b, lab in (("SwissTrees", "SwissTree"), ("TreeFam-A", "TreeFam-A"), ("VGNC", "VGNC")):
            for m in ("TPR", "PPV"):
                v, h = d.get((b, m), (math.nan, 0))
                bv, bh = base.get((b, m), (math.nan, 0))
                r[f"{lab} {m}"] = f"{v:.3f}±{h:.3f}"
                r[f"Δ{lab} {m}"] = "" if arm == "cand2_fast" else f"{v - bv:+.3f}"
        for b in ("GO", "EC"):
            v = d.get((b, "avg Schlicker"), (math.nan, 0))[0]
            r[f"{b} Schlicker"] = f"{v:.4f}"
            r[f"Δ{b}"] = "" if arm == "cand2_fast" else f"{v - base[(b, 'avg Schlicker')][0]:+.4f}"
            r[f"{b} n"] = f"{int(d.get((b, 'NR_ORTHOLOGS'), (0, 0))[0]):,}"
        arows.append(r)
    # context: published OMA resources (+ OrthoFinder) on the same tests (no Δ; reference only)
    for t in ["FastOMA_v0.1.4_nov"] + OMA + ["OrthoFinder"]:
        d = pub[t]
        r = {"arm": t, "description": "published reference (OpenEBench 2022)", "status": "reference"}
        for b, lab in (("SwissTrees", "SwissTree"), ("TreeFam-A", "TreeFam-A"), ("VGNC", "VGNC")):
            for m in ("TPR", "PPV"):
                v, h = d.get((b, m), (math.nan, 0))
                r[f"{lab} {m}"] = f"{v:.3f}±{h:.3f}"
                r[f"Δ{lab} {m}"] = f"{v - base[(b, m)][0]:+.3f}"
        for b in ("GO", "EC"):
            v = d.get((b, "avg Schlicker"), (math.nan, 0))[0]
            r[f"{b} Schlicker"] = f"{v:.4f}"
            r[f"Δ{b}"] = f"{v - base[(b, 'avg Schlicker')][0]:+.4f}"
            r[f"{b} n"] = f"{int(d.get((b, 'NR_ORTHOLOGS'), (0, 0))[0]):,}"
        arows.append(r)
    acols = []
    for r in arows:
        acols += [c for c in r if c not in acols]
    with open(os.path.join(OUT, "ablation.tsv"), "w") as fh:
        w = csv.DictWriter(fh, acols, delimiter="\t")
        w.writeheader()
        w.writerows(arows)

    # markdown
    L = ["# QfO 2022 — final evaluation (freeze-v1)", "",
         "Generated by `scripts/make_qfo_final_summary.py`. ± = reported stderr field (95% CI half-width for "
         "SwissTree/TreeFam-A/VGNC; SEM for GO/EC/FAS/STD). STD: lower is better; n = completed tree samplings.", "",
         "## Main table", ""]
    mcols = ["method", "SwissTree TPR", "SwissTree PPV", "TreeFam-A TPR", "TreeFam-A PPV", "VGNC TPR", "VGNC PPV",
             "GO avg Schlicker", "GO n", "EC avg Schlicker", "EC n", "FAS", "FAS n",
             "TreeFam-A Pareto", "SwissTree Pareto"]
    L += ["| " + " | ".join(mcols) + " |", "|" + "---|" * len(mcols)]
    L += ["| " + " | ".join(r.get(c, "") for c in mcols) + " |" for r in rows]
    L += ["", "## Species-tree discordance (frac incorrect trees ± SEM / RF ± SEM / n trees)", ""]
    scols = ["method"] + STD
    L += ["| " + " | ".join(scols) + " |", "|" + "---|" * len(scols)]
    for r in rows:
        cells = [r["method"]] + [f"{r[f'{b} frac incorrect']} / {r[f'{b} RF']} / {r[f'{b} n trees']}"
                                 if r[f"{b} frac incorrect"] else "" for b in STD]
        L.append("| " + " | ".join(cells) + " |")
    L += ["", "## Ablations (fast preset; Δ vs cand2_fast; last rows = published OMA resources + OrthoFinder "
          "for context, Δ = published − cand2_fast)", ""]
    ac = ["arm", "description", "SwissTree TPR", "ΔSwissTree TPR", "SwissTree PPV", "ΔSwissTree PPV",
          "TreeFam-A TPR", "ΔTreeFam-A TPR", "TreeFam-A PPV", "ΔTreeFam-A PPV", "VGNC PPV", "ΔVGNC PPV",
          "GO Schlicker", "ΔGO", "EC Schlicker", "ΔEC"]
    L += ["| " + " | ".join(ac) + " |", "|" + "---|" * len(ac)]
    L += ["| " + " | ".join(r.get(c, r.get("status", "")) if c in r else "" for c in ac) + " |" for r in arows]
    L += ["", "![Pareto](pareto.png)", ""]
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    plot(pub, our, os.path.join(OUT, "pareto.png"))
    print(f"[final] wrote {OUT}/summary.tsv, summary.md, ablation.tsv, pareto.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
