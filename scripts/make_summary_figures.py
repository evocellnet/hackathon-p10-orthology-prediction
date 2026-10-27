#!/usr/bin/env python
"""Build the two summary figures for docs/ONE_PAGE_SUMMARY.md from committed result files.

fig1_method.png     -- schematic of the Ancestra pipeline
fig2_benchmarks.png -- (a) OrthoBench test F vs CPU-h, (b) QfO TreeFam-A TPR/PPV, (c) runtime scaling

Inputs: results/qfo_final/pareto_data.tsv, results/scaling/scaling.tsv,
OrthoBench numbers from results/orthobench_final/summary.md and
bench_repo/results/orthobench_competitors/ (copied below with their source paths).
"""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
os.makedirs(OUT, exist_ok=True)

# Palette: validated with the dataviz validator (light mode; aqua needs direct labels -> all marks labelled)
ACC, FAST, OMA = "#2a78d6", "#eb6834", "#1baf7a"
GREY, INK, MUTED, GRID = "#8a8a85", "#1f1f1d", "#5f5e58", "#e6e5df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})


def fig_method():
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.set_xlim(0, 11); ax.set_ylim(0, 3.2); ax.axis("off")
    boxes = [
        (0.1, "Inputs", "proteomes (FASTA)\nNCBI guide tree\n(opt.) gene order GFF\n(opt.) 3Di structure"),
        (2.3, "S0 families", "DIAMOND cluster 30% →\nall genes vs reps →\nconnected comps (+Leiden)\nnear-linear, deterministic"),
        (4.5, "S2 in-family\nscoring", "DIAMOND blocks +\nMMseqs2 re-align\nonly within families\n(no all-vs-all)"),
        (6.7, "S3 level-synchronous\nmerge", "bottom-up along species tree,\nall families at once: mutual-best\n+ in-paralogs; FastTree local trees\n+ species overlap; synteny re-rank"),
        (9.15, "Outputs", "HOGs (OrthoXML)\northolog pairs\nroot OGs\nmanifest + checksums"),
    ]
    for i, (x, title, body) in enumerate(boxes):
        w = [1.9, 1.9, 1.9, 2.2, 1.75][i]
        col = ACC if i in (1, 2, 3) else GREY
        ax.add_patch(FancyBboxPatch((x, 0.35), w, 2.3, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    fc="white", ec=col, lw=1.6))
        ax.text(x + w / 2, 2.3, title, ha="center", va="center", fontsize=9.5, fontweight="bold", color=INK)
        ax.text(x + w / 2, 1.25, body, ha="center", va="center", fontsize=7.6, color=INK, linespacing=1.35)
        if i < 4:
            ax.annotate("", xy=(boxes[i + 1][0] - 0.02, 1.5), xytext=(x + w + 0.02, 1.5),
                        arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.4))
    ax.text(5.5, 3.02, "Ancestra: orthology as level-synchronous ancestral-genome reconstruction inside de novo families",
            ha="center", fontsize=11, fontweight="bold", color=INK)
    ax.text(5.5, 0.1, "Cost ≈ clustering(N) + Σ_f g_f² (within families) + merge; no reference orthology DB; presets: accurate (FastTree) / fast (NJ)",
            ha="center", fontsize=8.2, color=MUTED)
    fig.savefig(os.path.join(OUT, "fig1_method.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_bench():
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0))
    # (a) OrthoBench TEST (36 RefOGs): F vs CPU-h. Sources: results/orthobench_final/summary.md,
    #     bench_repo/results/orthobench_competitors/*/score.json + runtime.json
    ob = [("Ancestra accurate", 76.2, 3.91, ACC), ("Ancestra fast", 74.8, 2.88, FAST),
          ("OrthoFinder 3.2", 69.7, 8.2, GREY), ("SonicParanoid2", 70.4, 7.9, GREY),
          ("ProteinOrtho", 47.5, 1.6, GREY), ("FastOMA* (flagged)", 33.5, 12.7, OMA)]
    ax = axs[0]
    for name, f, cpu, c in ob:
        ax.scatter(cpu, f, s=55, color=c, edgecolor="white", linewidth=1.5, zorder=3)
        off = {"Ancestra accurate": (1.06, 3.0), "Ancestra fast": (0.62, -3.6), "SonicParanoid2": (1.06, -3.4)}
        mx, dy = off.get(name, (1.06, 1.8))
        ax.text(cpu * mx, f + dy, name, fontsize=8, color=INK, va="center")
    ax.set_xscale("log"); ax.set_xlim(1.2, 40); ax.set_ylim(25, 85)
    from matplotlib.ticker import FixedLocator, NullLocator, FixedFormatter
    ax.xaxis.set_major_locator(FixedLocator([2, 5, 10, 20])); ax.xaxis.set_major_formatter(FixedFormatter(["2", "5", "10", "20"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("CPU-hours (OrthoBench, 251k proteins)"); ax.set_ylabel("OrthoBench test F (36 RefOGs)")
    ax.set_title("a  OrthoBench test: accuracy vs cost", loc="left", fontsize=10, color=INK)
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.text(1.25, 27, "Δ vs OrthoFinder +6.5 F [95% CI −6.5, +18.1]\n*FastOMA run flagged (export/DB issues; re-run pending)",
            fontsize=7.2, color=MUTED)
    # (b) QfO TreeFam-A Pareto. Source: results/qfo_final/pareto_data.tsv
    ax = axs[1]
    rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "results/qfo_final/pareto_data.tsv")), delimiter="\t")
            if r["benchmark"] == "TreeFam-A"]
    keep = {"OrthoFinder": "OrthoFinder", "Hieranoid_2": "Hieranoid 2", "FastOMA": "FastOMA", "OMA_HOGs": "OMA HOGs",
            "OMA_Pairs": "OMA Pairs", "InParanoid": "InParanoid 5", "SonicParanoid2": "SonicParanoid2"}
    for r in rows:
        m = r["method"]; tpr, ppv = float(r["TPR"]), float(r["PPV"])
        if "accurate" in m.lower():
            c, lab = ACC, "Ancestra accurate"
        elif "fast" in m.lower() and "ancestra" in (m.lower() + r["group"].lower()) or m.lower().startswith("cand2_fast"):
            c, lab = FAST, "Ancestra fast"
        elif "OMA" in m:
            c, lab = OMA, ("FastOMA" if m.startswith("FastOMA") else m.replace("_", " "))
        else:
            c, lab = GREY, next((v for k, v in keep.items() if m.startswith(k)), None)
        ax.scatter(tpr, ppv, s=50 if c != GREY else 28, color=c, edgecolor="white", linewidth=1.2, zorder=3)
        offs = {"Ancestra accurate": (-0.012, -0.016), "Ancestra fast": (0.008, -0.002), "FastOMA": (-0.06, -0.009),
                "OMA HOGs": (-0.03, 0.012), "OMA Pairs": (-0.075, 0.004), "OMA Groups": (0.006, -0.003),
                "OrthoFinder": (-0.03, -0.008), "Hieranoid 2": (0.006, 0.004)}
        if lab in offs and not m.startswith("SonicParanoid2-sens") and not m.startswith("InParanoid"):
            dx, dy = offs[lab]
            ax.text(tpr + dx, ppv + dy, lab, fontsize=7.5, color=INK)
    ax.set_xlabel("TreeFam-A TPR (recall)"); ax.set_ylabel("TreeFam-A PPV (precision)")
    ax.set_title("b  QfO 2022 TreeFam-A (all published methods)", loc="left", fontsize=10, color=INK)
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    # (c) scaling. Source: results/scaling/scaling.tsv
    ax = axs[2]
    sc = [r for r in csv.DictReader(open(os.path.join(ROOT, "results/scaling/scaling.tsv")), delimiter="\t")
          if r["status"] == "ok" and not r["method"].startswith("Ancestra S0")]
    groups = {}
    for r in sc:
        groups.setdefault(r["method"], []).append((float(r["n_proteins"]), float(r["cpu_h"])))
    style = {"accurate": (ACC, "Ancestra accurate"), "fast": (FAST, "Ancestra fast"),
             "OrthoFinder": (GREY, "OrthoFinder 3.2"), "SonicParanoid": ("#b5b4ad", "SonicParanoid2"),
             "FastOMA": (OMA, "FastOMA")}
    for m, pts in groups.items():
        key = next((k for k in style if k in m), None)
        if key is None:
            continue
        c, lab = style[key]
        pts.sort()
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        ax.plot(xs, ys, "-o", color=c, lw=2, ms=5, mec="white", mew=1)
        if lab.startswith("Ancestra"):
            ax.text(xs[-1] * 1.1, ys[-1], lab, fontsize=7.5, color=INK, va="center")
        else:
            ax.text(xs[-1] / 1.15, ys[-1] * 1.18, lab, fontsize=7.5, color=INK, va="bottom", ha="right")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(4e4, 4e6)
    ax.set_xlabel("proteins (nested QfO subsets 8→78 species)"); ax.set_ylabel("CPU-hours (8–12 threads, shared node)")
    ax.set_title("c  Runtime scaling", loc="left", fontsize=10, color=INK)
    ax.grid(color=GRID, lw=0.6, which="major"); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_benchmarks.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_method()
    fig_bench()
    print("wrote", OUT)
