"""Schematic pictures and plots for the hackathon slide deck (technical_report/slides/).

Run from the repository root: python technical_report/slides/make_assets.py
All numbers are copied from REPORT.md (sections 3.1-3.9); the schematics are illustrations, not data.
"""
import os
import random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle, Ellipse

OUT = "technical_report/slides/assets"
os.makedirs(OUT, exist_ok=True)
NAVY, BLUE, GREEN, ORANGE, RED, GREY, LIGHT, INK = (
    "#1F3A5F", "#3E6FA3", "#4E8A5B", "#D9822B", "#B5473A", "#9AA1AB", "#E9EEF5", "#22303F")
PAL = [BLUE, ORANGE, GREEN, "#8E6BB0", RED, "#C9A227"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13, "axes.edgecolor": GREY})


def canvas(w=6.4, h=4.0):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, w * 10)
    ax.set_ylim(0, h * 10)
    ax.axis("off")
    return fig, ax


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=200, bbox_inches="tight", pad_inches=0.05, transparent=False,
                facecolor="white")
    plt.close(fig)


def box(ax, x, y, w, h, fc, ec=None, r=1.2, lw=0, alpha=1.0, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc,
                                ec=ec or fc, lw=lw, alpha=alpha, zorder=z))


def arrow(ax, x1, y1, x2, y2, c=GREY, lw=2.0, style="-|>", ms=14, ls="-", rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms, color=c, lw=lw,
                                 linestyle=ls, connectionstyle=f"arc3,rad={rad}", zorder=3))


def label(ax, x, y, s, c=INK, size=13, ha="center", va="center", weight="normal", style="normal"):
    ax.text(x, y, s, color=c, fontsize=size, ha=ha, va=va, fontweight=weight, fontstyle=style, zorder=5)


# ---------------------------------------------------------------- method: three steps
def method():
    fig, ax = canvas(13, 4.6)
    rnd = random.Random(11)
    # proteomes
    for i, c in enumerate([BLUE, ORANGE, GREEN, "#8E6BB0"]):
        y = 34 - i * 8.5
        box(ax, 2, y, 16, 6, LIGHT, r=1)
        for j in range(7):
            ax.add_patch(Circle((4 + j * 2.1, y + 3), 0.75, fc=PAL[rnd.randrange(4)], ec="none", zorder=2))
        label(ax, 10, y + 7.6 if i == 0 else y - 1.3, "", size=1)
    label(ax, 10, 42.5, "Proteomes", weight="bold", size=14)
    arrow(ax, 20, 22, 27, 22)
    # S0 families
    label(ax, 39, 42.5, "S0  Families", weight="bold", size=14)
    for k, (cx, cy, c) in enumerate([(33, 31, BLUE), (45, 31, ORANGE), (33, 13, GREEN), (45, 13, "#8E6BB0")]):
        ax.add_patch(Ellipse((cx, cy), 10.5, 13, fc=c, alpha=0.18, ec=c, lw=1.6))
        for j in range(9):
            ax.add_patch(Circle((cx - 3 + (j % 3) * 3, cy - 3.5 + (j // 3) * 3.5), 0.8, fc=c, ec="none"))
    label(ax, 39, 3.0, "no all-vs-all search", c=GREY, size=12)
    arrow(ax, 52, 22, 58, 22)
    # S2 scoring matrix
    label(ax, 66, 42.5, "S2  Scores", weight="bold", size=14)
    for i in range(6):
        for j in range(6):
            v = 0.25 + 0.75 * (1 - abs(i - j) / 5) * rnd.uniform(0.6, 1)
            ax.add_patch(Rectangle((60 + j * 2.2, 13 + i * 2.2), 2.1, 2.1, fc=BLUE, alpha=min(v, 1), ec="white"))
    label(ax, 66, 3.0, "inside families only", c=GREY, size=12)
    arrow(ax, 75, 22, 81, 22)
    # S3 merge along tree
    label(ax, 100, 42.5, "S3  Merge up the species tree", weight="bold", size=14)
    leaves = [84, 92, 100, 108, 116]
    ytop = 9
    for x, c in zip(leaves, [BLUE, ORANGE, GREEN, "#8E6BB0", RED]):
        ax.plot([x, x], [ytop, ytop + 6], color=INK, lw=2)
        ax.add_patch(Circle((x, ytop - 2), 1.2, fc=c, ec="none"))
    ax.plot([84, 92], [ytop + 6] * 2, color=INK, lw=2); ax.plot([88, 88], [ytop + 6, ytop + 13], color=INK, lw=2)
    ax.plot([108, 116], [ytop + 6] * 2, color=INK, lw=2); ax.plot([112, 112], [ytop + 6, ytop + 13], color=INK, lw=2)
    ax.plot([100, 100], [ytop + 6, ytop + 13], color=INK, lw=2)
    ax.plot([100, 112], [ytop + 13] * 2, color=INK, lw=2); ax.plot([106, 106], [ytop + 13, ytop + 20], color=INK, lw=2)
    ax.plot([88, 106], [ytop + 20] * 2, color=INK, lw=2); ax.plot([97, 97], [ytop + 20, ytop + 26], color=INK, lw=2)
    for (x, y) in [(88, ytop + 13), (106, ytop + 20), (97, ytop + 26)]:
        box(ax, x - 3.4, y - 1.6, 6.8, 3.2, ORANGE, r=0.8, z=4)
        label(ax, x, y, "HOG", c="white", size=9, weight="bold")
    label(ax, 100, 3.0, "best matches + small local trees", c=GREY, size=12)
    arrow(ax, 120, 22, 126, 22)
    label(ax, 128.5, 22, "OrthoXML\nHOGs", weight="bold", size=14, ha="left")
    ax.set_xlim(0, 140)
    save(fig, "method_schema")


# ---------------------------------------------------------------- QfO non-inferiority forest
def forest():
    rows = [("FastOMA", 0.043, -0.005, 0.093), ("OMA Pairs", 0.018, -0.028, 0.069),
            ("OMA HOGs", 0.002, -0.058, 0.067), ("Hieranoid 2", -0.036, -0.069, 0.004),
            ("InParanoid 5", -0.042, -0.081, -0.000), ("SonicParanoid2-sens-g", -0.087, -0.132, -0.038),
            ("OrthoFinder", -0.160, -0.210, -0.106)]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for i, (n, d, lo, hi) in enumerate(rows):
        y = len(rows) - i
        c = GREEN if lo >= -0.03 else (ORANGE if hi >= -0.03 else RED)
        ax.plot([lo, hi], [y, y], color=c, lw=3.2, solid_capstyle="round")
        ax.plot(d, y, "o", color=c, ms=10)
    ax.axvline(0, color=GREY, lw=1)
    ax.axvline(-0.03, color=RED, lw=1.4, ls="--")
    ax.text(-0.031, len(rows) + 0.75, "non-inferiority\nmargin −0.03", color=RED, fontsize=10, ha="right", va="center")
    ax.set_yticks(range(len(rows), 0, -1), [r[0] for r in rows])
    ax.set_xlabel("Δ TreeFam-A recall (Ancestra − method), 95% CI")
    ax.set_ylim(0.4, len(rows) + 1.3)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    save(fig, "qfo_forest")


# ---------------------------------------------------------------- layer pictures
def genes(ax, x0, y, cols, w=4.2, gap=0.8, h=3.2, labels=None):
    for i, c in enumerate(cols):
        box(ax, x0 + i * (w + gap), y, w, h, c, r=0.6, z=3)
        if labels and labels[i]:
            label(ax, x0 + i * (w + gap) + w / 2, y + h / 2, labels[i], c="white", size=10, weight="bold")
    return [x0 + i * (w + gap) + w / 2 for i in range(len(cols))]


def synteny():
    fig, ax = canvas()
    label(ax, 3, 33, "human", ha="left", c=GREY, size=12)
    label(ax, 3, 16, "mouse", ha="left", c=GREY, size=12)
    ax.plot([2, 62], [28, 28], color=GREY, lw=1.5, zorder=1)
    ax.plot([2, 62], [11, 11], color=GREY, lw=1.5, zorder=1)
    top = genes(ax, 8, 26.4, [GREEN, BLUE, ORANGE, BLUE, GREEN], labels=["", "", "A", "", ""])
    bot = genes(ax, 8, 9.4, [GREEN, BLUE, ORANGE, BLUE, GREEN], labels=["", "", "A1", "", ""])
    for a, b in zip(top, bot):
        ax.plot([a, a + (b - a)], [26.2, 12.8], color=GREEN if a != top[2] else ORANGE, lw=1.6, alpha=0.7)
    box(ax, 46, 3.2, 4.2, 3.2, ORANGE, r=0.6, z=3); label(ax, 48.1, 4.8, "A2", c="white", size=10, weight="bold")
    label(ax, 52, 4.8, "paralog: no shared neighbours", ha="left", c=RED, size=10)
    label(ax, 32, 37.5, "Neighbours agree → ortholog", weight="bold", size=13)
    save(fig, "layer_synteny")


def structure():
    fig, ax = canvas()
    import numpy as np
    t = np.linspace(0, 1, 300)
    for k, (y0, c, seq, di) in enumerate([(21, BLUE, "MKVLAEGTWRS", "DPVLLVCQVVD"), (4, ORANGE, "MRIFSDNSWKA", "DPVLLVCQVLD")]):
        ax.plot(4 + 22 * t, y0 + 5 + 3.2 * np.sin(t * 6 * np.pi) * (t < 0.55) + 2 * (t >= 0.55) * np.sin(t * 2 * np.pi),
                color=c, lw=3)
        label(ax, 30, y0 + 7.5, seq, ha="left", size=12, c=GREY)
        label(ax, 30, y0 + 3, di, ha="left", size=12, c=c, weight="bold")
    label(ax, 32, 35.5, "grey: amino acids · colour: 3Di structure letters", size=10, c=GREY)
    label(ax, 32, 39.5, "Similar fold, little shared sequence", weight="bold")
    label(ax, 32, 0.5, "finds distant homologs, not orthology", c=RED, size=11)
    save(fig, "layer_structure")


def expression():
    fig, ax = canvas()
    rnd = random.Random(3)
    rows = [("human A", [0.9, 0.8, 0.1, 0.2, 0.05]), ("mouse A1", [0.85, 0.7, 0.2, 0.1, 0.1]),
            ("mouse A2", [0.8, 0.75, 0.15, 0.25, 0.1])]
    tissues = ["brain", "liver", "heart", "lung", "skin"]
    for i, (n, vals) in enumerate(rows):
        y = 26 - i * 7
        label(ax, 13, y + 3, n, ha="right", size=12)
        for j, v in enumerate(vals):
            ax.add_patch(Rectangle((15 + j * 7, y), 6.6, 6.2, fc=GREEN, alpha=0.12 + 0.88 * v, ec="white"))
    for j, t in enumerate(tissues):
        label(ax, 18.3 + j * 7, 36, t, size=10, c=GREY)
    label(ax, 32, 6, "both copies look alike → no tie-break", c=RED, size=11)
    save(fig, "layer_expression")


def isoforms():
    fig, ax = canvas()
    exons = [(4, 6), (14, 5), (24, 7), (36, 5), (46, 6)]
    for y, use in [(30, [0, 1, 2, 3, 4]), (20, [0, 2, 3, 4]), (10, [0, 1, 3])]:
        ax.plot([4, 52], [y + 1.5, y + 1.5], color=GREY, lw=1.2)
        for k in use:
            x, w = exons[k]
            box(ax, x, y, w, 3, BLUE if y == 30 else ORANGE, r=0.4, z=3)
    label(ax, 56, 31.5, "canonical", ha="left", size=11, c=BLUE)
    label(ax, 56, 16, "alternative\nisoforms", ha="left", size=11, c=ORANGE)
    label(ax, 32, 38, "One gene, several proteins", weight="bold")
    label(ax, 32, 3, "2–5× more searching, no clear gain", c=RED, size=11)
    save(fig, "layer_isoforms")


def genomic():
    fig, ax = canvas()
    for y, name, mid in [(26, "dog", GREEN), (10, "human", None)]:
        label(ax, 3, y + 6.2, name, ha="left", c=GREY, size=12)
        ax.plot([2, 62], [y + 1.6, y + 1.6], color=GREY, lw=1.5)
        genes(ax, 8, y, [BLUE, ORANGE], w=7)
        genes(ax, 44, y, [BLUE, ORANGE], w=7)
        if mid:
            box(ax, 26, y, 10, 3.2, GREEN, r=0.6, z=3); label(ax, 31, y + 1.6, "UOX", c="white", size=10, weight="bold")
        else:
            box(ax, 26, y, 10, 3.2, "white", ec=GREEN, lw=1.8, r=0.6, z=3)
            label(ax, 31, y + 1.6, "UOX", c=GREEN, size=10, weight="bold")
            ax.plot([29.5], [y + 5.8], marker="*", color=RED, ms=16)
            label(ax, 31.5, y - 3.5, "stop codons → pseudogene", c=RED, size=11)
    arrow(ax, 31, 25.5, 31, 15, c=GREEN, ls="--")
    label(ax, 34, 20.3, "search the DNA\nat the syntenic locus", ha="left", size=10, c=GREEN)
    save(fig, "layer_genomic")


def tree(ax, x0, y0, w, h, order, c, title, conflict):
    xs = [x0 + i * w / 4 for i in range(5)]
    for x, k in zip(xs, order):
        ax.add_patch(Circle((x, y0 - 1.5), 1.0, fc=PAL[k], ec="none"))
        ax.plot([x, x], [y0, y0 + h * 0.3], color=c, lw=2)
    ax.plot([xs[0], xs[1]], [y0 + h * 0.3] * 2, color=c, lw=2)
    m1 = (xs[0] + xs[1]) / 2
    ax.plot([m1, m1], [y0 + h * 0.3, y0 + h * 0.65], color=c, lw=2)
    ax.plot([xs[3], xs[4]], [y0 + h * 0.3] * 2, color=c, lw=2)
    m2 = (xs[3] + xs[4]) / 2
    ax.plot([m2, m2], [y0 + h * 0.3, y0 + h * 0.5], color=c, lw=2)
    ax.plot([xs[2], xs[2]], [y0 + h * 0.3, y0 + h * 0.5], color=c, lw=2)
    ax.plot([xs[2], m2], [y0 + h * 0.5] * 2, color=c, lw=2)
    m3 = (xs[2] + m2) / 2
    ax.plot([m3, m3], [y0 + h * 0.5, y0 + h * 0.65], color=c, lw=2)
    ax.plot([m1, m3], [y0 + h * 0.65] * 2, color=c, lw=2)
    label(ax, x0 + w / 2, y0 + h + 2, title, weight="bold", size=12)
    label(ax, x0 + w / 2, y0 - 6, conflict, c=GREY, size=10)


def species_tree():
    fig, ax = canvas(7.4, 4.0)
    tree(ax, 3, 14, 18, 18, [0, 1, 2, 3, 4], GREEN, "NCBI taxonomy", "default")
    tree(ax, 27, 14, 18, 18, [1, 0, 2, 4, 3], ORANGE, "from the data", "12–17% conflict")
    tree(ax, 51, 14, 18, 18, [3, 0, 4, 1, 2], RED, "random", "negative control")
    save(fig, "layer_speciestree")


def network():
    fig, ax = canvas()
    import math
    for cx, name, col in [(16, "human", BLUE), (48, "mouse", ORANGE)]:
        pts = [(cx + 9 * math.cos(a), 20 + 9 * math.sin(a)) for a in [0.3, 1.5, 2.6, 3.8, 5.0]]
        for p in pts:
            ax.plot([cx, p[0]], [20, p[1]], color=GREY, lw=1.4, zorder=1)
            ax.add_patch(Circle(p, 1.6, fc=LIGHT, ec=GREY, zorder=2))
        ax.add_patch(Circle((cx, 20), 2.4, fc=col, ec="none", zorder=3))
        label(ax, cx, 34, name, c=GREY, size=12)
    ax.plot([18.5, 45.5], [20, 20], color=GREEN, lw=2, ls="--")
    label(ax, 32, 13, "shared partners?", c=GREEN, size=11)
    label(ax, 32, 4, "only 3.4% of candidates measured", c=RED, size=11)
    save(fig, "layer_network")


def abundance():
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    names = ["human A", "mouse A1", "mouse A2"]
    vals = [120, 95, 70]
    ax.bar(names, vals, color=[BLUE, ORANGE, ORANGE], width=0.6)
    ax.set_ylabel("ppm (PaxDb)")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("Which copy is as abundant as the human gene?", fontsize=12, color=INK)
    ax.text(1, 112, "right in 39 of 116\nhard decisions", ha="center", color=RED, fontsize=11)
    fig.tight_layout()
    save(fig, "layer_abundance")


def domains():
    fig, ax = canvas()
    arch = [("human A", [(4, 10, BLUE, "kinase"), (18, 8, GREEN, "SH2")]),
            ("mouse A1", [(4, 10, BLUE, "kinase"), (18, 8, GREEN, "SH2")]),
            ("mouse A2", [(4, 10, BLUE, "kinase"), (18, 8, GREEN, "SH2")]),
            ("mouse B", [(4, 10, BLUE, "kinase"), (18, 8, ORANGE, "PH")])]
    for i, (n, ds) in enumerate(arch):
        y = 31 - i * 8
        ax.plot([4, 34], [y + 1.6, y + 1.6], color=GREY, lw=1.5)
        for x, w, c, t in ds:
            box(ax, x, y, w, 3.2, c, r=0.8, z=3); label(ax, x + w / 2, y + 1.6, t, c="white", size=10, weight="bold")
        label(ax, 37, y + 1.6, n, ha="left", size=12)
    label(ax, 57, 21, "same\narchitecture", c=RED, size=11)
    label(ax, 57, 5, "other family", c=GREEN, size=11)
    save(fig, "layer_domains")


def plm():
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    rnd = random.Random(5)
    for cx, cy, c, n in [(0, 0, BLUE, "family 1"), (5, 3, ORANGE, "family 2"), (2, 6, GREEN, "family 3")]:
        xs = [cx + rnd.gauss(0, 0.8) for _ in range(40)]
        ys = [cy + rnd.gauss(0, 0.8) for _ in range(40)]
        ax.scatter(xs, ys, s=28, color=c, alpha=0.6, label=n)
    ax.scatter([0.3], [0.2], s=160, marker="*", color=NAVY, zorder=5)
    ax.scatter([0.6], [-0.1], s=110, marker="^", color=RED, zorder=5)
    ax.scatter([-0.1], [0.6], s=110, marker="^", color=GREEN, zorder=5)
    ax.annotate("query", (0.3, 0.2), (-1.0, 3.6), arrowprops=dict(arrowstyle="->", color=NAVY), color=NAVY)
    ax.annotate("ortholog and paralog\nequally close", (0.5, 0.0), (3.2, -1.2),
                arrowprops=dict(arrowstyle="->", color=RED), color=RED, fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("embedding dimension 1"); ax.set_ylabel("embedding dimension 2")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    fig.tight_layout()
    save(fig, "layer_plm")


def introns():
    fig, ax = canvas()
    for y, name, cuts in [(28, "human A", [14, 28, 41]), (17, "mouse A1", [14, 28, 41]), (6, "mouse A2", [14, 28, 41])]:
        box(ax, 4, y, 44, 3.4, LIGHT, ec=GREY, lw=1, r=0.6)
        for c in cuts:
            ax.plot([4 + c, 4 + c], [y - 0.8, y + 4.2], color=ORANGE, lw=3)
        label(ax, 51, y + 1.7, name, ha="left", size=12)
    label(ax, 26, 37.5, "intron positions on the protein", weight="bold", size=12)
    label(ax, 26, 1.5, "kept by both copies → weak (AUC 0.58)", c=RED, size=11)
    save(fig, "layer_introns")


def gpu():
    fig, ax = canvas()
    rnd = random.Random(7)
    sizes = [rnd.lognormvariate(0.6, 0.6) for _ in range(30)]
    x, y = 2, 2
    for s in sizes:
        w = min(2 + s * 1.4, 8)
        if x + w > 62:
            x, y = 2, y + 5
        if y > 14:
            break
        box(ax, x, y, w, 4, GREEN, r=0.5)
        x += w + 0.8
    for i, w in enumerate([28, 20, 12]):
        box(ax, 2 + sum([28, 20, 12][:i]) + i * 0.8, 20, w, 5, ORANGE, r=0.5)
    label(ax, 2, 29, "few large families → DIAMOND (CPU)", c=ORANGE, size=12, weight="bold", ha="left")
    label(ax, 2, 33.5, "many families ≤1,000 genes → Smith–Waterman on the GPU", c=GREEN, size=12, weight="bold", ha="left")
    save(fig, "layer_gpu")


def agents():
    fig, ax = canvas(9, 4.4)
    box(ax, 29, 32, 32, 9, NAVY, r=1.5)
    label(ax, 45, 36.5, "PI agent  (Opus 5.5)", c="white", size=13, weight="bold")
    groups = [("Opus 5.5", ["literature", "design reviews", "algorithms", "verification"], BLUE, 4),
              ("Sonnet 5.5", ["environment", "competitor runs", "scaling", "bookkeeping"], GREEN, 48)]
    for name, items, c, x0 in groups:
        label(ax, x0 + 19, 24, name + " subagents", c=c, size=12, weight="bold")
        for i, t in enumerate(items):
            box(ax, x0 + (i % 2) * 19.5, 13 - (i // 2) * 8.5, 18, 6.5, c, r=1.2, alpha=0.15)
            label(ax, x0 + (i % 2) * 19.5 + 9, 16.2 - (i // 2) * 8.5, t, c=INK, size=11)
        arrow(ax, 45, 31.5, x0 + 19, 26.5, c=c)
    label(ax, 45, 2, "git repository = shared memory: research tree · claims ledger · lessons · disclosures",
          c=GREY, size=11)
    save(fig, "agents")


for f in [method, forest, synteny, structure, expression, isoforms, genomic, species_tree, network, abundance,
          domains, plm, introns, gpu, agents]:
    f()
print("ok", sorted(os.listdir(OUT)))
