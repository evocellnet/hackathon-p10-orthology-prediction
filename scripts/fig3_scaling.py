"""Figure 3: CPU time vs number of proteins (measured points only, no extrapolation).

A: nested QfO subsets, 8 threads (results/scaling/scaling.tsv).
B: Ancestra on 125-1000 UniProt reference proteomes, 12 threads (results/scaling_large/scaling.tsv).
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

t2 = rows("results/scaling/scaling.tsv")
t3 = rows("results/scaling_large/scaling.tsv")
series = [
    ("Ancestra default (accurate preset, seq-only)", "Ancestra accurate", "#1b9e77", "o"),
    ("Ancestra fast preset (NJ merge, seq-only)", "Ancestra fast", "#66a61e", "s"),
    ("OrthoFinder 3.2", "OrthoFinder 3.2", "#d95f02", "^"),
    ("SonicParanoid2 (default)", "SonicParanoid2", "#7570b3", "D"),
    ("FastOMA 0.5.1 (OMAmer LUCA, DB build not counted)", "FastOMA 0.5.1", "#e7298a", "v"),
]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2))
for key, lab, col, mk in series:
    ok = [r for r in t2 if r["method"] == key and r["status"] == "ok"]
    a.plot([int(r["n_proteins"]) for r in ok], [float(r["cpu_h"]) for r in ok], marker=mk, color=col, label=lab)
    for r in t2:
        if r["method"] == key and r["status"] != "ok":
            x, y = int(r["n_proteins"]), float(r["cpu_h"])
            if r["status"].startswith(">"):
                a.annotate("", xy=(x, y * 1.6), xytext=(x, y), arrowprops=dict(arrowstyle="->", color=col))
                a.plot([x], [y], marker=mk, mfc="none", color=col, ls="none")
                a.text(x * 1.05, y * 0.9, "stopped at\n6 h cap", fontsize=7, color=col)
            else:
                a.plot([x], [y], marker="x", color=col, ls="none", ms=8)
                a.text(x * 1.05, y * 0.75, "failed", fontsize=7, color=col)
a.set(xscale="log", yscale="log", xlabel="proteins", ylabel="CPU-hours",
      title="A  Nested QfO subsets (8–78 species)")
a.legend(fontsize=7, loc="upper left")
fast = sorted((r for r in t3 if r["preset"] == "fast"), key=lambda r: int(r["n_proteins"]))
xs = [int(r["n_proteins"]) for r in fast]
b.plot(xs, [float(r["cpu_h"]) for r in fast], "o-", color="#66a61e", label="CPU-hours (fast)")
b.plot(xs, [float(r["max_rss_gb"]) for r in fast], "s--", color="#666666", label="peak memory, GB (fast)")
acc = [r for r in t3 if r["preset"] == "accurate"]
b.plot([int(r["n_proteins"]) for r in acc], [float(r["cpu_h"]) for r in acc], "o", color="#1b9e77",
       label="CPU-hours (accurate, 250 proteomes)")
for r in fast:
    b.annotate(f'{r["n_species"]} proteomes', (int(r["n_proteins"]), float(r["cpu_h"])),
               xytext=(4, -12), textcoords="offset points", fontsize=7)
b.set(xscale="log", yscale="log", xlabel="proteins", title="B  UniProt reference proteomes (125–1,000)")
b.legend(fontsize=7, loc="upper left")
b.set_xticks([5e5, 1e6, 2e6, 5e6], ["0.5 M", "1 M", "2 M", "5 M"])
b.minorticks_off()
b.set_xlim(4e5, 7e6)
b.set_ylim(2, 300)
for ax in (a, b):
    ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("docs/figures/fig3_scaling.png", dpi=200)
