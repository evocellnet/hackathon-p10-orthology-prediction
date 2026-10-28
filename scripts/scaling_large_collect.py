#!/usr/bin/env python3
"""Collect work/scale_large/<tag>/ -> results/scaling_large/{scaling.tsv,summary.md,scaling_large.png,fits.json}.
Fits: log10(y) ~ log10(proteins), OLS slope with t-interval (df=n-2); extrapolation clearly labelled (NOT measured)."""
import json, os, re, glob
import numpy as np
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT)
OUT = "results/scaling_large"
rows = []
for d in sorted(glob.glob("work/scale_large/*/")):
    tag = os.path.basename(d.rstrip("/")); p = d + "logs/"
    if not os.path.exists(p + "exit_code.txt") or open(p + "exit_code.txt").read().strip() != "0":
        continue
    m = json.load(open(d + "out/manifest.json"))["stages"]; t = open(p + "time_v.txt").read()
    g = lambda pat: re.search(pat, t).group(1)
    el = [float(x) for x in g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)").split(":")]
    wall = sum(v * 60 ** i for i, v in enumerate(reversed(el)))
    cpu = float(g(r"User time \(seconds\): (\S+)")) + float(g(r"System time \(seconds\): (\S+)"))
    loads = [float(x.replace(",", "")) for l in open(p + "monitor.log") for x in re.findall(r"load average: ([\d.]+)", l)]
    r = dict(tag=tag, preset=("accurate" if "accurate" in tag else "fast"), n_species=m["input"]["stats"]["n_species"], n_proteins=m["input"]["stats"]["n_genes"],
             wall_h=wall / 3600, cpu_h=cpu / 3600, max_rss_gb=int(g(r"Maximum resident set size \(kbytes\): (\d+)")) / 1048576,
             load1_mean=float(np.mean(loads)), load1_max=max(loads), threads=12,
             n_families=m["families"]["stats"]["n_families"], n_root_groups=m["merge"]["stats"].get("n_root_groups"))
    for s in ("input", "families", "scores", "merge", "output"):
        r[f"{s}_wall_h"] = m[s]["wall_s"] / 3600; r[f"{s}_rss_gb"] = m[s]["peak_rss_mb"] / 1024
    rows.append(r)
rows.sort(key=lambda r: (r["preset"], r["n_proteins"]))
cols = list(rows[0].keys())
with open(f"{OUT}/scaling.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows: f.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
fast = [r for r in rows if r["preset"] == "fast"]
fits = {}
def fit(key, pts):
    x = np.log10([r["n_proteins"] for r in pts]); y = np.log10([r[key] for r in pts]); n = len(x)
    if n < 3: return None
    res = stats.linregress(x, y); h = stats.t.ppf(0.975, n - 2) * res.stderr
    return dict(slope=res.slope, ci95=[res.slope - h, res.slope + h], intercept=res.intercept, n_points=n, r2=res.rvalue ** 2)
TARGETS = {"200M proteins (10k proteomes)": 2e8}
for key in ("wall_h", "cpu_h", "max_rss_gb", "families_wall_h", "scores_wall_h", "merge_wall_h"):
    f = fit(key, fast)
    if f:
        f["extrapolated_to_200M"] = 10 ** (f["intercept"] + f["slope"] * np.log10(2e8))
        f["extrapolated_to_200M_ci_slope_lo_hi"] = [10 ** (f["intercept"] + s * np.log10(2e8)) for s in f["ci95"]]  # crude: slope CI only
        fits[key] = f
json.dump(dict(fits=fits, note="fast preset only; extrapolation is NOT a measurement"), open(f"{OUT}/fits.json", "w"), indent=1)
L = ["# Large-scale demo (T3): Ancestra on UniProt proteome subsets", "",
     "12 threads, node shared with other agents (load logged per run). Measured values only in the first table.", "",
     "| run | preset | species | proteins | wall h | CPU-h | peak RSS GB | S0 h | S2 h | merge h | mean/max load1 | families | root groups |", "|" + "---|" * 13]
for r in rows:
    L.append(f"| {r['tag']} | {r['preset']} | {r['n_species']} | {r['n_proteins']:,} | {r['wall_h']:.2f} | {r['cpu_h']:.1f} | {r['max_rss_gb']:.1f} | {r['families_wall_h']:.2f} | {r['scores_wall_h']:.2f} | {r['merge_wall_h']:.2f} | {r['load1_mean']:.0f}/{r['load1_max']:.0f} | {r['n_families']:,} | {r['n_root_groups']:,} |")
L += ["", "## Log-log slopes vs number of proteins (fast preset, OLS, 95% t-interval)", "", "| quantity | slope | 95% CI | points | R2 |", "|---|---|---|---|---|"]
for k, f in fits.items():
    L.append(f"| {k} | {f['slope']:.2f} | [{f['ci95'][0]:.2f}, {f['ci95'][1]:.2f}] | {f['n_points']} | {f['r2']:.3f} |")
L += ["", "## EXTRAPOLATION (not measured) to 200 M proteins / ~10,000 proteomes", "",
      "Power-law extrapolation from the fits above over a ~40-400x range beyond the largest measured point; the interval reflects the slope CI only and ignores model misspecification (e.g. memory is not a power law once chunked). Treat as order of magnitude.", "",
      "| quantity | extrapolated | slope-CI range |", "|---|---|---|"]
for k, f in fits.items():
    lo, hi = sorted(f["extrapolated_to_200M_ci_slope_lo_hi"]); u = "GB" if "rss" in k else "h"
    L.append(f"| {k} | {f['extrapolated_to_200M']:,.0f} {u} | {lo:,.0f} - {hi:,.0f} {u} |")
if len(fast) >= 2:
    a, b = fast[-2], fast[-1]
    L += ["", f"## Local slopes between the two largest measured points ({a['tag']} -> {b['tag']}, {b['n_proteins']/a['n_proteins']:.2f}x proteins)", "",
          "Report these alongside the global fit: the cost is superlinear; the local exponent for S0 (families) is higher than the global one.", "", "| quantity | n_prev | n_last | ratio | local slope |", "|---|---|---|---|---|"]
    for k in ("wall_h", "cpu_h", "max_rss_gb", "families_wall_h", "scores_wall_h", "merge_wall_h", "output_wall_h"):
        rr = b[k] / a[k]; L.append(f"| {k} | {a[k]:.3g} | {b[k]:.3g} | {rr:.2f} | {np.log(rr)/np.log(b['n_proteins']/a['n_proteins']):.2f} |")
    L += ["", "Caveats: timings come from a shared node (load1 mean/max in the table); S0 (families) at n1000 was 3.6x the n500 time for 2.1x proteins, i.e. the global power-law slope understates the growth at the top end, so the extrapolation below is more likely an under- than an over-estimate for time. Peak RSS is the S2 (scores) stage; a 65 TB extrapolation is not meaningful, it only says S2 would need chunking/streaming at that scale."]
open(f"{OUT}/summary.md", "w").write("\n".join(L) + "\n")
# plot
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
xs = np.array([r["n_proteins"] for r in fast]); xx = np.logspace(np.log10(xs.min() * 0.8), np.log10(2e8), 50)
for a, key, lab, c in ((ax[0], "wall_h", "wall time (h)", "#2a78d6"), (ax[1], "cpu_h", "CPU time (h)", "#4a3aa7"), (ax[2], "max_rss_gb", "peak RSS (GB)", "#eb6834")):
    a.loglog(xs, [r[key] for r in fast], "o", color=c, label="measured (fast)")
    for r in rows:
        if r["preset"] != "fast": a.loglog(r["n_proteins"], r[key], "D", color="#1baf7a", label="measured (accurate)")
    f = fits.get(key)
    if f: a.loglog(xx, 10 ** (f["intercept"] + f["slope"] * np.log10(xx)), "--", color="gray", lw=1.5, label=f"extrapolation, slope {f['slope']:.2f}")
    a.axvline(2e8, color="lightgray", lw=1); a.set_xlabel("proteins"); a.set_ylabel(lab); a.grid(alpha=.3, which="both"); a.set_title(lab)
    h, l = a.get_legend_handles_labels(); d = dict(zip(l, h)); a.legend(d.values(), d.keys(), fontsize=7)
fig.suptitle("Ancestra at scale (12 threads, shared node); dashed = extrapolation, not measured", fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/scaling_large.png", dpi=150)
