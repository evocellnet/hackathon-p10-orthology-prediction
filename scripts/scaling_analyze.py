#!/usr/bin/env python3
"""Aggregate results/scaling/runs/*/runtime.json -> scaling.tsv, summary.md, scaling.png (DESIGN §5).
Fit: log10(CPU-h) ~ log10(x), x = total proteins or species, per method, uncensored points only.
CI: OLS t-interval (df=n-2) and percentile bootstrap over points (B=5000, seed 11; resamples with <3 distinct x skipped)."""
import glob, json, os, random, sys
import numpy as np
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
OUT = "results/scaling"
LAB = {"anc": "Ancestra (sequence-only)", "anc_go": "Ancestra (+gene order)", "of": "OrthoFinder 3.2", "sp": "SonicParanoid2 (default)"}
FITM = ["anc", "of", "sp"]


def load():
    rows = []
    for f in sorted(glob.glob(f"{OUT}/runs/*/runtime.json")):
        r = json.load(open(f))
        r["stages"] = {}
        sa = os.path.join(os.path.dirname(f), "stage_accounting.json")
        if os.path.exists(sa):
            acc = json.load(open(sa))
            # merge invocations (a fresh run has one)
            for inv in acc:
                for st, v in inv.get("stages", {}).items():
                    if not v.get("skipped"):
                        r["stages"][st] = v
        rows.append(r)
    rows.sort(key=lambda r: (list(LAB).index(r["kind"]), r["n_species"]))
    return rows


def ols(x, y):
    lx, ly = np.log10(x), np.log10(y)
    n = len(lx)
    res = stats.linregress(lx, ly)
    out = {"n": n, "slope": res.slope, "intercept": res.intercept, "r2": res.rvalue ** 2}
    if n >= 3:
        t = stats.t.ppf(0.975, n - 2)
        out["ols_lo"], out["ols_hi"] = res.slope - t * res.stderr, res.slope + t * res.stderr
        rng = random.Random(11)
        sl = []
        for _ in range(5000):
            idx = [rng.randrange(n) for _ in range(n)]
            if len({lx[i] for i in idx}) < 3:
                continue
            sl.append(np.polyfit(lx[idx], ly[idx], 1)[0])
        out["boot_lo"], out["boot_hi"] = np.percentile(sl, [2.5, 97.5])
        out["boot_n"] = len(sl)
    return out


def hms(h):
    return f"{h:.2f}"


def main():
    rows = load()
    if not rows:
        sys.exit("no runs")
    # ---- tsv
    stages = ["input", "families", "scores", "merge", "output"]
    cols = ["method", "label", "n_species", "n_proteins", "status", "wall_h", "cpu_h", "user_h", "sys_h", "max_rss_gb", "cap_h",
            "load1_mean", "load1_max", "start", "end"] + [f"{s}_{m}" for s in stages for m in ("wall_s", "cpu_s", "rss_mb")]
    with open(f"{OUT}/scaling.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            cens = r["censored"]
            st = f"> {r['cap_s']/3600:g} h (censored)" if cens else ("ok" if r.get("ok") else f"failed rc={r['returncode']}")
            L = r.get("node_load_1min", {})
            v = [LAB[r["kind"]], r["label"], r["n_species"], r["n_proteins"], st, f"{r.get('wall_s',0)/3600:.3f}", f"{r.get('cpu_h',float('nan')):.3f}",
                 f"{r.get('user_s',0)/3600:.3f}", f"{r.get('sys_s',0)/3600:.3f}", f"{r.get('max_rss_gb',float('nan')):.2f}", f"{r['cap_s']/3600:g}",
                 f"{L.get('load1_mean',float('nan')):.1f}", f"{L.get('load1_max',float('nan')):.1f}", r["start"], r["end"]]
            for s in stages:
                d = r["stages"].get(s, {})
                v += [d.get("wall_s", ""), d.get("cpu_total_s", ""), d.get("peak_rss_mb", "")]
            fh.write("\t".join(str(x) for x in v) + "\n")
    # ---- fits
    fits = {}
    for m in FITM:
        pts = [r for r in rows if r["kind"] == m and r.get("ok") and not r["censored"]]
        for xn, xk in (("proteins", "n_proteins"), ("species", "n_species")):
            if len(pts) >= 2:
                fits[(m, xn)] = ols([r[xk] for r in pts], [r["cpu_h"] for r in pts])
    json.dump({f"{m}|{x}": v for (m, x), v in fits.items()}, open(f"{OUT}/fits.json", "w"), indent=1, default=float)
    # ---- plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"anc": "#2a78d6", "anc_go": "#2a78d6", "of": "#eb6834", "sp": "#1baf7a"}
    mk = {"anc": "o", "anc_go": "D", "of": "s", "sp": "^"}
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for ax, (xn, xk, xl) in zip(axs, (("proteins", "n_proteins", "total proteins"), ("species", "n_species", "species"))):
        for k in LAB:
            rs = [r for r in rows if r["kind"] == k and "cpu_h" in r]
            ok = [r for r in rs if r.get("ok") and not r["censored"]]
            ce = [r for r in rs if r["censored"]]
            if ok:
                ax.plot([r[xk] for r in ok], [r["cpu_h"] for r in ok], mk[k], color=col[k], ms=7, mfc=col[k] if k != "anc_go" else "none", label=LAB[k], ls="none")
            for r in ce:
                ax.plot(r[xk], r["cpu_h"], mk[k], color=col[k], mfc="none", ms=7)
                ax.annotate("", xy=(r[xk], r["cpu_h"] * 2.2), xytext=(r[xk], r["cpu_h"]), arrowprops=dict(arrowstyle="->", color=col[k], lw=1.2))
            if k in FITM and (k, xn) in fits:
                f = fits[(k, xn)]
                xs = np.array([min(r[xk] for r in ok), max(max(r[xk] for r in rs), 1)])
                ax.plot(xs, 10 ** (f["intercept"] + f["slope"] * np.log10(xs)), "-", color=col[k], lw=1.4, alpha=.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(xl); ax.set_ylabel("CPU time, user+sys (h)")
        ax.grid(True, which="both", alpha=.2, lw=.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axs[0].legend(frameon=False, fontsize=8, loc="upper left")
    axs[1].set_xticks([8, 16, 32, 64]); axs[1].set_xticklabels(["8", "16", "32", "64"]); axs[1].minorticks_off()
    fig.suptitle("Runtime scaling on nested QfO 2022 subsets (8 threads, shared node); open markers with arrow = censored (> cap), plotted at CPU used so far", fontsize=9)
    fig.savefig(f"{OUT}/scaling.png", dpi=160)
    # ---- summary
    sizes = sorted({r["n_species"] for r in rows})
    npro = {r["n_species"]: r["n_proteins"] for r in rows}
    L = ["# Runtime scaling (DESIGN §5)", "",
         "8 threads per run, one run at a time, shared 24-core node (load logged per run in `runs/*/load.log`; wall times from a busy node are upper bounds, CPU-h = user+sys from `/usr/bin/time -v` is the primary metric). Fresh work dir per run. Censored runs: `> cap`.", "",
         "Subsets: nested, domain-stratified (random.Random(11)), `subsets/sub{n}.{txt,nwk,json}`; guide tree = NCBI tree pruned to the subset.", "",
         "| method | " + " | ".join(f"{n} sp ({npro[n]/1000:.0f}k)" for n in sizes) + " |", "|---|" + "---|" * len(sizes)]
    for k in LAB:
        cells = []
        for n in sizes:
            r = next((r for r in rows if r["kind"] == k and r["n_species"] == n), None)
            if r is None:
                cells.append("not run")
            elif r["censored"]:
                cells.append(f"> {r['cap_s']/3600:g} h wall ({r.get('cpu_h',float('nan')):.1f} CPU-h used, {r.get('max_rss_gb',float('nan')):.1f} GB)")
            elif not r.get("ok"):
                cells.append(f"FAILED rc={r['returncode']}")
            else:
                cells.append(f"{r['wall_h']:.2f} h / {r['cpu_h']:.2f} CPU-h / {r['max_rss_gb']:.1f} GB")
        L.append(f"| {LAB[k]} | " + " | ".join(cells) + " |")
    L += ["", "Cells: wall / CPU-h / peak RSS (max over processes as reported by `time -v`).", "",
          "## Log-log slope of CPU-h (uncensored points)", "",
          "| method | x | n pts | slope | OLS 95% CI | bootstrap 95% CI | R² |", "|---|---|---|---|---|---|---|"]
    for m in FITM:
        for xn in ("proteins", "species"):
            f = fits.get((m, xn))
            if not f:
                continue
            ci = f"[{f['ols_lo']:.2f}, {f['ols_hi']:.2f}]" if "ols_lo" in f else "n<3"
            bi = f"[{f['boot_lo']:.2f}, {f['boot_hi']:.2f}]" if "boot_lo" in f else "n<3"
            L.append(f"| {LAB[m]} | {xn} | {f['n']} | {f['slope']:.2f} | {ci} | {bi} | {f['r2']:.3f} |")
    L += ["", "Bootstrap: percentile, B=5000, seed 11, resampling the points (resamples with <3 distinct x skipped). With 3-5 points both intervals are rough.", ""]
    L += ["## Per-stage cost, Ancestra (CPU s / wall s / peak RSS MB)", "", "| run | " + " | ".join(stages) + " |", "|---|" + "---|" * len(stages)]
    for r in rows:
        if r["kind"] in ("anc", "anc_go"):
            L.append(f"| {r['label']} | " + " | ".join(
                f"{r['stages'][s].get('cpu_total_s',0):.0f} / {r['stages'][s].get('wall_s',0):.0f} / {r['stages'][s].get('peak_rss_mb',0):.0f}" if s in r["stages"] else "-" for s in stages) + " |")
    open(f"{OUT}/summary.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


main()
