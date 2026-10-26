#!/usr/bin/env python3
"""Run S0 family-partition configs on the OrthoBench 12 proteomes and score them on DEV RefOGs.

For each config name: runs `python -m ancestra.families cluster` under /usr/bin/time -v,
then scores families.tsv on the DEV RefOGs only, and writes
results/S0_families/<name>.json  (metrics + per-RefOG + family stats + runtime + manifest).

The TEST subset is scored only with --final <name> (once, for the chosen config); that
writes results/S0_families/<name>.TEST.json and appends a line to results/test_log.tsv.

Usage:
  python scripts/run_S0_grid.py --configs results/S0_families/configs.json [--only n1 n2] [--force]
  python scripts/run_S0_grid.py --final <name> --reason "..."
"""
import argparse, datetime, hashlib, json, os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra import families as F  # noqa: E402

INPUT_DIR = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/Input")
WORK = os.path.join(ROOT, "work/S0")
RES = os.path.join(ROOT, "results/S0_families")


def sha256(path, n=None):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def combined_fasta():
    out = os.path.join(WORK, "ob_all.fa")
    files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith(".fa"))
    if not os.path.exists(out):
        os.makedirs(WORK, exist_ok=True)
        with open(out + ".tmp", "w") as fo:
            for f in files:
                fo.write(open(os.path.join(INPUT_DIR, f)).read())
        os.replace(out + ".tmp", out)
    return out, {f: sha256(os.path.join(INPUT_DIR, f)) for f in files}


def git_commit():
    try:
        c = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "ancestra/families.py"],
                               capture_output=True, text=True).stdout.strip()
        return c + ("+dirty(families.py)" if dirty else "")
    except Exception:
        return "unknown"


def parse_time_v(text):
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", text)
    if m:
        parts = [float(x) for x in m.group(1).split(":")]
        s = 0.0
        for p in parts:
            s = s * 60 + p
        out["wall_s"] = s
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if m:
        out["peak_rss_gb"] = int(m.group(1)) / 1e6
    for key, pat in (("user_s", r"User time \(seconds\): (\S+)"), ("sys_s", r"System time \(seconds\): (\S+)")):
        m = re.search(pat, text)
        if m:
            out[key] = float(m.group(1))
    if "user_s" in out:
        out["cpu_h"] = (out["user_s"] + out.get("sys_s", 0)) / 3600
    return out


def run_one(name, cfg, fasta, input_sums, force=False):
    outdir = os.path.join(WORK, name)
    res_path = os.path.join(RES, f"{name}.json")
    if os.path.exists(res_path) and not force:
        print(f"[skip] {name} (exists)")
        return json.load(open(res_path))
    os.makedirs(outdir, exist_ok=True)
    cfg = dict(cfg, name=name)
    cfg_path = os.path.join(outdir, "config.json")
    json.dump(cfg, open(cfg_path, "w"), indent=1)
    started = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    load = open("/proc/loadavg").read().split()[:3]
    cmd = ["/usr/bin/time", "-v", sys.executable, "-m", "ancestra.families", "cluster",
           "--fasta", fasta, "--outdir", outdir, "--config", cfg_path]
    print(f"[run] {name}: {' '.join(cmd)}", flush=True)
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    open(os.path.join(outdir, "time_v.txt"), "w").write(p.stderr)
    if p.returncode != 0:
        print(p.stderr[-3000:])
        raise SystemExit(f"{name} failed")
    runtime = parse_time_v(p.stderr)
    stats = json.load(open(os.path.join(outdir, "stats.json")))
    g2f = F.read_families_tsv(os.path.join(outdir, "families.tsv"))
    split = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
    refdir = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs")
    dev = F.evaluate_refogs(g2f, F.read_refogs(refdir, split["dev"], remove_uncertain=True))
    dev_full = F.evaluate_refogs(g2f, F.read_refogs(refdir, split["dev"], remove_uncertain=False))
    out = {
        "name": name,
        "subset": "dev",
        "dev": dev,
        "dev_keep_uncertain_aggregate": dev_full["aggregate"],
        "family_stats": {k: v for k, v in stats.items() if k not in ("config", "timings")},
        "runtime": dict(runtime, stage_timings=stats["timings"], node_loadavg_at_start=load),
        "manifest": {
            "config": stats["config"], "config_hash": stats["config_hash"], "commit": git_commit(),
            "families_py_sha256": sha256(os.path.join(ROOT, "ancestra/families.py")),
            "inputs": {"orthobench_input_sha256": input_sums,
                       "split": "splits/orthobench_split.json",
                       "split_sha256": sha256(os.path.join(ROOT, "splits/orthobench_split.json"))},
            "families_tsv": os.path.relpath(os.path.join(outdir, "families.tsv"), ROOT),
            "families_tsv_sha256": sha256(os.path.join(outdir, "families.tsv")),
            "command": " ".join(cmd[2:]), "started": started, "host_threads": os.cpu_count(),
        },
    }
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(res_path, "w"), indent=1)
    a = dev["aggregate"]; fs = out["family_stats"]
    print(f"[done] {name}: compl={a['completeness']:.3f} fission={a['mean_n_fam']:.2f} "
          f"recall={a['pair_recall']:.3f} fusion={a['fusion_rate']:.3f} sizeratio_med={a['median_size_ratio']:.2f} "
          f"sumg2/N2={fs['sum_g2_over_N2']:.4f} max={fs['max_family']} wall={runtime.get('wall_s')}", flush=True)
    return out


def final_test(name, reason):
    res = json.load(open(os.path.join(RES, f"{name}.json")))
    fam = os.path.join(ROOT, res["manifest"]["families_tsv"])
    assert sha256(fam) == res["manifest"]["families_tsv_sha256"], "families.tsv changed since dev scoring"
    split = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
    refdir = os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs")
    g2f = F.read_families_tsv(fam)
    test = F.evaluate_refogs(g2f, F.read_refogs(refdir, split["test"], remove_uncertain=True))
    test_full = F.evaluate_refogs(g2f, F.read_refogs(refdir, split["test"], remove_uncertain=False))
    out = {"name": name, "subset": "test", "test": test, "test_keep_uncertain_aggregate": test_full["aggregate"],
           "manifest": res["manifest"], "reason": reason,
           "scored_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds")}
    p = os.path.join(RES, f"{name}.TEST.json")
    if os.path.exists(p):
        raise SystemExit(f"{p} exists: test set already scored once for {name}")
    json.dump(out, open(p, "w"), indent=1)
    log = os.path.join(ROOT, "results/test_log.tsv")
    new = not os.path.exists(log)
    with open(log, "a") as fh:
        if new:
            fh.write("time\tsubtask\tconfig_name\tconfig_hash\tbenchmark\treason\n")
        fh.write(f"{out['scored_at']}\tS0\t{name}\t{res['manifest']['config_hash']}\tOrthoBench-test-RefOGs\t{reason}\n")
    print(json.dumps(test["aggregate"], indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default=os.path.join(RES, "configs.json"))
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--final")
    ap.add_argument("--reason", default="S0 final chosen config (pre-stated rule)")
    a = ap.parse_args()
    if a.final:
        final_test(a.final, a.reason)
        sys.exit()
    fasta, sums = combined_fasta()
    cfgs = json.load(open(a.configs))
    for name, cfg in cfgs.items():
        if name.startswith("_") or (a.only and name not in a.only):
            continue
        run_one(name, cfg, fasta, sums, force=a.force)
