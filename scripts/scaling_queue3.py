#!/usr/bin/env python3
"""Sequential runner for the runtime-scaling benchmark (DESIGN §5). ONE run at a time, 8 threads.

Reads work/scaling/queue.txt (re-read before every job, so jobs can be appended/edited while it runs):
    label|kind|n_species|cap_h        kind in {anc, anc_go, of, sp}
A job is done when results/scaling/runs/<label>/runtime.json exists. Fresh work dir per run
(work/scaling/runs/<label>); refuses to reuse one. Each run: /usr/bin/time -v around `timeout`,
uptime logged every 60 s (load.log), session-wide cleanup after exit. Censored = timeout (exit 124).
Per-run logs: results/scaling/runs/<label>/{time.txt,load.log,cmd_tail.txt,runtime.json,[manifest.json,stage_accounting.json]}
"""
import json, os, re, shutil, signal, subprocess, sys, time, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
QUEUE = "work/scaling/queue3.txt"
THREADS = 8
GO_DIR = "data/inputs/qfo2022/gene_order"


SNAP = "work/scaling/code_snapshot"   # frozen copy of ancestra/ + driver (the code under development changes during the benchmark)


def cmd_for(kind, n, wd, extra):
    P = f"work/scaling/inputs/sub{n}"
    if kind in ("anc", "anc_go"):
        c = ["python", f"{SNAP}/scripts/qfo_pipeline_driver.py", "--proteomes", P, "--tree", f"results/scaling/subsets/sub{n}.nwk",
             "--out", wd, "--threads", str(THREADS), "--gene-order", GO_DIR if kind == "anc_go" else "none"]
        if "preset" in extra:
            c += ["--preset", extra["preset"]]
        if "reuse" in extra:
            c += ["--reuse-stages-from", f"work/scaling/runs/{extra['reuse']}"]
        return c
    if kind == "fo":
        return ["bash", "scripts/scaling_fastoma.sh", str(n), wd]
    if kind == "of":
        return ["orthofinder", "-f", P, "-t", str(THREADS), "-a", str(THREADS), "-mog", "-o", wd]
    if kind == "sp":
        return ["sonicparanoid", "-i", P, "-o", wd, "-p", f"scal{n}", "-t", str(THREADS), "-m", "default"]
    raise SystemExit(f"bad kind {kind}")


def parse_time(path):
    t = open(path).read()
    def g(pat):
        m = re.search(pat, t)
        return m.group(1) if m else None
    el = g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)")
    parts = [float(x) for x in el.split(":")]
    wall = sum(v * 60 ** i for i, v in enumerate(reversed(parts)))
    return {"wall_s": wall, "user_s": float(g(r"User time \(seconds\): (\S+)")), "sys_s": float(g(r"System time \(seconds\): (\S+)")),
            "max_rss_gb": int(g(r"Maximum resident set size \(kbytes\): (\d+)")) / 1048576,
            "exit_status": int(g(r"Exit status: (\d+)")), "major_faults": int(g(r"Major .*: (\d+)")),
            "cpu_pct": g(r"Percent of CPU this job got: (\S+)")}


def load_stats(path):
    v = []
    for l in open(path):
        m = re.search(r"load average: ([\d.]+)", l)
        if m:
            v.append(float(m.group(1)))
    return {"n": len(v), "load1_mean": sum(v) / len(v), "load1_max": max(v), "load1_min": min(v)} if v else {}


def next_job():
    if not os.path.exists(QUEUE):
        return None
    for l in open(QUEUE):
        l = l.split("#")[0].strip()
        if not l:
            continue
        p = l.split("|")
        label, kind, n, cap = p[:4]
        extra = dict(kv.split("=") for kv in p[4].split(",")) if len(p) > 4 and p[4] else {}
        if not os.path.exists(f"results/scaling/runs/{label}/runtime.json"):
            return label, kind, int(n), float(cap), extra
    return None


def run(label, kind, n, cap_h, extra):
    rd, wd = f"results/scaling/runs/{label}", f"work/scaling/runs/{label}"
    if os.path.exists(wd):
        raise SystemExit(f"work dir {wd} exists; fresh dir required (no cache reuse)")
    os.makedirs(rd, exist_ok=True)
    os.makedirs(os.path.dirname(wd), exist_ok=True)
    if kind in ("anc", "anc_go"):
        os.makedirs(wd)
    cmd = cmd_for(kind, n, wd, extra)
    cap_s = int(cap_h * 3600)
    full = ["/usr/bin/time", "-v", "-o", f"{rd}/time.txt", "timeout", "-s", "TERM", "-k", "120", str(cap_s)] + cmd
    cmdlog = f"work/scaling/runs/{label}.cmd.log"
    start = datetime.datetime.now().astimezone()
    print(f"[queue] {start.isoformat()} START {label}: {' '.join(cmd)} cap={cap_h}h", flush=True)
    open(f"{rd}/cmd.txt", "w").write(" ".join(cmd) + f"\ncap_s={cap_s}\nstart={start.isoformat()}\nload_at_start={open('/proc/loadavg').read()}")
    lf = open(f"{rd}/load.log", "w")
    with open(cmdlog, "w") as lg:
        p = subprocess.Popen(full, stdout=lg, stderr=subprocess.STDOUT, start_new_session=True)
        nxt = 0
        while p.poll() is None:
            if time.time() >= nxt:
                lf.write(f"{datetime.datetime.now().astimezone().isoformat()} {subprocess.run(['uptime'], capture_output=True, text=True).stdout.strip()}\n"); lf.flush()
                nxt = time.time() + 60
            time.sleep(2)
    end = datetime.datetime.now().astimezone()
    lf.close()
    try:                                  # leftover descendants of a killed/timed-out run
        os.killpg(p.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    rc = p.returncode
    rt = {"label": label, "kind": kind, "n_species": n, "n_proteins": json.load(open(f"results/scaling/subsets/sub{n}.json"))["n_proteins"],
          "threads": THREADS, "cap_s": cap_s, "start": start.isoformat(), "end": end.isoformat(), "node_cores": os.cpu_count(),
          "returncode": rc, "preset": extra.get("preset", "accurate(default)"), "reuse_from": extra.get("reuse"), "code": "snapshot" if kind.startswith("anc") else None, "censored": rc == 124 or rc == 137 or rc == -9 or rc == 143}
    try:
        rt.update(parse_time(f"{rd}/time.txt"))
        rt["cpu_s"] = rt["user_s"] + rt["sys_s"]; rt["cpu_h"] = rt["cpu_s"] / 3600; rt["wall_h"] = rt["wall_s"] / 3600
    except Exception as e:
        rt["time_parse_error"] = repr(e)
    rt["node_load_1min"] = load_stats(f"{rd}/load.log")
    if rc == 0:
        rt["ok"] = True
    for f in ("manifest.json", "stage_accounting.json"):
        if os.path.exists(f"{wd}/{f}"):
            shutil.copy(f"{wd}/{f}", f"{rd}/{f}")
    os.system(f"tail -n 40 {cmdlog} > {rd}/cmd_tail.txt")
    json.dump(rt, open(f"{rd}/runtime.json", "w"), indent=1)
    print(f"[queue] {end.isoformat()} END {label} rc={rc} wall={rt.get('wall_s')} cpu_h={rt.get('cpu_h')} rss={rt.get('max_rss_gb')}", flush=True)


def main():
    idle = 0
    while True:
        j = next_job()
        if j is None:
            idle += 1
            if idle > 30:      # 30 min without work
                print("[queue] queue empty, exiting", flush=True)
                return
            time.sleep(60)
            continue
        idle = 0
        run(*j)


main()
