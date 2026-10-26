#!/usr/bin/env python3
"""Parse `/usr/bin/time -v` output + node-load log of run_competitor_timed.sh into runtime.json.
Usage: make_runtime_json.py <logs dir> <out.json> --threads N [--note TEXT] [--include <logsdir2> ...]
(--include adds a further phase, e.g. several sequential steps, summed; wall times add, RSS = max)"""
import argparse, json, re, datetime as dt
def parse_time(fn):
    t = open(fn).read()
    def g(p): return re.search(p, t).group(1)
    w = g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)").split(":")
    w = [float(x) for x in w]; wall = w[-1] + 60 * w[-2] + (3600 * w[-3] if len(w) == 3 else 0)
    return dict(wall_s=wall, user_s=float(g(r"User time \(seconds\): ([\d.]+)")), sys_s=float(g(r"System time \(seconds\): ([\d.]+)")),
                max_rss_gb=int(g(r"Maximum resident set size \(kbytes\): (\d+)")) / 1048576,
                exit_status=int(g(r"Exit status: (\d+)")),
                major_faults=int(g(r"Major \(requiring I/O\) page faults: (\d+)")))
def parse_load(fn):
    v = []
    for l in open(fn):
        m = re.search(r"load average: ([\d.]+)", l)
        if m: v.append(float(m.group(1)))
    return dict(n=len(v), load1_mean=sum(v) / len(v), load1_max=max(v), load1_min=min(v)) if v else {}
p = argparse.ArgumentParser(); p.add_argument("logs"); p.add_argument("out"); p.add_argument("--threads", type=int, required=True)
p.add_argument("--note", default=""); p.add_argument("--include", nargs="*", default=[]); a = p.parse_args()
ph = {}
tot = dict(wall_s=0, user_s=0, sys_s=0, max_rss_gb=0)
for d in [a.logs] + a.include:
    r = parse_time(d + "/time.txt"); r["node_load_1min"] = parse_load(d + "/load.log")
    r["start"] = open(d + "/start.txt").read().strip(); r["end"] = open(d + "/end.txt").read().strip()
    ph[d.rstrip("/").split("/")[-2] + "/" + d.rstrip("/").split("/")[-1]] = r
    for k in ("wall_s", "user_s", "sys_s"): tot[k] += r[k]
    tot["max_rss_gb"] = max(tot["max_rss_gb"], r["max_rss_gb"])
tot["cpu_s"] = tot["user_s"] + tot["sys_s"]; tot["cpu_h"] = tot["cpu_s"] / 3600; tot["wall_h"] = tot["wall_s"] / 3600
tot["threads_requested"] = a.threads; tot["node_cores"] = 24; tot["node_shared_busy"] = True
tot["note"] = a.note; tot["phases"] = ph
json.dump(tot, open(a.out, "w"), indent=1)
print(json.dumps({k: v for k, v in tot.items() if k != "phases"}, indent=1))
