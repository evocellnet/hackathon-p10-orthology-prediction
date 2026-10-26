#!/usr/bin/env python3
"""Run ``ancestra.pipeline.run_pipeline`` with per-stage CPU / load accounting.

Thin wrapper (does not modify the pipeline): the pipeline's own stage cache
(``OUT/stages/<stage>/DONE``) makes every call resumable; this driver adds

* per-stage wall, user+sys CPU of this process and of finished children
  (``getrusage`` deltas around the pipeline's "stage X: running/done" log lines),
  node load average at stage start/end, peak RSS (own VmHWM from stage.json and
  max child RSS) -> ``OUT/stage_accounting.json`` (appended per invocation);
* ``--stop-after STAGE``: stop cleanly after STAGE has finished (or was found
  up to date), e.g. ``--stop-after scores`` to check cost before the merge.

    python scripts/qfo_pipeline_driver.py --proteomes DIR --tree T.nwk --out OUT \
        [--threads 12] [--stop-after families|scores|merge|output] [--config cfg.json]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import resource
import shutil
import sys
import time
from dataclasses import asdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

STAGES = ("input", "families", "scores", "merge", "output")


class StopAfter(Exception):
    pass


def _cpu():
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {"self": s.ru_utime + s.ru_stime, "children": c.ru_utime + c.ru_stime,
            "children_maxrss_mb": c.ru_maxrss / 1024}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--stop-after", choices=STAGES)
    ap.add_argument("--config", help="JSON {'scoring': {...}, 'merge': {...}} (default: pipeline defaults)")
    a = ap.parse_args(argv)

    from ancestra.merge import MergeConfig
    from ancestra.pipeline import run_pipeline
    from ancestra.scoring import ScoringConfig

    sc, mc = ScoringConfig(), MergeConfig()
    if a.config:
        c = json.load(open(a.config))
        for k, v in c.get("scoring", {}).items():
            setattr(sc, k, v)
        for k, v in c.get("merge", {}).items():
            setattr(mc, k, v)

    os.makedirs(a.out, exist_ok=True)
    # stale temp dirs of killed runs (the pipeline names them <stage>.tmp<pid>)
    sdir = os.path.join(a.out, "stages")
    if os.path.isdir(sdir):
        for f in os.listdir(sdir):
            if re.match(r"^\w+\.tmp\d+$", f):
                shutil.rmtree(os.path.join(sdir, f), ignore_errors=True)

    acc_path = os.path.join(a.out, "stage_accounting.json")
    runs = json.load(open(acc_path)) if os.path.exists(acc_path) else []
    run = {"started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "argv": sys.argv, "threads": a.threads,
           "stop_after": a.stop_after, "stages": {}}
    runs.append(run)
    cur: dict = {}

    def save():
        with open(acc_path + ".tmp", "w") as fh:
            json.dump(runs, fh, indent=1)
        os.replace(acc_path + ".tmp", acc_path)

    def log(msg: str):
        print(f"{time.strftime('%H:%M:%S')} {msg}", file=sys.stderr, flush=True)
        m = re.match(r"\[ancestra\] stage (\w+): (running|done|up to date)", msg)
        if not m:
            return
        name, what = m.groups()
        if what == "running":
            cur[name] = {"t0": time.time(), "cpu0": _cpu(), "load_start": os.getloadavg()}
            return
        e = {"skipped": what == "up to date"}
        if name in cur:
            s = cur.pop(name)
            c1 = _cpu()
            e.update(wall_s=round(time.time() - s["t0"], 2),
                     cpu_self_s=round(c1["self"] - s["cpu0"]["self"], 1),
                     cpu_children_s=round(c1["children"] - s["cpu0"]["children"], 1),
                     children_maxrss_mb_so_far=round(c1["children_maxrss_mb"], 1),
                     load_start=[round(x, 1) for x in s["load_start"]],
                     load_end=[round(x, 1) for x in os.getloadavg()])
            e["cpu_total_s"] = round(e["cpu_self_s"] + e["cpu_children_s"], 1)
        sj = os.path.join(a.out, "stages", name, "stage.json")
        if os.path.exists(sj):
            info = json.load(open(sj))
            e["pipeline_wall_s"] = info.get("wall_s")
            e["peak_rss_mb"] = info.get("peak_rss_mb")
            e["stats"] = info.get("stats")
        run["stages"][name] = e
        save()
        if a.stop_after == name:
            raise StopAfter(name)

    t0 = time.time()
    try:
        m = run_pipeline(a.proteomes, a.tree, a.out, threads=a.threads, scoring=sc, merge=mc, log=log)
        run["total_wall_s"] = m.get("total_wall_s")
        run["finished"] = True
    except StopAfter as s:
        run["finished"] = False
        run["stopped_after"] = str(s)
        log(f"[driver] stopped after stage {s} (resumable)")
    run["driver_wall_s"] = round(time.time() - t0, 1)
    c = _cpu()
    run["cpu_self_s"], run["cpu_children_s"] = round(c["self"], 1), round(c["children"], 1)
    save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
