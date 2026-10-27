#!/usr/bin/env python
"""L2-gpu: evaluate arm runs (OrthoBench dev via gpu_dev_eval; TreeFam-A dev via l2gpu_treefam_paired)
and write results/l2_gpu/<name>.json with a manifest (git commit, inputs, run configs, runtimes).

    python scripts/l2gpu_eval.py ob NAME REF_RUN RUN [RUN ...]
    python scripts/l2gpu_eval.py qfo NAME REF_RUN RUN [RUN ...]
"""
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))


def run_info(run):
    out = {"run": os.path.abspath(run)}
    m = os.path.join(run, "manifest.json")
    if os.path.exists(m):
        mm = json.load(open(m))
        out["config_scoring"] = mm.get("config", {}).get("scoring")
        out["merge_cfg_nondefault"] = {k: v for k, v in (mm.get("config", {}).get("merge") or {}).items()
                                       if k in ("link", "tree_ft_batch", "tree_ft_inproc", "tree_method")}
        out["git"] = mm.get("git")
        st = {}
        for s in ("scores", "rescue", "merge", "output"):
            p = os.path.join(run, "stages", s, "stage.json")
            if os.path.exists(p):
                sj = json.load(open(p))
                st[s] = {k: sj.get(k) for k in ("wall_s", "peak_rss_mb", "skipped")}
                if s == "merge":
                    st[s]["tree_s"] = sj.get("stats", {}).get("tree_s")
                if s == "rescue":
                    st[s]["stats"] = {k: v for k, v in sj.get("stats", {}).items() if not isinstance(v, (dict, list))}
        out["stages"] = st
    t = os.path.join(run, "time.txt")
    if os.path.exists(t):
        for line in open(t):
            line = line.strip()
            for key in ("User time (seconds)", "System time (seconds)", "Elapsed (wall clock) time",
                        "Maximum resident set size"):
                if line.startswith(key):
                    out.setdefault("time", {})[key] = line.split(": ", 1)[1]
    ox = os.path.join(run, "hogs.orthoxml")
    if os.path.exists(ox):
        out["hogs_md5"] = subprocess.run(["md5sum", ox], capture_output=True, text=True).stdout.split()[0]
    return out


def main():
    ds, name, ref, *runs = sys.argv[1:]
    t0 = time.time()
    if ds == "ob":
        from gpu_dev_eval import boot, one
        from orthobench_dev_pairs_eval import load_truth
        truth = load_truth()
        a = one(ref, truth)
        res = [{k: v for k, v in a.items() if not k.startswith("_")}]
        for r in runs:
            b = one(r, truth)
            d = {k: v for k, v in b.items() if not k.startswith("_")}
            d.update(boot(a, b))
            res.append(d)
            print(json.dumps(d))
    else:
        tmp = os.path.join(ROOT, "results/l2_gpu", f".{name}.tmp.json")
        cmd = [sys.executable, os.path.join(ROOT, "scripts/l2gpu_treefam_paired.py"), ref, *runs, "--json", tmp]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        res = json.load(open(tmp))
        os.remove(tmp)
        for d in res[1:]:
            print(json.dumps(d))
    git = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    man = {"script": "scripts/l2gpu_eval.py", "argv": sys.argv, "git_commit": git, "dataset": ds,
           "split": "dev only (OrthoBench 34 dev RefOGs / TreeFam-A 359 dev families)",
           "bootstrap": "paired, B=2000, seed 11", "eval_runtime_s": round(time.time() - t0, 1),
           "runs": [run_info(r) for r in [ref, *runs]]}
    out = os.path.join(ROOT, "results/l2_gpu", f"{name}.json")
    json.dump({"results": res, "manifest": man}, open(out, "w"), indent=1)
    print("->", out)


if __name__ == "__main__":
    main()
