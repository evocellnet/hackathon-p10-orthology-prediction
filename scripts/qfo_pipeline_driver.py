#!/usr/bin/env python3
"""Run ``ancestra.pipeline.run_pipeline`` with per-stage accounting and code-aware caching.

Thin wrapper (does not modify the pipeline). The pipeline's own stage cache
(``OUT/stages/<stage>/DONE``) makes every call resumable. This driver adds:

* **Per-stage accounting** -> ``OUT/stage_accounting.json`` (appended per invocation):
  wall time; user+sys CPU of this process and of finished children (``getrusage``
  deltas around the pipeline's "stage X: running/done" log lines); node load at
  stage start and end; peak RSS (own VmHWM from stage.json, max child RSS).
* **Code-aware invalidation.** The pipeline's stage keys hash configs and inputs but
  not code, so a cached stage survives e.g. a fix in ``ancestra/families.py``. The
  driver therefore keeps ``OUT/stages/code_keys.json``: for each finished stage, a
  sha256 over the source files that stage runs (STAGE_CODE; merge and output: all of
  ancestra/*.py), its *effective* config
  (for families: ``asdict(FamilyConfig(**DEFAULT_FAMILY_CONFIG))`` without threads),
  and the upstream stage's code key. Before running, any cached stage whose recorded
  key is missing or differs is deleted together with every downstream stage, so the
  pipeline recomputes it. ``pipeline.py`` is deliberately not hashed: its stage
  structure is covered by the pipeline's own keys, and hashing it would invalidate
  S0/S2 on every edit of the orchestration code.
* **Provenance.** ``OUT/manifest.json`` (written by the pipeline) gets a ``driver``
  block with the git commit and ``ancestra/`` dirty flag, per-file sha256, the stage
  code keys, the families config and its hash, and the gene-order directory.
* ``--gene-order DIR`` (S5 tables; enables the synteny merge) or ``--gene-order none``.
* ``--preset accurate|fast|custom`` passed to the pipeline (part of the merge code key).
* ``--reuse-stages-from RUN``: hard-link (``cp -al``) the cached input/families/scores
  stages of another run dir that lack here, e.g. the same S0+S2 merged with and
  without gene order. Code keys are checked as usual.
* ``--adopt-unrecorded``: accept cached stages that have no recorded code key (runs made
  before this driver existed) instead of recomputing them. For timing dry runs only; it is
  logged in stage_accounting.json and the manifest.
* ``--family-rescue MODE [--family-rescue-guard] [--structure-3di F]``: pipeline S0 rescue stage
  (code-keyed as stage ``rescue``: rescue.py + scoring.py + the rescue options).
* ``--transplant-from RUN``: guide-tree ablations without recomputing S0/S2 (see ``transplant``).
* ``--stop-after STAGE``: stop cleanly after STAGE has finished or was up to date.
* A lock file (``OUT/.driver.lock``) prevents two drivers on one run dir.

    python scripts/qfo_pipeline_driver.py --proteomes DIR --tree T.nwk --out OUT \
        [--threads 12] [--gene-order DIR|none] [--reuse-stages-from RUN] \
        [--stop-after families|scores|merge|output] [--config cfg.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time
from dataclasses import asdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

STAGES = ("input", "families", "scores", "merge", "output")
# source files each stage executes (relative to ancestra/); upstream keys are chained
# merge/output are cheap (minutes) and import many modules lazily (localtree, synteny, genome(s),
# merge_numba, ...): they are keyed on *every* ancestra/*.py so that no dependency is missed.
ALL_PY = sorted(f for f in os.listdir(os.path.join(ROOT, "ancestra")) if f.endswith(".py"))
STAGE_CODE = {
    "input": ["io.py", "newick.py"],
    "families": ["families.py", "io.py"],
    "scores": ["scoring.py"],
    "merge": ALL_PY,
    "output": ALL_PY,
}


RESCUE_CODE = ["rescue.py", "scoring.py"]


def stages_for(rescue: dict | None) -> tuple:
    return STAGES if not rescue else ("input", "families", "scores", "rescue", "merge", "output")


class StopAfter(Exception):
    pass


def _cpu():
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {"self": s.ru_utime + s.ru_stime, "children": c.ru_utime + c.ru_stime,
            "children_maxrss_mb": c.ru_maxrss / 1024}


def sha256_file(p: str) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def _h(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def git_info() -> dict:
    def g(*args):
        return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True).stdout.strip()
    try:
        return {"commit": g("rev-parse", "HEAD"),
                "ancestra_dirty": g("status", "--porcelain", "ancestra") != ""}
    except Exception:
        return {"commit": None}


def code_keys(family_cfg: dict, gene_order: str | None, preset: str | None = None,
              rescue: dict | None = None) -> tuple[dict, dict]:
    files = {f: sha256_file(os.path.join(ROOT, "ancestra", f))
             for fs in list(STAGE_CODE.values()) + [RESCUE_CODE] for f in fs
             if os.path.exists(os.path.join(ROOT, "ancestra", f))}
    keys, up = {}, None
    for st in stages_for(rescue):
        if st == "rescue":
            keys[st] = _h({"up": up, "files": {f: files.get(f) for f in RESCUE_CODE}, "rescue": rescue})
            up = keys[st]
            continue
        extra = {"families_config": family_cfg} if st == "families" else {}
        if st == "merge":
            extra = {"gene_order": gene_order}
            if preset is not None:
                extra["preset"] = preset
        keys[st] = _h({"up": up, "files": {f: files.get(f) for f in STAGE_CODE[st]}, **extra})
        up = keys[st]
    return keys, files


def invalidate(out: str, keys: dict, log, adopt: bool = False) -> tuple[list, list]:
    sdir = os.path.join(out, "stages")
    kp = os.path.join(sdir, "code_keys.json")
    rec = json.load(open(kp)) if os.path.exists(kp) else {}
    removed, adopted = [], []
    stale = False
    for st in keys:
        d = os.path.join(sdir, st)
        if adopt and not stale and st not in rec and os.path.exists(os.path.join(d, "DONE")):
            rec[st] = keys[st]      # explicit dry-run override: accept a pre-driver cached stage
            adopted.append(st)
            continue
        if stale or (os.path.exists(os.path.join(d, "DONE")) and rec.get(st) != keys[st]):
            stale = True
            if os.path.exists(d):
                shutil.rmtree(d)
                removed.append(st)
            rec.pop(st, None)
    if removed:
        log(f"[driver] code/config changed or unrecorded -> recomputing stages {removed}")
    if adopted:
        log(f"[driver] ADOPTED unrecorded cached stages {adopted} (--adopt-unrecorded; dry runs only)")
    os.makedirs(sdir, exist_ok=True)
    with open(kp + ".tmp", "w") as fh:
        json.dump(rec, fh, indent=1)
    os.replace(kp + ".tmp", kp)
    return removed, adopted


def record_key(out: str, stage: str, key: str) -> None:
    kp = os.path.join(out, "stages", "code_keys.json")
    rec = json.load(open(kp)) if os.path.exists(kp) else {}
    rec[stage] = key
    with open(kp + ".tmp", "w") as fh:
        json.dump(rec, fh, indent=1)
    os.replace(kp + ".tmp", kp)


def reuse_stages(src: str, out: str, log) -> list:
    """Hard-link input/families/scores of another run (+ their recorded code keys)."""
    s_sd, o_sd = os.path.join(src, "stages"), os.path.join(out, "stages")
    os.makedirs(o_sd, exist_ok=True)
    skp = os.path.join(s_sd, "code_keys.json")
    srec = json.load(open(skp)) if os.path.exists(skp) else {}
    done = []
    for st in ("input", "families", "scores"):
        sd, od = os.path.join(s_sd, st), os.path.join(o_sd, st)
        if os.path.exists(os.path.join(od, "DONE")) or not os.path.exists(os.path.join(sd, "DONE")):
            continue
        if st not in srec:
            break        # unrecorded upstream stage: do not reuse it or anything below it
        subprocess.run(["cp", "-al", sd, od], check=True)
        record_key(out, st, srec[st])
        done.append(st)
    if done:
        log(f"[driver] reused stages {done} from {src} (hard links)")
    return done


def transplant(src: str, out: str, proteomes: str, tree: str, keys: dict, log) -> list:
    """Hard-link the tree-independent ``families`` and ``scores`` stages of ``src`` into ``out`` and rewrite
    their DONE keys to the keys ``ancestra.pipeline.run_pipeline`` computes for (proteomes, tree).

    The pipeline keys its input stage on the guide-tree md5 and chains every later key to it, so a new guide
    tree would recompute S0 and S2 (~1.6 h on QfO) although neither reads the tree (``s_fam`` uses the gene
    table, ``s_scores`` the gene table + families). The keys are recomputed with the pipeline's own ``_hash``
    and ``md5_file``; the formula is first checked to reproduce ``src``'s own DONE keys exactly (else abort).
    The input stage itself is recomputed by the pipeline (it writes the pruned tree). DONE/stage.json are
    replaced (unlink + write), never modified through the hard link, so ``src`` is untouched.
    """
    from ancestra.pipeline import DEFAULT_FAMILY_CONFIG, _hash, md5_file
    from ancestra.scoring import ScoringConfig

    def pipeline_keys(prot, tr):
        inputs = {"proteomes": {f: md5_file(os.path.join(prot, f)) for f in sorted(os.listdir(prot))
                                if not f.startswith(".")}, "tree": md5_file(tr), "families": None}
        k_in = _hash({"inputs": inputs})
        k_fam = _hash({"up": k_in, "file": None, "cols": [0, 1], "cfg": DEFAULT_FAMILY_CONFIG})
        return {"input": k_in, "families": k_fam, "scores": _hash({"up": k_fam, "cfg": ScoringConfig().key()})}
    sm = json.load(open(os.path.join(src, "manifest.json")))
    src_tree = os.path.join(src, "stages", "input", "guide_tree.nwk")
    # src was run on the same proteome dir; its tree file path is in argv
    argv = sm.get("argv", [])
    src_tree_arg = argv[argv.index("--tree") + 1] if "--tree" in argv else None
    chk = pipeline_keys(proteomes, src_tree_arg or src_tree)
    for st in ("families", "scores"):
        have = open(os.path.join(src, "stages", st, "DONE")).read().strip()
        if have != chk[st]:
            raise SystemExit(f"transplant: key formula does not reproduce {src} {st} ({chk[st]} != {have})")
    new = pipeline_keys(proteomes, tree)
    skp = os.path.join(src, "stages", "code_keys.json")
    srec = json.load(open(skp)) if os.path.exists(skp) else {}
    o_sd = os.path.join(out, "stages")
    os.makedirs(o_sd, exist_ok=True)
    done = []
    for st in ("families", "scores"):
        if srec.get(st) != keys[st]:
            raise SystemExit(f"transplant: {src} stage {st} was built by different code (code key mismatch)")
        od = os.path.join(o_sd, st)
        if os.path.exists(os.path.join(od, "DONE")) and open(os.path.join(od, "DONE")).read().strip() == new[st]:
            continue
        shutil.rmtree(od, ignore_errors=True)
        subprocess.run(["cp", "-al", os.path.join(src, "stages", st), od], check=True)
        for fn, content in (("DONE", new[st] + "\n"),):
            os.remove(os.path.join(od, fn))
            open(os.path.join(od, fn), "w").write(content)
        sj = os.path.join(od, "stage.json")
        info = json.load(open(sj))
        info["key"] = new[st]
        info["transplanted_from"] = {"run": src, "key": chk[st]}
        os.remove(sj)
        json.dump(info, open(sj, "w"), indent=1, default=str)
        record_key(out, st, keys[st])
        done.append(st)
    log(f"[driver] transplanted tree-independent stages {done} from {src} (re-keyed for the new tree)")
    return done


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--gene-order", default="none", help="S5 gene-order table dir, or 'none'")
    ap.add_argument("--reuse-stages-from")
    ap.add_argument("--preset", default=None,
                    help="ancestra.pipeline preset (accurate|fast|custom); default: the pipeline's default")
    ap.add_argument("--stop-after", choices=STAGES + ("rescue",))
    ap.add_argument("--family-rescue", default="none", choices=["none", "3di", "profile", "both"],
                    help="ancestra.pipeline S0 family rescue on cached scores (docs/S0_rescue.md)")
    ap.add_argument("--family-rescue-guard", action="store_true")
    ap.add_argument("--structure-3di", help="3Di FASTA for --family-rescue 3di|both")
    ap.add_argument("--transplant-from", help="run dir whose tree-independent families+scores stages are "
                    "hard-linked and re-keyed for this run's inputs (guide-tree ablations; see docstring)")
    ap.add_argument("--adopt-unrecorded", action="store_true",
                    help="accept cached stages that predate code keys (dry runs only; logged)")
    ap.add_argument("--config", help="JSON {'scoring': {...}, 'merge': {...}} (default: pipeline defaults)")
    a = ap.parse_args(argv)

    from ancestra.families import FamilyConfig
    from ancestra.merge import MergeConfig
    from ancestra.pipeline import DEFAULT_FAMILY_CONFIG, run_pipeline
    from ancestra.scoring import ScoringConfig

    sc, mc = ScoringConfig(), MergeConfig()
    if a.config:
        c = json.load(open(a.config))
        for k, v in c.get("scoring", {}).items():
            setattr(sc, k, v)
        for k, v in c.get("merge", {}).items():
            setattr(mc, k, v)
    go = None if a.gene_order in (None, "", "none") else os.path.abspath(a.gene_order)

    def log(msg: str):
        print(f"{time.strftime('%H:%M:%S')} {msg}", file=sys.stderr, flush=True)

    os.makedirs(a.out, exist_ok=True)
    lock = os.path.join(a.out, ".driver.lock")
    if os.path.exists(lock):
        try:
            pid = int(open(lock).read().strip())
            os.kill(pid, 0)
            raise SystemExit(f"another driver (pid {pid}) is running on {a.out}")
        except (ValueError, ProcessLookupError):
            pass
    with open(lock, "w") as fh:
        fh.write(str(os.getpid()))
    try:
        return _run(a, sc, mc, go, log, FamilyConfig, DEFAULT_FAMILY_CONFIG, run_pipeline)
    finally:
        os.remove(lock)


def _run(a, sc, mc, go, log_print, FamilyConfig, DEFAULT_FAMILY_CONFIG, run_pipeline) -> int:
    fam_cfg = asdict(FamilyConfig(**DEFAULT_FAMILY_CONFIG))
    fam_cfg.pop("threads", None)
    rescue = None
    if a.family_rescue != "none" or a.family_rescue_guard:
        rescue = {"mode": a.family_rescue, "guard": a.family_rescue_guard,
                  "3di_sha256": sha256_file(a.structure_3di) if a.structure_3di else None}
    keys, files = code_keys(fam_cfg, go, a.preset, rescue)

    sdir = os.path.join(a.out, "stages")
    if os.path.isdir(sdir):    # stale temp dirs of killed runs (<stage>.tmp<pid>)
        for f in os.listdir(sdir):
            if re.match(r"^\w+\.tmp\d+$", f):
                shutil.rmtree(os.path.join(sdir, f), ignore_errors=True)
    if a.reuse_stages_from:
        reuse_stages(a.reuse_stages_from, a.out, log_print)
    transplanted = []
    if a.transplant_from:
        transplanted = transplant(a.transplant_from, a.out, a.proteomes, a.tree, keys, log_print)
    removed, adopted = invalidate(a.out, keys, log_print, adopt=a.adopt_unrecorded)

    acc_path = os.path.join(a.out, "stage_accounting.json")
    runs = json.load(open(acc_path)) if os.path.exists(acc_path) else []
    run = {"started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "argv": sys.argv, "threads": a.threads,
           "gene_order": go, "preset": a.preset, "stop_after": a.stop_after, "invalidated": removed, "adopted_unrecorded": adopted, "rescue": rescue,
           "code_keys": keys, "stages": {}}
    runs.append(run)
    cur: dict = {}

    def save():
        with open(acc_path + ".tmp", "w") as fh:
            json.dump(runs, fh, indent=1)
        os.replace(acc_path + ".tmp", acc_path)

    def log(msg: str):
        log_print(msg)
        m = re.match(r"\[ancestra\] stage (\w+): (running|done|up to date)", msg)
        if not m:
            return
        name, what = m.groups()
        if what == "running":
            cur[name] = {"t0": time.time(), "cpu0": _cpu(), "load_start": os.getloadavg()}
            return
        record_key(a.out, name, keys[name])
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
        kw = {"preset": a.preset} if a.preset is not None else {}
        if rescue:
            kw.update(family_rescue=a.family_rescue, rescue_guard=a.family_rescue_guard,
                      structure_3di=a.structure_3di)
        m = run_pipeline(a.proteomes, a.tree, a.out, threads=a.threads, scoring=sc, merge=mc, log=log,
                         gene_order_dir=go, **kw)
        run["total_wall_s"] = m.get("total_wall_s")
        run["finished"] = True
        mp = os.path.join(a.out, "manifest.json")
        man = json.load(open(mp))
        man["driver"] = {"git": git_info(), "code_sha256": files, "stage_code_keys": keys,
                         "stage_code_files": STAGE_CODE, "families_config_effective": fam_cfg,
                         "families_config_hash": _h(fam_cfg), "gene_order_dir": go, "preset": a.preset,
                         "invalidated_this_run": removed, "adopted_unrecorded": adopted,
                         "transplanted_from": a.transplant_from, "transplanted_stages": transplanted,
                         "rescue": rescue,
                         "threads": a.threads}
        with open(mp + ".tmp", "w") as fh:
            json.dump(man, fh, indent=1, default=str)
        os.replace(mp + ".tmp", mp)
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
