#!/usr/bin/env python
"""v2-gpu diagnostics: build a run directory whose cached stages are *given* arrays, so that a
following `python -m ancestra run` (same inputs, same --backend/--config) recomputes only what is
missing. Experiment tool only (docs/GPU_v2.md §2); never used by the method.

    python scripts/gpu_inject_stages.py --src RUN --out NEW [--families FAM.npy]
        [--scores SCORES.npz --families-final FAMF.npy] [--backend hybrid] [--config cfg.json]

- input stage: hard-linked from SRC.
- families stage: FAM.npy (default: SRC's) under the key of the default S0 config.
- scores stage (optional): SCORES.npz + FAMF.npy under the key of the given backend/config.
"""
import argparse
import json
import os
import shutil
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.pipeline import DEFAULT_FAMILY_CONFIG, _hash  # noqa: E402
from ancestra.scoring import ScoringConfig  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--families")
    ap.add_argument("--scores")
    ap.add_argument("--families-final")
    ap.add_argument("--backend", default="hybrid")
    ap.add_argument("--config")
    ap.add_argument("--proteomes", help="recompute the input key as `ancestra run` does (when SRC was made "
                                        "by another driver with different keys)")
    ap.add_argument("--tree")
    a = ap.parse_args()
    sd = os.path.join(a.out, "stages")
    os.makedirs(sd, exist_ok=True)
    src_in = os.path.join(a.src, "stages", "input")
    shutil.copytree(src_in, os.path.join(sd, "input"), copy_function=os.link, dirs_exist_ok=True)
    k_in = open(os.path.join(src_in, "DONE")).read().strip()
    if a.proteomes:
        from ancestra.pipeline import md5_file
        pf = sorted(os.listdir(a.proteomes))
        inputs = {"proteomes": {f: md5_file(os.path.join(a.proteomes, f)) for f in pf if not f.startswith(".")},
                  "tree": md5_file(a.tree), "families": None}
        k_in = _hash({"inputs": inputs})
        open(os.path.join(sd, "input", "DONE.tmp"), "w").write(k_in + "\n")
        os.remove(os.path.join(sd, "input", "DONE"))          # hard link: replace, do not write through
        os.replace(os.path.join(sd, "input", "DONE.tmp"), os.path.join(sd, "input", "DONE"))
        print("input key", k_in)
    fam_key = _hash({"up": k_in, "file": None, "cols": [0, 1], "cfg": DEFAULT_FAMILY_CONFIG})
    fd = os.path.join(sd, "families")
    os.makedirs(fd, exist_ok=True)
    fam_src = a.families or os.path.join(a.src, "stages", "families", "families.npy")
    fam = np.load(fam_src)
    np.save(os.path.join(fd, "families.npy"), fam)
    json.dump({"key": fam_key, "stats": {"injected_from": os.path.abspath(fam_src)}, "skipped": False},
              open(os.path.join(fd, "stage.json"), "w"))
    open(os.path.join(fd, "DONE"), "w").write(fam_key + "\n")
    print("families key", fam_key, "n_fam", int(fam.max()) + 1)
    if a.scores:
        sc = ScoringConfig(backend=a.backend)
        if a.config:
            for k, v in json.load(open(a.config)).get("scoring", {}).items():
                setattr(sc, k, v)
        k_sc = _hash({"up": fam_key, "cfg": sc.key()})
        scd = os.path.join(sd, "scores")
        os.makedirs(scd, exist_ok=True)
        os.link(a.scores, os.path.join(scd, "scores.npz"))
        np.save(os.path.join(scd, "families_final.npy"), np.load(a.families_final or fam_src))
        json.dump({"key": k_sc, "stats": {"injected_from": os.path.abspath(a.scores)}, "skipped": False},
                  open(os.path.join(scd, "stage.json"), "w"))
        open(os.path.join(scd, "DONE"), "w").write(k_sc + "\n")
        print("scores key", k_sc)


if __name__ == "__main__":
    main()
