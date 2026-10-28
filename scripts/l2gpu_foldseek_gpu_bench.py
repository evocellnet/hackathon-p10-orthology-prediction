#!/usr/bin/env python
"""L2-gpu task 3: Foldseek 3Di+AA search, CPU vs GPU (`--gpu 1`, padded target DB), on the QfO rescue units
(the R4 search: small-family query units vs target representatives), same parameters as ancestra.rescue
(-s 7.5, --max-seqs 100, -e 1e-3, --alignment-type 2). A random sample of query units is searched against ALL
target units; time and hit agreement are reported, and the full-search cost is extrapolated linearly in queries.

    python scripts/l2gpu_foldseek_gpu_bench.py RUN_BASE 3DI.fa OUTDIR N_QUERIES THREADS
"""
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra import structure as S  # noqa: E402
from ancestra.io import load_proteomes  # noqa: E402
from ancestra.rescue import RescueConfig, load_3di_any, select_units, stage1_rep_of  # noqa: E402

FS = os.environ.get("ANCESTRA_FOLDSEEK", "foldseek")


def run(cmd):
    t0 = time.time()
    subprocess.run([str(c) for c in cmd], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return time.time() - t0


def main():
    base, tdi, out, nq, th = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    os.makedirs(out, exist_ok=True)
    gt = load_proteomes(os.path.join(ROOT, "data/inputs/qfo2022/proteomes"))
    fam = np.load(os.path.join(base, "stages/scores/families_final.npy"))
    cfg = RescueConfig(mode="both")
    rep_of = stage1_rep_of(gt, os.path.join(base, "stages/families/s0"))
    queries, targets, _ = select_units(fam, rep_of, gt.lengths, cfg.small_max, cfg.reps_per_family)
    ss = load_3di_any(tdi, gt)
    rng = np.random.default_rng(11)
    qs = np.sort(rng.choice(queries, size=min(nq, len(queries)), replace=False))
    qk = S.build_db(gt, gt.seqs, ss, os.path.join(out, "qdb"), None, genes=qs)
    tk = S.build_db(gt, gt.seqs, ss, os.path.join(out, "tdb"), None, genes=targets)
    res = dict(n_query_units=int(len(queries)), n_target_units=int(len(targets)), sample=int(len(qk)),
               targets_with_3di=int(len(tk)), threads=th)
    common = ["--alignment-type", "2", "-e", cfg.evalue, "--max-seqs", cfg.fs_max_seqs, "-s", cfg.fs_sens,
              "--threads", th, "-v", "1"]
    q, t = os.path.join(out, "qdb"), os.path.join(out, "tdb")
    res["cpu_search_s"] = run([FS, "search", q, t, os.path.join(out, "res_cpu"), os.path.join(out, "tmp_c"), *common])
    res["pad_s"] = run([FS, "makepaddedseqdb", t, os.path.join(out, "tdb_pad")])
    res["gpu_search_s"] = run([FS, "search", q, os.path.join(out, "tdb_pad"), os.path.join(out, "res_gpu"),
                               os.path.join(out, "tmp_g"), *common, "--gpu", "1"])
    hits = {}
    for tag, tdb in (("cpu", t), ("gpu", os.path.join(out, "tdb_pad"))):
        m8 = os.path.join(out, f"{tag}.m8")
        run([FS, "convertalis", q, tdb, os.path.join(out, f"res_{tag}"), m8, "--format-output",
             "query,target,evalue", "--threads", th, "-v", "1"])
        hits[tag] = {(a, b) for a, b, e in (l.split("\t") for l in open(m8)) if float(e) <= cfg.evalue}
    c, g = hits["cpu"], hits["gpu"]
    res.update(cpu_hits=len(c), gpu_hits=len(g), shared=len(c & g),
               jaccard=round(len(c & g) / max(len(c | g), 1), 4))
    f = len(queries) / max(len(qk), 1)
    res["extrapolated_full_cpu_s"] = round(res["cpu_search_s"] * f, 1)
    res["extrapolated_full_gpu_s"] = round(res["gpu_search_s"] * f, 1)
    json.dump(res, open(os.path.join(out, "bench.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
