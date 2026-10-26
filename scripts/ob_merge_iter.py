"""Re-run only the S3 merge on a finished OrthoBench `ancestra run` and score DEV RefOGs.

    python scripts/ob_merge_iter.py RUN_DIR --arm default [--set rd_hidden=true rd_z=3] [--label X]

Loads the run's checkpointed gene table, guide tree and scores (no re-search), runs
run_merge with the arm (+ overrides), writes OUT/orthobench_ogs.txt and appends one line
to results/core_v0_orthobench/runs.tsv (every dev config is logged; test RefOGs never scored).
"""
import argparse
import hashlib
import json
import os
import sys
import time
from dataclasses import asdict, replace
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import load_guide_tree, read_gene_table  # noqa: E402
from ancestra.merge import ARMS, run_merge, to_orthoxml  # noqa: E402
from ancestra.orthoxml import to_orthobench_ogs  # noqa: E402
from ancestra.scoring import ScoreTable  # noqa: E402
from orthobench_dev_eval import dev_eval  # noqa: E402

LOG = os.path.join(ROOT, "results/core_v0_orthobench/runs.tsv")
COLS = ["time", "label", "run", "config_hash", "F", "precision", "recall", "exact", "n_refogs_split",
        "n_refogs_with_fp", "n_pred_groups", "merge_s", "config"]


def parse_val(v):
    try:
        return json.loads(v)
    except json.JSONDecodeError:
        return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--arm", default="default")
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--label")
    ap.add_argument("--out")
    ap.add_argument("--scores", default="scores.npz", help="file in RUN/stages/scores/ (e.g. scores_pass1.npz)")
    ap.add_argument("--synteny", default="off", choices=["off", "real", "shuffled"])
    ap.add_argument("--syn", nargs="*", default=[], help="SyntenyContext kwargs k=v (JSON values)")
    a = ap.parse_args()
    cfg = ARMS[a.arm]
    over = dict(kv.split("=", 1) for kv in a.set)
    cfg = replace(cfg, **{k: parse_val(v) for k, v in over.items()})
    ch = hashlib.sha256((json.dumps(asdict(cfg), sort_keys=True) + a.scores + a.synteny
                         + " ".join(a.syn)).encode()).hexdigest()[:10]
    out = a.out or os.path.join(a.run, "iter", ch)
    os.makedirs(out, exist_ok=True)
    gt = read_gene_table(os.path.join(a.run, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(a.run, "stages/input/guide_tree.nwk"), gt.species)
    st = ScoreTable.load(os.path.join(a.run, "stages/scores", a.scores))
    context = None
    syn_kw = {}
    if a.synteny != "off":
        import numpy as np
        from ancestra.synteny import SyntenyContext, gene_positions
        stems = [sp.split(".")[0] for sp in gt.species]
        cc, oo, cov = gene_positions(gt.names, gt.species_of, stems,
                                     os.path.join(ROOT, "data/genomes/orthobench"))
        famf = os.path.join(a.run, "stages/scores/families_final.npy")
        fam = np.load(famf if os.path.exists(famf) else os.path.join(a.run, "stages/families/families.npy"))
        syn_kw = dict(mode="multiplicative", lam=0.1, anchor="amb", w=5)
        syn_kw.update({k: parse_val(v) for k, v in (kv.split("=", 1) for kv in a.syn)})
        if a.synteny == "shuffled":
            syn_kw["shuffle_seed"] = 11
        context = SyntenyContext(cc, oo, family=fam, **syn_kw)
    t0 = time.time()
    res = run_merge(gt, tree, st, cfg, context=context)
    ms = time.time() - t0
    ox = to_orthoxml(res, gt, tree)
    ogs = os.path.join(out, "orthobench_ogs.txt")
    to_orthobench_ogs(ox, None, ogs)
    json.dump(res.node_stats, open(os.path.join(out, "node_stats.json"), "w"))
    ev = dev_eval(ogs)
    json.dump(ev, open(os.path.join(out, "dev_eval.json"), "w"), indent=1)
    row = dict(time=datetime.now().isoformat(timespec="seconds"), label=a.label or a.arm, run=a.run,
               config_hash=ch, F=round(ev["F"], 2), precision=round(ev["precision"], 2),
               recall=round(ev["recall"], 2), exact=ev["exact"], n_refogs_split=ev["n_refogs_split"],
               n_refogs_with_fp=ev["n_refogs_with_fp"], n_pred_groups=ev["n_pred_groups"],
               merge_s=round(ms, 1), config=json.dumps(dict(asdict(cfg), scores=a.scores, synteny=a.synteny,
                                                   synteny_kw=syn_kw), sort_keys=True))
    new = not os.path.exists(LOG)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a") as fh:
        if new:
            fh.write("\t".join(COLS) + "\n")
        fh.write("\t".join(str(row[c]) for c in COLS) + "\n")
    print(f"{row['label']} [{ch}] DEV F={row['F']} P={row['precision']} R={row['recall']} exact={row['exact']} "
          f"split={row['n_refogs_split']} withFP={row['n_refogs_with_fp']} merge={ms:.1f}s -> {out}")


if __name__ == "__main__":
    main()
