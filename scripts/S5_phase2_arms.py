#!/usr/bin/env python
"""S5 phase 2: synteny arms on OrthoBench DEV with the core merge, scored at the pair level.

Reuses a finished `ancestra run` (default work/core_v0_ob/run_v2: gene table, guide tree,
scores.npz, families_final.npy -> no re-search); each arm = run_merge + flatten.
Base merge config = the core owner's current dev reference ("E syn=off R1-only"):
ARMS['default'] + min_cov=0 (ctx_mode='r1', ctx_window=0.75 unless overridden).
Arms: off; real (SyntenyContext multiplicative, w x lam grid, + variants); shuffled
(the same configs with shuffle_seed 11/23/47).
Per arm: pair-level counts (scripts/orthobench_dev_pairs_eval.py), OG-level dev F
(scripts/orthobench_dev_eval.py), a line in results/core_v0_orthobench/runs.tsv (label "S5p2 ...").
Writes results/S5_synteny/phase2/arms/<arm>.json + .npy (counts [RefOG, class, TP/FN/FP/TN]).
"""
import argparse, hashlib, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from datetime import datetime

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
OUT = os.path.join(ROOT, "results/S5_synteny/phase2/arms")
LOG = os.path.join(ROOT, "results/core_v0_orthobench/runs.tsv")
COLS = ["time", "label", "run", "config_hash", "F", "precision", "recall", "exact", "n_refogs_split",
        "n_refogs_with_fp", "n_pred_groups", "merge_s", "config"]
G = {}


def arms():
    A = {"off": dict(syn=None)}
    for w in (5, 10, 20):
        for lam in (0.05, 0.1, 0.3):
            A[f"real_w{w}_l{lam}"] = dict(syn=dict(mode="multiplicative", lam=lam, w=w, anchor="amb"))
    A["real_w10_tiebreak"] = dict(syn=dict(mode="tiebreak", w=10, anchor="amb"))
    A["real_w10_l0.1_win0.9"] = dict(syn=dict(mode="multiplicative", lam=0.1, w=10, anchor="amb"),
                                     cfg=dict(ctx_window=0.9))
    A["real_w10_l0.1_win0.5"] = dict(syn=dict(mode="multiplicative", lam=0.1, w=10, anchor="amb"),
                                     cfg=dict(ctx_window=0.5))
    A["real_w10_l0.1_coreanch"] = dict(syn=dict(mode="multiplicative", lam=0.1, w=10, anchor="core"))
    A["real_w10_l0.1_ctxall"] = dict(syn=dict(mode="multiplicative", lam=0.1, w=10, anchor="amb"),
                                     cfg=dict(ctx_mode="all"))
    # grid 2 (added after grid 1 showed lam=0.3 best; logged in docs/S5_synteny.md §8)
    for w in (10, 20):
        for lam in (0.5, 1.0):
            A[f"real_w{w}_l{lam}"] = dict(syn=dict(mode="multiplicative", lam=lam, w=w, anchor="amb"))
        A[f"real_w{w}_l0.3_win0.5"] = dict(syn=dict(mode="multiplicative", lam=0.3, w=w, anchor="amb"),
                                          cfg=dict(ctx_window=0.5))
        A[f"real_w{w}_l0.5_win0.5"] = dict(syn=dict(mode="multiplicative", lam=0.5, w=w, anchor="amb"),
                                          cfg=dict(ctx_window=0.5))
        A[f"real_w{w}_l0.3_cap5"] = dict(syn=dict(mode="multiplicative", lam=0.3, w=w, cap=5, anchor="amb"))
        A[f"real_w{w}_l0.3_ctxall"] = dict(syn=dict(mode="multiplicative", lam=0.3, w=w, anchor="amb"),
                                          cfg=dict(ctx_mode="all"))
    for k in [k for k in A if k.startswith("real_")]:
        for sd in (11, 23, 47):
            v = json.loads(json.dumps(A[k]))
            v["syn"]["shuffle_seed"] = sd
            A[k.replace("real_", "shuf_") + f"_s{sd}"] = v
    return A


def init(run):
    from ancestra.io import load_guide_tree, read_gene_table
    from ancestra.scoring import ScoreTable
    from ancestra.synteny import gene_positions
    from orthobench_dev_pairs_eval import load_truth
    gt = read_gene_table(os.path.join(run, "stages/input/genes.tsv"))
    G.update(run=run, gt=gt, tree=load_guide_tree(os.path.join(run, "stages/input/guide_tree.nwk"), gt.species),
             st=ScoreTable.load(os.path.join(run, "stages/scores/scores.npz")))
    stems = [sp.split(".")[0] for sp in gt.species]
    G["cc"], G["oo"], G["cov"] = gene_positions(gt.names, gt.species_of, stems,
                                                os.path.join(ROOT, "data/genomes/orthobench"))
    famf = os.path.join(run, "stages/scores/families_final.npy")
    G["fam"] = np.load(famf if os.path.exists(famf) else os.path.join(run, "stages/families/families.npy"))
    G["truth"] = load_truth()
    G["genes"] = set(G["truth"].gene_a) | set(G["truth"].gene_b)


def run_arm(item):
    name, spec = item
    from ancestra.merge import ARMS, run_merge, to_orthoxml
    from ancestra.orthoxml import to_orthobench_ogs
    from ancestra.synteny import SyntenyContext
    from orthobench_dev_eval import dev_eval
    from orthobench_dev_pairs_eval import counts, metrics, read_pred
    cfg = replace(ARMS["default"], min_cov=0.0, **spec.get("cfg", {}))
    ctx = None
    if spec.get("syn"):
        kw = dict(spec["syn"])
        ctx = SyntenyContext(G["cc"], G["oo"], family=G["fam"], **kw)
    t0 = time.time()
    res = run_merge(G["gt"], G["tree"], G["st"], cfg, context=ctx)
    ms = time.time() - t0
    ox = to_orthoxml(res, G["gt"], G["tree"])
    pred = read_pred(ox, G["genes"])
    arr, refogs = counts(G["truth"], pred)
    m = metrics(arr)
    wd = os.path.join(ROOT, "work/S5/phase2", name)
    os.makedirs(wd, exist_ok=True)
    ogs = os.path.join(wd, "orthobench_ogs.txt")
    to_orthobench_ogs(ox, None, ogs)
    ev = dev_eval(ogs)
    ns = res.node_stats
    rec = dict(arm=name, spec=spec, cfg=asdict(cfg), pair=m, og=dict(F=ev["F"], P=ev["precision"], R=ev["recall"],
               exact=ev["exact"]), merge_s=round(ms, 1), refogs=refogs,
               n_ctx_bonus=int(sum(s.get("n_ctx_bonus", 0) for s in ns)),
               n_anchors=int(sum(s.get("n_anchors", 0) for s in ns)),
               ctx_stats=(dict(n_cand=sum(s["n_cand"] for s in ctx.stats), n_supported=sum(s["n_supported"] for s in ctx.stats))
                          if ctx is not None else None))
    os.makedirs(OUT, exist_ok=True)
    np.save(os.path.join(OUT, f"{name}.npy"), arr)
    json.dump(rec, open(os.path.join(OUT, f"{name}.json"), "w"), indent=1, default=str)
    ch = hashlib.sha256(json.dumps(rec["cfg"], sort_keys=True).encode() + json.dumps(spec, sort_keys=True).encode()).hexdigest()[:10]
    row = dict(time=datetime.now().isoformat(timespec="seconds"), label=f"S5p2 {name}", run=os.path.relpath(G["run"], ROOT),
               config_hash=ch, F=round(ev["F"], 2), precision=round(ev["precision"], 2), recall=round(ev["recall"], 2),
               exact=ev["exact"], n_refogs_split=ev["n_refogs_split"], n_refogs_with_fp=ev["n_refogs_with_fp"],
               n_pred_groups=ev["n_pred_groups"], merge_s=round(ms, 1),
               config=json.dumps(dict(rec["cfg"], synteny=spec.get("syn"), pairF1=m["all"]["F1"]), sort_keys=True))
    with open(LOG, "a") as fh:
        fh.write("\t".join(str(row[c]) for c in COLS) + "\n")
    print(f"{name}: pairF1={m['all']['F1']} P={m['all']['P']} R={m['all']['R']} | vt={m['vert-teleost']['F1']} "
          f"tt={m['teleost-teleost']['F1']} mm={m['mammal-mammal']['F1']} | OG F={ev['F']:.2f} merge={ms:.0f}s",
          flush=True)
    return name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=os.path.join(ROOT, "work/core_v0_ob/run_v2"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--skip-done", action="store_true")
    a = ap.parse_args()
    A = arms()
    items = [(k, v) for k, v in A.items() if (not a.only or k in a.only)
             and not (a.skip_done and os.path.exists(os.path.join(OUT, f"{k}.json")))]
    print(len(items), "arms", flush=True)
    os.environ.setdefault("NUMBA_NUM_THREADS", "1")
    with ProcessPoolExecutor(a.workers, initializer=init, initargs=(a.run,)) as ex:
        for _ in ex.map(run_arm, items):
            pass


if __name__ == "__main__":
    main()
