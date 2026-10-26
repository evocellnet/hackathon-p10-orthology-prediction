"""Core v0 toy benchmark grid (S2/S3/S3b on the in-house simulator).

    python scripts/run_core_v0_toy.py [--sims work/core_v0/sims] [--workers 6] [--out results/core_v0]

Conditions: every simulation directory (def/hard x 8/32 species x seeds 11/23/47)
x families {truth (simulator root families, isolates S3), s0 (default S0
clustering, ancestra.families DEFAULT_FAMILY_CONFIG)} x merge arms (ARMS below).
Scoring: DIAMOND blocks, e-value 10 (default); on truth families also MMseqs2
exhaustive (e 10) and DIAMOND e 1e-3 with the default arm (S2 arms).

Pre-stated selection rule (fixed before the grid ran, 2026-10-26 ~15:55 CET):
default merge arm = highest mean ortholog-pair F1 over the 12 truth-family
conditions (diamond, e 10); ties within 0.002 -> higher mean precision.

Writes one JSON line per (sim, families, backend, evalue, arm) to
OUT/toy_grid.jsonl and full evaluation JSONs under work/core_v0/runs/.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ancestra.evaluate_toy import evaluate, partition_metrics  # noqa: E402
from ancestra.io import families_from_mapping, load_families, load_guide_tree, load_proteomes, write_fasta  # noqa: E402
from ancestra.merge import ARMS, run_merge, to_orthoxml  # noqa: E402
from ancestra.orthoxml import read_orthoxml  # noqa: E402
from ancestra.scoring import ScoreTable, ScoringConfig, score_families, score_with_rescue  # noqa: E402

S2_ARM_MERGE = ["max+repair", "avg+repair"]


def summarize(r: dict, tree_depth: dict) -> dict:
    p = r["pairs"]
    h = r["hogs"]
    lv = {k: v for k, v in h.items() if k != "families"}
    deep = [v for v in lv.values() if v["depth"] <= 1]
    shallow = [v for v in lv.values() if v["depth"] >= 2]
    m = lambda xs, k: float(np.mean([x[k] for x in xs])) if xs else None  # noqa: E731
    return dict(P=p["precision"], R=p["recall"], F1=p["f1"], tp=p["tp"], n_pred=p["n_pred"],
                n_true=p["n_true"], fp_paralogs=p["fp_true_paralogs"],
                fp_cross_family=p["fp_different_families"], fp_other=p["fp_other"],
                root_exact=h["families"]["exact"], root_jacc=h["families"]["mean_best_jaccard"],
                root_fissions=h["families"]["fissions"], root_fusions=h["families"]["fusions"],
                root_pairP=h["families"]["pair_precision"], root_pairR=h["families"]["pair_recall"],
                root_pairF1=h["families"]["pair_f1"],
                lvl_exact=m(list(lv.values()), "exact"), lvl_jacc=m(list(lv.values()), "mean_best_jaccard"),
                lvl_pairF1=m(list(lv.values()), "pair_f1"),
                deep_exact=m(deep, "exact"), shallow_exact=m(shallow, "exact"),
                lvl_fissions=int(sum(v["fissions"] for v in lv.values())),
                lvl_fusions=int(sum(v["fusions"] for v in lv.values())),
                recall_by_level={k: v["recall"] for k, v in r["pair_recall_by_level"].items()})


def run_condition(sim_dir: str, famsrc: str, work: str, threads: int, scorings=None,
                  arm_names=None) -> list[dict]:
    name = os.path.basename(sim_dir.rstrip("/"))
    setting, rest = name.split("_s")[0], name.split("_s")[1]
    nsp = int("".join(c for c in setting if c.isdigit()))
    kind = "".join(c for c in setting if not c.isdigit())
    seed = int(rest)
    wd = os.path.join(work, name, famsrc)
    os.makedirs(wd, exist_ok=True)
    gt = load_proteomes(os.path.join(sim_dir, "proteomes"))
    tree = load_guide_tree(os.path.join(sim_dir, "species_tree.nwk"), gt.species)
    truth_path = os.path.join(sim_dir, "truth", "hogs.orthoxml")
    truth = read_orthoxml(truth_path)
    fam_truth = load_families(os.path.join(sim_dir, "truth", "gene_families.tsv"), gt, 0, 2)
    base = dict(sim=name, setting=kind, n_species=nsp, seed=seed, families=famsrc, n_genes=gt.n)
    s0_info = {}
    if famsrc == "truth":
        fam = fam_truth
    else:
        from ancestra.families import FamilyConfig, read_families_tsv, run_families
        from ancestra.pipeline import DEFAULT_FAMILY_CONFIG
        ftsv = os.path.join(wd, "s0", "families.tsv")
        t0 = time.time()
        if not os.path.exists(ftsv):
            cfg = FamilyConfig(**DEFAULT_FAMILY_CONFIG)
            cfg.threads = threads
            fa = os.path.join(wd, "genes.fa")
            write_fasta(gt, fa)
            run_families(fa, os.path.join(wd, "s0"), cfg)
            os.remove(fa)
        s0_info["s0_wall_s"] = round(time.time() - t0, 2)
        fam = families_from_mapping(gt, read_families_tsv(ftsv))
        tf = [set(np.flatnonzero(fam_truth == f)) for f in np.unique(fam_truth)]
        pf = [set(np.flatnonzero(fam == f)) for f in np.unique(fam)]
        pm = partition_metrics(tf, pf)
        s0_info.update(s0_exact=pm["exact"], s0_fissions=pm["fissions"], s0_fusions=pm["fusions"],
                       s0_n_families=len(pf), s0_pair_recall=pm["pair_recall"],
                       s0_pair_precision=pm["pair_precision"])
    depth = {}
    for nd in tree.preorder():
        depth[nd.name] = 0 if nd.parent is None else depth[nd.parent.name] + 1
    if scorings is None:
        scorings = [("diamond", 10.0)]
        if famsrc == "truth":
            scorings += [("mmseqs", 10.0), ("diamond", 1e-3)]
    out = []
    for backend, ev in scorings:
        sp = os.path.join(wd, f"scores_{backend}_e{ev:g}.npz")
        t0 = time.time()
        if os.path.exists(sp):
            st = ScoreTable.load(sp)
            sstats = json.load(open(sp + ".json"))
        else:
            if backend.endswith("+rescue"):
                st, _, sstats = score_with_rescue(gt, fam, ScoringConfig(backend=backend.split("+")[0],
                                                                         evalue=ev, threads=threads),
                                                  os.path.join(wd, f"tmp_{backend}_{ev:g}"))
            else:
                st, sstats = score_families(gt, fam, ScoringConfig(backend=backend, evalue=ev, threads=threads),
                                            os.path.join(wd, f"tmp_{backend}_{ev:g}"))
            st.save(sp)
            json.dump(sstats, open(sp + ".json", "w"), default=str)
        if arm_names is not None:
            arms = {k: ARMS[k] for k in arm_names}
        elif (backend, ev) == ("diamond", 10.0):
            arms = ARMS
        else:
            arms = {k: ARMS[k] for k in S2_ARM_MERGE}
        for arm, cfg in arms.items():
            t0 = time.time()
            res = run_merge(gt, tree, st, cfg)
            tm = time.time() - t0
            ox = to_orthoxml(res, gt, tree)
            r = evaluate(ox, truth)
            ej = os.path.join(wd, f"eval_{backend}_e{ev:g}_{arm.replace('/', '_')}.json")
            json.dump(r, open(ej, "w"))
            rec = dict(base, backend=backend, evalue=ev, arm=arm, merge_cfg=asdict(cfg),
                       merge_s=round(tm, 3), score_s=sstats.get("wall_s"),
                       n_scored_pairs=sstats.get("n_scored_pairs"), **s0_info)
            rec.update(summarize(r, depth))
            out.append(rec)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default=os.path.join(ROOT, "work/core_v0/sims"))
    ap.add_argument("--work", default=os.path.join(ROOT, "work/core_v0/runs"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results/core_v0"))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--only", nargs="*", help="sim names subset")
    ap.add_argument("--families", nargs="*", default=["truth", "s0"])
    ap.add_argument("--scorings", nargs="*", help="backend:evalue list (default: built-in set)")
    ap.add_argument("--arms", nargs="*", help="merge arm subset")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    sims = sorted(d for d in os.listdir(a.sims) if os.path.isdir(os.path.join(a.sims, d)))
    if a.only:
        sims = [s for s in sims if s in a.only]
    jobs = [(os.path.join(a.sims, s), f) for s in sims for f in a.families]
    # 32-species jobs first (longest)
    jobs.sort(key=lambda j: ("32" not in j[0], j))
    outp = os.path.join(a.out, "toy_grid.jsonl")
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex, open(outp, "a") as fh:
        sc = [(x.split(":")[0], float(x.split(":")[1])) for x in a.scorings] if a.scorings else None
        futs = {ex.submit(run_condition, s, f, a.work, a.threads, sc, a.arms): (s, f) for s, f in jobs}
        for fu in as_completed(futs):
            s, f = futs[fu]
            try:
                recs = fu.result()
            except Exception:
                print(f"FAILED {s} {f}\n{traceback.format_exc()}", flush=True)
                continue
            for r in recs:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            best = max(recs, key=lambda r: r["F1"])
            print(f"[{time.time() - t0:7.1f}s] {os.path.basename(s)} {f}: {len(recs)} records; "
                  f"best {best['arm']} F1={best['F1']:.4f}", flush=True)


if __name__ == "__main__":
    main()
