"""S6 exp 2/3: run the S3 merge (default config = hybrid local-tree mode) on run_v2's checkpointed
inputs with a structure-aware ScoreTable (ancestra.structure.combine_scores) and score DEV
RefOGs (OG-level F via orthobench_dev_eval, pair-level P/R/F1 via orthobench_dev_pairs_eval,
plus pair metrics split by 3Di coverage of the two genes: both / one / none).

Arms (pre-registered before the first merge run):
  scope  deep = structure only for pairs whose species LCA is Euteleostomi, Chordata,
               Ecdysozoa or Bilateria (the tree-mode nodes + root); all = every pair
  mode   avg (geometric mean of S_seq and calibrated S_3di) | max | struct (S_3di only)
  pairs  seq (pairs scored by sequence) [+miss: pairs sequence left unscored]
  3Di    real | shuf11 | shuf23 | shuf47 ; off = sequence scores only
    python scripts/S6_merge_eval.py --arms off deep:avg:seq all:avg:seq ... [--variants real shuf11 ...]
"""
import argparse, json, os, sys, time, hashlib
from datetime import datetime
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import load_guide_tree, read_gene_table
from ancestra.merge import ARMS, run_merge, to_orthoxml
from ancestra.orthoxml import to_orthobench_ogs
from ancestra.scoring import ScoreTable
from ancestra import structure as S
from orthobench_dev_eval import dev_eval
from orthobench_dev_pairs_eval import load_truth, read_pred, counts, metrics, paired_bootstrap

RUN = os.path.join(ROOT, "work/core_v0_ob/run_v2")
W = os.path.join(ROOT, "work/S6/merge")
OUT = os.path.join(ROOT, "results/S6_structure")
DEEP = {"Euteleostomi", "Chordata", "Ecdysozoa", "Bilateria"}
SHALLOW_CLS = ["mammal-mammal", "tetrapod-tetrapod", "teleost-teleost"]
DEEP_CLS = ["vert-teleost", "vert-invert", "invert-invert"]


def load_fs(variant, sets, n):
    qs, ts, bs, es = [], [], [], []
    for s in sets:
        z = np.load(os.path.join(ROOT, "work/S6/pairs", s, variant, "fs.npz"))
        qs.append(z["q"]); ts.append(z["t"]); bs.append(z["bits"]); es.append(z["evalue"])
    q, t, b, e = np.concatenate(qs), np.concatenate(ts), np.concatenate(bs), np.concatenate(es)
    self_ = np.full(n, np.nan)
    m = q == t
    self_[q[m]] = b[m]
    return dict(q=q[~m], t=t[~m], bits=b[~m], evalue=e[~m]), self_


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["off", "deep:avg:seq", "all:avg:seq", "deep:max:seq"])
    ap.add_argument("--variants", nargs="*", default=["real", "shuf11", "shuf23", "shuf47"])
    ap.add_argument("--tag", default="")
    ap.add_argument("--base-scores", default=os.path.join(RUN, "stages/scores/scores.npz"),
                    help="sequence ScoreTable (default run_v2; e.g. work/S6/s0_down/real/scores.npz)")
    a = ap.parse_args()
    os.makedirs(W, exist_ok=True)
    gt = read_gene_table(os.path.join(RUN, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(RUN, "stages/input/guide_tree.nwk"), gt.species)
    st0 = ScoreTable.load(a.base_scores)
    cfg = ARMS["default"]
    truth = load_truth()
    has = np.zeros(gt.n, bool); has[np.load(os.path.join(ROOT, "work/S6/db/real/genes.npy"))] = True
    hasn = {gt.names[i]: bool(has[i]) for i in range(gt.n)}
    cov = np.array([hasn[x] + hasn[y] for x, y in zip(truth.gene_a, truth.gene_b)])
    genes = set(truth.gene_a) | set(truth.gene_b)
    _, names = S.pair_lca_nodes(gt, tree, np.zeros(1, int), np.zeros(1, int))
    deep_ids = np.array([i for i, nm in enumerate(names) if nm in DEEP])

    def deep_mask(q, t):
        l, _ = S.pair_lca_nodes(gt, tree, q, t)
        return np.isin(l, deep_ids)

    off_arr = None
    resf = os.path.join(OUT, "merge_runs.jsonl")
    jobs = []
    for arm in a.arms:
        if arm == "off":
            jobs.append(("off", "-"))
        else:
            for v in a.variants:
                jobs.append((arm, v))
    for arm, v in jobs:
        t0 = time.time()
        if arm == "off":
            st, cst = st0, {}
        else:
            scope, mode, pairs = arm.split(":")
            fs, fs_self = load_fs(v, pairs.split("+"), gt.n)
            st, cst = S.combine_scores(gt, st0, fs, fs_self, mode=mode,
                                       node_mask=deep_mask if scope == "deep" else None,
                                       add_missing="miss" in pairs)
        tc = time.time() - t0
        t1 = time.time()
        res = run_merge(gt, tree, st, cfg)
        ms = time.time() - t1
        ox = to_orthoxml(res, gt, tree)
        h = hashlib.sha256(f"{arm}|{v}".encode()).hexdigest()[:8]
        out = os.path.join(W, f"{a.tag}{arm.replace(':', '_')}_{v}")
        os.makedirs(out, exist_ok=True)
        ogs = os.path.join(out, "orthobench_ogs.txt")
        to_orthobench_ogs(ox, None, ogs)
        ev = dev_eval(ogs)
        pred = read_pred(ox, genes)
        arr, refogs = counts(truth, pred)
        pm = metrics(arr)
        covm = {}
        for k, nm in ((2, "both_3di"), (1, "one_3di"), (0, "no_3di")):
            arr_k, _ = counts(truth[cov == k], pred)
            covm[nm] = metrics(arr_k)["all"]
        tot = lambda cls: metrics(arr)  # noqa
        grp = {}
        for gname, cl in (("shallow", SHALLOW_CLS), ("deep", DEEP_CLS)):
            c = sum(arr[:, [list(pm.keys()).index(x) - 1 for x in cl]].sum(axis=(0, 1)) for _ in [0])
            tp, fn, fp = c[0], c[1], c[2]
            P = tp / (tp + fp) if tp + fp else float("nan"); R = tp / (tp + fn) if tp + fn else float("nan")
            grp[gname] = dict(P=round(P, 4), R=round(R, 4), F1=round(2 * P * R / (P + R), 4), tp=int(tp), fn=int(fn), fp=int(fp))
        row = dict(time=datetime.now().isoformat(timespec="seconds"), arm=arm, variant=v, tag=a.tag,
                   F=round(ev["F"], 2), P_og=round(ev["precision"], 2), R_og=round(ev["recall"], 2),
                   exact=ev["exact"], n_refogs_split=ev["n_refogs_split"], n_refogs_with_fp=ev["n_refogs_with_fp"],
                   pairP=pm["all"]["P"], pairR=pm["all"]["R"], pairF1=pm["all"]["F1"],
                   by_class={k: pm[k]["F1"] for k in pm if k != "all"}, by_depth=grp, by_coverage=covm,
                   combine=cst, combine_s=round(tc, 1), merge_s=round(ms, 1),
                   node_stats_tree={k: {kk: vv for kk, vv in s.items() if "tree" in kk} for k, s in res.node_stats.items()}
                   if isinstance(res.node_stats, dict) else None)
        if arm == "off":
            off_arr = arr
        elif off_arr is not None:
            row["delta_vs_off"] = paired_bootstrap(off_arr, arr)["all"]
        np.save(os.path.join(out, "pair_counts.npy"), arr)
        with open(resf, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        print(f"{a.tag}{arm:16s} {v:7s} F={row['F']} P={row['P_og']} R={row['R_og']} ex={row['exact']} | pair P={row['pairP']} "
              f"R={row['pairR']} F1={row['pairF1']} | deep F1={grp['deep']['F1']} shallow F1={grp['shallow']['F1']} | "
              f"cov both={covm['both_3di']['F1']} one={covm['one_3di']['F1']} none={covm['no_3di']['F1']} | "
              f"fit={cst.get('fit_a', 0):.3f},{cst.get('fit_b', 1):.3f} changed={cst.get('n_pairs_changed')} "
              f"+{cst.get('n_struct_only_added')} | {tc:.0f}+{ms:.0f}s", flush=True)


if __name__ == "__main__":
    main()
