"""S6 exp 1 (link + DEV evaluation): join S0 families through Foldseek 3Di+AA hits between
stage-1 structure representatives (scripts/S6_s0_search.py), for the real and the shuffled-3Di
hits, over a small pre-registered grid; score DEV RefOGs only (completeness, fission, fusion).

Pre-registered (written before the first evaluation):
  grid: e-value <= {1e-3, 1e-5} x small_max {10, None (any size)}; aligned span >= 0.3 of the
        shorter protein (the S0 coverage rule); joins capped at 7,500 genes (S0 gmax).
  selection: max DEV completeness; ties -> lower sum g^2. Eligible: sum g^2 <= 0.05 N^2.
  keep rule (DESIGN §5): completeness(real) - mean(shuffled) > 0 AND
        completeness(real) - completeness(off) > seed spread (max - min over the 3 shuffles).
    python scripts/S6_s0_rescue_eval.py
"""
import json, os, sys, time, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import read_gene_table
from ancestra import structure as S
from ancestra_eval import families_eval as FE

S0 = os.path.join(ROOT, "work/S0/2sA_dc30_cc_sh30_RA")
W = os.path.join(ROOT, "work/S6/s0")
OUT = os.path.join(ROOT, "results/S6_structure")
GRID = [dict(evalue=e, small_max=s, cov_shorter=0.3) for e, s in itertools.product([1e-3, 1e-5], [10, None])]
# exploratory (added after seeing the real-3Di grid, flagged as post hoc in the doc): singletons only,
# and a bidirectional coverage guard against the fusion increase seen in the grid
EXPLORATORY = [dict(evalue=1e-3, small_max=1, cov_shorter=0.3),
               dict(evalue=1e-3, small_max=10, cov_shorter=0.3, cov_min=0.5),
               dict(evalue=1e-3, small_max=1, cov_shorter=0.3, cov_min=0.5),
               dict(evalue=1e-3, small_max=1, cov_shorter=0.3, mode="attach"),
               dict(evalue=1e-3, small_max=10, cov_shorter=0.3, mode="attach"),
               dict(evalue=1e-5, small_max=10, cov_shorter=0.3, mode="attach")]
VARIANTS = ["real", "shuf11", "shuf23", "shuf47"]
EXPORT = dict(evalue=1e-3, small_max=10, cov_shorter=0.3, mode="attach")   # for the downstream run


def fam_stats(lab):
    _, c = np.unique(lab, return_counts=True)
    n = len(lab)
    return dict(n_families=int(len(c)), n_singletons=int((c == 1).sum()), max_family=int(c.max()),
                sum_g2_over_N2=float((c.astype(float) ** 2).sum() / n / n))


def main():
    gt = read_gene_table(os.path.join(ROOT, "work/core_v0_ob/run_v2/stages/input/genes.tsv"))
    idx = {n: i for i, n in enumerate(gt.names)}
    fam = np.full(gt.n, -1, np.int64)
    fid = {}
    for line in open(os.path.join(S0, "families.tsv")):
        if line.startswith("#"):
            continue
        g, f = line.split()
        fam[idx[g]] = fid.setdefault(f, len(fid))
    miss = fam < 0
    fam[miss] = len(fid) + np.arange(miss.sum())
    split = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
    ref = FE.read_refogs(os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"),
                         split["dev"], remove_uncertain=True)
    has = np.zeros(gt.n, bool); has[np.load(os.path.join(ROOT, "work/S6/db/real/genes.npy"))] = True
    reps = np.load(os.path.join(W, "reps.npy"))
    isrep = np.zeros(gt.n, bool); isrep[reps] = True

    def ev(lab):
        g2f = {gt.names[i]: str(int(lab[i])) for i in range(gt.n)}
        r = FE.evaluate_refogs(g2f, ref)
        return r

    off = ev(fam)
    rows = [dict(variant="off", config="-", **off["aggregate"], **fam_stats(fam),
                 incomplete=sorted(k for k, v in off["per_refog"].items() if not v["complete"]))]
    # diagnostics: stray genes of incomplete dev RefOGs
    diag = {}
    for rid, v in off["per_refog"].items():
        if v["complete"]:
            continue
        R = [idx[g] for g in ref[rid]]
        fams, cnt = np.unique(fam[R], return_counts=True)
        main_f = fams[np.argmax(cnt)]
        diag[rid] = dict(size=len(R), main_family_has_3di_rep=bool((isrep & (fam == main_f)).any()),
                         strays=[dict(gene=gt.names[g], len=int(gt.lengths[g]), has_3di=bool(has[g]),
                                      is_struct_rep=bool(isrep[g]), family_size=int((fam == fam[g]).sum()))
                                 for g in R if fam[g] != main_f])
    hits = {v: dict(np.load(os.path.join(W, v, "hits.npz"))) for v in VARIANTS if
            os.path.exists(os.path.join(W, v, "hits.npz"))}
    # best structure hit between each stray's family reps and the main family's reps
    for rid, d in diag.items():
        R = [idx[g] for g in ref[rid]]
        fams, cnt = np.unique(fam[R], return_counts=True)
        main_f = fams[np.argmax(cnt)]
        for s in d["strays"]:
            g = idx[s["gene"]]
            for v, h in hits.items():
                m = ((fam[h["q"]] == fam[g]) & (fam[h["t"]] == main_f)) | ((fam[h["t"]] == fam[g]) & (fam[h["q"]] == main_f))
                s[f"best_e_{v}"] = float(h["evalue"][m].min()) if m.any() else None
    for v, h in hits.items():
        for cfg in GRID + EXPLORATORY:
            t0 = time.time()
            lab, st = S.link_families(fam, h, gt.lengths, max_size=7500, **cfg)
            r = ev(lab)
            if cfg == EXPORT:
                fo = os.path.join(W, v, "families_attach_e1e-3_sm10.tsv")
                with open(fo, "w") as fh:
                    fh.write("#gene_id\tfamily_id\n")
                    fh.write("".join(f"{gt.names[i]}\tS6F{int(lab[i])}\n" for i in range(gt.n)))
            rows.append(dict(variant=v, config=json.dumps(cfg), **r["aggregate"], **fam_stats(lab), link=st,
                             incomplete=sorted(k for k, x in r["per_refog"].items() if not x["complete"]),
                             s=round(time.time() - t0, 1)))
            print(v, cfg, round(r["aggregate"]["completeness"], 3), r["aggregate"]["fusion_rate"],
                  fam_stats(lab), st, flush=True)
    json.dump(dict(rows=rows, stray_diag=diag, grid=GRID), open(os.path.join(OUT, "s0_rescue.json"), "w"), indent=1)
    # table
    keys = ["variant", "config", "n_complete", "completeness", "mean_n_fam", "mean_best_frac", "fusion_rate",
            "mean_foreign_refog_frac", "pair_recall", "pair_precision", "n_families", "n_singletons", "max_family",
            "sum_g2_over_N2"]
    with open(os.path.join(OUT, "s0_rescue.tsv"), "w") as fh:
        fh.write("\t".join(keys + ["n_joins", "incomplete"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(round(r[k], 4) if isinstance(r[k], float) else r[k]) for k in keys)
                     + "\t" + str(r.get("link", {}).get("n_joins", 0)) + "\t" + ",".join(r["incomplete"]) + "\n")
    print(json.dumps(diag, indent=1))


if __name__ == "__main__":
    main()
