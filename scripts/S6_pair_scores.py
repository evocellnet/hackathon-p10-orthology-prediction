"""S6 exp 2 (scoring step): Foldseek 3Di+AA alignment scores for within-family gene pairs of the
families that hold DEV RefOG genes (evaluation scope only: the merge is family-local, so other
families cannot affect DEV metrics; the method would score every family -- cost extrapolated in
docs/S6_structure.md).
Pair set 'seq': every pair of the sequence ScoreTable (run_v2) inside those families with both
genes in the 3Di DB, plus self pairs. Pair set 'miss' (optional): in-family pairs (families <=
5,000 genes, the scoring cap_all) that sequence left unscored, both genes with 3Di.
    python scripts/S6_pair_scores.py --set seq --variants real shuf11 shuf23 shuf47
"""
import argparse, json, os, sys, time, resource
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import read_gene_table
from ancestra.scoring import ScoreTable, family_members
from ancestra import structure as S
from ancestra_eval import families_eval as FE

RUN = os.path.join(ROOT, "work/core_v0_ob/run_v2")
W = os.path.join(ROOT, "work/S6/pairs")


def cpu():
    a = resource.getrusage(resource.RUSAGE_SELF); b = resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def dev_families(gt, fam):
    split = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
    ref = FE.read_refogs(os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"),
                         split["dev"], remove_uncertain=True)
    idx = {n: i for i, n in enumerate(gt.names)}
    g = [idx[x] for r in ref.values() for x in r if x in idx]
    return np.unique(fam[g])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="seq", choices=["seq", "miss"])
    ap.add_argument("--variants", nargs="*", default=["real", "shuf11", "shuf23", "shuf47"])
    ap.add_argument("--threads", type=int, default=10)
    a = ap.parse_args()
    gt = read_gene_table(os.path.join(RUN, "stages/input/genes.tsv"))
    fam = np.load(os.path.join(RUN, "stages/scores/families_final.npy"))
    devf = dev_families(gt, fam)
    infam = np.isin(fam, devf)
    genes = np.load(os.path.join(ROOT, "work/S6/db/real/genes.npy"))
    has = np.zeros(gt.n, bool); has[genes] = True
    key = {int(g): k for k, g in enumerate(genes)}
    st = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
    q0, t0 = st.q.astype(np.int64), st.t.astype(np.int64)
    if a.set == "seq":
        m = infam[q0] & has[q0] & has[t0]
        q, t = q0[m], t0[m]
        selfg = np.flatnonzero(infam & has)
        q = np.concatenate([q, selfg]); t = np.concatenate([t, selfg])
    else:
        n = gt.n
        k0 = set((q0 * n + t0)[infam[q0]].tolist())
        qs, ts = [], []
        for mem in family_members(fam):
            if not infam[mem[0]] or len(mem) > 5000:
                continue
            mm = np.sort(mem[has[mem]])
            i, j = np.triu_indices(len(mm), 1)
            kk = mm[i] * n + mm[j]
            sel = np.fromiter((x not in k0 for x in kk.tolist()), bool, len(kk))
            qs.append(mm[i][sel]); ts.append(mm[j][sel])
        q, t = np.concatenate(qs), np.concatenate(ts)
    info = dict(set=a.set, n_dev_families=int(len(devf)), n_genes_dev_families=int(infam.sum()),
                n_genes_with_3di=int((infam & has).sum()), n_pairs=int(len(q)))
    print(info, flush=True)
    out_json = os.path.join(ROOT, f"results/S6_structure/pair_scores_{a.set}.json")
    for v in a.variants:
        d = os.path.join(W, a.set, v)
        os.makedirs(d, exist_ok=True)
        t1, c1 = time.time(), cpu()
        with open(os.path.join(d, "align.log"), "w") as log:
            h = S.align_pairs(os.path.join(ROOT, "work/S6/db", v, "ob"), key, q, t, d, threads=a.threads, log=log)
        np.savez(os.path.join(d, "fs.npz"), **h)
        for f in os.listdir(d):
            if f.startswith(("aln", "pref")):
                os.remove(os.path.join(d, f))
        info[v] = dict(n_aligned=int(len(h["q"])), wall_s=round(time.time() - t1, 1), cpu_s=round(cpu() - c1, 1))
        print(v, info[v], flush=True)
        json.dump(info, open(out_json, "w"), indent=1)


if __name__ == "__main__":
    main()
