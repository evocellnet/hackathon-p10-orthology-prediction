"""S6 exp 1 (search step): structure representatives of the S0 stage-1 clusters searched against
themselves with Foldseek 3Di+AA, for the real DB and the shuffled-3Di DBs.
    python scripts/S6_s0_search.py [--variants real shuf11 shuf23 shuf47] [--threads 10]
Structure rep of a stage-1 cluster = its DIAMOND rep if it has 3Di, else the longest member
with 3Di (clusters without any 3Di member are not represented)."""
import argparse, json, os, sys, time, resource
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import read_gene_table
from ancestra import structure as S

S0 = os.path.join(ROOT, "work/S0/2sA_dc30_cc_sh30_RA")
W = os.path.join(ROOT, "work/S6/s0")


def cpu():
    a = resource.getrusage(resource.RUSAGE_SELF); b = resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="*", default=["real", "shuf11", "shuf23", "shuf47"])
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--evalue", type=float, default=1e-2)
    a = ap.parse_args()
    os.makedirs(W, exist_ok=True)
    gt = read_gene_table(os.path.join(ROOT, "work/core_v0_ob/run_v2/stages/input/genes.tsv"))
    idx = {n: i for i, n in enumerate(gt.names)}
    genes = np.load(os.path.join(ROOT, "work/S6/db/real/genes.npy"))
    key = {int(g): k for k, g in enumerate(genes)}
    clus = {}
    for line in open(os.path.join(S0, "stage1_clusters.tsv")):
        r, m = line.split()
        clus.setdefault(idx[r], []).append(idx[m])
    reps = []
    for r, mem in clus.items():
        if r in key:
            reps.append(r)
            continue
        c = [m for m in mem if m in key]
        if c:
            reps.append(max(c, key=lambda m: (gt.lengths[m], -m)))
    reps = np.array(sorted(set(reps)), dtype=np.int64)
    np.save(os.path.join(W, "reps.npy"), reps)
    info = dict(n_clusters=len(clus), n_struct_reps=int(len(reps)))
    print(info, flush=True)
    for v in a.variants:
        d = os.path.join(W, v)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "keys.txt"), "w") as fh:
            fh.write("".join(f"{key[int(g)]}\n" for g in reps))
        db = os.path.join(ROOT, "work/S6/db", v, "ob")
        sub = os.path.join(d, "reps")
        t0, c0 = time.time(), cpu()
        with open(os.path.join(d, "search.log"), "w") as log:
            for suf in ("", "_ss"):
                S._run([S.FOLDSEEK, "createsubdb", os.path.join(d, "keys.txt"), db + suf, sub + suf, "-v", "1"], log)
            S._run([S.FOLDSEEK, "createsubdb", os.path.join(d, "keys.txt"), db + "_h", sub + "_h", "-v", "1"], log)
            h = S.search(sub, sub, os.path.join(d, "srch"), threads=a.threads, evalue=a.evalue, log=log)
        np.savez(os.path.join(d, "hits.npz"), **h)
        info[v] = dict(n_hits=int(len(h["q"])), wall_s=round(time.time() - t0, 1), cpu_s=round(cpu() - c0, 1))
        print(v, info[v], flush=True)
        json.dump(info, open(os.path.join(ROOT, "results/S6_structure/s0_search.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
