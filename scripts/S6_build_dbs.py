"""S6: build Foldseek DBs (AA + 3Di) for OrthoBench: real 3Di and shuffled 3Di (seeds 11/23/47).
    python scripts/S6_build_dbs.py
Writes work/S6/db/{real,shuf11,shuf23,shuf47}/ob (+ genes.npy = gene id per DB key) and
results/S6_structure/build_dbs.json (timings, coverage)."""
import json, os, sys, time, resource
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.io import load_proteomes
from ancestra import structure as S

def main():
    t0 = time.time()
    gt = load_proteomes(os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/Input"))
    ss = S.load_3di(os.path.join(ROOT, "data/structure/orthobench_3di_best.fa"), gt)
    info = dict(n_genes=gt.n, n_3di=len(ss), n_len_mismatch=S.load_3di.n_length_mismatch,
                load_s=round(time.time() - t0, 1))
    print(info, flush=True)
    variants = {"real": ss}
    for seed in (11, 23, 47):
        variants[f"shuf{seed}"] = S.shuffled_3di(ss, seed)
    # identity of shuffled vs real (per residue), sanity
    ids = list(ss)[:20000]
    for k, v in variants.items():
        if k == "real":
            continue
        same = sum(sum(a == b for a, b in zip(ss[i], v[i])) for i in ids) / sum(len(ss[i]) for i in ids)
        info[f"{k}_residue_identity_to_real"] = round(same, 4)
    for k, v in variants.items():
        d = os.path.join(ROOT, "work/S6/db", k)
        os.makedirs(d, exist_ok=True)
        t1 = time.time()
        with open(os.path.join(d, "build.log"), "w") as log:
            genes = S.build_db(gt, gt.seqs, v, os.path.join(d, "ob"), log=log)
        np.save(os.path.join(d, "genes.npy"), genes)
        info[f"build_{k}_s"] = round(time.time() - t1, 1)
        print(k, len(genes), info[f"build_{k}_s"], flush=True)
    info["wall_s"] = round(time.time() - t0, 1)
    ru = resource.getrusage(resource.RUSAGE_SELF); rc = resource.getrusage(resource.RUSAGE_CHILDREN)
    info["cpu_s"] = round(ru.ru_utime + ru.ru_stime + rc.ru_utime + rc.ru_stime, 1)
    json.dump(info, open(os.path.join(ROOT, "results/S6_structure/build_dbs.json"), "w"), indent=1)
    print(info)

if __name__ == "__main__":
    main()
