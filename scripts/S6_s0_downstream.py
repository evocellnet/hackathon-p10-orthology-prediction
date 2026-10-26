"""S6 exp 1 downstream: effect of the structure S0 rescue on the merge (DEV OG-F and pair F1).

Families = run_v2's final families (S0 + core cross-family rescue) joined by the 3Di attach rule
(e <= 1e-3, small families <= 10 genes attach to their best-hit family; ancestra.structure.
link_families mode="attach") with real or shuffled-3Di hits. New within-family pairs (genes of
different old families now in one family; only new families holding DEV RefOG genes = evaluation
scope) are aligned with the core's MMseqs2 explicit-pair aligner (scoring._mmseqs_align_pairs,
e <= 10, as the hybrid backend) and appended to run_v2's ScoreTable; then the default merge.
    python scripts/S6_s0_downstream.py [--variants real shuf11 shuf23 shuf47] [--threads 10]
"""
import argparse, json, os, sys, time, resource
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import load_proteomes
from ancestra.scoring import ScoreTable, ScoringConfig, _mmseqs_align_pairs, _reduce_records, family_members
from ancestra import structure as S
from S6_pair_scores import dev_families

RUN = os.path.join(ROOT, "work/core_v0_ob/run_v2")
W = os.path.join(ROOT, "work/S6/s0_down")
CFG = dict(evalue=1e-3, small_max=10, cov_shorter=0.3, mode="attach")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="*", default=["real", "shuf11", "shuf23", "shuf47"])
    ap.add_argument("--threads", type=int, default=10)
    a = ap.parse_args()
    gt = load_proteomes(os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/Input"))
    fam = np.load(os.path.join(RUN, "stages/scores/families_final.npy"))
    st0 = ScoreTable.load(os.path.join(RUN, "stages/scores/scores.npz"))
    info = {}
    for v in a.variants:
        t0 = time.time()
        h = dict(np.load(os.path.join(ROOT, "work/S6/s0", v, "hits.npz")))
        lab, lst = S.link_families(fam, h, gt.lengths, max_size=7500, **CFG)
        devf = dev_families(gt, lab)
        qs, ts = [], []
        for mem in family_members(lab):
            if lab[mem[0]] not in devf:
                continue
            old = fam[mem]
            if len(np.unique(old)) == 1:
                continue
            i, j = np.triu_indices(len(mem), 1)
            x = old[i] != old[j]
            qs.append(mem[i][x]); ts.append(mem[j][x])
        q = np.concatenate(qs) if qs else np.zeros(0, np.int64)
        t = np.concatenate(ts) if ts else np.zeros(0, np.int64)
        wd = os.path.join(W, v)
        os.makedirs(wd, exist_ok=True)
        stats = {}
        cfg = ScoringConfig(threads=a.threads)
        with open(os.path.join(wd, "align.log"), "w") as log:
            c = _mmseqs_align_pairs(gt, q.astype(np.int64), t.astype(np.int64), cfg, wd, log, stats)
        for f in os.listdir(wd):
            if f not in ("align.log",):
                p = os.path.join(wd, f)
                if os.path.isfile(p):
                    os.remove(p)
        new = _reduce_records(c[0].astype(np.int64), c[1].astype(np.int64), *c[2:], gt.n) if len(c[0]) else None
        if new is not None:
            st = ScoreTable(*(np.concatenate([getattr(st0, k), getattr(new, k)]) for k in
                              ("q", "t", "bits", "qcov", "tcov", "pident", "evalue")), self_bits=st0.self_bits)
        else:
            st = st0
        st.save(os.path.join(wd, "scores.npz"))
        np.save(os.path.join(wd, "families.npy"), lab)
        info[v] = dict(link=lst, n_new_pairs_aligned=int(len(q)), n_new_scored=int(0 if new is None else len(new)),
                       wall_s=round(time.time() - t0, 1))
        print(v, info[v], flush=True)
    json.dump(info, open(os.path.join(ROOT, "results/S6_structure/s0_downstream_scores.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
