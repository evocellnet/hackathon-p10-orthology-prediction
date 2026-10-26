"""Experiment: repair S0 fissions with the cross-family hits seen by the S2 DIAMOND pass.

For each toy sim: default S0 families -> hybrid scoring (keeps cross-family hits) ->
union two families when a cross-family hit has min(qcov, tcov) >= COV and evalue <= E
-> rescore -> default merge. Reports S0 partition quality and ortholog-pair P/R/F1.

    python scripts/exp_family_rescue.py > results/core_v0/exp_family_rescue.tsv
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.evaluate_toy import pair_metrics, partition_metrics  # noqa: E402
from ancestra.families import read_families_tsv  # noqa: E402
from ancestra.io import families_from_mapping, load_families, load_guide_tree, load_proteomes  # noqa: E402
from ancestra.merge import MergeConfig, run_merge, to_orthoxml  # noqa: E402
from ancestra.orthoxml import hogs_to_pairs, read_orthoxml  # noqa: E402
from ancestra.scoring import ScoringConfig, score_families  # noqa: E402

SIMS = os.path.join(ROOT, "work/core_v0/sims")
RUNS = os.path.join(ROOT, "work/core_v0/runs")
VARIANTS = [("none", None, None), ("cov50_e1e-5", 0.5, 1e-5), ("cov50_e1e-3", 0.5, 1e-3),
            ("cov30_e1e-3", 0.3, 1e-3)]


def union(fam, cross, cov, ev):
    parent = np.arange(fam.max() + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    m = (np.minimum(cross["qcov"], cross["tcov"]) >= cov) & (cross["evalue"] <= ev)
    n = 0
    for a, b in zip(fam[cross["q"][m]].tolist(), fam[cross["t"][m]].tolist()):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
            n += 1
    return np.array([find(f) for f in fam]), n


def one(sim):
    d = os.path.join(SIMS, sim)
    gt = load_proteomes(os.path.join(d, "proteomes"))
    tree = load_guide_tree(os.path.join(d, "species_tree.nwk"), gt.species)
    ftrue = load_families(os.path.join(d, "truth", "gene_families.tsv"), gt, 0, 2)
    fam0 = families_from_mapping(gt, read_families_tsv(os.path.join(RUNS, sim, "s0", "s0", "families.tsv")))
    truth = set(hogs_to_pairs(read_orthoxml(os.path.join(d, "truth", "hogs.orthoxml"))))
    tf = [set(np.flatnonzero(ftrue == f)) for f in np.unique(ftrue)]
    wd = f"/tmp/ancestra_rescue_{sim}"
    st0, _ = score_families(gt, fam0, ScoringConfig(threads=4), wd)
    out = []
    for name, cov, ev in VARIANTS:
        if cov is None:
            fam, nj, st = fam0, 0, st0
        else:
            fam, nj = union(fam0, st0.cross, cov, ev)
            st, _ = score_families(gt, fam, ScoringConfig(threads=4), wd)
        pf = [set(np.flatnonzero(fam == f)) for f in np.unique(fam)]
        pm = partition_metrics(tf, pf)
        res = run_merge(gt, tree, st, MergeConfig())
        m = pair_metrics(hogs_to_pairs(to_orthoxml(res, gt, tree)), truth)
        out.append((sim, name, nj, pm["exact"], pm["fissions"], pm["fusions"], m["precision"],
                    m["recall"], m["f1"]))
    return out


def main():
    sims = sorted(os.listdir(SIMS))
    rows = []
    with ProcessPoolExecutor(6) as ex:
        for r in ex.map(one, sims):
            rows.extend(r)
    print("sim\tvariant\tn_joins\tfam_exact\tfam_fissions\tfam_fusions\tP\tR\tF1")
    for r in rows:
        print("\t".join(map(str, r[:3])) + "\t" + "\t".join(f"{x:.4f}" if isinstance(x, float) else str(x)
                                                          for x in r[3:]))
    print("\n#variant\tsetting\tfam_exact\tfissions\tfusions\tP\tR\tF1")
    for name, _, _ in VARIANTS:
        for st in ("def8", "hard8", "def32", "hard32"):
            rs = [r for r in rows if r[1] == name and r[0].startswith(st + "_")]
            print(f"#{name}\t{st}\t" + "\t".join(f"{np.mean([r[i] for r in rs]):.4f}" for i in range(3, 9)))


if __name__ == "__main__":
    main()
