"""L2-domains: write the protein sets to annotate with Pfam (dev scope only).

  ob_dev_fam.fa   : every protein of the S0 families (run_det families_final) that hold a DEV RefOG gene
  qfo_tf_dev.fa   : the TreeFam-A DEV reference genes (standalone signal)
  qfo_tf_dev_fam.fa: every protein of the QfO S0 families (qfo_recall/base) that hold a TreeFam-A DEV gene

    python scripts/L2d_prepare_fasta.py OUTDIR
Benchmark-scoped helper (reads dev labels to choose *which proteins to annotate*; the annotation
itself is label-free and per protein).
"""
import os
import sys

import numpy as np

ROOT = os.environ.get("P10_ROOT", "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import load_proteomes, read_gene_table  # noqa: E402

OB_RUN = os.path.join(ROOT, "work/core_v0_ob/run_det")
QFO_RUN = os.path.join(ROOT, "work/qfo_recall/base")


def write(path, names, seqs):
    with open(path, "w") as fh:
        for n, s in zip(names, seqs):
            fh.write(f">{n}\n{s}\n")


def main(out):
    os.makedirs(out, exist_ok=True)
    from orthobench_dev_pairs_eval import load_truth
    t = load_truth()
    genes = set(t.gene_a) | set(t.gene_b)
    gt = load_proteomes(os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/Input"))
    g2 = read_gene_table(os.path.join(OB_RUN, "stages/input/genes.tsv"))
    assert gt.names == g2.names
    fam = np.load(os.path.join(OB_RUN, "stages/scores/families_final.npy"))
    idx = gt.index()
    fs = {fam[idx[g]] for g in genes}
    sel = [i for i in range(gt.n) if fam[i] in fs]
    write(os.path.join(out, "ob_dev_fam.fa"), [gt.names[i] for i in sel], [gt.seqs[i] for i in sel])
    print("ob", len(sel))
    from qfo_treefam_dev_fast import dev_refpairs
    qg = {x for r in dev_refpairs() for x in r[:2]}
    gq = load_proteomes(os.path.join(ROOT, "data/inputs/qfo2022/proteomes"))
    g3 = read_gene_table(os.path.join(QFO_RUN, "stages/input/genes.tsv"))
    assert gq.names == g3.names
    fq = np.load(os.path.join(QFO_RUN, "stages/scores/families_final.npy"))
    iq = gq.index()
    selq = sorted(iq[g] for g in qg)
    write(os.path.join(out, "qfo_tf_dev.fa"), [gq.names[i] for i in selq], [gq.seqs[i] for i in selq])
    fsq = {fq[i] for i in selq}
    selq2 = [i for i in range(gq.n) if fq[i] in fsq]
    write(os.path.join(out, "qfo_tf_dev_fam.fa"), [gq.names[i] for i in selq2], [gq.seqs[i] for i in selq2])
    print("qfo", len(selq), len(selq2))


if __name__ == "__main__":
    main(sys.argv[1])
