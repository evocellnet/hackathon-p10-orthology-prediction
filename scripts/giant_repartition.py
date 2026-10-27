"""Re-derive the S0 stage-2 partition from a finished S0 run's cached files (stage1_clusters.tsv +
reps_hits.tsv) with a different Leiden cap (s2_gmax), without re-searching. Replicates
ancestra.families.run_families after the stage-2 search exactly (verified: gmax = the run's own value
reproduces its families.tsv partition).

    python scripts/giant_repartition.py RUN_DIR GMAX OUT.npy [--threads]

RUN_DIR = pipeline run dir (stages/input/genes.tsv + stages/families/s0/). Writes the gene -> family int
array (gene-table order), plus OUT.json with stats (n families, split families, sizes).
"""
import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra.families import (FamilyConfig, best_hits, filter_hits, label_families, leiden_split,  # noqa
                               partition_from_pairs, read_pairs_tsv)
from ancestra.io import read_fasta  # noqa


def canonical_names(run):
    """canonical id -> gene index, rebuilt from the proteome sequences (families._canonical_fasta order)."""
    from ancestra.io import load_proteomes, read_gene_table
    man = json.load(open(os.path.join(run, "manifest.json"))) if os.path.exists(os.path.join(run, "manifest.json")) else {}
    gt = read_gene_table(os.path.join(run, "stages/input/genes.tsv"))
    return gt


def main():
    run, gmax, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    prot = sys.argv[4]
    from ancestra.io import load_proteomes, read_gene_table
    t0 = time.time()
    gt = read_gene_table(os.path.join(run, "stages/input/genes.tsv"))
    full = load_proteomes(prot)
    assert full.names == gt.names
    seqs = full.seqs
    order = sorted(range(gt.n), key=lambda g: (-len(seqs[g]), seqs[g], gt.names[g]))
    width = max(7, len(str(gt.n)))
    cid = [f"c{i:0{width}d}" for i in range(gt.n)]          # canonical id of rank i
    canon2gene = {cid[i]: order[i] for i in range(gt.n)}
    s0 = os.path.join(run, "stages/families/s0")
    cfg = FamilyConfig(**{k: v for k, v in json.load(open(os.path.join(s0, "stats.json")))["config"].items()})
    gene_ids = sorted(canon2gene)
    pairs = read_pairs_tsv(os.path.join(s0, "stage1_clusters.tsv"))
    g2r = partition_from_pairs(gene_ids, pairs)
    rep_of = {m: r for r, m in pairs}
    root2rep = {}
    for g in gene_ids:
        if rep_of.get(g) == g:
            root2rep[g2r[g]] = g
    for g in gene_ids:
        root2rep.setdefault(g2r[g], g2r[g])
    rep_set = set(root2rep.values())
    gene2rep = {g: root2rep[g2r[g]] for g in gene_ids}
    rep_size = Counter(gene2rep.values())
    hits = os.path.join(s0, "reps_hits.tsv")
    with open(hits) as fh:
        gene_edges = filter_hits((l.rstrip("\n").split("\t") for l in fh), cfg)
    best = {}
    for a, b, w in gene_edges:
        ra, rb = gene2rep[a], gene2rep[b]
        if ra == rb:
            continue
        k = (ra, rb) if ra < rb else (rb, ra)
        if w > best.get(k, -1.0):
            best[k] = w
    edges = sorted((a, b, w) for (a, b), w in best.items())
    rep2root = partition_from_pairs(sorted(rep_set), [(a, b) for a, b, _ in edges])
    comp = defaultdict(list)
    for r, root in rep2root.items():
        comp[root].append(r)
    big = {root for root, rs in comp.items() if sum(rep_size[r] for r in rs) > gmax}
    flagged = {}
    t1 = time.time()
    if big:
        by_comp = defaultdict(list)
        for a, b, w in edges:
            if rep2root[a] in big and rep2root[a] == rep2root[b]:
                by_comp[rep2root[a]].append((a, b, w))
        for root in sorted(big):
            parts = leiden_split(sorted(comp[root]), by_comp[root], rep_size, gmax, cfg.seed)
            for part in parts:
                proot = min(part)
                for r in part:
                    rep2root[r] = proot
                    flagged[r] = root
    leiden_s = time.time() - t1
    g2r = {g: rep2root[gene2rep[g]] for g in gene_ids}
    n_moved = 0
    if cfg.s2_reassign and flagged:
        with open(hits) as fh:
            bh = best_hits((l.rstrip("\n").split("\t") for l in fh), cfg)
        new = dict(g2r)
        for g, (s_, _) in bh.items():
            if gene2rep[g] in flagged and s_ in flagged and g2r[s_] != g2r[g]:
                new[g] = g2r[s_]
                n_moved += 1
        g2r = new
    lab = label_families(g2r)
    fam = np.empty(gt.n, np.int64)
    comp_of = np.full(gt.n, -1, np.int64)          # parent component id (flagged) per gene, -1 if not split
    croot = {r: i for i, r in enumerate(sorted(big))}
    for g, f in lab.items():
        fam[canon2gene[g]] = int(f[1:])
        r = gene2rep[g]
        if r in flagged:
            comp_of[canon2gene[g]] = croot[flagged[r]]
    np.save(out, fam)
    np.save(out.replace(".npy", "_component.npy"), comp_of)
    sizes = np.bincount(fam)
    st = dict(gmax=gmax, n_families=int((sizes > 0).sum()), max_family=int(sizes.max()),
              n_split_parts=int(len(np.unique(fam[comp_of >= 0]))), genes_in_split=int((comp_of >= 0).sum()),
              n_moved=n_moved, sum_g2=int((sizes.astype(np.int64) ** 2).sum()),
              n_fam_gt5000=int((sizes > 5000).sum()), leiden_s=round(leiden_s, 1), wall_s=round(time.time() - t0, 1))
    json.dump(st, open(out.replace(".npy", ".json"), "w"), indent=1)
    print(json.dumps(st))


if __name__ == "__main__":
    main()
