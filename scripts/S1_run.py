#!/usr/bin/env python
"""S1 driver: build guide-tree arms for a dataset (qfo | ob).

Stages (each idempotent; outputs in work/S1/<ds>/ and results/S1_speciestree/<ds>/):
  ncbi     NCBI taxdump tree
  search   arm A sample + per-species DIAMOND search
  armA     besthit distances -> BIONJ / FastME -> NCBI-clade rooting (+ midpoint variant)
  armB     near-single-copy families -> FAMSA/FastTree -> ASTRAL-Pro3 (+ concatenation)
  perturb  NNI-perturbed (k in K, seeds 11,23,47) of the NCBI tree and random trees
Only declared method inputs are read: proteomes, S0 families, the NCBI taxdump.
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402

from ancestra import speciestree as st  # noqa: E402
from ancestra.io import read_fasta, species_name  # noqa: E402
from ancestra.newick import read_newick  # noqa: E402

TAXDUMP = os.path.join(ROOT, "data", "taxonomy")
DS = {
    "qfo": {"proteomes": [os.path.join(ROOT, "data/qfo/proteomes_2022_02_updZEBRAFISH")],
            "naming": "qfo", "s0": None},
    "ob": {"proteomes": [os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/Input")],
           "naming": "ensembl",
           "s0": os.path.join(ROOT, "work/S0/2sA_dc30_cc_sh30_RA/families.tsv")},
}
SEEDS = [11, 23, 47]


def load_json(p, default=None):
    return json.load(open(p)) if os.path.exists(p) else default


def save_json(obj, p):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as fh:
        json.dump(obj, fh, indent=1)
    os.replace(p + ".tmp", p)


def timed(fn):
    t0, c0 = time.time(), os.times()
    r = fn()
    c1 = os.times()
    return r, {"wall_s": round(time.time() - t0, 2),
               "cpu_self_s": round((c1.user - c0.user) + (c1.system - c0.system), 2),
               "cpu_children_s": round((c1.children_user - c0.children_user) +
                                       (c1.children_system - c0.children_system), 2),
               "loadavg_start": os.getloadavg()[0]}


def domain_of(f):
    d = os.path.basename(os.path.dirname(f))
    return d if d in ("Archaea", "Bacteria", "Eukaryota") else "Eukaryota"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ds", choices=list(DS))
    ap.add_argument("stages", nargs="+")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--n-euk", type=int, default=500)
    ap.add_argument("--n-prok", type=int, default=1000)
    ap.add_argument("--sens", default="--mid-sensitive")
    ap.add_argument("--max-fams", type=int, default=400)
    a = ap.parse_args()
    cfg = DS[a.ds]
    work = os.path.join(ROOT, "work", "S1", a.ds)
    res = os.path.join(ROOT, "results", "S1_speciestree", a.ds)
    os.makedirs(work, exist_ok=True)
    os.makedirs(res, exist_ok=True)
    files = st.list_proteomes(cfg["proteomes"])
    species = [species_name(f) for f in files]
    rt_path = os.path.join(res, "runtime.json")
    rt = load_json(rt_path, {})
    ncbi_path = os.path.join(res, "ncbi.nwk")

    for stage in a.stages:
        print(f"[{time.strftime('%H:%M:%S')}] {a.ds} {stage}", flush=True)
        if stage == "ncbi":
            def run():
                if cfg["naming"] == "qfo":
                    lt = st.qfo_leaf_taxids(files)
                else:
                    nm = st.ensembl_leaf_names(files)
                    tx = st.lookup_taxids_by_name(TAXDUMP, list(nm.values()))
                    lt = {leaf: tx[n] for leaf, n in nm.items()}
                t = st.ncbi_tree(lt, TAXDUMP)
                st.write_tree(t, ncbi_path, lengths=False)
                with open(os.path.join(res, "ncbi_nodes.tsv"), "w") as fh:
                    fh.write("node\ttaxid\trank\tn_children\tn_leaves\n")
                    for nd in t.preorder():
                        fh.write(f"{nd.name}\t{nd.data.get('taxid')}\t{nd.data.get('rank','')}\t"
                                 f"{len(nd.children)}\t{len(nd.leaves())}\n")
                save_json({"leaf_taxid": lt, "domain": {species_name(f): domain_of(f) for f in files}},
                          os.path.join(res, "species.json"))
                return t
            t, r = timed(run)
            rt["ncbi"] = r
        elif stage == "search":
            q = os.path.join(work, "sample.fa")

            def run():
                if cfg["naming"] == "qfo":
                    n = {species_name(f): (a.n_euk if domain_of(f) == "Eukaryota" else a.n_prok)
                         for f in files}
                else:
                    n = a.n_euk
                info = st.sample_queries(files, n, 11, q)
                save_json(info, os.path.join(work, "sample_info.json"))
                dbs = st.build_species_dbs(files, os.path.join(work, "db"), threads=a.threads)
                return st.search_per_species(q, dbs, os.path.join(work, "hits"), a.threads, sens=a.sens)
            s, r = timed(run)
            r["search"] = s
            rt["search"] = r
        elif stage == "armA":
            def run():
                df = st.load_hits(os.path.join(work, "hits"), len(files))
                D, C = st.besthit_distances(df, len(files))
                np.save(os.path.join(work, "armA_dist.npy"), D)
                with open(os.path.join(res, "armA_distance.tsv"), "w") as fh:
                    fh.write("species\t" + "\t".join(species) + "\n")
                    for i, sp in enumerate(species):
                        fh.write(sp + "\t" + "\t".join(f"{x:.5f}" for x in D[i]) + "\n")
                ncbi = read_newick(ncbi_path)
                out = {}
                for nm, method, spr in (("bionj", "B", False), ("fastme", "B", True)):
                    t, cmd = st.fastme_tree(D, species, os.path.join(work, "fastme_" + nm), method, spr)
                    for rooting in ("ncbi", "midpoint"):
                        tt, info = st.finalize(st.copy_tree(t), ncbi, rooting)
                        tag = f"armA_{nm}_{'ncbiroot' if rooting == 'ncbi' else 'midroot'}"
                        st.write_tree(tt, os.path.join(res, tag + ".nwk"))
                        info["cmd"] = cmd
                        out[tag] = info
                save_json(out, os.path.join(res, "armA_info.json"))
                return out
            _, r = timed(run)
            rt["armA"] = r
        elif stage == "armB":
            def run():
                seqs = {}
                gene_sp = {}
                for f in files:
                    sp = species_name(f)
                    for pid, s in read_fasta(f):
                        seqs[pid] = s
                        gene_sp[pid] = sp
                if cfg["s0"]:
                    fams0 = st.families_from_s0(cfg["s0"], gene_sp, len(files))
                    keys = sorted(fams0)
                    import random
                    random.Random(11).shuffle(keys)
                    keys = keys[:a.max_fams]
                    fams = {k: [(gene_sp[g], g) for g in fams0[k]] for k in keys}
                    finfo = {"source": "S0 families " + cfg["s0"], "n_candidates": len(fams0),
                             "s0_sha256": st.sha256_file(cfg["s0"])}
                else:
                    df = st.load_hits(os.path.join(work, "hits"), len(files))
                    fams_raw, ncand = st.families_from_hits(df, species, max_fams=a.max_fams)
                    fams = {k: [(sp, pid.split("__", 1)[1]) for sp, pid in v] for k, v in fams_raw.items()}
                    finfo = {"source": "sampled best hits (arm A search)", "n_candidates": ncand}
                occ = [len({sp for sp, _ in v}) / len(files) for v in fams.values()]
                finfo.update({"n_families": len(fams), "occupancy_mean": float(np.mean(occ)),
                              "occupancy_min": float(np.min(occ)),
                              "genes_per_family_mean": float(np.mean([len(v) for v in fams.values()]))})
                with open(os.path.join(work, "armB_families.tsv"), "w") as fh:
                    for k, v in fams.items():
                        for sp, g in v:
                            fh.write(f"{k}\t{sp}\t{g}\n")
                gdir = os.path.join(work, "genetrees")
                (gmap, paths, cpu_gt), r_gt = timed(lambda: st.gene_trees(fams, seqs, gdir, a.threads))
                r_gt["cpu_children_s"] = cpu_gt
                ncbi = read_newick(ncbi_path)
                t, ainfo = st.run_astral_pro(paths, gmap, gdir, a.threads)
                tt, info = st.finalize(t, ncbi, "ncbi")
                st.write_tree(tt, os.path.join(res, "armB_astralpro_ncbiroot.nwk"))
                out = {"families": finfo, "genetrees_runtime": r_gt, "astral": ainfo, "astral_root": info}
                try:
                    ct, cinfo = st.concat_tree(fams, gdir, gmap, a.threads, species)
                    ct2, info2 = st.finalize(ct, ncbi, "ncbi")
                    st.write_tree(ct2, os.path.join(res, "armB_concat_fasttree_ncbiroot.nwk"))
                    out["concat"] = cinfo
                    out["concat_root"] = info2
                except Exception as e:  # noqa: BLE001
                    out["concat_error"] = repr(e)
                save_json(out, os.path.join(res, "armB_info.json"))
                return out
            _, r = timed(run)
            rt["armB"] = r
        elif stage == "perturb":
            def run():
                ncbi = read_newick(ncbi_path)
                n_int = sum(1 for n in ncbi.preorder() if n.children)
                ks = sorted({max(1, round(f * n_int)) for f in (0.05, 0.10, 0.25)})
                out = {"k_values": ks, "n_internal_ncbi": n_int, "trees": {}}
                for k in ks:
                    for s in SEEDS:
                        t, done = st.nni_perturb(ncbi, k, s)
                        st.name_internal_nodes(t, ncbi)
                        p = os.path.join(res, f"perturbed_ncbi_nni{k}_s{s}.nwk")
                        st.write_tree(t, p, lengths=False)
                        out["trees"][os.path.basename(p)] = {"k": k, "seed": s, "nni_done": done}
                for s in SEEDS:
                    t = st.random_topology(sorted(ncbi.leaf_names()), s)
                    st.name_internal_nodes(t, ncbi)
                    p = os.path.join(res, f"random_s{s}.nwk")
                    st.write_tree(t, p, lengths=False)
                    out["trees"][os.path.basename(p)] = {"seed": s}
                save_json(out, os.path.join(res, "perturb_info.json"))
                return out
            _, r = timed(run)
            rt["perturb"] = r
        save_json(rt, rt_path)
        print(json.dumps(rt.get(stage), default=str)[:600], flush=True)


if __name__ == "__main__":
    main()
