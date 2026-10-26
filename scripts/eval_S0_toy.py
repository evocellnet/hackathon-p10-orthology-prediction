#!/usr/bin/env python3
"""S0 on toy simulations (ancestra.simulate): exact-truth fission / fusion check.

For each simulation dir and config: concatenate proteomes, run ancestra.families, score vs the
true root families (truth/gene_families.tsv `family`) with the same metrics as OrthoBench.
Root families are independent random sequences, so ANY fusion is a false merge (negative control),
and completeness measures fission at the simulated divergence (positive control).
Writes results/S0_families/toy/<sim>__<config>.json.

  python scripts/eval_S0_toy.py --sims work/S0/toy/h1.0 work/S0/toy/h2.0 --configs NAME [NAME ...] [--ceiling]

--ceiling also scores the full all-vs-all single-linkage reference (DIAMOND -k 0, connected components)
at (--sensitive, e<=1e-5) and (--ultra-sensitive, e<=1e-3): the best any hit-graph partition can do.
"""
import subprocess
import argparse, json, os, sys, time
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra import families as F  # noqa: E402
from ancestra_eval import families_eval as FE  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", nargs="+", required=True)
    ap.add_argument("--configs", nargs="+", required=True)
    ap.add_argument("--config-file", default=os.path.join(ROOT, "results/S0_families/configs.json"))
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--ceiling", action="store_true")
    a = ap.parse_args()
    allcfg = json.load(open(a.config_file))
    outdir = os.path.join(ROOT, "results/S0_families/toy")
    os.makedirs(outdir, exist_ok=True)
    for sim in a.sims:
        simname = os.path.basename(sim.rstrip("/"))
        fa = os.path.join(sim, "all_proteins.fa")
        if not os.path.exists(fa):
            with open(fa, "w") as fo:
                for f in sorted(os.listdir(os.path.join(sim, "proteomes"))):
                    fo.write(open(os.path.join(sim, "proteomes", f)).read())
        truth = {}
        for line in open(os.path.join(sim, "truth/gene_families.tsv")):
            t = line.rstrip("\n").split("\t")
            if t[0] == "protein":
                continue
            truth.setdefault(t[2], set()).add(t[0])
        params = json.load(open(os.path.join(sim, "params.json")))
        if a.ceiling:
            wd = os.path.join(ROOT, "work/S0/toy_runs", f"{simname}__ceiling")
            os.makedirs(wd, exist_ok=True)
            subprocess.run(["diamond", "makedb", "--in", fa, "-d", os.path.join(wd, "all"), "--quiet"], check=True)
            for sens, ev in (("--sensitive", "1e-5"), ("--ultra-sensitive", "1e-3")):
                hits = os.path.join(wd, f"ava{sens}.tsv")
                subprocess.run(["diamond", "blastp", "-q", fa, "-d", os.path.join(wd, "all"), sens, "-e", ev, "-k", "0",
                                "--threads", str(a.threads), "--quiet", "--outfmt", "6", "qseqid", "sseqid", "-o", hits],
                               check=True)
                g2f = F.label_families(F.partition_from_pairs(F.read_fasta_ids(fa), F.read_pairs_tsv(hits)))
                ag = FE.evaluate_refogs(g2f, truth)["aggregate"]
                name = f"ceiling_allvsall_cc{sens}_e{ev}"
                json.dump({"sim": sim, "config": name, "aggregate": ag, "family_stats": F.family_stats(g2f)},
                          open(os.path.join(outdir, f"{simname}__{name}.json"), "w"), indent=1)
                print(f"{simname} {name}: complete={ag['completeness']:.3f} pairR={ag['pair_recall']:.3f} "
                      f"pairP={ag['pair_precision']:.3f}", flush=True)
        for name in a.configs:
            cfg = dict(allcfg[name], name=name, threads=a.threads)
            cfg.pop("s1_reuse", None)
            wd = os.path.join(ROOT, "work/S0/toy_runs", f"{simname}__{name}")
            t0 = time.time()
            st = F.run_families(fa, wd, F.FamilyConfig(**cfg))
            g2f = F.read_families_tsv(os.path.join(wd, "families.tsv"))
            ev = FE.evaluate_refogs(g2f, truth)
            gene2true = {g: f for f, gs in truth.items() for g in gs}
            fam_true = {}
            for g, f in g2f.items():
                fam_true.setdefault(f, set()).add(gene2true[g])
            n_fused = sum(1 for s in fam_true.values() if len(s) > 1)
            res = {"sim": sim, "sim_params": {k: params.get(k) for k in ("seed", "params") if k in params},
                   "config": name, "aggregate": ev["aggregate"],
                   "n_pred_families_with_gt1_true_family": n_fused,
                   "family_stats": {k: v for k, v in st.items() if k not in ("config", "timings")},
                   "timings": st["timings"], "wall_s": time.time() - t0}
            json.dump(res, open(os.path.join(outdir, f"{simname}__{name}.json"), "w"), indent=1)
            ag = ev["aggregate"]
            print(f"{simname} {name}: true_fams={ag['n_refogs']} complete={ag['completeness']:.3f} "
                  f"mean_n_fam={ag['mean_n_fam']:.2f} fused_pred_fams={n_fused} pairR={ag['pair_recall']:.3f} "
                  f"pairP={ag['pair_precision']:.3f}", flush=True)


if __name__ == "__main__":
    main()
