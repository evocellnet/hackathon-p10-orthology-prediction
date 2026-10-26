#!/usr/bin/env python
"""results/S5_synteny/manifest.json: inputs (md5), code (md5), outputs, commands."""
import glob, hashlib, json, os, subprocess, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel(ps):
    return {os.path.relpath(p, ROOT): md5(p) for p in sorted(ps) if os.path.isfile(p)}


try:
    commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
except Exception:
    commit = None
m = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), git_commit=commit or "not a git repo",
         code=rel([f"{ROOT}/ancestra/genomes.py", f"{ROOT}/ancestra/synteny.py", f"{ROOT}/ancestra/genomes_qfo.py"]
                  + glob.glob(f"{ROOT}/scripts/S5_*.py") + glob.glob(f"{ROOT}/scripts/S5_*.sh")
                  + [f"{ROOT}/scripts/build_gene_order.sh", f"{ROOT}/scripts/build_gene_order_orthobench.py",
                     f"{ROOT}/scripts/summarize_S5.py", f"{ROOT}/tests/test_synteny.py"]),
         inputs=rel([f"{ROOT}/splits/orthobench_split.json", f"{ROOT}/splits/orthobench_refog_components.tsv",
                     f"{ROOT}/work/S0/2sA_dc30_cc_sh30_RA/families.tsv"]
                    + glob.glob(f"{ROOT}/data/orthobench/ensembl_99/gff3/*.gff3.gz")),
         derived=rel(glob.glob(f"{ROOT}/data/genomes/orthobench/*") + glob.glob(f"{ROOT}/data/genomes/qfo2022/coverage.tsv")
                     + [f"{ROOT}/work/S5/refog_pair_labels.parquet"] + glob.glob(f"{ROOT}/work/S5/diamond/*.tsv")),
         outputs=rel(glob.glob(f"{ROOT}/results/S5_synteny/*.json") + glob.glob(f"{ROOT}/results/S5_synteny/*.md") + glob.glob(f"{ROOT}/results/S5_synteny/phase2/*.json") + glob.glob(f"{ROOT}/results/S5_synteny/phase2/*.md") + glob.glob(f"{ROOT}/results/S5_synteny/phase2/arms/*")),
         commands=["scripts/build_gene_order.sh", "scripts/S5_diamond_pairs.sh (THREADS=8)",
                   "python scripts/S5_refog_labels.py", "python scripts/S5_pair_benchmark.py",
                   "python scripts/S5_pair_benchmark.py --toy --threads 4",
                   "python scripts/S5_ancestral_toy.py  (wgd8_s11)",
                   "python scripts/S5_ancestral_toy.py --name def8_s23 --seed 23 --wgd none --dup 0.2",
                   "python scripts/S5_merge_preview.py --workers 4", "python scripts/summarize_S5.py", "python scripts/orthobench_dev_pairs_eval.py truth work/core_v0_ob/run_v2/pairs.tsv", "python scripts/S5_phase2_arms.py --workers 4 --skip-done", "python scripts/S5_phase2_analyze.py"],
         diamond="diamond blastp --sensitive -k 100 --max-hsps 1 -e 1e-3 (whole proteomes, both directions)")
m["outputs"].pop("results/S5_synteny/manifest.json", None)
json.dump(m, open(f"{ROOT}/results/S5_synteny/manifest.json", "w"), indent=1)
print(len(m["code"]), len(m["inputs"]), len(m["derived"]), len(m["outputs"]))
