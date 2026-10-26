#!/usr/bin/env python
"""results/S7_isoforms/manifest.json and results/S9_expression/manifest.json (inputs, code, outputs md5)."""
import glob, hashlib, json, os, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel(ps):
    return {os.path.relpath(p, ROOT): md5(p) for p in sorted(ps) if os.path.isfile(p)}


common_code = [f"{ROOT}/scripts/S5_pair_benchmark.py", f"{ROOT}/scripts/S5_refog_labels.py",
               f"{ROOT}/tests/test_expression_isoforms.py", f"{ROOT}/scripts/S7_S9_manifest.py"]
common_in = [f"{ROOT}/splits/orthobench_split.json", f"{ROOT}/splits/orthobench_refog_components.tsv",
             f"{ROOT}/work/S5/refog_pair_labels.parquet"]
S7 = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), git_commit="not a git repo",
          code=rel([f"{ROOT}/ancestra/isoforms.py", f"{ROOT}/scripts/S7_isoform_prep.py",
                    f"{ROOT}/scripts/S7_isoform_benchmark.py"] + common_code),
          inputs=rel(common_in + glob.glob(f"{ROOT}/data/orthobench/ensembl_99/pep/*.gz")
                     + glob.glob(f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/Input/*.fa")
                     + glob.glob(f"{ROOT}/work/S5/pair_scores_*.parquet")),
          derived=rel([f"{ROOT}/work/S7/isoforms.parquet", f"{ROOT}/work/S7/isoforms.fa", f"{ROOT}/work/S7/hits.tsv",
                       f"{ROOT}/work/S7/prep_stats.json", f"{ROOT}/work/S7/pair_scores.parquet",
                       f"{ROOT}/work/S7/query_rows.parquet"]),
          outputs=rel(glob.glob(f"{ROOT}/results/S7_isoforms/*.json")),
          commands=["python scripts/S7_isoform_prep.py   (DIAMOND 6 threads)", "python scripts/S7_isoform_benchmark.py",
                    "python scripts/S7_S9_manifest.py"],
          diamond="diamond blastp --sensitive -k 0 --max-hsps 1 -e 1e-3 -p 6, all-vs-all over 1,993 isoforms of 957 DEV-RefOG genes")
S9 = dict(created=S7["created"], git_commit="not a git repo",
          code=rel([f"{ROOT}/ancestra/expression.py", f"{ROOT}/scripts/S9_expression_benchmark.py"] + common_code),
          inputs=rel(common_in + glob.glob(f"{ROOT}/data/extra/bgee_v15_2/*.tsv.gz")
                     + glob.glob(f"{ROOT}/work/S5/pair_scores_Homo_sapiens__*.parquet") + [f"{ROOT}/work/S7/hits.tsv"]
                     + glob.glob(f"{ROOT}/data/genomes/orthobench/*.parquet")),
          derived=rel(glob.glob(f"{ROOT}/work/S9/pair_expr_*.parquet")),
          outputs=rel(glob.glob(f"{ROOT}/results/S9_expression/*.json")),
          commands=["python scripts/S9_expression_benchmark.py  (needs work/S7/hits.tsv)", "python scripts/S7_S9_manifest.py"])
for d, m in (("S7_isoforms", S7), ("S9_expression", S9)):
    m["outputs"].pop(f"results/{d}/manifest.json", None)
    json.dump(m, open(f"{ROOT}/results/{d}/manifest.json", "w"), indent=1)
    print(d, len(m["code"]), len(m["inputs"]), len(m["derived"]), len(m["outputs"]))
