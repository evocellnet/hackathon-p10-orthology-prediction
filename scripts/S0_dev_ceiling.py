#!/usr/bin/env python3
"""Upper bound for hit-graph S0 partitions on the OrthoBench DEV RefOGs.

1. Single-linkage ceiling: connected components of stage-1 clusters + all N x R hits (e <= 1e-5,
   no coverage filter, no size cap) from an existing 2sA run directory.
2. Stray probe: every DEV RefOG gene that is outside its RefOG's majority component is searched
   against ALL 251,378 proteins with DIAMOND --ultra-sensitive (e <= 1, -k 2000); we record the
   best hit to any other member of the same RefOG. A RefOG is "unreachable" when some stray has no
   hit with e <= 1e-3 to a member of its majority component.
Writes results/S0_families/dev_ceiling.json.

  python scripts/S0_dev_ceiling.py --run work/S0/2sA_dc30_cc_nocov
"""
import argparse, json, os, subprocess, sys, tempfile
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ancestra import families as F  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--run", default=os.path.join(ROOT, "work/S0/2sA_dc30_cc_nocov"))
ap.add_argument("--fasta", default=os.path.join(ROOT, "work/S0/ob_all.fa"))
ap.add_argument("--threads", type=int, default=24)
a = ap.parse_args()

split = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
R = F.read_refogs(os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"), split["dev"])
ids = F.read_fasta_ids(a.fasta)
edges = F.read_pairs_tsv(os.path.join(a.run, "stage1_clusters.tsv"))
for line in open(os.path.join(a.run, "reps_hits.tsv")):
    t = line.split("\t")
    if float(t[2]) <= 1e-5:
        edges.append((t[0], t[1]))
g2c = F.label_families(F.partition_from_pairs(ids, edges))
ev = F.evaluate_refogs(g2c, R)
stats = F.family_stats(g2c)
strays = {}
for rid, genes in R.items():
    best = Counter(g2c[g] for g in genes).most_common(1)[0][0]
    for g in genes:
        if g2c[g] != best:
            strays[g] = (rid, best)
probe = {}
with tempfile.TemporaryDirectory(dir=os.path.join(ROOT, "work/S0")) as td:
    F._write_subset_fasta(a.fasta, set(strays), os.path.join(td, "q.fa"))
    subprocess.run(["diamond", "makedb", "--in", a.fasta, "-d", os.path.join(td, "all"), "--quiet"], check=True)
    subprocess.run(["diamond", "blastp", "-q", os.path.join(td, "q.fa"), "-d", os.path.join(td, "all"),
                    "--ultra-sensitive", "-e", "1", "-k", "2000", "--threads", str(a.threads), "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "evalue", "bitscore", "-o", os.path.join(td, "h.tsv")],
                   check=True)
    hits = defaultdict(list)
    for line in open(os.path.join(td, "h.tsv")):
        q, s, e, b = line.rstrip("\n").split("\t")
        if q != s:
            hits[q].append((float(e), s))
for g, (rid, best) in strays.items():
    to_main = sorted(e for e, s in hits[g] if g2c.get(s) == best and s in R[rid])
    probe[g] = {"refog": rid, "best_evalue_to_majority_refog_members": to_main[0] if to_main else None}
unreachable = sorted({v["refog"] for v in probe.values()
                      if v["best_evalue_to_majority_refog_members"] is None
                      or v["best_evalue_to_majority_refog_members"] > 1e-3})
out = {"run": os.path.relpath(a.run, ROOT),
       "single_linkage_ceiling": {"completeness": ev["aggregate"]["completeness"],
                                  "n_complete": ev["aggregate"]["n_complete"], "n_refogs": ev["aggregate"]["n_refogs"],
                                  "max_component": stats["max_family"], "sum_g2_over_N2": stats["sum_g2_over_N2"],
                                  "incomplete": {k: {"size": v["size"], "n_fam": v["n_fam"]}
                                                 for k, v in ev["per_refog"].items() if not v["complete"]}},
       "stray_probe": probe, "unreachable_refogs_at_e1e-3": unreachable,
       "max_reachable_completeness": 1 - len(unreachable) / len(R)}
json.dump(out, open(os.path.join(ROOT, "results/S0_families/dev_ceiling.json"), "w"), indent=1)
print(json.dumps({k: out[k] for k in ("single_linkage_ceiling", "unreachable_refogs_at_e1e-3",
                                      "max_reachable_completeness")}, indent=1))
