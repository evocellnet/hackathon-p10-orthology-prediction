"""Verifier: S0 completeness = share of RefOGs (low-certainty genes removed) whose genes all lie in one family;
genes absent from families.tsv count as singletons. Independent of ancestra_eval."""
import json, os
ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
B = f"{ROOT}/data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"
split = json.load(open(f"{ROOT}/splits/orthobench_split.json"))
g2f = {}
for l in open(f"{ROOT}/work/S0/2sA_dc30_cc_sh30_RA/families.tsv"):
    if l.startswith("#"): continue
    g, f = l.split()[:2]; g2f[g] = f
out = {}
for sub in ("dev", "test"):
    for keep_unc in (False, True):
        c = 0
        for r in split[sub]:
            R = {x.strip() for x in open(f"{B}/{r}.txt") if x.strip()}
            u = f"{B}/low_certainty_assignments/{r}.txt"
            if not keep_unc and os.path.exists(u):
                R -= {x.strip() for x in open(u) if x.strip()}
            c += len({g2f.get(g, "single:" + g) for g in R}) == 1
        out[f"{sub}{'_keepunc' if keep_unc else ''}"] = [c, len(split[sub]), round(c / len(split[sub]), 4)]
print(out)
json.dump(out, open(f"{ROOT}/results/verification/s0_completeness.json", "w"), indent=1)
