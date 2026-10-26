#!/usr/bin/env python
"""Markdown tables for docs/S5_synteny.md, generated from results/S5_synteny/*.json."""
import glob, json, os, sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = f"{ROOT}/results/S5_synteny"
out = []


def md(df):
    return df.to_markdown(index=False)


def pair_tables(path, title):
    d = json.load(open(path))
    rows = []
    for p in d["pairs"]:
        s = p["summary"]
        if s["seq"]["auc_all"] is None:
            continue
        r = dict(pair=p["pair"].replace("__", " / "), orth=p["n_orth"], par=p["n_par"], hard=p["n_hard_pairs"],
                 anchors=p["n_anchors"], seq_auc=s["seq"]["auc_all"], seq_auc_hard=s["seq"]["auc_hard"])
        for w in d["windows"]:
            x = s[f"w{w}"]
            r[f"w{w} AUC"] = x["auc_all"]; r[f"w{w} AUC hard"] = x["auc_hard"]
            r[f"w{w} shuf"] = x["shuf_auc_all"]; r[f"w{w} shuf hard"] = x["shuf_auc_hard"]
        rows.append(r)
    out.append(f"### {title}: ROC-AUC per species pair (synteny score vs ortholog/paralog label)\n")
    out.append(md(pd.DataFrame(rows)))
    pool = d["pooled"]
    rows = [dict(score="sequence (bits)", auc_all=pool["seq"]["auc_all"], ci_all=pool["seq"].get("auc_all_ci95"),
                 auc_hard=pool["seq"]["auc_hard"], ci_hard=pool["seq"].get("auc_hard_ci95"))]
    for w in d["windows"]:
        x = pool[f"syn_w{w}"]
        rows.append(dict(score=f"synteny w={w}", auc_all=x["auc_all"], ci_all=x.get("auc_all_ci95"),
                         auc_hard=x["auc_hard"], ci_hard=x.get("auc_hard_ci95"), ap_all=x["ap_all"],
                         ap_hard=x["ap_hard"], shuf_auc_all=x["shuf_auc_all"], shuf_auc_hard=x["shuf_auc_hard"],
                         shuf_ap_hard=x["shuf_ap_hard"]))
    out.append(f"\n\n### {title}: pooled over pairs (n={pool['n']}, orthologs {pool['n_orth']}, hard {pool['n_hard']}; "
               "95% CI = bootstrap over RefOGs/families, B=1000)\n")
    out.append(md(pd.DataFrame(rows)))
    for key, nm in (("pooled_hard_queries", "hard queries"), ("pooled_all_queries", "all queries with >=1 ortholog and >=1 paralog candidate")):
        q = d.get(key, {})
        rows = []
        for col in sorted(q):
            if col.startswith("shuf") and not col.startswith("shuf11"):
                continue
            v = q[col]
            rows.append(dict(score=col, n=v["n"], seq_top1=v["seq_acc"], syn_correct=v["syn_correct"],
                             syn_wrong=v["syn_wrong"], syn_abstain=v["syn_abstain"],
                             syn_precision=v["syn_precision_when_decides"],
                             **{"comb@0.9": v["comb0.9_acc"], "comb@0.9-seq CI": v["comb0.9_minus_seq_ci95"],
                                "comb@0.75": v["comb0.75_acc"], "comb@0.75-seq CI": v["comb0.75_minus_seq_ci95"],
                                "ovr2@0.9": v.get("ovr2@0.9_acc"), "ovr2@0.75": v.get("ovr2@0.75_acc")}))
        out.append(f"\n\n### {title}: top-1 partner choice, {nm}\n")
        out.append(md(pd.DataFrame(rows)))
    rows = []
    for p in d["pairs"]:
        q = p["queries"].get("syn_w10", {})
        if not q.get("n"):
            continue
        rows.append(dict(pair=p["pair"].replace("__", " / "), hard_queries=q["n"], seq_correct=q["seq_correct"],
                         syn_correct=q["syn_correct"], syn_wrong=q["syn_wrong"], syn_abstain=q["syn_abstain"],
                         tandem_competitor=q["n_tandem_competitor"], abstain_tandem=q["syn_abstain_tandem"],
                         **{"comb@0.75 correct": q["comb0.75_correct"]}))
    out.append(f"\n\n### {title}: hard queries per pair (w=10)\n")
    out.append(md(pd.DataFrame(rows)))
    out.append("\n")


def ancestral(path):
    d = json.load(open(path))
    rows = []
    for r in d["nodes"]:
        x = dict(node=r["node"], depth=r["depth"], cand=r["n_cand"], orth=r["n_orth"], hard=r["n_hard"],
                 anchors=r["n_anchors"], anchor_prec=round(r["anchor_precision"], 3) if r["anchor_precision"] is not None else None,
                 seq_auc=r["seq_auc"])
        for w in d["windows"]:
            y = r[f"w{w}"]
            x.update({f"w{w} AUC": y["auc"], f"w{w} AUC hard": y["auc_hard"], f"w{w} oracle": y["auc_oracle"],
                      f"w{w} shuf": y["shuf_auc"], f"w{w} orth sup>0": y["frac_orth_sup"],
                      f"w{w} par sup>0": y["frac_par_sup"]})
        rows.append(x)
    out.append(f"### Ancestral context_support on toy `{os.path.basename(d['sim'])}` "
               f"(wgd={d['params']['params'].get('wgd_branch') if 'params' in d['params'] else ''}), per node\n")
    out.append(md(pd.DataFrame(rows)))
    out.append("\n")


def preview(path):
    d = json.load(open(path))
    rows = [dict(arm=k, **{kk: (vv if not isinstance(vv, list) else ", ".join(f"{z:+.4f}" for z in vv))
                         for kk, vv in v.items()}) for k, v in d["summary"].items()]
    out.append("### Merge integration preview on toy sims (unmodified core merge + SyntenyContext)\n")
    out.append(md(pd.DataFrame(rows)))
    out.append("\n")


if os.path.exists(f"{R}/pair_benchmark_orthobench_dev.json"):
    pair_tables(f"{R}/pair_benchmark_orthobench_dev.json", "OrthoBench dev")
for p in sorted(glob.glob(f"{R}/pair_benchmark_toy_*.json")):
    pair_tables(p, f"Toy ({os.path.basename(p)[19:-5]})")
for p in sorted(glob.glob(f"{R}/ancestral_toy_*.json")):
    ancestral(p)
if os.path.exists(f"{R}/merge_preview_toy.json"):
    preview(f"{R}/merge_preview_toy.json")
open(f"{R}/tables.md", "w").write("\n".join(out))
print("\n".join(out)[:3000])
