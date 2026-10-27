"""TreeFam-A DEV eval (fast in-memory scorer; dev families only) of one QfO rescue arm run dir,
plus the missed-pair breakdown (S0 split / different root HOG / paralog call) by LCA level.
Appends to results/rescue/qfo_runs.tsv; writes RUN/tf_dev_rescue.json.

    python scripts/rescue_eval_qfo.py work/rescue/qfo/R1 --label R1
"""
import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from ancestra.io import load_guide_tree, read_gene_table  # noqa: E402
from ancestra.orthoxml import read_orthoxml  # noqa: E402
from qfo_merge_iter import fn_breakdown  # noqa: E402
from qfo_treefam_dev_fast import score as tf_score  # noqa: E402

LOG = os.path.join(ROOT, "results/rescue/qfo_runs.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    gt = read_gene_table(os.path.join(a.run, "stages/input/genes.tsv"))
    tree = load_guide_tree(os.path.join(a.run, "stages/input/guide_tree.nwk"), gt.species)
    famf = os.path.join(a.run, "stages/rescue/families_final.npy")
    fam = np.load(famf if os.path.exists(famf) else os.path.join(a.run, "stages/scores/families_final.npy"))
    ox = read_orthoxml(os.path.join(a.run, "hogs.orthoxml"))
    tf = tf_score(ox)
    fnb = fn_breakdown(ox, gt, fam, tree)
    man = json.load(open(os.path.join(a.run, "manifest.json")))
    st = man["stages"]
    rs = st.get("rescue", {}).get("stats", {})
    res = dict(label=a.label, TPR=round(tf["TPR"], 4), PPV=round(tf["PPV"], 4),
               **{k: tf[k] for k in ("n_TP", "n_FP", "n_FN", "n_TN")},
               fn_classes={k: sum(v.values()) for k, v in fnb.items()},
               s0_fn_opisthokonta=fnb.get("i_diff_family", {}).get("Opisthokonta"),
               rescue_s=st.get("rescue", {}).get("wall_s"), merge_s=st["merge"]["wall_s"],
               rescue_attach=rs.get("attach"), rescue_rescore=rs.get("rescore"),
               rescue_search_s={k: rs.get(k) for k in ("search_3di_s", "search_profile_s")},
               guard=st["merge"].get("stats", {}).get("rescue_guard"))
    json.dump(dict(res, fn_breakdown=fnb, tf=tf), open(os.path.join(a.run, "tf_dev_rescue.json"), "w"),
              indent=1, default=str)
    cols = ["time", "label", "TPR", "PPV", "n_TP", "n_FP", "n_FN", "n_TN", "fn_classes", "s0_fn_opisthokonta",
            "rescue_s", "merge_s", "rescue_attach", "rescue_rescore", "rescue_search_s", "guard"]
    row = dict(time=datetime.now().isoformat(timespec="seconds"), **res)
    new = not os.path.exists(LOG)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a") as fh:
        if new:
            fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(json.dumps(row[c], default=str) if isinstance(row[c], dict) else str(row[c])
                           for c in cols) + "\n")
    print(json.dumps(res, default=str))


if __name__ == "__main__":
    main()
