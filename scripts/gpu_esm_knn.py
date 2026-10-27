#!/usr/bin/env python
"""v2-gpu task (e): ESM-2 mean embeddings + FAISS-GPU kNN as an S0 candidate generator (DEV only).

    python scripts/gpu_esm_knn.py --model facebook/esm2_t12_35M_UR50D --out work/gpuv2/esm_knn [--k 50]

Runs in conda_envs/gpu. Embeds every OrthoBench protein (fp16, sequences > 1022 aa as mean of 1022-aa windows),
L2-normalises, exact inner-product kNN with faiss-gpu, then measures on the 34 DEV RefOGs against the frozen S0
partition (work/ob_final/main_accurate_go):
  - stray genes = DEV RefOG genes outside the family holding most of their RefOG;
  - reach@k = strays with >= 1 of their k nearest neighbours in that majority family
    (sequence-search reference: S0 ceiling 27-28/34, docs/S0_families.md);
  - specificity: share of kNN neighbours inside the query's own S0 family (all genes), and for DEV RefOG genes the
    share of neighbours in a *different* DEV RefOG (fusion risk).
Test RefOGs are never read.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def embed(seqs, model_name, max_tokens=60000, win=1022):
    import torch
    from transformers import AutoTokenizer, EsmModel
    tok = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name, add_pooling_layer=False).half().cuda().eval()
    # windows of <= win residues; a protein's embedding = length-weighted mean of its windows
    pieces, owner, wlen = [], [], []
    for i, s in enumerate(seqs):
        for a in range(0, len(s), win):
            pieces.append(s[a:a + win])
            owner.append(i)
            wlen.append(len(s[a:a + win]))
    order = np.argsort([-len(p) for p in pieces], kind="stable")
    d = model.config.hidden_size
    out = np.zeros((len(pieces), d), np.float32)
    t0 = time.time()
    k = 0
    with torch.inference_mode():
        while k < len(order):
            L = len(pieces[order[k]]) + 2
            bs = max(1, max_tokens // L)
            idx = order[k:k + bs]
            enc = tok([pieces[j] for j in idx], return_tensors="pt", padding=True)
            enc = {x: y.cuda() for x, y in enc.items()}
            h = model(**enc).last_hidden_state.float()
            m = enc["attention_mask"].clone()
            m[:, 0] = 0                                            # drop <cls>
            lens = enc["attention_mask"].sum(1) - 1                # position of <eos>
            m[torch.arange(len(idx)), lens] = 0                    # drop <eos>
            e = (h * m[..., None]).sum(1) / m.sum(1, keepdim=True)
            out[idx] = e.cpu().numpy()
            k += bs
    torch.cuda.synchronize()
    t_emb = time.time() - t0
    owner = np.array(owner)
    wl = np.array(wlen, np.float32)
    E = np.zeros((len(seqs), d), np.float32)
    np.add.at(E, owner, out * wl[:, None])
    E /= np.bincount(owner, weights=wl, minlength=len(seqs))[:, None]
    E /= np.linalg.norm(E, axis=1, keepdims=True)
    return E, t_emb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="facebook/esm2_t12_35M_UR50D")
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=50)
    ap.add_argument("--run", default=os.path.join(ROOT, "../ob_final/main_accurate_go"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    import pandas as pd
    genes = pd.read_csv(os.path.join(a.run, "stages/input/genes.tsv"), sep="\t")
    names = genes["protein"].tolist()
    seqs = {}
    cur = None
    for line in open(os.path.join(ROOT, "../S0/ob_all.fa")):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
    seqs = ["".join(seqs[n]).rstrip("*") for n in names]
    res = {"model": a.model, "n": len(seqs), "residues": int(sum(map(len, seqs)))}
    emb_path = os.path.join(a.out, "emb_" + a.model.split("/")[-1] + ".npy")
    if os.path.exists(emb_path):
        E = np.load(emb_path)
    else:
        E, t = embed(seqs, a.model)
        res["embed_s"] = round(t, 1)
        res["residues_per_s"] = round(res["residues"] / t)
        np.save(emb_path, E)
    import faiss
    t0 = time.time()
    resg = faiss.StandardGpuResources()
    index = faiss.GpuIndexFlatIP(resg, E.shape[1])
    index.add(E)
    D, I = index.search(E, a.k + 1)
    res["knn_s"] = round(time.time() - t0, 1)
    # determinism of the kNN: rerun and compare
    D2, I2 = index.search(E, a.k + 1)
    res["knn_repeat_identical"] = bool((I == I2).all() and (D == D2).all())
    fam = np.load(os.path.join(a.run, "stages/families/families.npy"))
    # drop self (first hit unless a tie)
    nb = np.where(I[:, :1] == np.arange(len(E))[:, None], I[:, 1:], I[:, :-1])
    from ancestra_eval import families_eval as FE
    split = json.load(open(os.path.join(ROOT, "splits/orthobench_split.json")))
    R = FE.read_refogs(os.path.join(ROOT, "data/orthobench/Open_Orthobench/BENCHMARKS/RefOGs"), split["dev"])
    gid = {n: i for i, n in enumerate(names)}
    g2ref = {}
    for rid, gs in R.items():
        for g in gs:
            if g in gid:
                g2ref[gid[g]] = rid
    rows = {}
    for kk in (1, 5, 10, 25, a.k):
        n_stray = n_reach = 0
        complete_after = 0
        for rid, gs in R.items():
            ii = [gid[g] for g in gs if g in gid]
            fams, cnt = np.unique(fam[ii], return_counts=True)
            main = fams[np.argmax(cnt)]
            strays = [i for i in ii if fam[i] != main]
            reach = [i for i in strays if (fam[nb[i, :kk]] == main).any()]
            n_stray += len(strays)
            n_reach += len(reach)
            complete_after += int(len(reach) == len(strays))
        same_fam = float((fam[nb[:, :kk]] == fam[:, None]).mean())
        dev = np.array(sorted(g2ref))
        other_ref = np.mean([[g2ref.get(int(j), None) not in (None, g2ref[i]) for j in nb[i, :kk]] for i in dev])
        rows[kk] = dict(strays=n_stray, strays_reached=n_reach,
                        refogs_complete_if_reached_joined=complete_after,
                        completeness=round(complete_after / len(R), 3),
                        nb_same_S0_family=round(same_fam, 3), dev_nb_other_dev_refog=round(float(other_ref), 4))
    res["base_completeness"] = round(np.mean([len(np.unique(fam[[gid[g] for g in gs if g in gid]])) == 1
                                              for gs in R.values()]), 3)
    res["by_k"] = rows
    json.dump(res, open(os.path.join(a.out, "summary_" + a.model.split("/")[-1] + ".json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
