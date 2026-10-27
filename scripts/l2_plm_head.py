#!/usr/bin/env python
"""L2-plm idea 2: contrastive projection head on frozen mean-pooled pLM embeddings, trained on the self-derived
labels of scripts/l2_plm_selflabels.py (synteny-supported anchors vs same-family paralogs; no benchmark labels).

Head: d -> 512 -> 256 (GELU, dropout 0.1), unit-normalised output; score = cosine. Loss: InfoNCE over
(positive, up to K=16 sampled negatives) per query, temperature 0.1. AdamW lr 1e-3 wd 1e-4, batch 256 queries,
<= 30 epochs, the epoch with the best val triplet accuracy is kept. Inputs standardised with train-gene stats.
  --neg both|i      negative types: (i)+(ii) (pre-registered) or (i) only (logged deviation)
  --labels real|shuffled   shuffled = control C2: positive replaced by a random same-family gene of another
                    species (seed 11)
Writes data/extra/l2_plm/heads/<name>.pt, work/l2_plm/head_<name>_refog.npz (gids, z of the DEV RefOG genes),
results/l2_plm/head_<name>.json (curves, manifest).  Runs in conda_envs/gpu.
"""
import argparse, json, os, sys, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2_plm_common import DATA, WORK, RES, manifest, md5

K = 16


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--neg", default="both", choices=["both", "i"])
    ap.add_argument("--labels", default="real", choices=["real", "shuffled"])
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--label-set", default="prereg", choices=["prereg", "hard"],
                    help="hard = deviation D2 labels (work/l2_plm/selflabels_hard.npz)")
    a = ap.parse_args()
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    t0 = time.time()
    name = f"{a.model}_{a.neg}_{a.labels}_s{a.seed}" + ("_D2hard" if a.label_set == "hard" else "")
    torch.manual_seed(a.seed); np.random.seed(a.seed)
    rng = np.random.default_rng(a.seed)
    z = np.load(os.path.join(WORK, "selflabels.npz" if a.label_set == "prereg" else "selflabels_hard.npz"))
    q, p, n_, ty, fam, val = z["q"], z["p"], z["n"], z["ty"], z["fam"], z["val"]
    if a.neg == "i":
        m = ty == 1
        q, p, n_, ty, fam, val = q[m], p[m], n_[m], ty[m], fam[m], val[m]
    import glob
    Es, Gs = [], []
    pat = f"emb_{a.model}_train.mean.npy" if a.label_set == "prereg" else f"emb_{a.model}_train*.mean.npy"
    for f in sorted(glob.glob(os.path.join(WORK, pat))):
        Es.append(np.load(f)); Gs.append(np.load(f.replace(".mean.npy", ".gids.npy")))
    E = np.concatenate(Es); gids = np.concatenate(Gs)
    gids, first = np.unique(gids, return_index=True); E = E[first]
    row = {int(g): i for i, g in enumerate(gids)}
    # queries: (q, p) -> list of negatives
    groups = defaultdict(list)
    for qq, pp, nn_, vv in zip(q.tolist(), p.tolist(), n_.tolist(), val.tolist()):
        groups[(qq, pp, vv)].append(nn_)
    keys = list(groups)
    if a.labels == "shuffled" and a.label_set == "hard":
        # D2 control: the positive is a random member of the near-tie set instead of the synteny winner
        r11 = np.random.default_rng(11); ng = {}
        for (qq, pp, vv), negs in groups.items():
            c = [pp] + negs; w = int(r11.choice(c))
            ng[(qq, w, vv)] = [x for x in c if x != w]
        groups = ng; keys = list(groups)
    elif a.labels == "shuffled":
        from l2_plm_common import gene_table
        gt = gene_table(); sp = np.asarray(gt.species_of)
        famall = np.load(os.path.join(DATA, "work/core_v0_ob/run_v2/stages/scores/families_final.npy"))
        mem = defaultdict(list)
        for g in gids.tolist():
            mem[famall[g]].append(g)
        r11 = np.random.default_rng(11)
        newk = {}
        for (qq, pp, vv) in keys:
            cand = [g for g in mem[famall[qq]] if sp[g] != sp[qq]]
            newk[(qq, pp, vv)] = (qq, int(r11.choice(cand)), vv)
        groups = {newk[k_]: [x for x in v if x != newk[k_][1]] for k_, v in groups.items()}
        groups = {k_: v for k_, v in groups.items() if v}
        keys = list(groups)
    tr = [k_ for k_ in keys if not k_[2]]; va = [k_ for k_ in keys if k_[2]]
    mu = E[[row[g] for g in {k_[0] for k_ in tr} | {k_[1] for k_ in tr}]].mean(0)
    sd = E.std(0) + 1e-6
    X = torch.tensor((E - mu) / sd, dtype=torch.float32, device="cuda")
    d = X.shape[1]
    head = nn.Sequential(nn.Linear(d, 512), nn.GELU(), nn.Dropout(0.1), nn.Linear(512, 256)).cuda()
    opt = torch.optim.AdamW(head.parameters(), lr=1e-3, weight_decay=1e-4)
    # val triplets (all negatives, no sampling)
    vq = np.array([row[k_[0]] for k_ in va for _ in groups[k_]]); vp = np.array([row[k_[1]] for k_ in va for _ in groups[k_]])
    vn = np.array([row[x] for k_ in va for x in groups[k_]])
    vfam = np.array([k_[0] for k_ in va for _ in groups[k_]])

    def val_acc(raw=False):
        head.eval()
        with torch.no_grad():
            Z = F.normalize(X if raw else head(X), dim=1)
            s1 = (Z[vq] * Z[vp]).sum(1); s2 = (Z[vq] * Z[vn]).sum(1)
            c = torch.where(s1 > s2, 1.0, torch.where(s1 == s2, 0.5, 0.0))
        head.train()
        return float(c.mean())

    def train_acc_sample(Zf, idxs):
        return None

    curves = dict(epoch=[], train_loss=[], val_acc=[], train_acc=[])
    base_val = val_acc(raw=True)
    best = (-1, None, -1)
    trq = np.array([row[k_[0]] for k_ in tr]); trp = np.array([row[k_[1]] for k_ in tr])
    trn = [np.array([row[x] for x in groups[k_]]) for k_ in tr]
    for ep in range(a.epochs):
        perm = rng.permutation(len(tr))
        tot, nb = 0.0, 0
        for s in range(0, len(perm), 256):
            ii = perm[s:s + 256]
            neg = np.full((len(ii), K), -1, np.int64)
            for r_, j in enumerate(ii):
                c = trn[j]
                pick = c if len(c) <= K else rng.choice(c, K, replace=False)
                neg[r_, :len(pick)] = pick
            nm = torch.tensor(neg >= 0, device="cuda")
            Zq = F.normalize(head(X[trq[ii]]), dim=1); Zp = F.normalize(head(X[trp[ii]]), dim=1)
            Zn = F.normalize(head(X[torch.tensor(np.maximum(neg, 0), device="cuda")]), dim=2)
            lp = (Zq * Zp).sum(1, keepdim=True)
            ln = torch.einsum("bd,bkd->bk", Zq, Zn).masked_fill(~nm, -1e4)
            logits = torch.cat([lp, ln], 1) / 0.1
            loss = F.cross_entropy(logits, torch.zeros(len(ii), dtype=torch.long, device="cuda"))
            opt.zero_grad(); loss.backward(); opt.step()
            tot += float(loss) * len(ii); nb += len(ii)
        # train triplet accuracy on a fixed sample of 20k train triplets
        head.eval()
        with torch.no_grad():
            Z = F.normalize(head(X), dim=1)
            rs = np.random.default_rng(0).choice(len(tr), min(5000, len(tr)), replace=False)
            tq = np.concatenate([np.full(len(trn[j]), trq[j]) for j in rs]); tp = np.concatenate([np.full(len(trn[j]), trp[j]) for j in rs])
            tn = np.concatenate([trn[j] for j in rs])
            s1 = (Z[tq] * Z[tp]).sum(1); s2 = (Z[tq] * Z[tn]).sum(1)
            tacc = float(torch.where(s1 > s2, 1.0, torch.where(s1 == s2, 0.5, 0.0)).mean())
        head.train()
        va_ = val_acc()
        curves["epoch"].append(ep + 1); curves["train_loss"].append(round(tot / nb, 4))
        curves["val_acc"].append(round(va_, 4)); curves["train_acc"].append(round(tacc, 4))
        print(f"ep {ep + 1} loss {tot / nb:.4f} train_acc {tacc:.4f} val_acc {va_:.4f}", flush=True)
        if va_ > best[0]:
            best = (va_, {k_: v.detach().clone() for k_, v in head.state_dict().items()}, ep + 1)
    head.load_state_dict(best[1]); head.eval()
    os.makedirs(os.path.join(DATA, "data/extra/l2_plm/heads"), exist_ok=True)
    wf = os.path.join(DATA, "data/extra/l2_plm/heads", name + ".pt")
    torch.save(dict(state=best[1], mu=mu, sd=sd, model=a.model, d=d), wf)
    # project DEV RefOG genes (mean embeddings from the refog run of the same backbone)
    Er = np.load(os.path.join(WORK, f"emb_{a.model}_refog.mean.npy")); gr = np.load(os.path.join(WORK, f"emb_{a.model}_refog.gids.npy"))
    with torch.no_grad():
        Zr = F.normalize(head(torch.tensor((Er - mu) / sd, dtype=torch.float32, device="cuda")), dim=1).cpu().numpy()
    np.savez(os.path.join(WORK, f"head_{name}_refog.npz"), gids=gr, z=Zr)
    out = dict(name=name, n_train_queries=len(tr), n_val_queries=len(va), n_val_triplets=int(len(vq)),
               val_acc_raw_cosine=round(base_val, 4), best_val_acc=round(best[0], 4), best_epoch=best[2], curves=curves,
               weights=wf, weights_md5=md5(wf))
    out["manifest"] = manifest(vars(a) | dict(K=K, temp=0.1, lr=1e-3, wd=1e-4, batch=256, hidden=512, out=256),
                               dict(labels="work/l2_plm/selflabels.npz" if a.label_set == "prereg" else "work/l2_plm/selflabels_hard.npz",
                                     emb=f"work/l2_plm/emb_{a.model}_train*"), t0)
    json.dump(out, open(os.path.join(RES, f"head_{name}.json"), "w"), indent=1)
    print(json.dumps({k_: v for k_, v in out.items() if k_ not in ("curves", "manifest")}))


if __name__ == "__main__":
    main()
