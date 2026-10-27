#!/usr/bin/env python
"""L2-plm: embed OrthoBench proteins with a pLM (run in conda_envs/gpu).

    python scripts/l2_plm_embed.py --model e650 --genes GIDS.npy --out work/l2_plm/emb_e650_refog [--residue]

Writes OUT.mean.npy (n x d float32, mean over residues, L2-normalised NOT applied), OUT.gids.npy and, with
--residue, OUT.res.npy (concatenated per-residue fp16) + OUT.off.npy (offsets). Long proteins are windowed
(ESM-2: 1022 aa; ProtT5: 4000 aa); the mean is length-weighted over windows, per-residue is concatenated.
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from l2_plm_common import HF, DATA, gene_table, read_seqs, manifest

MODELS = {"e650": ("facebook/esm2_t33_650M_UR50D", 1022), "e3b": ("facebook/esm2_t36_3B_UR50D", 1022),
          "e35": ("facebook/esm2_t12_35M_UR50D", 1022), "pt5": ("Rostlab/prot_t5_xl_half_uniref50-enc", 4000)}


def load(key):
    os.environ["HF_HOME"] = HF
    import torch
    name, win = MODELS[key]
    if key == "pt5":
        from transformers import T5EncoderModel, T5Tokenizer
        tok = T5Tokenizer.from_pretrained(name, do_lower_case=False, legacy=True)
        model = T5EncoderModel.from_pretrained(name, torch_dtype=torch.float16)
    else:
        from transformers import AutoTokenizer, EsmModel
        tok = AutoTokenizer.from_pretrained(name)
        model = EsmModel.from_pretrained(name, add_pooling_layer=False, torch_dtype=torch.float16)
    return tok, model.cuda().eval(), win


def prep(key, s):
    s = s.upper()
    if key == "pt5":
        s = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in s)
        return " ".join(s)
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=list(MODELS))
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--residue", action="store_true")
    ap.add_argument("--max-tokens", type=int, default=12000)
    a = ap.parse_args()
    import torch
    torch.manual_seed(11)
    t0 = time.time()
    gt = gene_table()
    gids = np.load(a.genes).astype(np.int64)
    names = [gt.names[g] for g in gids]
    seqs = read_seqs(names)
    tok, model, win = load(a.model)
    pieces, owner = [], []
    for i, s in enumerate(seqs):
        for st in range(0, len(s), win):
            pieces.append(s[st:st + win]); owner.append(i)
    order = np.argsort([-len(p) for p in pieces], kind="stable")
    d = model.config.d_model if a.model == "pt5" else model.config.hidden_size
    pm = np.zeros((len(pieces), d), np.float32)
    pres = [None] * len(pieces) if a.residue else None
    k = 0
    torch.cuda.reset_peak_memory_stats()
    with torch.inference_mode():
        while k < len(order):
            L = len(pieces[order[k]]) + 2
            bs = max(1, a.max_tokens // L)
            idx = order[k:k + bs]
            enc = tok([prep(a.model, pieces[j]) for j in idx], return_tensors="pt", padding=True)
            enc = {x: y.cuda() for x, y in enc.items() if x in ("input_ids", "attention_mask")}
            h = model(**enc).last_hidden_state.float()
            for r, j in enumerate(idx):
                n = len(pieces[j])
                e = h[r, 1:n + 1] if a.model != "pt5" else h[r, :n]     # ESM: drop <cls>/<eos>; T5: drop </s>
                pm[j] = e.mean(0).cpu().numpy()
                if a.residue:
                    pres[j] = e.half().cpu().numpy()
            k += len(idx)
    torch.cuda.synchronize()
    owner = np.array(owner)
    wl = np.array([len(p) for p in pieces], np.float32)
    E = np.zeros((len(seqs), d), np.float32)
    np.add.at(E, owner, pm * wl[:, None])
    E /= np.bincount(owner, weights=wl, minlength=len(seqs))[:, None]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.save(a.out + ".mean.npy", E)
    np.save(a.out + ".gids.npy", gids)
    if a.residue:
        per = [[] for _ in seqs]
        for j in range(len(pieces)):
            per[owner[j]].append(pres[j])
        per = [np.concatenate(x) for x in per]
        off = np.zeros(len(per) + 1, np.int64); off[1:] = np.cumsum([len(x) for x in per])
        assert all(len(x) == len(s) for x, s in zip(per, seqs))
        np.save(a.out + ".res.npy", np.concatenate(per)); np.save(a.out + ".off.npy", off)
    res = sum(map(len, seqs))
    info = manifest(dict(model=MODELS[a.model][0], window=win, dtype="fp16", residue=a.residue,
                         max_tokens=a.max_tokens),
                    dict(genes=a.genes, fasta="work/S0/ob_all.fa"), t0,
                    dict(n=len(seqs), residues=int(res), residues_per_s=round(res / (time.time() - t0)),
                         peak_gpu_gb=round(torch.cuda.max_memory_allocated() / 2**30, 2)))
    json.dump(info, open(a.out + ".json", "w"), indent=1)
    print(json.dumps(info))


if __name__ == "__main__":
    main()
