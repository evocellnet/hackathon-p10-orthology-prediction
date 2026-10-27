#!/usr/bin/env python
"""v2-gpu task (d): batched ProstT5 encoder + CNN 3Di head (HF transformers, fp16) on the GPU.

    python scripts/gpu_prostt5_batched.py IN.fa OUT_3di.fa [--max-tokens 8000]

Weights: Rostlab/ProstT5_fp16 (HF hub, cached under data/models) and the encoder-only CNN head
cnn_chkpnt/model.pt from github.com/mheinzinger/ProstT5 (data/models/prostt5_cnn/model.pt). This
follows the authors' predict_3Di_encoderOnly.py (prefix <AA2fold>, U/Z/O/B -> X, argmax of the CNN).
Prints residues/s. Runs in conda_envs/gpu.
"""
import argparse
import os
import re
import sys
import time

import numpy as np
import torch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SS = "ACDEFGHIKLMNPQRSTVWY"


class CNN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.classifier = torch.nn.Sequential(
            torch.nn.Conv2d(1024, 32, kernel_size=(7, 1), padding=(3, 0)), torch.nn.ReLU(),
            torch.nn.Dropout(0.0), torch.nn.Conv2d(32, 20, kernel_size=(7, 1), padding=(3, 0)))

    def forward(self, x):
        return self.classifier(x.permute(0, 2, 1).unsqueeze(-1)).squeeze(-1)


def read_fasta(p):
    out, h = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            h = line[1:].split()[0]
            out[h] = []
        elif h:
            out[h].append(line)
    return {k: "".join(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--max-tokens", type=int, default=8000)
    ap.add_argument("--max-batch", type=int, default=256)
    a = ap.parse_args()
    from transformers import T5EncoderModel, T5Tokenizer
    tok = T5Tokenizer.from_pretrained("Rostlab/ProstT5_fp16", do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/ProstT5_fp16", torch_dtype=torch.float16).cuda().eval()
    cnn = CNN()
    ck = torch.load(os.path.join(ROOT, "data/models/prostt5_cnn/model.pt"), map_location="cpu")
    cnn.load_state_dict(ck.get("state_dict", ck))
    cnn = cnn.cuda().eval()
    seqs = read_fasta(a.inp)
    names = sorted(seqs, key=lambda k: -len(seqs[k]))
    res = {}
    torch.cuda.synchronize()
    t0 = time.time()
    k = 0
    with torch.inference_mode():
        while k < len(names):
            L = len(seqs[names[k]]) + 2
            bs = max(1, min(a.max_batch, a.max_tokens // L))
            batch = names[k:k + bs]
            txt = ["<AA2fold> " + " ".join(re.sub(r"[UZOB]", "X", seqs[n].upper())) for n in batch]
            enc = tok(txt, add_special_tokens=True, padding="longest", return_tensors="pt")
            ids, mask = enc.input_ids.cuda(), enc.attention_mask.cuda()
            h = model(ids, attention_mask=mask).last_hidden_state
            pred = cnn(h.float()).argmax(dim=1).cpu().numpy()          # B x Ltok
            for i, n in enumerate(batch):
                res[n] = "".join(SS[x] for x in pred[i, 1:len(seqs[n]) + 1])
            k += bs
    torch.cuda.synchronize()
    dt = time.time() - t0
    nres = sum(len(seqs[n]) for n in names)
    with open(a.out, "w") as fh:
        for n in seqs:
            fh.write(f">{n}\n{res[n]}\n")
    print(f"proteins {len(names)} residues {nres} seconds {dt:.1f} residues/s {nres / dt:.0f}", file=sys.stderr)


if __name__ == "__main__":
    main()
