#!/usr/bin/env python
"""L2-plm: embedding-alignment score (EBA / pLM-BLAST spirit, own simplified variant) on the GPU.

For a pair (x, y) of per-residue embeddings: S = cosine similarity matrix (Lx x Ly); EBA-style enhancement
S' = 0.5 * (row-z(S) + col-z(S)); local Smith-Waterman on S' with a linear gap GAP and zero floor; raw score
= best cell. Normalised score = raw / min(self(x), self(y)). Pre-registered constants: GAP = 1.0.

    python scripts/l2_plm_alnscore.py --emb work/l2_plm/emb_e650_refog --pairs P.npy --out O.npy
P.npy: (m, 2) int64 of *row indices into the embedding file* (self pairs allowed). Output: raw scores (m,).
Runs in conda_envs/gpu (torch + cupy).
"""
import argparse, json, os, sys, time
import numpy as np

GAP = 1.0
KERNEL = r"""
extern "C" __global__ void swd(const float* S, const long long* off, const int* L1, const int* L2,
                               float* scratch, const long long* soff, float* out, int n, float gap) {
  // one block per pair; anti-diagonal wavefront; buffers indexed by row i (0..a)
  int p = blockIdx.x;
  if (p >= n) return;
  const float* M = S + off[p];
  int a = L1[p], b = L2[p];
  float* B0 = scratch + soff[p];          // diagonal d-2
  float* B1 = B0 + (a + 1);               // diagonal d-1
  float* B2 = B1 + (a + 1);               // diagonal d
  float best = 0.f;
  for (int d = 2; d <= a + b; d++) {
    int lo = max(1, d - b), hi = min(a, d - 1);
    int lo1 = max(1, d - 1 - b), hi1 = min(a, d - 2);
    int lo2 = max(1, d - 2 - b), hi2 = min(a, d - 3);
    for (int i = lo + threadIdx.x; i <= hi; i += blockDim.x) {
      int j = d - i;
      float dg = (i - 1 >= lo2 && i - 1 <= hi2) ? B0[i - 1] : 0.f;
      float up = (i - 1 >= lo1 && i - 1 <= hi1) ? B1[i - 1] : 0.f;
      float lf = (i >= lo1 && i <= hi1) ? B1[i] : 0.f;
      float h = fmaxf(0.f, fmaxf(dg + M[(long long)(i - 1) * b + (j - 1)], fmaxf(up - gap, lf - gap)));
      B2[i] = h;
      best = fmaxf(best, h);
    }
    __syncthreads();
    float* t = B0; B0 = B1; B1 = B2; B2 = t;
  }
  __shared__ float red[256];
  red[threadIdx.x] = best;
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) red[threadIdx.x] = fmaxf(red[threadIdx.x], red[threadIdx.x + s]);
    __syncthreads();
  }
  if (threadIdx.x == 0) out[p] = red[0];
}
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-cells", type=float, default=1.5e8)
    a = ap.parse_args()
    import torch, cupy as cp
    t0 = time.time()
    R = np.load(a.emb + ".res.npy", mmap_mode="r"); off = np.load(a.emb + ".off.npy")
    P = np.load(a.pairs).astype(np.int64)
    lens = np.diff(off)
    # per-gene GPU tensors (unit-normalised residues), loaded once
    Rt = torch.empty(R.shape, dtype=torch.float16, device="cuda")
    for c in range(0, len(R), 50000):
        x = torch.from_numpy(np.array(R[c:c + 50000])).cuda().float()
        Rt[c:c + 50000] = torch.nn.functional.normalize(x, dim=1).half()
        del x
    k = cp.RawKernel(KERNEL, "swd")
    # shorter sequence on the rows (fewer serial outer iterations do not matter; inner loop is per row)
    cells = lens[P[:, 0]] * lens[P[:, 1]]
    order = np.argsort(-cells, kind="stable")
    out = np.zeros(len(P), np.float32)
    i, done, last = 0, 0, time.time()
    while i < len(order):
        tot, j = 0, i
        while j < len(order) and (tot == 0 or tot + cells[order[j]] <= a.max_cells):
            tot += cells[order[j]]; j += 1
        idx = order[i:j]
        mats, L1, L2 = [], [], []
        for q in idx:
            x, y = P[q]
            A = Rt[off[x]:off[x + 1]]; B = Rt[off[y]:off[y + 1]]
            S = (A @ B.T).float()
            z1 = (S - S.mean(1, keepdim=True)) / (S.std(1, keepdim=True) + 1e-6)
            z2 = (S - S.mean(0, keepdim=True)) / (S.std(0, keepdim=True) + 1e-6)
            mats.append((0.5 * (z1 + z2)).reshape(-1)); L1.append(S.shape[0]); L2.append(S.shape[1])
        flat = torch.cat(mats).contiguous()
        sizes = np.array([m.numel() for m in mats], np.int64)
        o = np.zeros(len(idx), np.int64); o[1:] = np.cumsum(sizes)[:-1]
        L2a = np.array(L2, np.int32)
        so = np.zeros(len(idx), np.int64); so[1:] = np.cumsum(3 * (np.array(L1, np.int64) + 1))[:-1]
        scratch = cp.zeros(int(so[-1] + 3 * (L1[-1] + 1)), cp.float32)
        res = cp.zeros(len(idx), cp.float32)
        Sc = cp.from_dlpack(flat)
        k((len(idx),), (256,), (Sc, cp.asarray(o), cp.asarray(np.array(L1, np.int32)),
                                               cp.asarray(L2a), scratch, cp.asarray(so), res,
                                               np.int32(len(idx)), np.float32(GAP)))
        out[idx] = cp.asnumpy(res)
        del flat, Sc, scratch, mats
        i = j
        done += tot
        if time.time() - last > 60:
            last = time.time()
            print(f"progress {i}/{len(order)} pairs, {done / cells.sum():.3f} of cells, {time.time() - t0:.0f}s, "
                  f"{done / (time.time() - t0) / 1e9:.2f} Gcells/s", file=sys.stderr, flush=True)
    np.save(a.out, out)
    print(json.dumps(dict(pairs=len(P), cells=float(cells.sum()), wall_s=round(time.time() - t0, 1),
                          peak_gpu_gb=round(torch.cuda.max_memory_allocated() / 2**30, 2))))


if __name__ == "__main__":
    main()
