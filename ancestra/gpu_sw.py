"""v2-gpu: batched Smith-Waterman (Gotoh, BLOSUM62 11/1) on the GPU for explicit pair lists.

Post-freeze v2 work (branch ``v2-gpu``). Not used by any default path: the S2 backend
``gpusw`` (``ancestra.scoring``) calls this module in a *subprocess* running the GPU python
(``$ANCESTRA_GPU_PYTHON``, default ``conda_envs/gpu/bin/python``), so the CPU environment
needs no CUDA packages.

Kernel: one CUDA thread aligns one pair (inter-sequence parallelism, as in CUDASW++'s
inter-task kernel). The shorter sequence runs along the rows, the longer one is processed in
column strips of ``W`` residues kept in registers; the strip boundary (H and E of the last
column) lives in a warp-interleaved global scratch (coalesced). Pairs are sorted by (number
of strips, rows) so the lanes of a warp do the same amount of work. Pure int32 arithmetic and
one thread per pair: results are bit-identical across runs, batch sizes and GPUs.

Outputs per pair: raw SW score, end positions (first maximum in the kernel's scan order),
start positions (second, reversed pass on the pairs that are kept), and from these the bit
score (BLAST lambda/K for BLOSUM62 11/1), an e-value against a fixed database size, and the
query/target coverages. No traceback, so no percent identity (NaN).

CLI (used by ``scoring._backend_gpusw``)::

    python -m ancestra.gpu_sw families IN.npz OUT.npz [--device 0]   # or: pairs IN OUT

IN.npz: ``seq`` (uint8 codes), ``off`` (int64), ``len`` (int32), ``q``/``t`` (int64 pair
indices into the sequences), optional ``db_residues`` and ``evalue``.
OUT.npz: ``keep`` (index into the input pairs of the reported pairs), ``score``, ``bits``,
``evalue``, ``qcov``, ``tcov``, ``qstart``/``qend``/``tstart``/``tend``, ``self_score`` of every
input pair (raw), plus ``stats`` (json).
"""
from __future__ import annotations

import json
import math
import sys
import time

import numpy as np

ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"
_B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""
GAP_OPEN, GAP_EXTEND = 11, 1          # BLAST/DIAMOND convention: a gap of k costs 11 + k
LAMBDA, K = 0.267, 0.041              # BLAST gapped Karlin-Altschul for BLOSUM62 11/1
PAD = 31                              # padding code (score -10000 against everything)
W = 64                                # strip width (columns kept in registers; 64 fastest on A100)
BLOCK = 64                            # threads per CUDA block
SCRATCH_BYTES = 1 << 30               # global scratch per kernel launch
N_STREAMS = 4                         # concurrent launches (hide the tail of long pairs)
WW = 32                               # strip width per lane of the warp kernel
LONG_CELLS = 4_000_000                # pairs with rows x columns above this use the warp kernel


def blosum62_32() -> np.ndarray:
    m = np.array(_B62.split(), dtype=np.int32).reshape(24, 24)
    out = np.full((32, 32), -10000, dtype=np.int32)
    out[:24, :24] = m
    return out


def encode(seqs, lower_as_x: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sequences -> (codes uint8 concatenated, offsets int64, lengths int32). Letters outside
    the BLOSUM alphabet (U, O, J) map to X; lower-case letters (soft-masked low-complexity
    regions) map to X with ``lower_as_x``, else to their upper-case residue."""
    lut = np.full(256, ALPHABET.index("X"), dtype=np.uint8)
    for i, a in enumerate(ALPHABET):
        lut[ord(a)] = i
        if not lower_as_x:
            lut[ord(a.lower())] = i
    lens = np.fromiter((len(s) for s in seqs), dtype=np.int64, count=len(seqs))
    off = np.zeros(len(seqs) + 1, dtype=np.int64)
    np.cumsum(lens, out=off[1:])
    buf = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8)
    return lut[buf], off[:-1].copy(), lens.astype(np.int32)


def bits_of(score) -> np.ndarray:
    return (LAMBDA * np.asarray(score, np.float64) - math.log(K)) / math.log(2)


_SRC = r"""
#define W %(W)d
#define GOE %(GOE)d
#define GE %(GE)d
#define NEG (-1000000)
extern "C" __global__ void sw(const unsigned char* __restrict__ seq, const long long* __restrict__ off,
        const int* __restrict__ rid, const int* __restrict__ cid, const int* __restrict__ rlen,
        const int* __restrict__ clen, const int npairs, const int* __restrict__ mat,
        int* __restrict__ scratch, const long long* __restrict__ woff, const int* __restrict__ wm,
        const int reverse, int* __restrict__ o_score, int* __restrict__ o_i, int* __restrict__ o_j)
{
    __shared__ int smat[1024];
    for (int k = threadIdx.x; k < 1024; k += blockDim.x) smat[k] = mat[k];
    __syncthreads();
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= npairs) return;
    const int lane = p & 31, warp = p >> 5;
    const int m = rlen[p], n = clen[p];
    const unsigned char* R = seq + off[rid[p]];
    const unsigned char* C = seq + off[cid[p]];
    int* bH = scratch + woff[warp];
    int* bE = bH + (long long)(wm[warp] + 1) * 32;
    for (int i = 0; i <= m; i++) { bH[i * 32 + lane] = 0; bE[i * 32 + lane] = NEG; }
    int best = 0, bi = -1, bj = -1;
    for (int j0 = 0; j0 < n; j0 += W) {
        int cr[W], Hr[W], Fr[W];
#pragma unroll
        for (int w = 0; w < W; w++) {
            const int j = j0 + w;
            cr[w] = (j < n) ? (int)(reverse ? C[n - 1 - j] : C[j]) : %(PAD)d;
            Hr[w] = 0; Fr[w] = NEG;
        }
        int hdiag = 0;
        const int lim = n - j0;
        for (int i = 0; i < m; i++) {
            const int a = reverse ? R[m - 1 - i] : R[i];
            const int* srow = smat + a * 32;
            int hl = bH[(i + 1) * 32 + lane];
            int e = bE[(i + 1) * 32 + lane];
            int hd = hdiag;
            hdiag = hl;
            int rbest = 0;
#pragma unroll
            for (int w = 0; w < W; w++) {
                const int f = max(Fr[w] - GE, Hr[w] - GOE);
                e = max(e - GE, hl - GOE);
                const int h = max(max(hd + srow[cr[w]], 0), max(e, f));
                hd = Hr[w]; Hr[w] = h; Fr[w] = f; hl = h;
                rbest = max(rbest, h);
            }
            if (rbest >= best && rbest > 0) {
                /* first valid column holding the row maximum (padding columns excluded);
                   ties between cells are broken by (row, column) so that the reported end
                   does not depend on the strip width or scan order */
                int vb = 0, vw = -1;
#pragma unroll
                for (int w = 0; w < W; w++) if (w < lim && Hr[w] > vb) { vb = Hr[w]; vw = w; }
                if (vb > best || (vb == best && vb > 0 && (i < bi || (i == bi && j0 + vw < bj)))) {
                    best = vb; bi = i; bj = j0 + vw;
                }
            }
            bH[(i + 1) * 32 + lane] = hl;
            bE[(i + 1) * 32 + lane] = e;
        }
    }
    o_score[p] = best; o_i[p] = bi; o_j[p] = bj;
}
"""


_SRC_WARP = r"""
#define WW %(WW)d
#define GOE %(GOE)d
#define GE %(GE)d
#define NEG (-1000000)
/* one warp per pair: lane k owns columns [c0 + k*WW, c0 + (k+1)*WW) of a 32*WW-column pass and
   runs k rows behind lane k-1 (wavefront); the strip boundary (H, E of the last column) moves
   lane -> lane by shuffles, and pass -> pass through bH/bE in global memory (lane 31 writes row i
   at step i+31, lane 0 reads it at step i of the next pass). Same recurrences and the same
   (score, row, column) tie rule as the per-thread kernel -> identical results. */
extern "C" __global__ void sww(const unsigned char* __restrict__ seq, const long long* __restrict__ off,
        const int* __restrict__ rid, const int* __restrict__ cid, const int* __restrict__ rlen,
        const int* __restrict__ clen, const int npairs, const int* __restrict__ mat,
        int* __restrict__ scratch, const long long* __restrict__ soff,
        const int reverse, int* __restrict__ o_score, int* __restrict__ o_i, int* __restrict__ o_j)
{
    __shared__ int smat[1024];
    for (int k = threadIdx.x; k < 1024; k += blockDim.x) smat[k] = mat[k];
    __syncthreads();
    const int lane = threadIdx.x & 31;
    const int p = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    if (p >= npairs) return;
    const int m = rlen[p], n = clen[p];
    const unsigned char* R = seq + off[rid[p]];
    const unsigned char* C = seq + off[cid[p]];
    int* bH = scratch + soff[p];
    int* bE = bH + (m + 1);
    for (int i = lane; i <= m; i += 32) { bH[i] = 0; bE[i] = NEG; }
    __syncwarp();
    int best = 0, bi = -1, bj = -1;
    for (int c0 = 0; c0 < n; c0 += 32 * WW) {
        const int j0 = c0 + lane * WW;
        const int lim = n - j0;
        int cr[WW], Hr[WW], Fr[WW];
#pragma unroll
        for (int w = 0; w < WW; w++) {
            const int j = j0 + w;
            cr[w] = (j < n) ? (int)(reverse ? C[n - 1 - j] : C[j]) : %(PAD)d;
            Hr[w] = 0; Fr[w] = NEG;
        }
        int hdiag = 0;            /* H[i-1][j0-1] */
        int out_h = 0, out_e = NEG;
        for (int step = 0; step < m + 31; step++) {
            const int i = step - lane;
            int in_h = __shfl_up_sync(0xffffffffu, out_h, 1);
            int in_e = __shfl_up_sync(0xffffffffu, out_e, 1);
            if (i >= 0 && i < m) {
                if (lane == 0) { in_h = bH[i + 1]; in_e = bE[i + 1]; }
                if (j0 >= n) { in_h = 0; in_e = NEG; }
                const int a = reverse ? R[m - 1 - i] : R[i];
                const int* srow = smat + a * 32;
                int hl = in_h, e = in_e, hd = hdiag;
                hdiag = in_h;
                int rbest = 0;
#pragma unroll
                for (int w = 0; w < WW; w++) {
                    const int f = max(Fr[w] - GE, Hr[w] - GOE);
                    e = max(e - GE, hl - GOE);
                    const int h = max(max(hd + srow[cr[w]], 0), max(e, f));
                    hd = Hr[w]; Hr[w] = h; Fr[w] = f; hl = h;
                    rbest = max(rbest, h);
                }
                if (rbest >= best && rbest > 0) {
                    int vb = 0, vw = -1;
#pragma unroll
                    for (int w = 0; w < WW; w++) if (w < lim && Hr[w] > vb) { vb = Hr[w]; vw = w; }
                    if (vb > best || (vb == best && vb > 0 && (i < bi || (i == bi && j0 + vw < bj)))) {
                        best = vb; bi = i; bj = j0 + vw;
                    }
                }
                out_h = hl; out_e = e;
                if (lane == 31) { bH[i + 1] = hl; bE[i + 1] = e; }
            }
        }
        __syncwarp();
    }
    /* warp reduction with the tie rule: max score, then smallest row, then smallest column */
    for (int d = 16; d > 0; d >>= 1) {
        const int ob = __shfl_down_sync(0xffffffffu, best, d);
        const int oi = __shfl_down_sync(0xffffffffu, bi, d);
        const int oj = __shfl_down_sync(0xffffffffu, bj, d);
        if (ob > best || (ob == best && ob > 0 && (oi < bi || (oi == bi && oj < bj)))) {
            best = ob; bi = oi; bj = oj;
        }
    }
    if (lane == 0) { o_score[p] = best; o_i[p] = best > 0 ? bi : -1; o_j[p] = best > 0 ? bj : -1; }
}
"""


class GpuSW:
    """Batched SW on one GPU. Short pairs: one thread per pair (``sw``). Long pairs
    (rows x columns > ``LONG_CELLS``): one warp per pair (``sww``), launched first on its own
    stream so that it overlaps the short-pair launches. Both kernels give identical results."""

    def __init__(self, device: int = 0, w: int | None = None, ww: int | None = None):
        import cupy as cp
        self.cp = cp
        cp.cuda.Device(device).use()
        self.w = w or W
        self.ww = ww or WW
        self.kern = cp.RawKernel(_SRC % dict(W=self.w, GOE=GAP_OPEN + GAP_EXTEND, GE=GAP_EXTEND,
                                             PAD=PAD), "sw")
        self.kern_w = cp.RawKernel(_SRC_WARP % dict(WW=self.ww, GOE=GAP_OPEN + GAP_EXTEND,
                                                    GE=GAP_EXTEND, PAD=PAD), "sww")
        self.mat = cp.asarray(blosum62_32().ravel())
        self.main = cp.cuda.Stream(non_blocking=True)
        self.long_s = cp.cuda.Stream(non_blocking=True)
        self.workers = [cp.cuda.Stream(non_blocking=True) for _ in range(N_STREAMS)]

    def set_db(self, seq, off):
        """Upload the sequence buffer once (re-used by every run on the same arrays)."""
        if getattr(self, "_db_id", None) != (id(seq), id(off), len(seq)):
            cp = self.cp
            cp.cuda.Device().synchronize()
            self.d_seq = cp.asarray(seq)
            self.d_off = cp.asarray(off)
            cp.cuda.Device().synchronize()
            self._db_id = (id(seq), id(off), len(seq))

    def _launch_long(self, rid, cid, rl, cl, rev):
        """Warp-per-pair kernel on the long stream; returns device outputs (not synchronised)."""
        cp = self.cp
        n = len(rid)
        with self.long_s:
            soff = np.zeros(n, np.int64)
            soff[1:] = np.cumsum(2 * (rl[:-1].astype(np.int64) + 1))
            scratch = cp.empty(int(soff[-1] + 2 * (int(rl[-1]) + 1)), dtype=cp.int32)
            o = [cp.empty(n, cp.int32) for _ in range(3)]
            blk = 128
            self.kern_w(((n * 32 + blk - 1) // blk,), (blk,),
                        (self.d_seq, self.d_off, cp.asarray(rid.astype(np.int32)),
                         cp.asarray(cid.astype(np.int32)), cp.asarray(rl.astype(np.int32)),
                         cp.asarray(cl.astype(np.int32)), np.int32(n), self.mat, scratch,
                         cp.asarray(soff), rev, o[0], o[1], o[2]))
        return o, scratch

    def _launch_short(self, rid, cid, rl, cl, rev, block):
        """Thread-per-pair kernel over the worker streams; returns device outputs in the input
        order (on the main stream, after waiting for the workers)."""
        cp = self.cp
        npairs = len(rid)
        Wk = self.w
        with self.main:
            d_rl = cp.asarray(rl.astype(np.int32))
            d_cl = cp.asarray(cl.astype(np.int32))
            # balance warps: same number of strips, similar rows (descending work); ties do not
            # matter for the results (every pair is computed independently)
            key = ((d_cl.astype(cp.int64) + Wk - 1) // Wk) * 65536 + d_rl.astype(cp.int64)
            order = cp.argsort(-key)
            rl_s, cl_s = d_rl[order], d_cl[order]
            rid_s = cp.asarray(rid.astype(np.int32))[order]
            cid_s = cp.asarray(cid.astype(np.int32))[order]
            nw = (npairs + 31) // 32
            pad = nw * 32 - npairs
            wmax = cp.concatenate([rl_s, cp.zeros(pad, cp.int32)]).reshape(nw, 32).max(axis=1).get()
            per_w = 2 * (wmax.astype(np.int64) + 1) * 32          # ints of scratch per warp
            cum = np.concatenate([[0], np.cumsum(per_w)])
            budget = SCRATCH_BYTES // 4
            bounds = [0]
            while bounds[-1] < nw:
                s0 = bounds[-1]
                e = int(np.searchsorted(cum, cum[s0] + budget, side="right")) - 1
                bounds.append(max(s0 + 1, min(e, s0 + 65536, nw)))
            need = max(int(cum[b1] - cum[b0]) for b0, b1 in zip(bounds[:-1], bounds[1:]))
            ns = max(1, min(N_STREAMS, len(bounds) - 1))
            scratches = [cp.empty(need, dtype=cp.int32) for _ in range(ns)]
            woffs = [cp.asarray(cum[b0:b1] - cum[b0]) for b0, b1 in zip(bounds[:-1], bounds[1:])]
            wms = [cp.asarray(wmax[b0:b1].astype(np.int32)) for b0, b1 in zip(bounds[:-1], bounds[1:])]
            d_s = cp.empty(npairs, cp.int32)
            d_i = cp.empty(npairs, cp.int32)
            d_j = cp.empty(npairs, cp.int32)
            ready = self.main.record()
        for st in self.workers[:ns]:
            st.wait_event(ready)
        for li, (b0, b1) in enumerate(zip(bounds[:-1], bounds[1:])):
            p0, p1 = b0 * 32, min(b1 * 32, npairs)
            k = p1 - p0
            with self.workers[li % ns]:
                self.kern(((k + block - 1) // block,), (block,),
                          (self.d_seq, self.d_off, rid_s[p0:p1], cid_s[p0:p1], rl_s[p0:p1],
                           cl_s[p0:p1], np.int32(k), self.mat, scratches[li % ns], woffs[li], wms[li],
                           rev, d_s[p0:p1], d_i[p0:p1], d_j[p0:p1]))
        for st in self.workers[:ns]:
            self.main.wait_event(st.record())
        with self.main:
            inv = cp.empty_like(order)
            inv[order] = cp.arange(npairs, dtype=order.dtype)
            out = (d_s[inv], d_i[inv], d_j[inv])
        return out, (scratches, woffs, wms, rid_s, cid_s, rl_s, cl_s)

    def run(self, seq, off, rid, cid, reverse=False, rlen=None, clen=None, block=None):
        """Score pairs (rid[k] on the rows, cid[k] on the columns; rows should be the shorter
        sequence). Returns (score, end_row, end_col) in the input order (numpy int32)."""
        npairs = len(rid)
        out = [np.zeros(npairs, np.int32), np.full(npairs, -1, np.int32), np.full(npairs, -1, np.int32)]
        if npairs == 0:
            return tuple(out)
        self.set_db(seq, off)
        rid = np.asarray(rid, np.int64)
        cid = np.asarray(cid, np.int64)
        rl = np.asarray(rlen, np.int64)
        cl = np.asarray(clen, np.int64)
        rev = np.int32(1 if reverse else 0)
        is_long = rl * cl > LONG_CELLS
        iL, iS = np.flatnonzero(is_long), np.flatnonzero(~is_long)
        dl = self._launch_long(rid[iL], cid[iL], rl[iL], cl[iL], rev) if len(iL) else None
        ds = self._launch_short(rid[iS], cid[iS], rl[iS], cl[iS], rev, block or BLOCK) if len(iS) else None
        if ds is not None:
            with self.main:
                for x, d in zip(out, ds[0]):
                    x[iS] = d.get()
        if dl is not None:
            with self.long_s:
                for x, d in zip(out, dl[0]):
                    x[iL] = d.get()
        self.cp.cuda.Device().synchronize()
        return tuple(out)


def align_pairs(seq, off, lens, q, t, *, db_residues: float = 5e7, evalue: float = 10.0,
                device: int = 0, sw: GpuSW | None = None) -> dict:
    """Full computation for pairs (q[k], t[k]): forward pass (score, ends), threshold on the
    e-value, reversed pass for the starts of the kept pairs. Returns a dict of numpy arrays."""
    sw = sw or GpuSW(device)
    q = np.asarray(q, np.int64)
    t = np.asarray(t, np.int64)
    lq, lt = lens[q].astype(np.int64), lens[t].astype(np.int64)
    qrow = lq <= lt                              # shorter sequence on the rows
    rid = np.where(qrow, q, t)
    cid = np.where(qrow, t, q)
    rl, cl = np.minimum(lq, lt), np.maximum(lq, lt)
    t0 = time.time()
    s, ei, ej = sw.run(seq, off, rid, cid, rlen=rl, clen=cl)
    t_fwd = time.time() - t0
    bits = bits_of(s)
    ev = lq.astype(np.float64) * db_residues * np.power(2.0, -bits)
    keep = np.flatnonzero((s > 0) & (ev <= evalue))
    # reversed pass on prefixes ending at the forward maximum -> start positions
    t0 = time.time()
    s2, ri, rj = sw.run(seq, off, rid[keep], cid[keep], reverse=True,
                        rlen=ei[keep].astype(np.int64) + 1, clen=ej[keep].astype(np.int64) + 1)
    t_rev = time.time() - t0
    rs = ei[keep] - ri               # start on rows = (m_prefix - 1) - rev_i
    cs = ej[keep] - rj
    qrk = qrow[keep]
    qstart = np.where(qrk, rs, cs)
    qend = np.where(qrk, ei[keep], ej[keep])
    tstart = np.where(qrk, cs, rs)
    tend = np.where(qrk, ej[keep], ei[keep])
    cells_fwd = float((rl * cl).sum())
    cells_rev = float(((ei[keep].astype(np.float64) + 1) * (ej[keep] + 1)).sum())
    return dict(keep=keep, score=s[keep], bits=bits[keep].astype(np.float32),
                evalue=ev[keep].astype(np.float32),
                qcov=((qend - qstart + 1) / lq[keep]).astype(np.float32),
                tcov=((tend - tstart + 1) / lt[keep]).astype(np.float32),
                qstart=qstart.astype(np.int32), qend=qend.astype(np.int32),
                tstart=tstart.astype(np.int32), tend=tend.astype(np.int32),
                rev_mismatch=int((s2 != s[keep]).sum()),
                stats=dict(n_pairs=int(len(q)), n_kept=int(len(keep)), fwd_s=round(t_fwd, 3),
                           rev_s=round(t_rev, 3), cells_fwd=cells_fwd, cells_rev=cells_rev,
                           gcups_fwd=round(cells_fwd / max(t_fwd, 1e-9) / 1e9, 1),
                           gcups_rev=round(cells_rev / max(t_rev, 1e-9) / 1e9, 1)))


def family_pairs_chunks(fam: np.ndarray, chunk_pairs: int = 30_000_000):
    """All unordered within-family pairs q <= t (self pairs included) of families with >= 2
    members, in family-id order, yielded in chunks of about ``chunk_pairs`` pairs."""
    order = np.argsort(fam, kind="stable")
    sf = fam[order]
    cuts = np.flatnonzero(np.diff(sf)) + 1
    qs, ts, n = [], [], 0
    for m in np.split(order, cuts):
        if len(m) < 2:
            continue
        m = np.sort(m)
        i, j = np.triu_indices(len(m))
        qs.append(m[i])
        ts.append(m[j])
        n += len(i)
        if n >= chunk_pairs:
            yield np.concatenate(qs), np.concatenate(ts)
            qs, ts, n = [], [], 0
    if qs:
        yield np.concatenate(qs), np.concatenate(ts)


def score_families_gpu(seq, off, lens, fam, *, db_residues: float, evalue: float = 10.0,
                       device: int = 0, chunk_pairs: int = 30_000_000) -> dict:
    """Exhaustive within-family SW for every family with >= 2 members (self pairs give the
    self-scores). Returns the kept records (e-value <= evalue) and stats."""
    import queue
    import threading
    sw = GpuSW(device)
    out = {k: [] for k in ("q", "t", "bits", "qcov", "tcov", "evalue")}
    st = dict(n_pairs=0, n_kept=0, fwd_s=0.0, rev_s=0.0, cells_fwd=0.0, cells_rev=0.0, n_chunks=0)
    # pair lists are generated by a producer thread while the GPU works on the previous chunk
    qq: queue.Queue = queue.Queue(maxsize=1)

    def produce():
        for c in family_pairs_chunks(fam, chunk_pairs):
            qq.put(c)
        qq.put(None)
    th = threading.Thread(target=produce, daemon=True)
    th.start()
    while True:
        c = qq.get()
        if c is None:
            break
        q, t = c
        r = align_pairs(seq, off, lens, q, t, db_residues=db_residues, evalue=evalue, sw=sw)
        k = r["keep"]
        out["q"].append(q[k].astype(np.int32))
        out["t"].append(t[k].astype(np.int32))
        for c in ("bits", "qcov", "tcov", "evalue"):
            out[c].append(r[c])
        for c in ("n_pairs", "n_kept", "fwd_s", "rev_s", "cells_fwd", "cells_rev"):
            st[c] += r["stats"][c]
        st["n_chunks"] += 1
        if r["rev_mismatch"]:
            raise RuntimeError(f"reverse pass disagrees with forward scores ({r['rev_mismatch']} pairs)")
    res = {k: (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in out.items()}
    st["gcups_fwd"] = round(st["cells_fwd"] / max(st["fwd_s"], 1e-9) / 1e9, 1)
    st["gcups_rev"] = round(st["cells_rev"] / max(st["rev_s"], 1e-9) / 1e9, 1)
    st["fwd_s"], st["rev_s"] = round(st["fwd_s"], 2), round(st["rev_s"], 2)
    res["stats"] = st
    return res


def sw_reference(a: str, b: str) -> int:
    """Pure-python Gotoh local alignment score (tests)."""
    m = blosum62_32()
    ea, _, _ = encode([a])
    eb, _, _ = encode([b])
    go, ge = GAP_OPEN + GAP_EXTEND, GAP_EXTEND
    n = len(b)
    H = [0] * (n + 1)
    F = [-10**9] * (n + 1)
    best = 0
    for i in range(len(a)):
        diag, e, Hn = 0, -10**9, [0] * (n + 1)
        for j in range(n):
            F[j + 1] = max(F[j + 1] - ge, H[j + 1] - go)
            e = max(e - ge, Hn[j] - go)
            h = max(0, H[j] + int(m[ea[i], eb[j]]), e, F[j + 1])
            Hn[j + 1] = h
            best = max(best, h)
        H = Hn
    return best


def main(argv=None):
    """``python -m ancestra.gpu_sw [pairs|families] IN.npz OUT.npz [--device N]``"""
    import argparse
    import os
    ap = argparse.ArgumentParser(prog="python -m ancestra.gpu_sw")
    ap.add_argument("mode", choices=["pairs", "families"])
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--device", type=int, default=0)
    a = ap.parse_args(argv)
    z = np.load(a.inp)
    t0 = time.time()
    dbr = float(z["db_residues"]) if "db_residues" in z else 5e7
    ev = float(z["evalue"]) if "evalue" in z else 10.0
    if a.mode == "pairs":
        res = align_pairs(z["seq"], z["off"], z["len"], z["q"], z["t"], db_residues=dbr, evalue=ev,
                          device=a.device)
    else:
        res = score_families_gpu(z["seq"], z["off"], z["len"], z["fam"], db_residues=dbr, evalue=ev,
                                 device=a.device)
    res["stats"]["total_s"] = round(time.time() - t0, 3)
    stats = res.pop("stats")
    tmp = a.out + ".tmp.npz"
    np.savez(tmp, stats=np.array(json.dumps(stats)), **res)
    os.replace(tmp, a.out)
    print(json.dumps(stats), file=sys.stderr)


if __name__ == "__main__":
    main()
