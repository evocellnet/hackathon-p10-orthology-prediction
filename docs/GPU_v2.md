# GPU speed-ups for Ancestra (post-freeze v2 work, branch `v2-gpu`)

Status: in progress (interim notes; started 2026-10-27 11:11 CET, T+19h34 by `scripts/clock.sh`).
Node: 48 CPUs, 500 GB, 1× A100-PCIE-40GB (driver 580, CUDA 13.0). This work used the GPU plus ≤16 CPU threads;
other agents ran on the rest of the node (load 25–45), so CPU wall times are upper bounds.

Rules followed: `ancestra/` at tag `freeze-v1` is not changed in behaviour. All GPU code is opt-in and lives on branch
`v2-gpu` (git worktree `work/v2gpu_wt`). Accuracy is checked only on OrthoBench DEV RefOGs
(`splits/orthobench_split.json`) and TreeFam-A dev, or by equivalence to the CPU path.

## 1. MMseqs2-GPU and DIAMOND sanity check (task 1)

Setup: `work/gputest/`, human UP000005640 (20,652 proteins, 11.5 M residues) against itself; MMseqs2 18-8cc5c from
`conda_envs/orth` (built with CUDA; the binary carries SASS for sm_75/80/86/89/90, so there is **no PTX JIT** on the A100).

| run | wall | notes |
|---|---|---|
| CPU `search -s 7.5`, 48 threads (earlier) | 48 s (88 s at 24 thr) | prefilter 35 s + align 13 s |
| GPU `search --gpu 1` (first run, earlier) | 252 s | GPU ungapped prefilter 197 s + CPU align 54 s (8 threads) |
| GPU `ungappedprefilter --gpu 1`, warm rerun | **198 s** | identical to the first run → not a warm-up/JIT effect |
| same, 1,033 queries vs human | 10.3 s | **≈ 10 ms per query** |
| same, 1,033 queries vs all OrthoBench (251k proteins, 128 M residues) | 23.6 s | 11× the residues for 2.3× the time |
| 4 concurrent processes × 258 queries vs human | 12.1 s | = sequential; GPU at 96–99 % SM busy throughout |
| `--prefilter-mode 3` (ungapped + gapped on GPU), 1,033 queries | 30 s | ≈ 29 ms per query |

**Diagnosis.** The flags were right (`--gpu 1` switches the prefilter to `--prefilter-mode 1`, the GPU ungapped
filter, automatically). The cost model fitted from the runs above is ≈ **8.8 ms fixed GPU time per query** plus
0.11 ms per million target residues. The fixed per-query cost dominates whenever the target database is small.
MMseqs2-GPU processes one query at a time against the whole (padded) database. That design suits a few queries
against a huge database (UniRef, colabfold), not many-vs-many searches of ~10⁴–10⁶ proteins. Effective throughput
was 0.6 TCUPS on human-vs-human and 3.6 TCUPS against OrthoBench. Concurrent processes do not help, because the
GPU is already saturated. The 251 s was not a warm-up effect.

**Sensitivity** (all non-self pairs, e ≤ 1e-3, `--max-seqs 300`): the GPU gave 946,183 pairs and CPU `-s 7.5` gave
899,109. 849,951 pairs are shared; 96,232 are GPU-only (median e 1e-8) and 49,158 CPU-only (median e 2e-17). Bit
scores of shared pairs are identical (same CPU `align`). The GPU ungapped filter is ~5 % more sensitive. It misses
some strong hits because the ungapped top-300 cut is applied before gapped alignment (large families).

**Implication for Ancestra.** On OrthoBench, S2 searches 251k queries in 100k-protein blocks and S0 searches 251k
genes vs representatives. At ≥ 9–15 ms per query, the MMseqs2-GPU prefilter alone needs ≈ 40–60 min of GPU time
per pass. DIAMOND needs 13 min (S2, 10 threads) and 2 min (S0 stage 2), so MMseqs2-GPU is **not** a useful
S0/S2 backend here. It would pay off only for a large target DB, e.g. the 1000-proteome demo with few queries.

**DIAMOND:** 2.2.8 has no GPU code path (`diamond help` and the binary have no CUDA/GPU strings). The GPU DIAMOND
announced by its authors is not in the open-source release.

## 2. Where the time goes (OrthoBench, CPU path, frozen `accurate` preset)

From the frozen main run (`work/ob_final/main_accurate_go`, 10 threads, shared node at load 23–31). Its manifest
holds the per-stage numbers:

| stage | wall | CPU | share of wall | GPU candidate |
|---|---|---|---|---|
| S0 families (DIAMOND cluster + genes vs reps + CC/Leiden) | 259 s | 1,798 s | 14 % | (a) MMseqs2-GPU: no, see §1; (e) pLM kNN: not equivalent |
| **S2 scores** (DIAMOND 3 blocks 797 s, MMseqs2 rescue of 4.4 M missed pairs 29 s, cross-family rescue re-score 248 s) | **1,129 s** | 8,516 s | **60 %** | **(b) GPU Smith-Waterman of the explicit within-family pairs** |
| merge (two-pass synteny; FAMSA+FastTree on 15,240 clusters, `tree_s` 441 s) | 491 s | 3,745 s | 26 % | (c) GPU distances/NJ instead of FastTree |
| output | 11 s | 10 s | 1 % | – |
| **total** | **1,893 s** | 14,074 s (3.9 CPU-h) | | |

With the `fast` preset the merge takes 56 s, so S2 is ≈ 75 % of the wall time. On QfO the shares are similar:
S2 8.9 CPU-h, S0 3.0 CPU-h, and the accurate merge 7 CPU-h.

**Decision.** S2 is the target with the largest end-to-end effect for both presets. S2 is "score the pairs inside
each family". That workload fits a GPU directly: Σ_f g_f(g_f+1)/2 = 207 M pairs and 7.5·10¹³ DP cells on the
OrthoBench S0 partition. Half of the cells are in the 6 families larger than 5,000 genes. The cell count follows
from the families and lengths, so the cost is predictable. ProstT5 (d) is a separate capability (complete 3Di
coverage) rather than a speed-up of the frozen method. (c) only affects the `accurate` preset. (e) changes S0
itself.

## 3. GPU Smith-Waterman S2 backend (`--backend gpusw`, opt-in)

**Code.**
- `ancestra/gpu_sw.py`: a CUDA kernel via CuPy NVRTC. It runs in a subprocess with `conda_envs/gpu` (or
  `$ANCESTRA_GPU_PYTHON`), so the CPU env needs no CUDA packages.
- `ancestra/scoring.py`: `_backend_gpusw` and the `gpu_*` config fields. `ScoringConfig.key()` is unchanged for
  the CPU backends; this was checked against `freeze-v1`, so existing caches stay valid.
- `ancestra/pipeline.py`: the `--backend gpusw` choice. Default behaviour is unchanged.

**Algorithm.**
- Full Smith-Waterman/Gotoh with BLOSUM62 and gap 11/1 (DIAMOND's defaults), in int32.
- Every unordered within-family pair is aligned, including self pairs, which give the self-scores. This is the
  exhaustive reference that the `hybrid` backend approximates: "DIAMOND + re-align the pairs DIAMOND missed".
- The forward pass gives the score and the end cell. The tie rule is (max score, then smallest row, then smallest
  column), so the result does not depend on scan order.
- A reversed pass on the kept pairs gives the start cell, and from it the coverages. Bit score uses
  λ = 0.267, K = 0.041. The e-value uses a fixed DB size of one DIAMOND block (100k × mean length), and the
  threshold is e ≤ 10 as in the CPU path.
- The cap rule for families over 5,000 genes and the ScoreTable format are the same code as on the CPU path.

**Kernels.**
- (i) One thread per pair: the shorter sequence on the rows, the longer one in register strips of 64 columns, and
  a warp-interleaved scratch for the strip boundary. Pairs are sorted by (strips, rows) so that warps are balanced.
- (ii) One warp per pair for long pairs (rows × cols > 4·10⁶): a wavefront over 32 lanes × 32 columns with shuffles.
  It is launched first, on its own stream.
- Short-pair launches are spread round-robin over 4 streams. Without (ii) and the streams, the ~40 titin-class
  proteins (> 8,000 aa) made a single-thread tail: a 30 M-pair chunk took 88.6 s, against 33.7 s with them.
- Throughput on A100-40GB: 600–640 GCUPS forward and 540–610 GCUPS reverse on real OrthoBench chunks, with pair
  generation overlapped in a producer thread.

**Exactness tests.**
- Scores equal a pure-Python Gotoh reference on 683 random and mutated pairs, including non-standard letters.
- Re-scoring the reported [start, end] substrings gives the same score.
- The warp kernel equals the thread kernel on 40k pairs, for scores and all coordinates.
- Strip widths 16, 32, 48, 64, 80, 96 and 128 give identical outputs.

## 4. ProstT5 on the GPU (task d): complete 3Di coverage becomes feasible

Tool: `foldseek createdb IN.fa DB --prostt5-model data/structure/prostt5_weights --gpu 1` (foldseek 10.941cd33,
the same fp16 GGUF ProstT5 as the CPU runs in `docs/structure_feasibility.md`). Work dir: `work/gpuv2/prostt5/`.

| run | proteins | residues | wall | residues/s |
|---|---|---|---|---|
| GPU, 500 random OrthoBench proteins with AFDB 3Di (50–1000 aa) | 500 | 188,947 | 28.7 s (incl. ~3 s model load) | 6,600–7,300 |
| GPU, random 6 % of the OrthoBench proteins without 3Di | 2,957 | 2,412,447 | 396 s | **6,100** |
| CPU (same tool, 8 threads, node load ~25) | 50 | 21,788 | 217 s | 100 |
| CPU (structure_feasibility.md, 12–24 threads, loaded) | 50–67 | 26–32k | ~200 s | 121–161 |

- The GPU run is ≈ 40–60× faster than the CPU runs on this node. It is sequence-by-sequence: GPU ~87 % busy and
  one CPU core used, with no batching in foldseek.
- **Determinism:** two GPU runs give byte-identical `_ss` DBs.
- **GPU vs CPU equivalence** (50 proteins): 99.91 % of residues agree and 34/50 sequences are identical (fp16 GPU
  arithmetic vs CPU). Accuracy against the true AFDB 3Di is the same (0.7031 vs 0.7029).
- **Sanity check of predicted vs true 3Di** (500 proteins with exact-sequence AFDB structures): per-residue
  identity **0.672**, against 0.161 for shuffled pairs and 0.278 for the majority class (D). This is in the
  range ProstT5 reports.

**Cost of completing the 3Di layer** at the measured 6,100 res/s:

| set | proteins | residues | GPU time | CPU estimate (100–250 res/s) |
|---|---|---|---|---|
| OrthoBench proteins without AFDB 3Di (19.5 %) | 48,903 | 38.9 M | **1.8 h** | 43–108 h |
| QfO 2022 proteins without AFDB 3Di (6.4 % of the 988,778 run genes) | 63,237 | 79.7 M | **3.6 h** | 89–221 h |
| all of OrthoBench (ProstT5 for every protein) | 251,378 | 127.8 M | 5.8 h | 142–355 h |
| all of QfO | 988,778 | 445.3 M | 20 h | 495–1,237 h |

A single A100 thus fills the 3Di gaps of either benchmark in a few hours, where the CPU needed days. The
structure layer itself (3Di family rescue) gave no merge gain and is off by default (ONE_PAGE_SUMMARY). So this
widens what can be tested rather than speeding up the frozen method.
