# GPU speed-ups for Ancestra (post-freeze v2 work, branch `v2-gpu`)

> Copy of `docs/GPU_v2.md` from branch `v2-gpu` (commit e3e8ffd). Code lives on that branch only; `main` keeps the frozen method.

Status: done for this time-box (2026-10-27, 11:11–15:45 CET, T+19h34 → T+25h05 by `scripts/clock.sh`).
- Node: 48 CPUs, 500 GB, 1× A100-PCIE-40GB (driver 580, CUDA 13.0). This work used the GPU plus ≤ 16 CPU threads.
  Other agents ran on the rest of the node (load 20–30), so all wall times are from a shared node.
- Rules followed: `ancestra/` at tag `freeze-v1` is unchanged in behaviour. All GPU code is opt-in and lives on
  branch `v2-gpu` (git worktree `work/v2gpu_wt`).
- Accuracy is checked only on OrthoBench DEV RefOGs (`splits/orthobench_split.json`) and TreeFam-A dev, plus
  byte-equivalence of the CPU path.
- Run directories: `work/gpuv2/` (gitignored).

## Summary

| GPU option | speed effect | accuracy vs CPU path (DEV) | verdict |
|---|---|---|---|
| (a) MMseqs2-GPU (`--gpu 1`) for S0/S2 | ≈ 9 ms fixed GPU time per query, so slower than DIAMOND for many-vs-many searches (§1) | ~5 % more sensitive than CPU `-s 7.5` | **not useful** at these sizes |
| (b) **`--backend gpusw`** (default split: families ≤ 1000 genes exhaustive SW on GPU, larger ones DIAMOND concurrently) | S2 1.32× / e2e 1.18× (OrthoBench accurate, 16 thr); S2 1.21× / e2e 1.18× (QfO fast) | OrthoBench OG-F **80.2 vs 79.7** (Δ +0.5 [−3.5, +3.3]); macro pair F1 −0.003 [−0.018, +0.010]; QfO TreeFam-A dev TPR/PPV **0.614/0.929 vs 0.615/0.922** | **equivalent within noise, deterministic**; modest gain |
| (b') `gpusw` all-GPU (`gpu_max_family=0`) | S2 2.97× / e2e 1.57× (OrthoBench accurate); QfO 16× fewer CPU-h | OrthoBench OG-F −5.3 [−12.0, −0.7]; QfO PPV −0.016 | **not equivalent**; opt-in for speed only |
| (c) GPU instead of FAMSA+FastTree local trees | cost is per-process fixed time (0.18 s per call), which a GPU cannot remove (§8) | – | not pursued; CPU batching recommended |
| (d) ProstT5 on GPU (foldseek `--gpu 1`) | 6,100–7,300 res/s vs 100–250 on CPU (40–60×); the OrthoBench 3Di gap was **done in 1 h 56 min** (100 % coverage); the QfO gap ≈ 3.6 h | 99.91 % residue agreement with CPU ProstT5; 0.672 identity to true 3Di; deterministic | **feasible**: complete 3Di coverage in hours (§4) |
| (e) ESM-2 + FAISS-GPU kNN S0 candidates | 592 s embedding + 4 s kNN (OrthoBench) | DEV completeness upper bound 0.853 at k = 25 (sequence ceiling 0.82), but 42 % of neighbours are outside the own family | research lead; not integrated (§5) |

**Main lesson.** S2 dominates the CPU wall time (60–75 %), and its GPU part is cheap: all families ≤ 1000 genes take
12 s on OrthoBench and 32 s on QfO. But **78–85 % of DIAMOND's time is the few families with more than 1,000 genes**
(30 on OrthoBench, 145 on QfO). There the frozen merge depends on DIAMOND's *sparse* hit graph.
- Exhaustive GPU scoring of those families adds ~10 M weak paralog edges. These dilute average linkage and cost
  4–8 OG-F (recall).
- No score threshold, masking, cap rule or composition correction reproduces DIAMOND's seed-based sparsity
  (§6, §7). The best pair-set Jaccard is 0.67.
- A large *equivalent* GPU speed-up therefore needs either a merge that is robust to graph density (a v2 method
  change, to be developed on dev) or a GPU seed-and-extend search that matches DIAMOND. Neither fits a 12-h box.

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

**Usage.** `python -m ancestra run ... --backend gpusw`. The options go in `--config '{"scoring": {...}}'`:

| option | default | meaning |
|---|---|---|
| `gpu_max_family` | 1000 | families with ≤ this many genes are scored on the GPU; larger ones go through DIAMOND on the CPU, concurrently; 0 = all on the GPU |
| `gpu_mask` | `none` | `tantan` soft-masks low-complexity regions to X |
| `gpu_comp_bias` | false | Hauser-2016 composition-bias correction, both directions |
| `gpu_db_residues` | 0 | e-value DB size; 0 = one DIAMOND block |
| `gpu_device` | 0 | CUDA device |

The defaults are the configuration that §6–7 found accuracy-equivalent on OrthoBench dev.

**Algorithm (GPU part).**
- Full Smith-Waterman/Gotoh with BLOSUM62 and gap 11/1 (DIAMOND's defaults), in int32.
- Every unordered pair within a GPU-scored family is aligned, including self pairs, which give the self-scores.
  This is the exhaustive reference that `hybrid` approximates with "DIAMOND + re-align the pairs DIAMOND missed".
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

**Done: the OrthoBench 3Di gap is filled.** All 48,903 OrthoBench proteins without an AFDB 3Di were predicted
with foldseek `--gpu 1`.
- 38.9 M residues in **1 h 56 min** wall (5,580 res/s). For ~25 min of that the GPU was shared with other GPU jobs.
- All lengths match the input.
- Output: `work/gpuv2/prostt5/orthobench_gap_prostt5_3di.fa` (md5 9205cad1…).
- Together with `data/structure/orthobench_3di_best.fa` (AFDB, 80.5 %) this gives **100 % 3Di coverage of
  OrthoBench**. The predicted part is ProstT5 3Di, about 67 % identical to structure-derived 3Di. It is not copied
  into `data/`; that is for the PI to decide (DATA.md entry needed).

**Batched alternative** (`scripts/gpu_prostt5_batched.py`): HF `Rostlab/ProstT5_fp16` encoder plus the authors'
encoder-only CNN 3Di head (`cnn_chkpnt/model.pt` from github.com/mheinzinger/ProstT5). Both were downloaded on
2026-10-27 to `data/models/`. Batches are length-sorted, with ≤ 8,000 tokens per batch.
- Speed: **11,900 res/s** on the 500-protein validation set, 1.7× foldseek's per-sequence GPU path. On the
  ≤ 4,000-aa gap sample (2,911 proteins, 2.1 M residues) it reaches 7,900 res/s, 1.3×.
- Its output agrees with foldseek-GPU on 99.64 % of residues, with the same accuracy against the true 3Di (0.671 vs
  0.672).
- Limitations:
  - Proteins longer than about 4,000 aa run out of memory: full attention needs ~31 GB at the longest gap protein.
    They are 46 of 2,957 proteins and 11 % of the residues of the gap sample, and must go through foldseek.
  - fp16 outputs depend on batch composition: max-tokens 8,000 vs 16,000 differ at a few residues. A fixed batch
    plan is reproducible; foldseek's per-sequence path is batch-free.
- With foldseek for the long proteins, the OrthoBench gap would take ≈ 1.4 h instead of 1.8 h. The gain is real but
  small; foldseek `--gpu 1` is the simpler tool.

A single A100 thus fills the 3Di gaps of either benchmark in a few hours, where the CPU needed days. The
structure layer itself (3Di family rescue) gave no merge gain and is off by default (ONE_PAGE_SUMMARY). So this
widens what can be tested rather than speeding up the frozen method.

## 5. pLM embeddings + FAISS-GPU kNN as an S0 candidate generator (task e; DEV only)

Script: `scripts/gpu_esm_knn.py` (conda_envs/gpu). It computes ESM-2 (`esm2_t12_35M_UR50D`, fp16) mean-pooled
embeddings of all 251,378 OrthoBench proteins (> 1,022 aa as length-weighted means of 1,022-aa windows), L2-normalises
them, and runs exact inner-product kNN with `faiss.GpuIndexFlatIP`. Its output is evaluated against the frozen S0
partition on the 34 DEV RefOGs. Output: `work/gpuv2/esm_knn/summary_esm2_t12_35M_UR50D.json`.

- **Cost:** embedding took 592 s (216k residues/s; the GPU was shared with an S2 job for the first ~3 min). The
  all-vs-all kNN (k = 51) took **4.4 s**, and a repeat search was identical.
- **Reach of the strays.** The S0 partition puts 26/34 DEV RefOGs in one family; the other 8 RefOGs leave 23
  "stray" genes outside their RefOG's main family. "Strays reached" counts strays with ≥ 1 kNN neighbour in that
  main family. The completeness columns count RefOGs whose strays are *all* reached; this is an upper bound,
  because every reached stray is assumed to be joined:

| k | strays reached | completeness upper bound | neighbours in the query's own S0 family | DEV-gene neighbours in another DEV RefOG |
|---|---|---|---|---|
| 1 | 5/23 | 28/34 = 0.824 | 0.835 | 0.0000 |
| 10 | 7/23 | 28/34 = 0.824 | 0.715 | 0.0014 |
| 25 | 12/23 | **29/34 = 0.853** | 0.584 | 0.0018 |
| 50 | 13/23 | 29/34 = 0.853 | 0.485 | 0.0036 |

- **Reading.** Embedding neighbours reach the sequence-search ceiling (0.82, docs/S0_families.md) at k = 1, and
  the 3Di level (0.853) at k = 25. At k = 25, however, 42 % of all neighbours fall outside the query's own family.
  A usable S0 rescue would therefore need a verification step, so that it does not fuse families.
- The S0 ceiling analysis showed that 6 RefOGs have strays with no e ≤ 1e-3 sequence hit to their main family.
  A sequence-based verification (e.g. the GPU SW kernel of §3 on the cross-family kNN pairs) would therefore remove
  most of the gain. An embedding-similarity verification would need a shuffled control and a keep rule (LESSONS #9).
- **Verdict:** feasible and cheap on the GPU (10 min + 5 s). The reachable completeness is the same as the 3Di
  layer's, which also gave no merge gain. Not integrated into the pipeline; it is a v2 research lead, not a speed-up.

## 6. Accuracy of GPU scoring on OrthoBench DEV: an exhaustive score graph is *not* equivalent

The DEV metrics come from `scripts/gpu_dev_eval.py`:
- OG-F/P/R: official `benchmark.py` on the 34 DEV RefOGs.
- Pair F1 micro and macro (mean over RefOGs) from `hogs.orthoxml`.
- Paired bootstrap over RefOGs (B = 2000, seed 11) for ΔOG-F and Δmacro pair F1.

All runs use the frozen S0 partition and the `accurate` preset with gene order. Outputs: `work/gpuv2/eval_*.json`.

| run | S2 scoring | OG-F | OG-P | OG-R | pair F1 micro / macro | ΔOG-F [95 % CI] | Δmacro [95 % CI] |
|---|---|---|---|---|---|---|---|
| CPU frozen (`main_accurate_go`) | hybrid + cross-family rescue | **79.69** | 75.7 | 84.1 | 0.804 / 0.803 | ref | ref |
| A: CPU, rescue off (pass-1 scores) | hybrid | 78.74 | 73.7 | 84.6 | 0.791 / 0.796 | −0.95 [−3.7, +0.6] | −0.007 [−0.022, +0.002] |
| GPU all families, unmasked | exhaustive SW | 74.36 | 73.5 | 75.2 | 0.784 / 0.790 | −5.32 [−12.0, −0.7] | −0.013 [−0.033, +0.005] |
| B: GPU all, on CPU-joined families | exhaustive SW | 73.93 | 68.4 | 80.4 | 0.782 / 0.786 | −5.76 [−11.1, −1.8] | −0.017 [−0.038, +0.003] |
| C: GPU all, tantan-masked | exhaustive SW | 72.02 | 69.0 | 75.3 | 0.777 / 0.782 | −7.66 [−13.3, −1.8] | −0.021 [−0.040, −0.005] |
| D: C + e ≤ 1e-3 in families > 1000 | exhaustive SW | 73.74 | 72.0 | 75.5 | 0.791 / 0.790 | −5.95 [−12.1, −1.2] | −0.012 [−0.031, +0.002] |
| E: unmasked + e ≤ 1e-3 in families > 1000 | exhaustive SW | 73.92 | 69.5 | 78.9 | 0.784 / 0.792 | −5.77 [−11.8, −0.9] | −0.011 [−0.030, +0.006] |
| H: unmasked + cap rule from 1000 genes (`cap_all` 1000) | exhaustive SW | 71.37 | 64.6 | 79.8 | 0.781 / 0.787 | −8.31 [−15.1, −2.8] | −0.015 [−0.037, +0.004] |
| I1: unmasked + composition-bias correction (both directions, best kept) | exhaustive SW | 74.46 | 70.6 | 78.8 | 0.791 / 0.792 | −5.23 [−11.1, −0.2] | −0.011 [−0.030, +0.007] |

**Where the loss comes from.** Score tables were spliced by family size, with no rescue and the reference = A:

| spliced table | OG-F | OG-R | ΔOG-F vs A [95 % CI] | Δmacro vs A [95 % CI] |
|---|---|---|---|---|
| F: GPU for families ≤ 1000, CPU (DIAMOND) for larger | 77.22 | 84.7 | −1.52 [−3.5, −0.2] | −0.004 [−0.019, +0.009] |
| G: CPU for families ≤ 1000, GPU for larger | 75.71 | 75.1 | −3.03 [−9.4, +0.7] | −0.005 [−0.021, +0.011] |
| GPU everywhere | 74.36 | 75.2 | −4.38 [−11.0, −0.2] | −0.006 [−0.027, +0.013] |

**Interpretation.**
1. **Small families (≤ 1000) look the same.** There the CPU `hybrid` is already exhaustive (DIAMOND hits plus
   MMseqs2 SW of every missed pair). The GPU reproduces its pair set: 5.59 M of 5.62 M CPU pairs, plus 0.19 M weak
   GPU-only pairs (median 34 bits).
2. **Big families (> 1000) are where the GPU differs.** There the CPU graph is *DIAMOND-only*. Exhaustive SW adds
   10.3 M pairs in families of 1,001–5,000 genes (+80 %). Bit thresholds cannot reproduce DIAMOND's seeding: the best
   Jaccard is 0.675, at 45 bits. The extra edges cost recall (OG-R 84.6 → 75.1). Average linkage over *existing*
   edges (merge.py, `es/ec`) is diluted by many weak paralog edges, so true orthologous units stop merging.
3. **Neither masking nor score corrections help.** Tantan masking, an e-value cut, the cap rule from 1000 genes and
   a composition-bias correction all fail to fix it (C, D, E, H, I1). I1 uses the Hauser-2016 local
   composition-bias correction (MMseqs2's `calcLocalAaBiasCorrection`, the method DIAMOND cites for its default
   composition statistics). It is implemented in both kernels, plus an exact GPU port of the per-residue
   corrections. It brings GPU bit scores much closer to DIAMOND's on 1 M sampled DIAMOND pairs:
   - big families: median GPU/DIAMOND ratio 1.049 → 1.000, mean |log ratio| 0.115 → 0.071;
   - small families: mean |log ratio| 0.081 → 0.046.

   Dev accuracy does not change (I1 ≈ unbiased). The loss therefore comes from the *pair set* (graph density), not
   from the bit scores. The frozen merge was
   tuned on DIAMOND-sampled graphs, and the big families are the repeat- and composition-rich ones (ZNF-like).
   There DIAMOND's composition-based statistics also change the bit scores.
4. **The cross-family rescue matters little:** −0.95 OG-F when it is off.

**Design decision (made before any final-config result was seen).** `gpusw` reproduces `hybrid`'s semantics per
family-size class:
- families ≤ `gpu_max_family` = 1000 (= `rescue_max_family`) are exhaustive on the GPU;
- larger families go through exactly the DIAMOND block path, on the CPU, concurrently with the GPU. DIAMOND
  also supplies their cross-family hits, so the rescue is kept for them.
- `gpu_max_family=0` gives the fully exhaustive GPU variant above (faster, not accuracy-equivalent).

## 7. End-to-end on OrthoBench: speed and DEV accuracy of the final `gpusw` (split) backend

The two runs below are fresh e2e runs on the same node, minutes apart: 16 threads, `--preset accurate`, gene order on,
node load 20–30 from other agents.
- CPU: `work/gpuv2/e2e_cpu_acc`.
- GPU: `work/gpuv2/e2e_gpu_acc`, with `--backend gpusw` = families ≤ 1000 genes on the GPU, larger ones via DIAMOND
  concurrently.

The `fast` rows are merge-only reruns on the same S0/S2 stages. Their S0/S2 times are taken from the accurate runs.

| stage | CPU (hybrid) | GPU (gpusw split) | speed-up | all-GPU (`gpu_max_family=0`)† |
|---|---|---|---|---|
| S0 families | 150 s | 143 s | (same code) | 143 s |
| **S2 scores** | **721 s** (DIAMOND 3 blocks 492 s, MMseqs2 17 s, rescue 159 s / 101 joins) | **545 s** (GPU 12 s for 9.1 M pairs / 138k genes, *concurrently* DIAMOND on the 30 families > 1000 genes 382 s, rescue 115 s / 10 joins) | **1.32×** | 243 s (2.97×) |
| merge `accurate` (FastTree trees 273–276 s) | 323 s | 327 s | – | 370 s |
| merge `fast` | 56 s | 60 s | – | ~60 s |
| output | 11 s | 11 s | – | 11 s |
| **total `accurate`** | **1,207 s** (3.85 CPU-h) | **1,027 s** (3.19 CPU-h) | **1.18×** | ≈ 767 s (1.57×) |
| **total `fast`** | ≈ 938 s | ≈ 761 s | **1.23×** | ≈ 458 s (2.05×) |

† All-GPU timings come from `det_nomask_rep` and `ob_gpusw_acc`. That mode is **not** accuracy-equivalent (§6).

**DEV accuracy vs the frozen CPU outputs** (`gpu_dev_eval.py`; frozen = `work/ob_final/{main_accurate_go,fast_go}`):

| preset | CPU OG-F / macro pair F1 | GPU split OG-F / macro pair F1 | ΔOG-F [95 % CI] | Δmacro [95 % CI] |
|---|---|---|---|---|
| accurate | 79.69 / 0.803 | **80.19** / 0.800 | +0.51 [−3.5, +3.3] | −0.003 [−0.018, +0.010] |
| fast | 75.71 / 0.773 | 73.95 / 0.797 | −1.76 [−4.4, +1.0] | +0.023 [+0.000, +0.056] |

- Both are inside the project's noise band (±2–4 OG-F, LESSONS #10), and every OG-F CI contains 0. On `fast`,
  macro pair F1 is better while OG-F is lower: a different trade-off, not a uniform shift.
- **This is equivalence within noise, not identity.** The GPU path aligns small families with plain SW instead of
  DIAMOND+MMseqs2 bit scores. It also gets only 10 cross-family rescue joins instead of 101, because only the
  big-family DIAMOND block sees cross-family hits.

**CPU path unchanged.** The v2-gpu branch's CPU e2e run reproduces the frozen run byte for byte:
- `hogs.orthoxml` md5 cc261e73…;
- `scores.npz` d37a0c13…;
- S0 `families.npy` 407f2ec4… for both the CPU and GPU runs.

**Determinism.**
- The final split config, re-run end-to-end from S0 (`e2e_gpu_acc_rep`, while another GPU job shared the device),
  gives byte-identical `scores.npz` (cafa3114…), `hogs.orthoxml` (f96e1f0d…) and `pairs.tsv` (9caa604e…).
- Exhaustive GPU S2 (`gpu_max_family=0`) run twice with different kernel implementations (the first single-stream
  thread kernel vs the final warp + 4-stream version) gives byte-identical `scores.npz` (md5 1f98bba7…).
- The kernels use int32 only and have a fixed tie rule. Batch size, strip width, streams and the thread/warp kernel
  choice do not change outputs (tests/test_gpu_sw.py).
- The DIAMOND part is the CPU path's code with pinned `-b 4 -c 1`.

**Why the speed-up is modest.** 78 % of the CPU path's DIAMOND time (382 of 492 s) is the block of the 30 families
over 1,000 genes. Those families (97.5k genes, 3.6·10¹³ of the 7.5·10¹³ DP cells) are where exhaustive GPU scoring
is not equivalent (§6), so they stay on DIAMOND. The GPU finishes its share (all families ≤ 1000) in 12 s, so the
GPU is idle > 95 % of S2.

**Can a GPU filter reproduce DIAMOND's sparser graph in big families?** No, at least not with score thresholds. On the
93 M pairs of families with 1,001–5,000 genes (2 min on the GPU), we computed both the gapped and the ungapped
(diagonal) SW score:
- The best Jaccard similarity to DIAMOND's 12.8 M-pair set is 0.67, for an ungapped score ≥ 30–35 bits or a gapped
  score ≥ 45 bits.
- DIAMOND's selection comes from seeding in a reduced alphabet, plus masking and composition statistics; a score
  cut-off does not reproduce it.
- Equivalence in big families would need either DIAMOND itself, or a merge that is robust to the density of the
  score graph. The latter is a method change, i.e. a v2 research question.

## 8. Accurate-preset local trees (task c): not a GPU problem

In the e2e runs above, the `accurate` merge spends 273–276 s of its 323–327 s in FAMSA + FastTree. There is one
job per ambiguous cluster: 15,240 clusters on OrthoBench and 61,216 jobs on QfO, with a median of ~5 units each.
Measured per call on 5 proteins:

| step | time per call |
|---|---|
| `famsa -t 1` | 0.073 s |
| `FastTree -lg -nosupport` | 0.105 s |
| empty process | 0.0008 s |

So the cost is ~0.18 CPU-s per cluster, almost all of it fixed per-call work: model set-up and initialisation
for alignments of 5 sequences. This agrees with S4 (0.24 CPU-s at 5 units).
- A GPU cannot remove a per-process fixed cost.
- A GPU replacement of the trees themselves (batched distances + NJ, or pLM-embedding distances) would be a
  *different* tree method. S4 already showed that NJ-type trees on score distances lose the label-accuracy gain
  (+0.027) that FastTree buys, and that complete-but-rescaled distances (e_fsbits) do not recover it.

**Not implemented.** The cheap, equivalence-preserving fix is on the CPU, and was left as a recommendation:
- run FastTree once per batch of clusters (`FastTree -n N` reads N alignments from one file);
- align in-process (pyfamsa);
- keep `-t 1` and the submission order.
- Byte-identity with the current outputs would still have to be checked: FastTree's RNG state may carry over
  between the alignments of one `-n` run.

## 9. QfO 2022: speed and TreeFam-A DEV equivalence

Setup:
- cached deterministic S0 of `work/qfo2022_det/det_nogo`;
- `--preset fast`, gene order on (= the frozen `cand2_fast` configuration), 16 threads, same node, run back-to-back;
- TreeFam-A DEV pairs scored with `scripts/qfo_treefam_dev_fast.py` (official formula, dev families only).

| | CPU `hybrid` (`work/gpuv2/qfo_cpu_fast`) | GPU `gpusw` split (`work/gpuv2/qfo_gpusw_fast`) | speed-up |
|---|---|---|---|
| S2 scores | 2,790 s (DIAMOND 9 blocks 1,865 s, MMseqs2 77 s, rescue 611 s / 358 joins) | 2,303 s (GPU 32 s for 36.5 M pairs / 384k genes in parallel with DIAMOND 5 blocks of the 145 families > 1000 genes 1,579 s; rescue 520 s / 40 joins) | **1.21×** |
| merge (fast) + output | 389 + 40 s | 385 + 41 s | – |
| S2 + merge + output | 3,228 s, 9.05 CPU-h, 30 GB peak RSS | 2,737 s, 7.33 CPU-h, 42 GB peak RSS | **1.18×** wall, 0.81× CPU-h |
| TreeFam-A DEV TPR / PPV | 0.615 / 0.922 (± 0.011 / 0.008) | **0.614 / 0.929** (± 0.011 / 0.007) | Δ −0.001 / +0.006 |

- The branch's CPU run is byte-identical to the frozen `cand2_fast` (`hogs.orthoxml` e4f6d7ba…), and its
  `scores.npz` is byte-identical to `det_nogo` (65ad834c…, computed at 12 threads). The default CPU path is
  unchanged and thread-independent.
- **TreeFam-A DEV:** the GPU split is equivalent. TPR is equal, and PPV is +0.006, inside the binomial CIs.
- **Speed:** small, as on OrthoBench. 85 % of QfO's DIAMOND time (1,579 of 1,865 s) is the blocks of the 145
  families over 1,000 genes; they hold 472k of 856k family genes and 3.0·10¹⁴ of the 3.1·10¹⁴ DP cells.
- The GPU path peaks at 42 GB, against 30 GB on the CPU. The cause was float64 concatenation of the DIAMOND and GPU
  records; it is fixed to float32 in commit "float32 merge of records" (output identical, verified by the OrthoBench
  rerun above). The QfO run predates the fix.

**All-GPU on QfO** (`gpu_max_family=0`, `work/gpuv2/qfo_allgpu_fast`):
- **Workload:** 1.12·10⁹ within-family pairs and 3.1·10¹⁴ DP cells.
- **GPU time:** 1,202 s (forward 942 s, reverse 187 s). The rate was only 328 GCUPS because a foldseek ProstT5 job
  shared the GPU, so this wall time is an upper bound; at the uncontended 600 GCUPS it would be ≈ 650 s.
- **Totals:** S2 took 1,541 s. S2 + merge + output took 2,083 s wall and **0.55 CPU-h** (user 1,993 s), against
  3,228 s and 9.05 CPU-h on the CPU path: 16× fewer CPU-hours.
- **TreeFam-A DEV:** TPR 0.608 / PPV **0.906**. Against the CPU path that is −0.007 TPR and **−0.016 PPV**, outside
  the CIs.
- **Verdict:** as on OrthoBench, the fully exhaustive score graph is *not* accuracy-equivalent; here it costs
  precision. It remains an opt-in mode for very large inputs, where CPU-hours dominate.
