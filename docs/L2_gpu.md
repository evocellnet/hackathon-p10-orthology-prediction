# L2-gpu: the GPU ideas dropped from GPU_v2 (density-robust all-GPU S2, batched local trees, GPU structure rescue)

Agent `gpu`, branch `l2-gpu` (from `v2-gpu` 1f9a343), worktree `work/l2_gpu_wt`, run dirs `work/l2gpu/` (gitignored).
Started 2026-10-27 18:00 CET (`scripts/clock.sh`). Rules: docs/L2_protocol.md (dev only, pre-registration, ≤ 8 threads,
`nice -n 5`, GPU ≤ 18 GB, shared with agent `plm`). The frozen presets are unchanged; every new behaviour is opt-in.

## 0. Pre-registration (written 18:10 CET, before any new arm was run)

### 0.0 Label-free diagnosis done before writing this section
`scripts/l2gpu_pairset_diag.py` → `results/l2_gpu/diag_pairset_ob.json` (OrthoBench, families > 1000 genes, no labels).
DIAMOND (CPU `hybrid` pass 1) keeps 13.69 M pairs there, exhaustive GPU SW 24.07 M; 13.48 M shared. DIAMOND's set
is **not** a top-k set: it contains only 88.5 % of each query's best GPU hit per target species (rank 0), 79–82 %
of ranks 1–2, and 38 % of ranks 100–200. Its inclusion probability rises with the self-normalised score (28 % at
norm_max < 0.02, 98 % at > 0.5). Best Jaccard of a rule with DIAMOND's set: bits ≥ 45 (0.67, GPU_v2), norm_max ≥ 0.02
(0.57), top-100 per species (0.55), ratio-to-best ≥ 0.2 (0.49). So reproducing DIAMOND's set is not the goal; the
goal is a GPU graph on which the frozen merge is as accurate as on the CPU graph.

Existing per-RefOG view of the loss (GPU_v2 runs, already reported there): spliced "CPU small + GPU big" (G) vs
CPU A: RefOG021 (125 genes, largest weight in OG-F) loses TP 7,261 → 2,485 (split), RefOG033 FP 324 → 1,368,
RefOG006 FN 0 → 223, RefOG052 FN 0 → 26, RefOG008 FP 0 → 28; RefOG002/004/005/023/065 gain FP. OG-F is dominated by
RefOG021; macro pair F1 moves little (−0.005, CI includes 0).

### 0.1 Task 1: density-robust input for the all-GPU S2 (equivalence task)
**Hypothesis.** The all-GPU loss comes from avg-link dilution: the merge averages over *present* leaf-pair edges
(`es/ec`), so exhaustive scoring adds many weak paralog edges inside large families that lower avg-link scores
of multi-gene units unevenly and change mutual-best / in-paralog / repair decisions. A sparsification of the GPU graph
that keeps the edges the merge actually relies on (a query's strongest partners per species) should restore the
CPU-path accuracy without DIAMOND.

**Mechanism diagnostic (not a candidate):** M-max — merge with `link=max` (insensitive to added weak edges) on
the CPU pass-1 table (A) and on the all-GPU table. If Δ(GPU − CPU) shrinks under max-link vs avg-link, dilution
is the mechanism. Plus per-node merge statistics (R1 unions/blocked, R2 absorbed, R3 repaired, tree clusters)
of all-GPU vs A, and a per-RefOG trace of RefOG021.

**Arms.** Each is a deterministic transform of the all-GPU score table (`gpu_max_family=0`), applied only to the
non-self pairs of families with > 1000 genes (= `rescue_max_family`; below that the CPU `hybrid` graph is also
exhaustive). A pair is kept if *either* direction keeps it (as in DIAMOND's two-direction union and `apply_cap`).
Everything downstream is the frozen preset (OrthoBench `accurate` + gene order; QfO `fast` + gene order).
- **K-dm**: top-k per (query, target species) by bit score (same-species included as its own "species"), k = the
  smallest k whose kept count on OrthoBench ≥ DIAMOND's count in the same families (density-matched, label-free).
- **K10**: same with k = 10 (sparser point for the Pareto).
- **R-dm**: bits ≥ α · best(query, target species), α density-matched as above (largest α with kept ≥ DIAMOND's count).
- **F-dm**: self-normalised floor norm_max = bits / max(self_a, self_b) ≥ x, x density-matched as above ("weak =
  missing in aggregation").
- **Control C-rand**: random subset of the same big-family pairs at the K-dm density (seed 11; seeds 23/47 if
  C-rand were to pass) — tests whether *which* edges are dropped matters, not just density.
- Reference arms (existing runs, not rerun): CPU frozen path (OrthoBench `work/ob_final/main_accurate_go`; QfO
  `work/gpuv2/qfo_cpu_fast`, byte-identical to frozen `cand2_fast`), all-GPU unsparsified (OrthoBench
  `work/gpuv2/ob_gpusw_acc`; QfO `work/gpuv2/qfo_allgpu_fast`).
The k, α, x values are fixed on OrthoBench and applied unchanged to QfO. "Composition-corrected GPU e-values" is
not re-run: GPU_v2 I1 (Hauser-2016 bias correction) already showed no dev change.

**Equivalence rule (vs the CPU frozen path; all must hold):**
- E1 OrthoBench dev OG-F: Δ ≥ −1.0 and paired-bootstrap 95 % CI lower bound ≥ −4.0 (B = 2000, seed 11, over RefOGs).
- E2 OrthoBench dev macro pair F1: Δ ≥ −0.005 and CI lower bound ≥ −0.02.
- E3 TreeFam-A dev: ΔTPR ≥ −0.01 and ΔPPV ≥ −0.01 (points), paired bootstrap over dev families (B = 2000, seed 11)
  CI lower bounds ≥ −0.02.
- E4 `verify_orthoxml` valid.
- E5 cost: QfO S2+merge+output CPU-h ≤ 0.5 × the CPU path's (9.05 CPU-h) — i.e. most of the 16× saving is kept.
If several arms pass: the one with the highest TreeFam-A dev TPR + PPV; ties → higher OrthoBench macro pair F1.
C-rand passing would mean density alone explains the loss (reported, not adopted). An arm that passes gets a pipeline
option (`ScoringConfig.gpu_sparsify`, default off) and a frozen config; the in-pipeline run must reproduce the
post-hoc table byte-for-byte.

### 0.2 Task 2: batched FastTree / in-process FAMSA (equivalence task, CPU)
**Arms.** (a) one FastTree process per batch of clusters (`FastTree -n N`, alignments concatenated in submission
order) with FAMSA per cluster as now; (b) (a) + in-process FAMSA (pyfamsa, if installable) instead of a subprocess.
Same `tree_ft_workers`, `-t 1`.
**Rule.** Primary: OrthoBench `accurate` (cached S0/S2 of the frozen run) `hogs.orthoxml` byte-identical to the
frozen output. If not byte-identical: Newick per cluster compared; a non-identical variant is acceptable only with
E1/E2 above (Δ within noise) and is then reported as "equivalent, not identical". Speed: merge-stage wall and
CPU time vs the current path, same threads, back-to-back.

### 0.3 Task 3: structure rescue (R4) at full 3Di coverage with GPU ProstT5 (improvement task)
**Arms** (all on the frozen S0/S2 caches, R4 code of docs/S0_rescue.md unchanged: 3Di ∪ profile rescue + RBH merge
guard): R0 = off; R4-afdb = AFDB 3Di only (as in S0_rescue); **R4-full** = AFDB 3Di + ProstT5-GPU 3Di for every
protein without an AFDB structure; R4-shuf = R4-full with each protein's 3Di string shuffled (seeds 11/23/47).
Foldseek search on GPU if the installed foldseek supports `--gpu 1` for search (checked and reported), else CPU.
**Primary dev metric:** OrthoBench dev macro pair F1 (the most stable OrthoBench metric). **Keep rule (protocol):**
Δ(R4-full − R0) paired-bootstrap 95 % CI lower bound > 0; R4-full beats R4-shuf by more than the shuffled seed
spread; guards: TreeFam-A dev PPV not worse by > 0.01, OrthoBench dev OG-F not worse by > 1.0, OrthoXML valid.
Secondary: R4-full vs R4-afdb (does completing coverage add anything?), family completeness on dev RefOGs, cost
(GPU-h for ProstT5, search CPU-h/GPU-h, guard merge).
QfO arm only if the ProstT5 gap (≈ 63k proteins, ≈ 3.6 GPU-h) fits next to the shared GPU; otherwise QfO is run with
AFDB 3Di (93.6 % coverage) and reported as such.

**Order:** task 1 first (no GPU needed: post-hoc transforms of existing GPU tables), task 2 in parallel on CPU,
task 3 after.

### 0.4 Deviation log
- **D1 (19:25 CET, after the first five OrthoBench task-1 arms and the RefOG021 trace were seen; logged before running
  any of the arms below).** All five pre-registered OrthoBench arms fail E1 (§2). The trace (§2.2) shows the
  loss is not "weak edges diluting strong decisions" in general (top-39 removes 10.2 M of 24.1 M big-family pairs and
  changes nothing) but the uneven *presence* of pairs between multi-gene units. Round 2 adds, under the same
  equivalence rule and selection rule as §0.1:
  - **MK39**: mutual top-39 per (query, target species) (both directions must keep the pair);
  - **Pareto splices SP3000 / SP5000**: GPU table for families ≤ T genes, DIAMOND (CPU pass-1) table for larger
    families, no cross-family rescue (a conservative stand-in for `gpusw` with `gpu_max_family = T`; cost estimated
    from measured DIAMOND block times, stated as an estimate).
  QfO runs of the round-1 arms continue unchanged (they were queued before any result was seen).
- **D2 (19:40 CET, after X1/X2 were seen).** X1 (GPU scores on DIAMOND's exact pair set) is as far from the CPU
  path as the full GPU graph (OG-F −9.1 vs A), so OrthoBench OG-F reacts strongly to *any* score perturbation in
  large families. To calibrate what "within noise" means for this kind of change, a **score-jitter noise floor** is
  added: the CPU pass-1 table (A) with every non-self bit score multiplied by exp(0.01·N(0,1)) (1 % noise, far below
  the GPU/DIAMOND score differences), seeds 11/23/47 on OrthoBench (accurate) and seed 11 on QfO (fast, CPU
  pass-1 table). Reported as a control next to every arm; the pre-registered rule E1–E5 is not changed.
