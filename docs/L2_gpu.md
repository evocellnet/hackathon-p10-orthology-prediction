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

## 1. Data
No new external data. Used: the GPU_v2 score tables (`work/gpuv2/ob_gpusw_acc`, `work/gpuv2/qfo_allgpu_fast`,
all-GPU exhaustive SW), the CPU tables (`work/gpuv2/diagA_cpu_norescue` = OrthoBench `hybrid` pass 1;
`work/gpuv2/qfo_cpu_fast/stages/scores/scores_pass1.npz`), the frozen runs (`work/ob_final/main_accurate_go`,
`work/gpuv2/qfo_cpu_fast` = frozen `cand2_fast` byte for byte), AFDB 3Di (`data/structure/*_3di_best.fa`, DATA.md) and the
OrthoBench ProstT5-GPU 3Di gap file of GPU_v2 (`work/gpuv2/prostt5/orthobench_gap_prostt5_3di.fa`, md5 9205cad1…).
New: QfO ProstT5-GPU 3Di for the 63,237 run proteins without AFDB 3Di (`work/l2gpu/prostt5_qfo/`, same weights
`data/structure/prostt5_weights`, foldseek 10.941cd33 `createdb --gpu 1`). Software: pyfamsa 0.7.0 + scoring-matrices
from PyPI (MIT), installed with `pip --target work/l2gpu/pylib` (not into the shared env), retrieved 2026-10-27.

## 2. Task 1 — density-robust input for all-GPU S2

### 2.1 Arms and results (dev only)
Reference = CPU frozen path. `scripts/l2gpu_summarize_t1.py` → `results/l2_gpu/task1_summary.{json,md}`; per-arm runs
`work/l2gpu/{ob,qfo}/<arm>/`; transforms `scripts/l2gpu_sparsify_posthoc.py` (= `ancestra.scoring.sparsify_table`,
the pipeline option `scoring.gpu_sparsify`); density matching `results/l2_gpu/density_match_ob.json` (k = 39,
α = 0.207, x = 0.0528 reproduce DIAMOND's 13.69 M big-family pairs on OrthoBench; on QfO the same values keep
84.5 M / 94.2 M / 68.5 M of 109.1 M pairs vs DIAMOND's 77.5 M — `density_match_qfo_reportonly.json`).
OrthoBench: `accurate` + gene order, Δ vs frozen (OG-F, macro pair F1; paired bootstrap over 34 dev RefOGs).
QfO: `fast` + gene order, TreeFam-A dev (paired bootstrap over 359 dev families). All B = 2000, seed 11.

| arm | OB OG-F | ΔOG-F [95 % CI] | OB macro pair F1 | Δmacro [95 % CI] | TF-A dev TPR / PPV | ΔTPR [CI] | ΔPPV [CI] | E1 | E2 | E3 |
|---|---|---|---|---|---|---|---|---|---|---|
| CPU frozen (reference) | 79.69 | – | 0.8026 | – | 0.6151 / 0.9224 | – | – | | | |
| CPU, no cross-family rescue (A) | 78.74 | −0.95 [−3.7, +0.6] | 0.7960 | −0.0066 [−0.022, +0.002] | 0.6158 / 0.9216 | +0.0007 [−0.006, +0.007] | −0.0008 [−0.003, +0.001] | pass | fail | pass |
| GPU split (`gpu_max_family` 1000, GPU_v2) | 80.19 | +0.51 [−3.5, +3.3] | 0.7996 | −0.0030 [−0.018, +0.010] | 0.6138 / 0.9287 | −0.0013 [−0.007, +0.006] | +0.0062 [−0.003, +0.024] | pass | pass | pass |
| all-GPU, unsparsified (GPU_v2) | 74.36 | −5.32 [−11.9, −0.7] | 0.7896 | −0.0130 [−0.033, +0.005] | 0.6081 / 0.9061 | −0.0069 [−0.025, +0.014] | −0.0163 [−0.052, +0.018] | fail | fail | fail |
| **K39** top-39 / species | 74.38 | −5.31 [−11.9, −0.6] | 0.7896 | −0.0130 [−0.033, +0.005] | 0.6061 / 0.9178 | −0.0090 [−0.028, +0.012] | −0.0047 [−0.036, +0.025] | fail | fail | fail |
| **K10** top-10 / species | 66.03 | −13.65 [−23.8, −2.9] | 0.7913 | −0.0113 [−0.030, +0.006] | 0.6131 / 0.9229 | −0.0019 [−0.033, +0.036] | +0.0005 [−0.026, +0.029] | fail | fail | fail |
| **R207** bits ≥ 0.207·best / species | 71.29 | −8.40 [−14.5, −2.4] | 0.7879 | −0.0146 [−0.035, +0.004] | **0.6366 / 0.9300** | **+0.0215** [−0.003, +0.050] | +0.0075 [−0.012, +0.030] | fail | fail | pass |
| **F0528** norm_max ≥ 0.0528 | 72.07 | −7.61 [−13.8, −2.3] | 0.7916 | −0.0109 [−0.032, +0.007] | 0.6241 / 0.9086 | +0.0090 [−0.018, +0.040] | −0.0139 [−0.040, +0.012] | fail | fail | fail |
| MK39 mutual top-39 (D1) | 63.48 | −16.20 [−26.4, −4.4] | 0.7793 | −0.0233 [−0.055, +0.001] | 0.6104 / 0.9073 | −0.0046 [−0.021, +0.015] | −0.0151 [−0.069, +0.026] | fail | fail | fail |
| SP3000 splice (D1) | 76.91 | −2.78 [−5.9, −0.4] | 0.7876 | −0.0150 [−0.034, +0.002] | 0.6162 / 0.9213 | +0.0011 [−0.010, +0.013] | −0.0011 [−0.027, +0.024] | fail | fail | fail |
| SP5000 splice (D1) | 77.48 | −2.21 [−5.7, +0.1] | 0.7901 | −0.0125 [−0.032, +0.005] | 0.6224 / 0.9216 | +0.0074 [−0.006, +0.023] | −0.0008 [−0.026, +0.024] | fail | fail | fail |
| *control* C-rand, random 57.5 % | 70.97 | −8.72 [−15.7, −2.8] | 0.7283 | −0.0743 [−0.113, −0.041] | 0.5321 / 0.8478 | −0.0829 [−0.107, −0.058] | −0.0746 [−0.189, +0.010] | fail | fail | fail |
| *noise floor* J11: CPU A, bits × e^{0.01·N(0,1)} (D2) | 75.52 | −4.17 [−9.5, −0.1] | 0.7854 | −0.0172 [−0.035, −0.004] | 0.6070 / 0.9226 | −0.0081 [−0.034, +0.012] | +0.0001 [−0.010, +0.009] | fail | fail | fail |
| *noise floor* J23 (D2) | 71.44 | −8.24 [−15.7, −0.7] | 0.7825 | −0.0201 [−0.040, −0.004] | – | – | – | fail | fail | – |
| *noise floor* J47 (D2) | 73.76 | −5.92 [−12.7, +0.4] | 0.7862 | −0.0164 [−0.037, −0.001] | – | – | – | fail | fail | – |

Same arms against the same-rescue CPU reference A (no rescue in any all-GPU arm), OrthoBench:
all-GPU −4.38 OG-F / −0.0064 macro; K39 −4.36 / −0.0064; K10 −12.71 / −0.0047; R207 −7.45 / −0.0081; F0528 −6.67 /
−0.0044; SP3000 −1.83 / −0.0085; SP5000 −1.26 / −0.0060; **1 % jitter of A itself: −3.22 / −0.0106, −7.30 / −0.0136,
−4.98 / −0.0098** (mean −5.2 OG-F, −0.011 macro). QfO vs A (`t1_qfo_vs_cpuA.json`): jitter J11 ΔTPR −0.0088
[−0.035, +0.013], ΔPPV +0.0010; R207 +0.0208 / +0.0083; K10 −0.0026 / +0.0013; all-GPU −0.0076 / −0.0155.

**Verdict (pre-registered rule): no arm passes.** Every all-GPU arm fails E1 (OrthoBench OG-F) and E2 (macro pair F1);
only R207 passes E3 (TreeFam-A dev). The controls show why E1/E2 cannot be met by *any* perturbation of the S2 scores:
**multiplying the CPU table's bit scores by 1 % log-normal noise costs −3 to −7 OG-F and −0.010 to −0.014 macro pair
F1** (three seeds), as much as or more than switching to the full GPU graph (−4.4 / −0.006 vs A). The frozen merge's
OrthoBench dev numbers sit on a tuned optimum of the exact DIAMOND table; every perturbation regresses. In other
words, the all-GPU graph is within the score-noise floor on OrthoBench, but so is a 1 % jitter, and the pre-registered
rule (written before this was known) is stricter than that floor. On QfO (359 families, the larger dev set) the picture
is different: **R207 is the best arm of the whole grid (TPR +0.021, PPV +0.008 vs CPU; above the jitter's −0.009 / +0.001)**,
K10 and the splices are equivalent at point level, unsparsified all-GPU and F0528/MK39 lose PPV (−0.014 to −0.016).
The random-sparsification control (C-rand) is far worse than every rank/ratio rule at the same density (QfO TPR
−0.083, OB macro −0.074): *which* edges are dropped matters, not just how many.

### 2.2 Mechanism (how extra edges change merge decisions)
Diagnostics, all OrthoBench accurate, no rescue (`results/l2_gpu/t1_ob_diag_*.json`, `trace_fam77_*.json`):
1. **Not "weak edges diluting strong decisions" in general.** Top-39 per species removes 10.2 M of the 24.1 M
   big-family GPU pairs and changes the result by +0.02 OG-F / 0.0000 macro: the dropped low-rank edges are never
   decisive. What differs between GPU and DIAMOND are the *mid-rank* paralog edges: DIAMOND misses 11.5 % of a
   query's best hit per species (median 86 bits, norm 0.078, min coverage 0.37 — the query has no ortholog in that
   species, so its best is a distant paralog) and ~20–30 % of ranks 1–20 (`diag_pairset_ob.json`).
2. **Score values matter as much as the pair set.** X1 (GPU scores on DIAMOND's exact pair set) is −9.1 OG-F vs A;
   swapping only the self-scores (X3: CPU table with GPU self-scores; 5 % of big-family genes have GPU self-scores
   ≥ 1.67× DIAMOND's composition-adjusted ones) is −17 OG-F. Max-link (insensitive to added weak edges) does not make
   the two tables agree either (Dmax: CPU 67.1, GPU 71.9). OG-F on 34 RefOGs is dominated by four RefOGs
   (021, 023, 007, 033) whose outcome flips under any of these perturbations.
3. **Trace of RefOG021** (125 genes, the largest weight in OG-F; S0 family 77, 4,056 genes;
   `scripts/l2gpu_trace_family.py`, `scripts/l2gpu_trace_clusters.py`, one-family reruns reproduce the full-run
   outcome: CPU 1 root group of 121 genes, GPU 3 groups [91, 28, 2]). The split happens at *Chordata*, and the rule
   steps (R1–R3) are not what unites RefOG021 on the CPU path: without local trees the CPU table also leaves it in 27
   pieces. It is the **ambiguity clustering of the local-tree step** that differs: on the CPU table all 20 RefOG021
   units at Chordata fall into one ambiguous cluster (28 units) that FastTree resolves into one group. On the GPU table
   one Euteleostomi unit (34 genes, 28 of RefOG021) has **degree 281 instead of ~50–60**, its best avg-link partner is a
   non-RefOG021 singleton (0.142) while its best avg-link score to any RefOG021 unit drops to 0.056 (CPU-path analogue
   0.109): the many extra member pairs to its true partners are weaker than the few DIAMOND kept (avg over *present*
   pairs is diluted), while a sparsely connected singleton keeps a high average. So the unit leaves the RefOG021
   cluster and lands in a 57-unit cluster of other paralogs; a third piece ends in a 61-unit cluster above
   `tree_max_units` = 60 (raising the cap to 150 does not reunite it).
   **Mechanism in one sentence:** average linkage over *present* leaf pairs rewards partners connected by few strong
   pairs; exhaustive scoring makes presence uniform for close homologs but leaves it sparse for distant ones, which
   reorders best partners of multi-gene units and therefore the near-best (τ = 0.9) ambiguity clusters that the
   local-tree step resolves.
4. Implication: a fix on the input side must reproduce DIAMOND's *joint* pattern of presence and score scale (seeding,
   masking, composition-based statistics), not a density. Rank/ratio rules help on QfO (R207, K10) because they cut the
   mid-rank paralog edges per species, but cannot reproduce the OrthoBench optimum. A density-robust *merge*
   (e.g. avg-link with missing pairs imputed or weighted, or best-k-mean link) is the principled fix, but it is a
   change of the frozen method and would need its own dev cycle.

### 2.3 Speed / accuracy Pareto (QfO 2022, `fast`, S2 + merge + output; shared node)
| configuration | wall | CPU-h | GPU-h | TF-A dev TPR / PPV | OB dev OG-F / macro |
|---|---|---|---|---|---|
| CPU `hybrid` (frozen) | 3,228 s | 9.05 | – | 0.615 / 0.922 | 79.69 / 0.803 |
| GPU split `gpu_max_family` 1000 | 2,737 s | 7.33 | ≤ 0.33 | 0.614 / 0.929 | 80.19 / 0.800 |
| splice T = 3000 (DIAMOND only for families > 3000)* | ≈ 2,600 s | ≈ 5.4–6.5 | ≤ 0.33 | 0.616 / 0.921 | 76.91 / 0.788 |
| splice T = 5000* | ≈ 2,400 s | ≈ 3.8–5.2 | ≤ 0.33 | 0.622 / 0.922 | 77.48 / 0.790 |
| all-GPU + **R207** | ≈ 2,220 s | **0.59** | 0.18–0.33 | **0.637 / 0.930** | 71.29 / 0.788 |
| all-GPU + K10 | ≈ 2,220 s | 0.59 | 0.18–0.33 | 0.613 / 0.923 | 66.03 / 0.791 |
| all-GPU unsparsified | 2,083 s | 0.55 | 0.18–0.33 | 0.608 / 0.906 | 74.36 / 0.790 |

Measured rows from GPU_v2 §9 (16 threads). The sparsification adds 134–142 s single-thread on QfO (24 s on OrthoBench)
in the post-hoc script (`tables/*.npz.json`). GPU-h: 1,202 s on a shared A100 (≈ 650 s uncontended). *Splice rows are
accuracy from post-hoc splices without cross-family rescue; their cost is an estimate: all-GPU + the DIAMOND share of
the families above T, bounded by their share of DIAMOND rows (0.61 / 0.16) and of Σg² (0.88 / 0.68) of the 6.8
DIAMOND CPU-h of the split run.
