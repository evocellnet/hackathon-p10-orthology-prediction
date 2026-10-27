# L2-plm: protein language models as an extra layer

Agent `plm`, branch `l2-plm` (worktree `work/l2_plm_wt`), started 2026-10-27 18:01 CET. Post-freeze, DEV only
(docs/L2_protocol.md). Nothing here changes `freeze-v1` or any reported number.

## 0. Pre-registration (written and committed before any dev result was computed)

**Hypothesis.** Protein language model (pLM) representations carry a signal that ranks a gene's ortholog above
its same-species paralog better than the frozen sequence score does, above all on the hard cases where
sequence nearly ties (teleost ohnologs, 2R paralogs, tandem copies), and at deep LCAs. The structure layer failed
because paralogs share folds (S6: deep triplet accuracy 3Di+AA 0.694 vs sequence 0.751); a pLM is trained on
sequence alone and could keep the post-duplication divergence that 3Di loses — or it could behave like 3Di. A
small head trained on *self-derived* ortholog/paralog labels (synteny-supported 1:1 anchors from the project's
gene-order data) may sharpen the pLM signal towards orthology.

### Step (a): standalone discrimination (idea 1, idea 2 evaluation)

*Evaluation set.* Triplets (a, o, p) from the DEV pair labels (`scripts/orthobench_dev_pairs_eval.load_truth`,
34 DEV RefOGs, within-RefOG pairs only): o is an ortholog and p a paralog of a, and o and p are in the same
species. Correct if S(a,o) > S(a,p); ties count 0.5 (= query-conditioned AUC). Main set `T_seq`: both (a,o) and
(a,p) have a sequence score in the frozen ScoreTable (`work/core_v0_ob/run_v2`, norm_max), as in S6. Secondary
set `T_all`: all triplets, a missing sequence score = 0.
Groups: all; **deep** (species-pair class vert-teleost, vert-invert, invert-invert; as S6); shallow (the other
classes); **hard** = triplets whose query (a, species of o) satisfies the S5 rule on the frozen sequence bits:
best paralog bits ≥ 0.9 × best ortholog bits. Caveat fixed in advance: the hard set is selected on the sequence
score, so sequence accuracy on it is depressed by construction (regression to the mean favours any other
score); the shuffled control and the deep group (selected without sequence) are the guards against this.

*Arms (scores S(x,y)).*
- `seq`: frozen sequence score (reference).
- `e35_mean`: cosine of mean-pooled ESM-2 35M (existing `work/gpuv2/esm_knn` embeddings).
- `e650_mean`: cosine of mean-pooled ESM-2 650M (`facebook/esm2_t33_650M_UR50D`), fp16 inference.
- `e3b_mean`: ESM-2 3B (`facebook/esm2_t36_3B_UR50D`), only if it fits the 18 GB GPU budget.
- `pt5_mean`: ProtT5-XL-U50 encoder (`Rostlab/prot_t5_xl_half_uniref50-enc`).
- `e650_aln`, `pt5_aln`: embedding-alignment score in the spirit of EBA / pLM-BLAST: per-residue cosine
  similarity matrix, EBA-style row/column z-score enhancement, local Smith-Waterman on it, score normalised by the
  smaller of the two self-scores. Our own simplified settings (not EBA's or pLM-BLAST's exact algorithm), fixed
  before any dev look: **local SW, linear gap 1.0 in enhanced-score units, zero floor** (one fixed value; no
  tuning on dev).
- `comb_<arm>@0.9`: sequence decides unless the two sequence scores are within 0.9× of each other (the S5 near-tie
  rule); then the pLM arm decides. This is the standalone form of the tie-breaker in idea 3.
- ESM-C is used only if installable without a licence question; otherwise recorded as not pursued.
- Long proteins: > 1,022 aa (ESM-2) are embedded as 1,022-aa windows (mean arm: length-weighted mean; per-residue
  arm: concatenation of windows); ProtT5 takes full length up to 4,000 aa, longer as windows.

*Controls.* Shuffled layer: embeddings permuted among genes of the same species (fixed seeds 11, 23, 47), applied
to every pLM arm (and the per-residue arms) identically.

*Statistics.* Paired bootstrap over RefOGs (B = 2000, seed 11) for accuracy(arm) − accuracy(seq) and
accuracy(arm) − mean accuracy(shuffled 11/23/47). Also pooled pair ROC-AUC with RefOG bootstrap CI, for
comparability with S5, and coverage (fraction of triplets scorable).

*Step-(a) pass rule (decides whether idea 3 is attempted).* For at least one arm: on the **hard** or the
**deep** group, accuracy(arm or comb_arm) − accuracy(seq) has a 95% CI lower bound > 0 **and**
accuracy(arm) − shuffled has CI lower bound > 0. If no arm passes, the layer stops here as a negative result.

### Idea 2: a head trained on self-derived labels ("a new pLM on the extra data")

*Labels (no benchmark labels).*
- Training families: S0 families of the frozen partition (`run_v2 families_final.npy`) that contain **no RefOG
  gene of any split** (dev and test RefOG membership is read only to exclude these families entirely).
- Positives: synteny-supported unambiguous reciprocal best hits between two species inside a family: frozen
  sequence bits, runner-up < 0.9 × best on both sides (S5 anchor rule), and synteny context support ≥ 1 at w = 10
  (S5 `context` definition, genome-wide anchors outside the pair's family; gene order from
  `data/genomes/orthobench`).
- Negatives, per positive (a, b): (i) out-paralog candidates (a, b′), b′ ≠ b in the species of b and the same family,
  with no synteny support; (ii) within-genome paralogs (a, a′), a′ in the species of a and the same family.
  Training objective is ranking: S(a,b) > S(a,b′) — the same form as the DEV triplets.
- Split: training families 90 / 10 train / val by family (seed 11). Validation metric: val triplet accuracy on
  the self-derived labels; early stopping on it. Train/val curves are reported.
*Model.* Frozen mean-pooled embeddings of the best backbone **by val accuracy on self-derived labels** (not by DEV),
among those embedded; a 2-layer MLP projection head (d → 512 → 256, GELU, dropout 0.1) trained with a triplet
InfoNCE loss (positive vs its negatives, temperature 0.1), AdamW lr 1e-3, wd 1e-4, ≤ 30 epochs, batch 256.
Seeds 11 / 23 / 47 (three heads → seed spread). Optional extra input (3Di tokens) only if time allows.
*Controls.* (C1) shuffled embeddings within species at evaluation, as above; (C2) **label-shuffled head**: the same
head trained with positives replaced by random same-family cross-species pairs (seed 11) — does the synteny label
add anything beyond family membership?
*Evaluation.* DEV triplets as in step (a); the same pass rule.

### Step (b): integration (idea 3), only for an arm that passes step (a)

Harness: the S6 merge evaluation (`scripts/S6_merge_eval.py` pattern): run the frozen S3 merge (default arm) on
run_v2's checkpointed inputs with a modified ScoreTable, then pair-level DEV metrics.
Arms (fixed now):
- `tie`: for every scored pair whose score is within 0.9× of a competing partner of the same query in the same
  species, multiply S by (1 + 0.1·z), z = the pLM score standardised within the family (clipped to ±2).
  Pairs outside such near-ties keep S_seq.
- `avg`: S = geometric mean of S_seq and the pLM score calibrated to S_seq (linear fit of −log S_seq on −log S_plm,
  as S6 `combine_scores`).
Each with real and shuffled (11/23/47) pLM scores; `off` = frozen scores.
*Primary metric:* OrthoBench DEV **macro pair F1** (mean over RefOGs). *Guards:* pair PPV (−0.01 max), DEV OG-F
(−1.0 max), OrthoXML validity (`verify_orthoxml`). TreeFam-A DEV TPR/PPV (`fast` preset) as the second keep
metric if an arm passes on OrthoBench.
*Keep rule (protocol §6b):* Δ(real − off) paired-bootstrap 95% CI (B = 2000, seed 11) lower bound > 0; real beats
shuffled by more than the seed spread; no guard worse. A passing arm gets a frozen config for the test agent.

### Idea 4 (optional): embedding kNN with a family guard for S0 completeness
Only if time remains. Guard fixed in advance: a cross-family kNN link (k = 25, ESM-2) is accepted only if it is
mutual and its cosine is above the 99th percentile of shuffled-within-species cosines (seed 11), and a family join
needs ≥ 2 such links. Metric: DEV RefOG completeness and RefOG fusions (RefOG010 must not fuse).

### Fixed in advance
Seeds 11/23/47; B = 2000; no test RefOGs, TreeFam-A test, SwissTree or QfO test; GPU ≤ 18 GB, 4 threads, nice 5.

## 1. Data

All downloads in `data/extra/l2_plm/hf` (Hugging Face hub cache, gitignored), retrieved 2026-10-27; md5 of the
weight files in `data/extra/l2_plm/weights_md5.txt`.

| model | HF repo @ revision | file actually loaded (md5) | licence |
|---|---|---|---|
| ESM-2 650M | `facebook/esm2_t33_650M_UR50D` @ 08e4846e | model.safetensors (9595441349dfbef0f579bd2710f86df3) | MIT |
| ESM-2 3B | `facebook/esm2_t36_3B_UR50D` @ 476b6399 | pytorch_model-0000{1,2}-of-00002.bin (377bf54b…, 333bfff5…) | MIT |
| ProtT5-XL-U50 encoder (fp16) | `Rostlab/prot_t5_xl_half_uniref50-enc` @ 94a6abc0; transformers 5 loaded the safetensors of the hub's auto-conversion `refs/pr/1` (d53b8ce6) | model.safetensors (b1be00dd3b3580fdfbe074a44476b487) | **licence-unclear**: no licence on the HF model card (checked via the hub API); used locally only, flagged |
| ESM-2 35M (reused) | `facebook/esm2_t12_35M_UR50D`, embeddings from `work/gpuv2/esm_knn` (GPU_v2 §5) | – | MIT |

ESM-C (EvolutionaryScale) was not used: it needs the `esm>=3` package (the envs carry fair-esm 2.0) and the
hub lists a custom licence (`license:other` next to `mit` for `esmc-600m-2024-12`); not pursued for time and
licence reasons. Licences were read from the hub API tags on 2026-10-27 (ESM-2: `license:mit`). No other external data:
gene order (`data/genomes/orthobench`, S5), the frozen S0 partition and sequence ScoreTable (`work/core_v0_ob/run_v2`)
and the DEV pair labels (`work/S5/refog_pair_labels.parquet`) are project data.

## 2. Methods (as run)

Code on branch `l2-plm`: `scripts/l2_plm_{common,embed,alnscore,discrim,selflabels,head,knn_guard,complement,summarize}.py`.
Benchmark/exploration code only; nothing in `ancestra/` was changed, so the frozen presets are byte-identical.

- **Embeddings** (`l2_plm_embed.py`, conda_envs/gpu, fp16, last layer). ESM-2: 1,022-aa windows; ProtT5: whole
  protein up to 4,000 aa, longer proteins as windows; the protein vector is the length-weighted mean of residue
  vectors; per-residue matrices are kept for the 957 DEV RefOG genes. Internal `*` are passed as unknown tokens
  (ESM `<unk>`, ProtT5 `X`) so residue counts stay aligned.
- **Embedding alignment** (`l2_plm_alnscore.py`). Cosine matrix of unit residue vectors (fp16 matmul on GPU),
  EBA-style enhancement 0.5·(row z-score + column z-score), local SW with linear gap 1.0 and zero floor, on a CuPy
  kernel (one 256-thread block per pair, anti-diagonal wavefront). Checked against a NumPy reference on 36 pairs
  (max relative difference 4·10⁻⁴, fp16 matmul) and between the earlier thread-per-pair kernel and the wavefront
  kernel on 1,000 pairs (identical). Throughput ≈ 1·10⁹ cells/s on the shared A100. Normalised score =
  raw / min(self(x), self(y)).
- **Triplet discrimination** (`l2_plm_discrim.py`): as pre-registered (§0). 25,562 triplets from 34 DEV RefOGs;
  T_seq 17,750; deep 10,429 (T_all) / 7,110 (T_seq); hard 6,416 / 5,125. Shuffled controls permute embeddings
  among the RefOG genes of the same species (seeds 11/23/47); the alignment arms re-align the permuted pairs.
  Paired bootstrap over RefOGs, B = 2,000, seed 11. Pooled pair ROC-AUC uses 200 RefOG resamples (descriptive).
- **Self-derived labels** (`l2_plm_selflabels.py`, pre-registered mode): 508,666 genome-wide anchors; 78 S0
  families hold a RefOG gene (dev or test) and are excluded entirely; families > 2,000 genes are skipped (cost);
  798 families sampled (seed 11) until 25,011 genes, 88,257 queries, 3.94 M (positive, negative) pairs; val = 79
  families. Overlap of training genes with any RefOG gene: 0.
- **Heads** (`l2_plm_head.py`): as pre-registered; per epoch, each query sees ≤ 16 sampled negatives (all
  negatives are used for val accuracy). Weights `data/extra/l2_plm/heads/*.pt` (md5 in the head JSONs).
- **kNN guard** (`l2_plm_knn_guard.py`): as pre-registered; the null of "shuffled-within-species cosines" was
  operationalised as the cosine of (i, j′), j′ a random gene of j's species, for every mutual cross-family link.

## 3. Deviations from the pre-registration (logged)

- **D1 (sensitivity, e650 only):** a head trained with negatives of type (i) only (`e650_i_real_s11`), because
  type (ii) (same-genome paralogs) does not match the DEV triplet form. Same verdict as the pre-registered head.
- **D2 (post hoc, after seeing that the pre-registered heads do not help on hard triplets):** the pre-registered
  label rule cannot contain sequence-hard examples. The S5 anchor rule demands runner-up < 0.9 × best in the other
  species, so **0 of the 1.96 M type-(i) training triplets have a negative within 0.9× of the positive's bits**
  (12,210 within 0.75×). D2 builds labels where synteny decides a sequence near-tie instead. For query a and
  species s, the near-tie set C holds the same-family genes of s with bits ≥ 0.9 × a's best bits in s. The positive
  is the unique member of C with synteny support ≥ 2 while all others have 0, and the negatives are the rest of C.
  This gave 28,516 queries and 45,815 triplets from 1,786 training families (every training family ≤ 2,000 genes;
  31,231 genes). The 23,756 genes not already embedded were embedded with e3b, the backbone with the best
  pre-registered val accuracy. Heads use seeds 11/23/47. The control replaces the synteny winner by a random
  member of C (seed 11). D2 is exploratory: it was designed after a DEV look, so even a pass would only have
  justified a new pre-registered test.
- Families > 2,000 genes were skipped when building training labels (cost; not stated in §0).
- The kNN-guard null had to be operationalised (§2); the choice was made before the run.

## 4. Results

All tables: `results/l2_plm/tables.md` (generated from `discrimination_final.json`, which holds every arm and
control of this round in one run; the per-step JSONs `discrimination_{mean,aln,head_*}.json` are earlier
subsets with identical numbers). Triplet accuracy on T_seq. Sequence: all 0.866, deep 0.757, shallow 0.939,
hard 0.617, hard∩deep 0.607. Deep sequence accuracy matches S6's 0.751 on its slightly different set.

### 4.1 Step (a), idea 1: pretrained pLMs

| arm (T_seq) | all | deep | hard | hard∩deep | Δ vs seq, deep [95% CI] | Δ vs seq, hard [95% CI] |
|---|---|---|---|---|---|---|
| seq | 0.866 | 0.757 | 0.617 | 0.607 | – | – |
| ESM-2 35M mean | 0.754 | 0.608 | 0.528 | 0.485 | −0.149 [−0.200, −0.099] | −0.089 [−0.143, −0.010] |
| ESM-2 650M mean | 0.765 | 0.632 | 0.558 | 0.532 | −0.125 [−0.184, −0.040] | −0.059 [−0.123, +0.046] |
| ESM-2 3B mean | 0.782 | 0.656 | 0.584 | 0.564 | −0.101 [−0.140, −0.037] | −0.033 [−0.084, +0.053] |
| ProtT5 mean | 0.761 | 0.624 | 0.529 | 0.489 | −0.133 [−0.172, −0.077] | −0.088 [−0.147, −0.007] |
| ESM-2 650M aln | 0.774 | 0.720 | 0.618 | 0.629 | −0.037 [−0.097, +0.044] | +0.001 [−0.062, +0.111] |
| ESM-2 3B aln | 0.827 | 0.715 | 0.608 | 0.604 | −0.042 [−0.113, +0.045] | −0.009 [−0.092, +0.116] |
| ProtT5 aln | 0.731 | 0.715 | 0.620 | 0.643 | −0.042 [−0.082, +0.007] | +0.003 [−0.032, +0.064] |
| seq, ProtT5-aln breaks near-ties (comb@0.9) | 0.868 | 0.766 | 0.625 | 0.630 | +0.010 [−0.006, +0.022] | +0.008 [−0.013, +0.028] |

- **Signal is real.** Every pLM arm is far above its shuffled control. The shuffled means are 0.49–0.52, and
  Δ(real − shuffled) on all triplets is +0.23 to +0.33 with CI lower bounds ≥ +0.18.
- **The signal is weaker than sequence everywhere.** No arm or combination has a CI lower bound > 0 against
  sequence on deep, hard or hard∩deep (`PASSED: []`).
  - The closest is ProtT5-alignment as a near-tie breaker on hard∩deep: +0.022 [−0.001, +0.046]. It beats its
    own shuffled control (+0.048 [+0.007, +0.070]).
  - About 40 arm × group comparisons were run, so one CI touching 0 is what chance would give.
- **What helps.** Scale helps mean-pooled ESM-2 (35M < 650M < 3B). Embedding alignment helps more than scale on
  deep pairs (650M: 0.632 → 0.720). ProtT5 is the worst pooled model but the best aligner on hard∩deep.
- **On the hard set the arms look closer to sequence than they are.** Hard triplets are selected on the sequence
  score, so sequence accuracy there is depressed by construction (pre-registered caveat). The pLM arms reach
  sequence level on hard, yet stay below it on deep, which is selected without sequence.
- **Pooled pair ROC-AUC** (descriptive; raw bits are not comparable across families) agrees: deep AUC is 0.634
  for sequence and 0.49–0.64 for the pLM arms.
- **Why it fails** (`complement_final.json`): pLM errors are correlated with sequence errors. On the 1,944
  hard triplets that sequence gets wrong, the pLM arms are *below* chance: mean-pooled 0.24–0.37, alignment
  0.37–0.41, ProtT5-alignment 0.53 (shuffled ≈ 0.48–0.53). On the triplets sequence gets right they are at
  0.67–0.76.
  A pLM trained on single sequences sees the same post-duplication divergence that the alignment score sees, only
  blurred. Where an out-paralog happens to look more like the query than the true ortholog does, the pLM agrees
  with sequence. This is the sequence analogue of S6's "paralogs share folds".

### 4.2 Idea 2: heads trained on self-derived synteny labels

| head (pre-registered labels) | val acc (raw cosine → head) | DEV all | DEV deep | DEV hard |
|---|---|---|---|---|
| ESM-2 650M, s11 | 0.955 → 0.966 | 0.804 | 0.678 | 0.583 |
| ProtT5, s11 | 0.941 → 0.964 | 0.792 | 0.652 | 0.556 |
| **ESM-2 3B, s11 / s23 / s47** (selected by val) | 0.961 → 0.973 / – / – | 0.822 / 0.814 / 0.818 | 0.680 / 0.667 / 0.690 | 0.577 / 0.559 / 0.588 |
| ESM-2 3B, C2 label-shuffled control | – | 0.789 | 0.631 | 0.552 |
| D2 (deviation), 3B, s11 / s23 / s47 | 0.560 → 0.596 / 0.618 / 0.605 | 0.769 / 0.774 / 0.779 | 0.607 / 0.634 / 0.629 | 0.564 / 0.599 / 0.584 |
| D2 control (random near-tie winner) | 0.484 → 0.512 | 0.757 | 0.616 | 0.554 |

- **The heads improve the frozen embeddings on easy triplets.** e3b head − e3b mean: all +0.040
  [+0.010, +0.056], shallow +0.051 [+0.026, +0.074].
  - For e3b, the synteny label adds information beyond family membership. Head s11 − C2: all +0.032
  [+0.008, +0.050], deep +0.049 [+0.012, +0.080]. All three seeds are > C2 on deep.
  - For e650, the head and C2 are equal (+0.007 [−0.025, +0.033]).
- **They do not help on hard triplets.** e3b head − e3b mean on hard is −0.007 [−0.065, +0.056]. Every head
  stays below sequence on all groups (deep −0.077 [−0.122, −0.031]). Every comb@0.9 variant is ≤ sequence.
- **Why.** The self-derived labels are easy by construction (§3 D2). The training task saturates: train accuracy
  is 0.999 after 1 epoch, and val accuracy with the raw cosine is already 0.955.
- **D2 does not fix it.**
  - Synteny-decided near-ties are hard for a pLM: the raw cosine gets 0.56 on val, the head 0.60–0.62, against
    0.51 for the random-winner control.
  - On DEV, D2 heads are below sequence on all, deep and hard (hard −0.018 to −0.053, CIs include 0).
  - They are not separable from their control (D2 − control on hard: +0.010 to +0.045, all CIs include 0).
  - D2 heads de-correlate pLM and sequence errors: accuracy on hard triplets where sequence is wrong rises
    from 0.37 (e3b mean) to 0.43–0.49. But it does not get above chance.

### 4.3 Idea 4: kNN with a family guard (S0 completeness)

The pre-registered guard is far too permissive.
- Mean-pooled ESM-2 35M cosines are compressed near 1: the null's 99th percentile is τ = 0.973, and 90% of
  the 376k mutual cross-family links pass it.
- 38,027 family pairs are joined, and union-find chains them into **one family of 229,289 of 251,378 genes**.
- DEV completeness rises from 26/34 to 32/34, but with 457 new DEV RefOG–RefOG fusions.
- The within-species-shuffled controls do the same (largest family about 239k, completeness 31–33/34), so the
  completeness gain is not specific.
- Killed. Any stricter guard would be tuned on DEV after this look and would need its own pre-registration.
  GPU_v2 §5 already noted that sequence verification removes most of the reach.

### 4.4 Step (b), idea 3: integration

**Not attempted.** No arm passed the pre-registered step-(a) rule (beat sequence on hard or deep, with a CI lower
bound > 0, and beat shuffled). There is therefore no frozen config for a test agent.

## 5. Cost

| item | wall (shared A100 / node) |
|---|---|
| embeddings, 957 DEV RefOG genes (e650 / e3b / ProtT5, per-residue) | 28 s / 132 s / 146 s; peak 1.9 / 6.1 / 9.9 GB |
| embeddings, 25,011 training genes (e650 / e3b / ProtT5) | 13 / 28 / 52 min (17k / 7.9k / 4.3k res/s, GPU shared with the `gpu` agent) |
| embeddings, 23,756 D2 genes (e3b) | 22 min |
| embedding alignment, 120,730 pairs × 3 models (1.4·10¹¹ cells each) | 2.3–2.7 min each; peak 6.7–9.3 GB |
| 12 heads | 37 min total |
| kNN guard (4 variants, FAISS-GPU) | 2 min |
| **total GPU** | **≈ 3 GPU-h** (wall, shared device; incl. ~15 min of killed/restarted runs) |
| CPU (labels, triplets, bootstrap) | < 1 CPU-h (≤ 4 threads) |

Throughput for comparison: embedding the full OrthoBench (128 M residues) would take about 2 h with ESM-2 650M
alone, and about 4.5 h with 3B, on a dedicated A100.

## 6. Verdict

**Negative.**
- Protein language models carry a strong homology signal: they are far above the shuffled controls.
- They do not separate orthologs from same-species paralogs better than the frozen sequence score. This holds
  for every arm on deep pairs and on the sequence-hard cases (ESM-2 35M/650M/3B and ProtT5; mean-pooled cosine
  or EBA-style embedding alignment).
- Used only to break sequence near-ties, they are neutral (best +0.022 [−0.001, +0.046] on hard∩deep).
- Their errors are positively correlated with sequence errors.
- A contrastive head trained on the project's own synteny-derived labels sharpens the embeddings on easy cases,
  but not where it matters. The label rule that makes anchors trustworthy (unambiguous RBH) is exactly the rule
  that removes the hard cases. Labels built from synteny-decided near-ties (D2) are too noisy or too hard to
  learn from 28k examples.
- Synteny itself (S5: +0.16 top-1 on hard queries) remains the only extra layer that wins on these cases.
- No integration, no frozen config.

## 7. Literature credit (all checked via PubMed or Crossref; agent-verified 2026-10-27)

- Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language
  model. *Science* 2023;379(6637):1123–1130. doi:10.1126/science.ade2574 (ESM-2).
- Elnaggar A, Heinzinger M, Dallago C, et al. ProtTrans: toward understanding the language of life through
  self-supervised learning. *IEEE TPAMI* 2022;44(10):7112–7127. doi:10.1109/TPAMI.2021.3095381 (ProtT5).
- Kaminski K, Ludwiczak J, Pawlicki K, et al. pLM-BLAST: distant homology detection based on direct comparison of
  sequence representations from protein language models. *Bioinformatics* 2023;39(10):btad579.
  doi:10.1093/bioinformatics/btad579.
- Pantolini L, Studer G, Pereira J, et al. Embedding-based alignment: combining protein language models with
  dynamic programming alignment to detect structural similarities in the twilight-zone. *Bioinformatics*
  2024;40(1):btad786. doi:10.1093/bioinformatics/btad786 (EBA; our enhancement step follows its idea; our DP
  settings are our own).
- Heinzinger M, Littmann M, Sillitoe I, et al. Contrastive learning on protein embeddings enlightens midnight zone.
  *NAR Genom Bioinform* 2022;4(2):lqac043. doi:10.1093/nargab/lqac043 (contrastive heads on frozen pLM embeddings).
- Hong L, Hu Z, Sun S, et al. Fast, sensitive detection of protein homologs using deep dense retrieval. *Nat
  Biotechnol* 2025;43(6):983–995. doi:10.1038/s41587-024-02353-6 (embedding kNN homolog retrieval).
- Liu W, Wang Z, You R, et al. PLMSearch: protein language model powers accurate and fast sequence search for
  remote homology. *Nat Commun* 2024;15:2775. doi:10.1038/s41467-024-46808-5.
- Hamamsy T, Morton JT, Blackwell R, et al. Protein remote homology detection and structural alignment using deep
  learning. *Nat Biotechnol* 2024;42(6):975–985. doi:10.1038/s41587-023-01917-2.
- Johnson J, Douze M, Jégou H. Billion-scale similarity search with GPUs. *IEEE Trans Big Data*
  2021;7(3):535–547. doi:10.1109/TBDATA.2019.2921572 (FAISS).
- van den Oord A, Li Y, Vinyals O. Representation learning with contrastive predictive coding. arXiv:1807.03748
  (InfoNCE).
- Prior work on pLMs *for orthology* is thin. A PubMed search found no paper using pLM embeddings for orthology
  inference or ortholog/paralog discrimination. The closest are:
  - Endoh et al., *Proteomes* 2026 (doi:10.3390/proteomes14030036): ESM-2 divergence within RBH orthologs;
  - Shaw, Love & McWhite, *J Mol Evol* 2025;93(6):765–776 (doi:10.1007/s00239-025-10282-4): embeddings vs
    functional similarity and the ortholog conjecture.
  arXiv and bioRxiv were not searched.

## 8. Failed or not pursued (for S1 Text)

- Pretrained pLM similarity (ESM-2 35M/650M/3B, ProtT5; mean-pooled cosine or EBA-style embedding alignment) ranks
  orthologs above same-species paralogs worse than sequence on DEV. Deep: 0.61–0.72 vs 0.757. On hard cases
  it is at most equal to sequence. Its errors are correlated with sequence errors.
- As a near-tie breaker it is neutral: best +0.022 [−0.001, +0.046]. Step (a) failed, so no merge integration.
- Contrastive heads on frozen embeddings, trained on synteny-supported anchors (no benchmark labels), improve
  easy cases (+0.04) but not hard ones. Unambiguous-anchor labels contain no sequence near-ties by construction.
  Near-tie labels decided by synteny (post hoc) did not separate from their control.
- pLM kNN family rescue with a mutual-link + null-percentile guard chained 91% of genes into one family. Killed.
- Not pursued: ESM-C (package and licence), LoRA fine-tuning (the head result gave no reason to expect a gain on
  hard cases), 3Di tokens as extra head input (S6 already showed folds are shared by paralogs), and pLM distances
  for `fast` NJ trees (gated on step (a)).
