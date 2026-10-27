# S0 family rescue after the stronger merge (re-test of S6's 3Di rescue + MMseqs2 profiles)

Status: done (grid complete 03:42 CET, T+13h02). Author: method-engineer agent, started 2026-10-27 00:21 CET (T+9h40, `scripts/clock.sh`).
Code: `ancestra/rescue.py` (new), additive option in `ancestra/pipeline.py`
(`--family-rescue {none,3di,profile,both}`, `--structure-3di FASTA`, `--family-rescue-guard`; default `none`
= no new stage, identical stage keys). Tests: `tests/test_rescue.py`. Drivers: `scripts/rescue_ob_arms.sh`,
`scripts/rescue_qfo_arms.sh`, `scripts/rescue_eval_ob.py`, `scripts/rescue_eval_qfo.py`.
Results: `results/rescue/` (`ob_runs.tsv`, `qfo_runs.tsv`, per-arm JSON copies). Work: `work/rescue/` (hard-linked
copies of the cached deterministic S0/S2 stages; `work/qfo2022_det` and `work/core_v0_ob/run_det` are never written).

## 0. Pre-registration (written 00:40 CET, before any arm result was looked at)
**Base.** QfO: cached det S0 + S2 (`work/qfo_recall/base` = `det_go` stages), merge with the `fast` preset + gene
order (= A8check, TreeFam-A dev TPR 0.6163 / PPV 0.9198). OrthoBench: cached `run_det` S0 + S2, `accurate` preset +
gene order (OG-F 79.69). Everything after S2 is rerun through `python -m ancestra run` (so the pipeline option is
what is measured). R0 is rerun through the same path as a same-code baseline.

**Rescue (one fixed configuration, `RescueConfig` defaults; no tuning):** small family = ≤ 10 genes; query units =
one gene per S0 stage-1 cluster inside small families; target units = representatives of the 20 largest stage-1
clusters of every family > 10 genes; link = e ≤ 1e-3 and aligned span ≥ 0.3 of the shorter protein; **attach rule**:
each small family joins only the family of its single best link (lowest e), joins applied in increasing e, a join
that would exceed 7,500 genes is skipped; target families are never joined to each other and small families never
chain. Only new pairs are scored (DIAMOND `--sensitive -k0` + MMseqs2 SW of missed pairs in families ≤ 1,000 genes,
cap rule for families > 5,000), so S0/S2 caches are reused.

| arm | evidence | merge |
|---|---|---|
| R0 | none (current preset) | – |
| R1 | Foldseek 3Di+AA (`--alignment-type 2`, -s 7.5, max-seqs 100), AFDB exact-sequence 3Di (`data/structure/{qfo,orthobench}_3di_best.fa`) | – |
| R2 | MMseqs2 profiles: one profile per target representative, `result2profile` of its alignment (e ≤ 1e-3) to ≤ 100 family members (other reps first, then its own cluster); profile query vs small-family sequences, -s 7.5 | – |
| R3 | R1 ∪ R2 hits, same attach rule (lowest e wins) | – |
| R4 | R3 | **merge-side guard**: after a first merge, an attached small family is kept only if ≥ 1 of its genes is predicted orthologous to an original member of the target family *and* the pair is a reciprocal best hit (norm_max) between the two species; otherwise it gets its old family back (its cross-family pairs are dropped) and the merge is run again |

**Adoption rule (PI):** vs R0, TreeFam-A dev TPR up with PPV drop ≤ 0.01, AND OrthoBench dev OG-F ≥ −1, AND OrthoBench
dev macro pair F1 ≥ −0.005. If several arms pass: highest TreeFam-A dev TPR, ties → cheaper. Presets are not changed
by this agent; a passing arm is reported to the PI first.

**Pre-arm diagnosis (dev only, `results/rescue/qfo_dev_s0fn_pairs.json`).** Of the 1,286 S0-split TreeFam-A dev FN
pairs (S2 final families), only 198 involve a singleton (168) or a ≤ 10-gene family (30); 377 involve families of
11–100 genes and **711 have both genes in families > 100 genes** (mostly parts of the Leiden-split giant component,
e.g. 5,115 vs 6,554 genes). 1,148 genes are involved, almost all pairs are *S. pombe* (UP000002485) vs another
eukaryote. So an attach rule for small families can reach at most ~200 FN records (≤ +0.014 TPR) on QfO; the Leiden
cut is out of reach for every arm in this grid.

## 1. Implementation (`ancestra/rescue.py`, pipeline option)
- `stage1_rep_of` rebuilds S0's canonical order (length desc, sequence, id) to map `stage1_clusters.tsv`
  (`c<rank>` ids) back to genes; `select_units` picks query and target units (§0). QfO: 172,165 query units
  (167,021 with 3Di) vs 41,320 target units (37,787 with 3Di) in 5,259 target families; OrthoBench: 20,639 vs 12,264
  (17,857 / 9,710 with 3Di).
- `load_3di_any` accepts both header styles (OrthoBench `protein|uniprot`, QfO `proteome|protein`); 3Di search
  reuses `structure.build_db/search`.
- `search_profiles`: `mmseqs createdb --shuffle 0` → explicit centre→member prefilter (`tsv2db`) → `align -e 1e-3`
  → `result2profile` → `search prof qdb -s 7.5 -e 1e-3`.
- `attach_small` (rule §0), `rescore_new_pairs` (only new pairs; S2 hybrid logic), `run_rescue` (driver, optional
  content-keyed search-hit cache in `<run>/stages/rescue_hits/`), `write_families` (updated `families.tsv`),
  `guard_detach` + `revert_families` (R4).
- `pipeline.py` (additive): with `--family-rescue ≠ none` a stage `rescue` runs after `scores`
  (key = scores key + `RescueConfig.key()` + 3Di md5); it writes `scores.npz`, `families_final.npy`,
  `families_before.npy`, `attach.npz`, `families.tsv`. The merge then reads these; with `--family-rescue-guard`
  the merge stage runs merge → guard → merge. With the default `none` there is no new stage and every stage key is
  unchanged (`test_pipeline_family_rescue_stage`).
- Deterministic: no random step; DIAMOND `-b4 -c1`, MMseqs2/Foldseek ordered inputs, stable sorts.
- Tests (`tests/test_rescue.py`, 9 tests): attach rules (best link, no target–target or small–small joins, coverage,
  size cap, two-source union, determinism), unit selection, canonical stage-1 mapping, revert, guard on a hand-built
  merge result, a toy end-to-end profile rescue (remote homolog attached, random singleton not; cache reuse), and the
  pipeline stage/key behaviour. `tests/test_no_leak.py` still passes.

## 2. Results (dev only; `results/rescue/summary.{md,json}`, `qfo_runs.tsv`, `ob_runs.tsv`)
**Code change during the grid.** Another agent edited `ancestra/merge.py` at 02:08 (orthologGroup nesting fix,
`results/core_v0_qfo/runs.tsv` B0→B1). QfO R1–R4 (started 00:36–01:50) ran on the old merge; QfO R0 (02:20) on the new
one. The old-code baseline is therefore B0 = A8check (TPR 0.6163 / PPV 0.9198, same code as R1–R4), and the two
passing candidates were re-merged on the new code (R2n, R4n: same rescue stage, hard-linked; only the merge reruns,
`scripts/rescue_rerun_newcode.sh`). On OrthoBench the new code changes nothing: R0 = 79.69 / macro 0.8064 with
per-RefOG counts identical to `run_det` accurate; OB R1–R4 ran on the old code, R4n/R2n on the new one.

| arm | QfO code | TF-A dev TPR | PPV | S0-split FN | OB OG-F (P/R) | OB macro pair F1 | ΔTPR | ΔPPV | ΔOG-F | Δmacro | rule |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 | current merge.py | 0.6151 | 0.9224 | 1286 | 79.69 (75.7/84.1) | 0.8064 | – | – | – | – | – |
| R1 | old merge.py (pre 02:08) | 0.6135 | 0.9204 | 1276 | 68.66 (62.8/75.8) | 0.8099 | -0.0028 | +0.0006 | -11.03 | +0.0035 | fail |
| R2 | old merge.py (pre 02:08) | 0.6200 | 0.9165 | 1283 | 79.79 (75.7/84.4) | 0.8063 | +0.0037 | -0.0033 | +0.10 | -0.0001 | pass |
| R3 | old merge.py (pre 02:08) | 0.6126 | 0.9187 | 1276 | 68.78 (62.8/76.0) | 0.8122 | -0.0037 | -0.0011 | -10.91 | +0.0058 | fail |
| R4 | old merge.py (pre 02:08) | 0.6183 | 0.9214 | 1276 | 79.93 (75.8/84.5) | 0.8109 | +0.0020 | +0.0016 | +0.24 | +0.0045 | pass |
| R2n | current merge.py | 0.6163 | 0.9230 | 1283 | 79.79 (75.7/84.4) | 0.8063 | +0.0012 | +0.0006 | +0.10 | -0.0001 | pass |
| R4n | current merge.py | 0.6170 | 0.9248 | 1276 | 79.93 (75.8/84.5) | 0.8109 | +0.0019 | +0.0024 | +0.24 | +0.0045 | pass |

Δ vs R0 (OrthoBench, all rows) and vs the same-code QfO baseline (R1–R4: B0 = 0.6163 / 0.9198; R2n/R4n: R0). Rule: ΔTPR > 0, ΔPPV ≥ −0.01, ΔOG-F ≥ −1, Δmacro ≥ −0.005.


**Per-RefOG (OrthoBench dev, vs R0).**
- R1 (3Di): RefOG057 and RefOG028 become exact, RefOG013 +10 TP — but RefOG010 (3 genes) absorbs 102 foreign genes
  (FP 0 → 306) and RefOG021 splits (5 → 13 groups, TP 6,793 → 1,635). Same failure as S6 §2: RefOG010's root HOG is
  99 genes of its own S0 family plus 6 attached genes from two 3-gene families (Ciona/worm outgroups) that bridge the
  family's Bilateria-level groups. OG-F −11.0.
- R2 (profiles): RefOG066 exact (the 92-aa Ciona fragment, one of the six "unreachable by sequence" RefOGs of
  docs/S0_families.md), RefOG013 +10 TP; nothing else changes. OG-F +0.10, macro −0.0001.
- R3 = R1's damage + R2's gain (−10.9).
- **R4 (guard): exactly RefOG057 and RefOG066 become exact; every other RefOG is unchanged** (RefOG010/021 damage
  gone; RefOG013/028 gains also gone). The guard keeps 500 of 3,192 attachments on OrthoBench and 6,338 of 28,148 on
  QfO. OG-F +0.24, macro pair F1 +0.0045, micro pair F1 +0.0036 [−0.0002, +0.0097] (paired bootstrap over RefOGs).
  Family-level S0 completeness (deterministic partition): 26/34 → 28/34 with 3Di (013, 057) or profiles (013, 066),
  29/34 = 0.853 with both (fusion proxy unchanged, 0.353); the dev ceiling of sequence-only S0 was ~0.82.

**QfO TreeFam-A dev.** The S0-split FN class barely moves (1,286 → 1,276–1,283): as diagnosed in §0, only ~200 of
these records involve a small family at all. The TPR changes come mostly from pairs *inside receiving families*
(merge perturbation), not from the attached genes: R4n vs R0 changes 43 TP records (10 touch an attached gene, 33 sit
in receiving families) and 16 FN / 6 FP / 28 TN records (`results/rescue/qfo_diff_R4n_vs_R0.json`); R2n vs R0:
64 TP (3 attached). Unguarded 3Di (R1, R3) lowers TPR (paralog calls 1,776 → 2,104–2,134: attached remote homologs
are mostly paralogs, cf. S6 §3.4). Every QfO Δ is far inside the TPR 95% CI half-width (±0.011).

## 3. Cost
Shared node (load 40–60), 5 threads per arm; QfO `fast`, OrthoBench `accurate`.

| step | QfO (988k proteins) | OrthoBench (251k) |
|---|---|---|
| unit selection + canonical stage-1 map | ~20 s | 1 s |
| 3Di search (172k × 41k units, -s 7.5) | 533 s wall, ≈0.45 CPU-h | 17 s wall (10 thr) |
| profile build + search (41k profiles) | 1,233 s wall, ≈0.95 CPU-h | 96 s wall (10 thr) |
| re-scoring new pairs only (DIAMOND + MMseqs2) | 173–434 s (11k–48k attached genes, 0.3–0.4 M new pairs) | 16–27 s |
| rescue stage total, R3/R4 (both searches) | ≈ 30 min wall, ≈1.7 CPU-h | ≈2 min |
| R4 guard: guard check + second merge | +7 min (fast; QfO `accurate` would add a full FastTree pass ≈ +1 h wall) | +~20 min wall, +1.0 CPU-h (accurate) |
Full S2 re-scoring was **not** needed: only new pairs are aligned (vs 8.9 CPU-h for S2 on QfO).

## 4. Verdict and recommendation (for the PI; presets NOT changed)
- **R4 (3Di ∪ profiles + merge-side RBH guard) passes the pre-registered rule on the current code**: TreeFam-A dev TPR
  +0.0019, PPV +0.0024; OrthoBench dev OG-F +0.24, macro pair F1 +0.0045. It also passes on the old code
  (+0.0020 / +0.0016). R2 (profiles only) passes on the old code (+0.0037 / −0.0033 / +0.10 / −0.0001); on the
  current code it still passes (QfO +0.0012 / +0.0006; OB R2n = R2). By the pre-fixed tie rule (highest TF-A dev TPR
  on the current code) R4 is the candidate. R1 and R3 (unguarded 3Di) fail (OG-F −11).
- **The effect size is small and mostly not the targeted mechanism on QfO.** Every QfO Δ is inside ±0.011 (TPR CI
  half-width); only 10 S0-split records are repaired; most changed pairs are merge perturbations in receiving
  families. On OrthoBench the gain is real and mechanistic (two RefOGs, 057 and 066, become exact; one through 3Di,
  one through profiles) but it is 2 of 34 RefOGs.
- **The guard is what makes structure usable**: it removes the S6 failure mode (RefOG010 fusion) on both benchmarks
  by undoing ~80% of the attachments, keeping only those that the merge confirms with a reciprocal best ortholog.
- **Cost**: QfO rescue stage ≈ 30 min wall / ≈1.7 CPU-h at 5 threads (+ the guard's second merge: +7 min for
  `fast`, ≈ +1 h for `accurate` on QfO); OrthoBench ≈ 2 min + ≈ 20 min (second accurate merge).
- Recommendation: R4 meets the letter of the rule, but its QfO benefit is within noise while the guard doubles the
  merge cost of `accurate`. Options for the PI: (a) adopt R4 into `accurate` only (OrthoBench-type gains, 2 RefOGs);
  (b) keep it as a documented optional mode (`--family-rescue both --structure-3di … --family-rescue-guard`) and
  report it as a positive sub-task result with a null QfO effect. I lean to (b) unless the frozen-config QfO run can
  absorb +1 h. The Leiden cut of the giant S0 component (711 of 1,286 S0-split FN records, both genes in > 100-gene
  families, mostly *S. pombe*) is the larger lever and needs a different fix (e.g. gene-level reassignment across
  Leiden parts by best-hit family), not small-family attachment.

## 5. Caveats
- Pre-arm diagnosis is dev only; no test family or test RefOG was read. TreeFam-A was scored only with
  `qfo_treefam_dev_fast.py` (dev), never `--allow-test`.
- No shuffled-3Di control was run in this grid (S6 §2 already showed 0 rescued RefOGs under shuffled 3Di); the guard
  arm was pre-registered before results but its RBH criterion is one choice among several.
- `merge.py` changed mid-grid (§2); old/new baselines are kept apart in the table.
- QfO R2n/R4n reused the old runs' rescue stages (rescue does not import merge code; stage keys identical).
