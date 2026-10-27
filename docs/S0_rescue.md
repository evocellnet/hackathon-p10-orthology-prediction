# S0 family rescue after the stronger merge (re-test of S6's 3Di rescue + MMseqs2 profiles)

Status: in progress. Author: method-engineer agent, started 2026-10-27 00:21 CET (T+9h40, `scripts/clock.sh`).
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
