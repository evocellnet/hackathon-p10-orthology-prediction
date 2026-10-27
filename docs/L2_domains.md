# L2-domains: protein domain architecture as an extra data layer

Agent `domains` (L2 round, docs/L2_protocol.md). Branch `l2-domains`, worktree `work/l2_domains_wt`.
Post-freeze exploration on DEV data only. Nothing here changes the frozen presets (`freeze-v1`).
Started 2026-10-27 18:00 CET (`scripts/clock.sh`).

## 0. Pre-registration (written 18:10 CET, committed before any domain-layer result was looked at)

At the time of writing, the only output that exists is a timing test of hmmsearch on 500 proteins
(throughput only, no annotation was inspected). No similarity, AUC or merge result has been computed.

### 0.1 Hypotheses
- **H1 (standalone).** Agreement of Pfam domain architectures separates orthologs from paralogs among homologous
  pairs, including the hard cases where sequence similarity is ambiguous, better than a shuffled-architecture
  control. Prior work: orthologs conserve domain architecture more than paralogs at equal sequence divergence
  (Forslund et al. 2011), and feature-architecture similarity is used as an orthology quality signal (FAS; QfO
  FAS benchmark).
- **H2 (recall at the merge, rescue).** 78% of missed TreeFam-A dev pairs are merge failures at the pair's LCA
  (RESEARCH_TREE 8b3). Relaxing the R1 mutual-best tolerance (0.95 → 0.8) recovered recall before, but cost
  precision (global tol 0.8: OrthoBench dev OG-F 71.2; deep_tol: PPV −0.033). If architecture agreement is
  informative, relaxing the tolerance **only for candidates with identical domain architecture** keeps part of
  the recall gain without the precision loss.
- **H3 (re-ranking).** An architecture bonus in R1 partner choice (same mechanics as the adopted synteny bonus:
  multiplicative, only among candidates within 0.75 of the best raw score) improves pair-level accuracy.
- **H4 (precision guard).** Vetoing R1 mutual-best seeds whose two units share no Pfam clan/family (both
  annotated) removes false orthology calls without recall loss.

### 0.2 Data and annotation (fixed)
- Pfam-A 38.2 HMMs (InterPro FTP), `hmmsearch --cut_ga` (HMMER 3.4), per-domain envelope coordinates from
  `--domtblout`. Pfam families are mapped to their clan (`Pfam-A.clans.tsv`) when they have one; otherwise the
  family accession is used.
- Overlap resolution per protein: domains sorted by independent E-value, a domain is kept unless it overlaps an
  already kept domain by > 50% of the shorter envelope. Architecture = kept domains ordered by envelope start,
  consecutive identical units collapsed (repeat counts are not compared). Unannotated = no kept domain.
- Similarities between two proteins: `J` = Jaccard of the clan sets (primary); `ident` = identical collapsed
  ordered architecture (secondary). Undefined (NaN) if either protein is unannotated.
- Proteins annotated (scope chosen without looking at any result): (i) all 58,130 proteins of the 37 S0 families
  (`work/core_v0_ob/run_det`) that hold a DEV RefOG gene; (ii) the 3,581 TreeFam-A DEV reference genes (QfO 2022);
  (iii) only if the standalone gate (0.5) passes: all 446,439 proteins of the 276 QfO S0 families
  (`work/qfo_recall/base`) that hold a TreeFam-A DEV gene (estimated ~26 CPU-h, needed for the QfO merge arms).

### 0.3 Step (a): standalone signal (dev only)
- **OrthoBench dev:** labelled pairs from `work/S5/refog_pair_labels.parquet` (19,146 within-RefOG cross-species
  pairs, 4,485 orthologs; plus 7,058 "deep" paralog pairs between RefOGs of one homology component).
  Sequence score = the pipeline's leaf score (`norm_max` of `run_det` `scores.npz`; 0 if unscored).
  **Hard queries** (S5 definition, with this score): a gene q of species A that has ≥1 ortholog and ≥1 paralog
  candidate in species B, and whose best paralog score ≥ 0.9 × its best ortholog score; hard pairs = all candidates
  of hard queries.
- **TreeFam-A dev:** reference pairs (S = speciation/ortholog, D = duplication/paralog) of the 359 dev families
  (`scripts/qfo_treefam_dev_fast.dev_refpairs`). Sequence score = `norm_max` of `work/qfo_recall/base`. Hard
  queries defined in the same way (per gene and other species).
- **Metrics:** ROC-AUC of `J` (ties count 0.5) for ortholog vs paralog, on all pairs with both genes annotated and
  on hard pairs; 95% CI by bootstrap over RefOGs / TreeFam families (B=1000, seed 11); coverage (fraction of pairs
  with both annotated). Same for `ident`. Hard-query top-1 partner accuracy: sequence top-1 vs re-ranking by
  `score × (1 + 0.3·J)` (the H3 rule at leaf level); Δ with bootstrap CI.
- **Controls:** architectures permuted among *annotated* proteins of the same species (seeds 11, 23, 47;
  annotation status is kept, so the control isolates the content of the architecture). Stricter secondary control:
  permutation within species *and* S0 family.
- **Gate to step (b):** on OrthoBench dev hard pairs **or** TreeFam-A dev hard pairs, the lower bound of the 95% CI
  of the real `J` AUC is above the largest of the three species-scope shuffled AUCs. Otherwise stop and record a
  negative result (no QfO family annotation, no merge arms).

### 0.4 Step (b): integration arms (flags off by default; frozen presets unchanged)
Implementation: `MergeConfig.arch_mode` (`off` default) and an `arch` object passed to `run_merge` (both passes of
the constrained two-pass synteny merge receive it). Unit architecture at a node = the most frequent architecture
among the unit's annotated member genes (ties → smaller code); unit `J` = Jaccard of the clan sets of the two
unit architectures; NaN if either unit has no annotated member (then the arm does nothing for that candidate).

| arm | rule | H |
|---|---|---|
| `off` | current preset (sequence + synteny), arch off | – |
| `bonus` | R1 partner choice on `score + 0.3·score·J`, only for candidates within 0.75 × best raw score of both endpoints (same window as synteny) | H3 |
| `gate` | R1 seeds = mutual best at tol 0.95 (as now) ∪ mutual best at tol 0.8 **with identical unit architecture** | H2 |
| `gate_deep` | as `gate`, only at deep nodes (median seed −ln S ≥ 1.5; QfO only, OrthoBench has no deep node) | H2 |
| `veto` | R1 seeds with unit `J` = 0 (both annotated) are removed | H4 |

Every real arm is run with shuffled architectures (seeds 11/23/47, species scope as in 0.3). One fixed setting per
arm (λ=0.3, tol 0.8, J=0); no grid. Bases: QfO = `work/qfo_recall/base` cached stages, `fast` preset + gene order
(`scripts/qfo_merge_iter.py` logic); OrthoBench = `work/core_v0_ob/run_det` cached stages, `accurate` preset + gene
order (OG-F 79.69 / macro pair F1 0.8064 reference).

- **Primary dev metric:** TreeFam-A dev pair F1 = 2TP/(2TP+FP+FN) over the dev reference pairs (raw counts),
  Δ(real − off) with a paired bootstrap over dev families (B=2000, seed 11).
- **Keep rule (protocol 6b):** (1) the Δ(real − off) CI lower bound > 0; (2) real − mean(shuffled) > the seed spread
  (max − min of the three shuffled values); (3) guards: TreeFam-A dev PPV not lower than off by > 0.01; OrthoBench dev
  OG-F not lower by > 1.0; OrthoXML of the QfO output passes `scripts/verify_orthoxml.py`.
- **Secondary (reported, not decisive):** TreeFam-A dev TPR/PPV; OrthoBench dev micro and macro pair F1 with the paired
  RefOG bootstrap; FN breakdown (S0 split / different root HOG / paralog call).
- If more than one arm passes, the one with the largest primary Δ is the candidate; combinations are post hoc.
  A passing arm is re-checked on the `accurate` preset on QfO (guards only) if time permits; otherwise this is
  stated as not done. The frozen config of a passing arm is written to `results/l2_domains/frozen_config.json`.

### 0.5 Not pre-registered / out of scope
Domain-level re-scoring of pairs on shared domains only (idea 4) and architecture-aware S0 splitting (idea 3, S0
side) are not run in this round (cost; listed in §8). TED/CATH domains are not used (Pfam covers the needed
proteins at lower cost; listed in §8).
