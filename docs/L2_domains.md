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

## 1. Data

| source | URL | version | retrieved | licence | md5 |
|---|---|---|---|---|---|
| Pfam-A HMMs | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz | Pfam 38.2 (06/2026; 30,134 families; built on UniProtKB reference proteomes 2025_03) | 2026-10-27 | CC0 1.0 | 7ab3c4e215d0daaea3004e37c4e24f8a (gz; matches the release `md5_checksums`); uncompressed 589ccc8370eaf079022db0dbb25c8aa2 |
| Pfam-A clans | .../current_release/Pfam-A.clans.tsv.gz | 38.2 | 2026-10-27 | CC0 1.0 | a70eeffcd8105ff094ea5834456d8bab |
| Pfam-A.hmm.dat | .../current_release/Pfam-A.hmm.dat.gz | 38.2 | 2026-10-27 | CC0 1.0 | aab364689372b1846e8e880976e7f750 |

Local copies in `data/extra/l2_domains/` (gitignored). Annotation: `scripts/L2d_hmmsearch.sh` (hmmsearch 3.4,
`--cut_ga --noali --domtblout`, single-thread processes over sequence chunks).
- Stage 1 (OrthoBench dev families 58,130 + TreeFam-A dev genes 3,581 = 61,711 proteins, 41 M residues): 30 min wall,
  4.9 CPU-h (8 processes). 59,655 proteins (96.7%) have ≥ 1 GA hit; 56,108 of the 58,130 OrthoBench family
  proteins (96.5%) and 3,547 of 3,581 TreeFam-A dev genes (99.1%) have an architecture after overlap resolution.
- Stage 2 (QfO S0 families holding a TreeFam-A dev gene, 442,858 further proteins, 228 M residues): run after the
  standalone gate passed (§2), 6 processes; cost in §6.
- Not used: the SonicParanoid2-bundled MMseqs2 Pfam 35 profile DB (a 500-protein test gave implausible hit lists;
  not debugged), the InterPro REST API (no per-proteome location download; QfO 2022 sequences may differ from current
  UniProt), TED/CATH (not needed once Pfam covered > 96% of the proteins).

## 2. Step (a): standalone signal (`results/l2_domains/standalone.json`, `scripts/L2d_standalone.py`)

J = Jaccard of clan sets; ident = identical collapsed ordered architecture. AUC CI = bootstrap over RefOGs
(OrthoBench) or TreeFam families (B=1000, seed 11). Shuffled = architectures permuted among annotated genes of the
same species (seeds 11/23/47); "fam-shuffled" = within species × S0 family (a weak null: species × family groups
with one annotated gene are not permuted, so part of the real signal survives it).

| set | pairs (orth) | coverage | AUC J [95% CI] | AUC ident | shuffled J (11/23/47) | fam-shuffled J | AUC seq (same pairs) |
|---|---|---|---|---|---|---|---|
| OB dev, all labelled pairs (within-RefOG + deep paralogs) | 25,712 (4,393) | 0.981 | 0.688 [0.599, 0.814] | 0.638 | 0.497 / 0.496 / 0.498 | 0.614 / 0.626 / 0.626 | 0.844 |
| OB dev, all, **hard** | 6,990 (917) | | **0.637 [0.547, 0.738]** | 0.592 | 0.501 / 0.496 / 0.490 | 0.556 / 0.572 / 0.581 | 0.743 |
| OB dev, within-RefOG only | 18,807 (4,393) | 0.982 | 0.542 [0.465, 0.704] | 0.530 | 0.497 / 0.496 / 0.498 | 0.568 / 0.579 / 0.584 | 0.801 |
| OB dev, within-RefOG only, hard | 5,175 (887) | | 0.488 [0.446, 0.582] | 0.490 | 0.500 / 0.496 / 0.488 | 0.499 / 0.514 / 0.528 | 0.681 |
| TreeFam-A dev, all reference pairs | 25,023 (13,729) | 0.985 | 0.493 [0.447, 0.552] | 0.497 | 0.501 / 0.500 / 0.498 | 0.449 / 0.442 / 0.447 | 0.616 |
| TreeFam-A dev, hard | 4,658 (1,455) | | 0.498 [0.434, 0.572] | 0.499 | 0.505 / 0.497 / 0.495 | 0.480 / 0.476 / 0.474 | 0.444 |

Mean J, orthologs vs paralogs: OB within-RefOG 0.900 vs 0.846 (hard 0.854 vs 0.858); OB deep paralogs pull the
"all" paralog mean down to 0.600. TreeFam-A dev 0.876 vs 0.904 (paralogs slightly *more* similar).

**Hard-query partner choice** (sequence top-1 vs `score × (1 + 0.3·J)` within 0.75 × best; Δ with bootstrap CI):

| set | hard queries | seq top-1 | real Δ [95% CI] | shuffled Δ (11/23/47) |
|---|---|---|---|---|
| OB dev, all | 1,012 | 0.415 | +0.001 [−0.008, +0.011] | −0.014 / −0.006 / −0.008 |
| OB dev, within-RefOG | 971 | 0.390 | +0.001 [−0.009, +0.012] | −0.014 / −0.006 / −0.008 |
| TreeFam-A dev | 1,141 | 0.202 | +0.004 [−0.003, +0.016] | +0.004 / −0.004 / +0.001 |

**Family-level signal (idea 3).** Same RefOG vs different RefOG of one homology component (OB dev, all labelled pairs):
AUC J **0.979 [0.839, 0.995]** vs sequence 0.708; shuffled 0.501 / 0.498 / 0.499.

**Reading.**
1. The pre-registered gate (§0.3) passes on OrthoBench dev hard pairs (CI lower bound 0.547 > max shuffled 0.501),
   so step (b) was run as pre-registered.
2. But the signal is almost entirely **between RefOGs**: architecture separates distant paralogs from different
   orthogroups (AUC 0.98), not orthologs from paralogs *inside* an orthogroup (within-RefOG AUC 0.54, hard 0.49;
   TreeFam-A dev 0.49). Inside an orthogroup, orthologs and in-paralogs share the same architecture ~70–85% of the
   time. This agrees with Forslund et al. (2011): orthologs and their closest paralogs have very similar architectures.
3. Architecture adds nothing to sequence for partner choice (Δ top-1 +0.001, CI includes 0). Shuffled architectures
   hurt (−0.006 to −0.014), so the real layer is neutral rather than harmful.
4. Expectation for step (b), stated before any merge arm was run: H3 (re-ranking) should be null; H2/H4 can only act
   through deep/between-orthogroup decisions; the primary metric (TreeFam-A dev pairs, all *within* TreeFam families)
   is mostly blind to between-family fusions, so a null primary result is the most likely outcome.

## 3. Methods (step b)

- **Code** (branch `l2-domains`): `ancestra/domains.py` (domtbl parsing with overlap resolution and clan mapping;
  `DomainArch`: gene → architecture code, shuffled control, per-node unit architecture = most frequent architecture of
  the unit's annotated members, vectorised over the merge's member CSR); `ancestra/merge.py`: `MergeConfig.arch_mode`
  (`off` default, omitted from `key()` at its default so stage hashes are unchanged), `arch_lam`, `arch_tol`, and an
  `arch=` argument of `run_merge` → `decide_node` (both passes of `run_merge_constrained` receive it through `**kw`).
  `bonus` adds `arch_lam·score·J` to R1's selection score inside the synteny window; `gate` adds mutual-best seeds at
  `arch_tol` whose two units have the same architecture; `gate_deep` does so only at nodes flagged deep; `veto` drops
  mutual-best seeds with unit J = 0. With `arch_mode="off"` (or `arch=None`) the code path is unchanged:
  `tests/test_domains.py::test_run_merge_arch_off_identical`, and the OrthoBench/QfO `off` arms reproduce the frozen
  dev numbers exactly (OrthoBench OG-F 79.69 = rescue R0; QfO TP/FP/FN/TN 8,607/722/5,386/10,698 = runs.tsv B3).
  `pytest -q`: 157 passed, 5 skipped.
- **Drivers:** `scripts/L2d_merge_arm.py` (one arm on cached stages; constrained two-pass synteny merge + veto + R3
  choice, exactly as `pipeline.s_merge`; QfO `fast`, OrthoBench `accurate` with 2 FastTree workers),
  `scripts/L2d_queue.sh`, `scripts/L2d_analyze.py` (paired bootstraps, keep rule). Each arm JSON in
  `results/l2_domains/arms/` carries its manifest (git commit, config key, inputs, wall/CPU time).
- **Annotated scope in the merge:** QfO: all proteins of the 276 S0 families that hold a TreeFam-A dev gene (446,439;
  Pfam architecture for 411,523); OrthoBench: all proteins of the 37 dev-holding S0 families (56,108 with architecture).
  Proteins outside the scope have no architecture; the arms then do nothing for their units, and the merge is
  family-local, so dev metrics depend only on the annotated families.

## 4. Results, step (b) (`results/l2_domains/integration.json`, `integration_table.md`, `arms/*.json`)

Primary: TreeFam-A dev pair F1 (raw pair records), Δ(real − off) with a paired bootstrap over the 359 dev families
(B=2000, seed 11). Shuffled = architectures permuted within species among annotated genes, seeds 11/23/47.
QfO base `fast` (off: TPR 0.6151 / PPV 0.9226 / F1 0.7381); OrthoBench base `accurate` (off: OG-F 79.69, micro pair
F1 0.8035, macro 0.8026).

| arm | ΔF1 [95% CI] (primary) | F1 shuffled 11/23/47 | real − shuf (spread) | ΔTPR [CI] | ΔPPV [CI] | OB OG-F (shuffled) | OB pair F1 Δ [CI] | events (real) | keep? |
|---|---|---|---|---|---|---|---|---|---|
| bonus | −0.0017 [−0.0188, +0.0134] | 0.7376 / 0.7291 / 0.7310 | +0.0038 (0.0085) | −0.0032 [−0.024, +0.016] | +0.0020 [−0.004, +0.009] | **76.67** (79.35 / 79.45 / 81.43) | −0.0085 [−0.020, +0.001] | 602,249 bonused edges | no (CI, shuffle, OG-F) |
| gate | **+0.0182 [+0.0032, +0.0336]** | 0.7421 / 0.7321 / 0.7366 | +0.0194 (0.0100) | +0.0383 [+0.014, +0.064] | **−0.0250 [−0.046, −0.006]** | **76.82** (80.12 / 78.75 / 78.75) | −0.0139 [−0.027, +0.000] | 116,609 seeds added | **no** (PPV guard, OG-F guard) |
| gate_deep | +0.0103 [−0.0147, +0.0329] | 0.7410 / 0.7294 / 0.7376 | +0.0124 (0.0116) | +0.0302 [+0.013, +0.050] | −0.0319 [−0.085, +0.002] | – (no deep node) | – | 30,297 seeds added | no (CI, PPV guard) |
| veto | +0.0013 [+0.0000, +0.0044] | 0.6791 / 0.6822 / 0.6766 | +0.0601 (0.0057) | +0.0017 [0, +0.006] | +0.0002 [0, +0.001] | 79.69 (71.59 / 72.42 / 68.79) | +0.0000 | 1,923 seeds vetoed (OB 61) | no (CI lower bound = 0; one family) |

Official-formula TreeFam-A dev values (pseudocount, /2; comparable to `results/core_v0_qfo/runs.tsv`): off .6151/.9224,
bonus .6118/.9244, gate .6534/.8974, gate_deep .6452/.8905, veto .6168/.9226.

**Post hoc (labelled; not pre-registered, run after the real arms to interpret the gate):** the same relaxation
*without* the architecture condition (`results/l2_domains/posthoc/`).

| arm | TF-A dev TPR / PPV / F1 | OB OG-F / pair F1 |
|---|---|---|
| off | .6151 / .9226 / .7381 | 79.69 / .8035 |
| ungated tol 0.8 (all nodes) | .6477 / .8869 / .7486 | 71.25 / .7959 |
| **gate** (tol 0.8 only if same architecture) | .6534 / .8976 / .7563 | 76.82 / .7896 |
| ungated deep_tol 0.8 | .6499 / .8754 / .7460 | – |
| gate_deep | .6453 / .8907 / .7484 | – |

**Findings.**
1. **No arm passes the keep rule.** The only arm with a primary CI above 0 is `gate` (+0.018 F1, TPR +0.038), and it
   breaks two guards: TreeFam-A dev PPV −0.025 (rule: ≥ −0.01) and OrthoBench dev OG-F −2.87 (rule: ≥ −1.0). On the QfO
   Pareto plot it is dominated by an existing sequence-only arm (C2, accurate + so_struct + deep absorption: .6548/.9136 official formula,
   `results/core_v0_qfo/runs.tsv`), although that arm is on the accurate base.
2. **Architecture carries some information at the merge, but not enough.** Gated relaxation dominates ungated
   relaxation on QfO (TPR +0.006 *and* PPV +0.011) and costs far less OG-F on OrthoBench (76.8 vs 71.2). So requiring
   the same architecture removes part of the false seeds that a looser tolerance lets in. It does not remove enough of
   them: within an orthogroup, paralogs share the architecture as often as orthologs (§2), so the gate admits paralog
   seeds too.
3. **Re-ranking (`bonus`, H3) is null on QfO and harmful on OrthoBench** (OG-F −3.0, pair F1 −0.0085 with CI touching 0),
   as predicted by the standalone null within orthogroups. Shuffled bonuses are mostly harmless (OG-F 79.4–81.4), so the
   harm comes from real architectures pulling partner choice towards same-architecture paralogs.
4. **Veto (H4) almost never fires with real architectures** (1,923 of ~596k QfO seeds, 61 of ~169k OrthoBench seeds; one
   TreeFam family changes, +24 TP). Its shuffled control is not a useful null: random architectures produce J = 0 for
   most seeds and veto 333k of them (TPR −0.08), so "real beats shuffled" is trivially met. Pre-registered as is; noted
   as a design weakness of the control for veto-type rules.
5. The OrthoBench and QfO `off` arms reproduce the frozen dev numbers exactly, and every arm leaves the frozen presets
   untouched (`arch_mode="off"` default).

## 5. Verdict

**Negative result for the merge; keep as a diagnostic.** Pfam domain architecture separates orthogroups from each
other (same vs different RefOG: AUC 0.98 vs sequence 0.71), but it does not separate orthologs from in-paralogs inside
an orthogroup (OrthoBench within-RefOG AUC 0.54, hard 0.49; TreeFam-A dev 0.49), which is the decision the level-
synchronous merge actually makes. None of the four pre-registered integrations (R1 bonus, architecture-gated relaxed
tolerance at all or deep nodes, R1 veto) passes the keep rule; the gate buys recall (+0.038 TPR, CI > 0) at a precision
cost that breaks both guards. No frozen config is written (no arm passed), so there is nothing for a test agent to run.

Possible later use (not tested here; would need a new pre-registration): the strong between-orthogroup signal fits the
S0 side (family fusion guard, e.g. the 3Di rescue that fused RefOG010, or Leiden splitting of giant components), not the
merge.

## 6. Cost

| item | CPU-h | wall | notes |
|---|---|---|---|
| Pfam stage 1 (61,711 proteins, 41 M aa) | 4.9 | 30 min | 8 × hmmsearch --cpu 1 |
| Pfam stage 2 (442,858 proteins, 228 M aa) | 26.4 | 3.2 h | 6 processes |
| merge arms (30 pre-registered + 3 post hoc + 1 duplicate killed) | ≈ 16.5 | – | QfO fast 0.11 CPU-h/run; OrthoBench accurate ≈ 1.0 CPU-h/run at 2 FastTree workers |
| standalone benchmark | 0.02 | 55 s | |
| total | ≈ 48 | ~6 h compute wall (overlapping) | no GPU |

Annotation dominates: ≈ 0.22 CPU-s per protein, i.e. ≈ 60 CPU-h for all 988,778 QfO 2022 proteins and ≈ 15 CPU-h for the
251,378 OrthoBench proteins, which is ~4× the whole accurate OrthoBench pipeline (3.9 CPU-h, RESEARCH_TREE 8.3). The merge
overhead of the layer itself is small (QfO merge 334 s off vs 336–381 s with architectures, ≤ 14%).

## 7. Literature credit (each checked on PubMed; DOIs)

- Forslund K, Pekkari I, Sonnhammer ELL. Domain architecture conservation in orthologs. *BMC Bioinformatics* 2011;12:326.
  doi:10.1186/1471-2105-12-326 (orthologs conserve architecture more than paralogs at equal divergence, but the closest
  paralogs have very similar architectures, which is the limitation we hit).
- Dosch J, Bergmann H, Tran V, Ebersberger I. FAS: assessing the similarity between proteins using multi-layered feature
  architectures. *Bioinformatics* 2023;39(5):btad226. doi:10.1093/bioinformatics/btad226 (feature-architecture similarity;
  basis of the QfO FAS benchmark).
- Persson E, Kaduk M, Forslund SK, Sonnhammer ELL. Domainoid: domain-oriented orthology inference. *BMC Bioinformatics*
  2019;20:523. doi:10.1186/s12859-019-3137-2 (domain-level orthology; idea 4, not pursued here).
- Bitard-Feildel T, Kemena C, Greenwood JM, Bornberg-Bauer E. Domain similarity based orthology detection. *BMC
  Bioinformatics* 2015;16:154. doi:10.1186/s12859-015-0570-8 (porthoDom: domain-content similarity to group proteins
  before ortholog search, i.e. the family-level use that our between-orthogroup signal supports).
- Cosentino S, Sriswasdi S, Iwasaki W. SonicParanoid2: fast, accurate, and comprehensive orthology inference with machine
  learning and language models. *Genome Biol* 2024;25:195. doi:10.1186/s13059-024-03298-4 (domain-architecture-aware
  orthology at scale).
- Paysan-Lafosse T, et al. The Pfam protein families database: embracing AI/ML. *Nucleic Acids Res* 2025;53(D1):D523–D534.
  doi:10.1093/nar/gkae997 (Pfam; release 38.2 used).
- Eddy SR. Accelerated profile HMM searches. *PLoS Comput Biol* 2011;7(10):e1002195. doi:10.1371/journal.pcbi.1002195
  (HMMER3, hmmsearch).

## 8. Failed or not pursued (for S1 Text)

- Pfam domain architecture (Pfam 38.2, hmmsearch GA): separates orthogroups (same vs different RefOG AUC 0.98) but not
  orthologs from in-paralogs (OrthoBench dev within-RefOG AUC 0.54, hard 0.49; TreeFam-A dev 0.49; shuffled 0.50).
- Four merge integrations failed the pre-registered keep rule on dev: R1 architecture bonus (QfO ΔF1 −0.002; OrthoBench
  OG-F −3.0), architecture-gated relaxed tolerance (ΔF1 +0.018 [+0.003, +0.034], TPR +0.038, but PPV −0.025 and
  OrthoBench OG-F −2.9), the same at deep nodes only (PPV −0.032), and an R1 veto of seeds sharing no domain (fires on
  0.3% of seeds; one family changes).
- Post hoc: the architecture gate dominates the same relaxation without it (TPR +0.006, PPV +0.011, OrthoBench OG-F
  +5.6), so architecture is informative but too weak to license extra merges.
- Not pursued: domain-level re-scoring of pairs (Domainoid-style), architecture-aware S0 splitting/fusion guard, TED/CATH
  domains, InterPro precomputed matches. Annotation alone (~60 CPU-h for QfO) costs more than the whole pipeline.

## 9. Deviations and incidents

- No change to the pre-registration. Post hoc arms (§4) are labelled and were run after the real arms.
- Scheduling: one duplicate OrthoBench run (`veto`, seed 47) was started by two queues at once; the duplicate was killed
  before writing, and the kept output is from the first run (merges are deterministic).
- A first test run of `pytest` failed only because `diamond` was not on PATH; rerun with the conda env on PATH: 157
  passed, 5 skipped.
