# L2-network: interaction networks, co-fractionation and protein abundance as extra data layers

**Verdict: NEGATIVE, stopped at step (a).** Pre-registered gate failed for the primary network layer on
OrthoBench dev and (extension) on TreeFam-A dev, and for abundance. No integration arm was run, so there is
no frozen config. Summary:
- **Real signal, wrong place.** Interolog support separates orthologs from paralogs overall: OrthoBench dev
  AUC 0.630 [0.532, 0.699], TreeFam-A dev 0.622 [0.559, 0.663], shuffled 0.49–0.50. On the sequence hard
  cases it decides no better than a coin: precision when it decides is 0.44 (OrthoBench, n = 9 decisions)
  and 0.52 (TreeFam-A, n = 52). comb@0.9 − seq is −0.18 [−0.24, 0] and −0.025 [−0.075, 0].
- **Coverage is the binding limit on OrthoBench.** Only 3.4% of dev candidate pairs and 1.0% of dev hard
  pairs are evaluable. Teleosts, chicken, dog, opossum and Ciona have almost no own-species interaction data.
  Oracle bound: fixing every evaluable error of the frozen run gives at most +0.0016 micro pair F1.
- **TreeFam-A has better coverage but no hard-case signal.** 26% of candidates and 19% of hard pairs are
  evaluable. Used as a rescue of unpredicted pairs, its precision is only 0.48–0.54 at every support
  threshold, against the 0.94 PPV the merge works at.
- **Abundance (PaxDb) fails like expression did.** AUC 0.676 [0.604, 0.732], but on hard queries it is
  right only 34% of the time when it decides; comb@0.9 − seq is −0.131 [−0.246, −0.046].

Agent `network` (4 threads, no GPU), branch `l2-network`, worktree `work/l2_network_wt`. Post-freeze
exploration under docs/L2_protocol.md: dev only, nothing here changes `freeze-v1` or any reported number.

## 1. Pre-registration (committed 2026-10-27 ~18:15 CET, before any score or dev result was computed)

What was known when this was written: the protocol, the S5/S9 docs, and the per-species *data
coverage* of STRING and PaxDb after ID mapping (§2; no labels used). No interolog or abundance score
had been computed.

### 1.1 Hypotheses
- **H1 (network neighbourhood conservation, interologs).** An ortholog pair (a, b) shares more
  interaction partners, connected through anchor orthologs, than a paralog pair does. This holds above
  all for the sequence hard cases, where a paralog scores ≥ 0.9 × the ortholog. Prior work:
  interologs (Walhout et al. 2000; Matthews et al. 2001), network-based functional orthologs
  (Bandyopadhyay et al. 2006), IsoRank (Singh et al. 2008).
- **H2 (abundance conservation).** Orthologs have more similar protein abundance (PaxDb) than
  paralogs. Expected outcome: it fails like expression (S9), unless it does not.

### 1.2 Data and circularity rules
- STRING v12.0, per organism, OrthoBench's 12 species. Only the organism's **own** channels count:
  `experiments` and `coexpression`. The following are never used: every `*_transferred` column,
  `textmining`, `homology`, `neighborhood`, `fusion` and `cooccurence`. They are either transferred by
  orthology or are genomic-context or phylogenetic-profile evidence.
- `database` (curated pathways and complexes) is direct in STRING's sense. But KEGG and Reactome
  membership in non-model species is itself assigned by orthology. It is therefore run **only as a
  flagged circular reference** (`db400_CIRCULAR`). It is never eligible for integration and is never
  reported as evidence.
- PaxDb v5.0: the integrated whole-organism dataset per species, or the only dataset for chimp and dog.
- IDs are mapped by Ensembl protein ID, then by sequence md5, then by gene. Coverage is reported (§2).

### 1.3 Score definitions (fixed)
- **Anchors** are data-derived and never come from labels.
  - Standalone test: the 1:1 ortholog pairs of the frozen Ancestra accurate OrthoBench run
    (`work/ob_final/main_accurate_go/pairs.tsv`). Within the species pair, each gene must have exactly
    one predicted ortholog.
  - Integration: the merge's own strict 1:1 mutual-best units at each node (`anchor="core"`), the same
    anchors synteny uses.
  - Anchors with either gene in the S0 family of a or b are ignored. This also excludes the pair
    being scored.
- **Interolog count** `count(a,b)` = the number of anchors (x, y) with x a neighbour of a in A's network
  and y a neighbour of b in B's network. This is the **primary** score, and the same kernel as
  synteny's `context_support`, with the network in place of the gene-order window.
- **Secondary score:** `cos` = count / sqrt(n_a · n_b), where n_a and n_b are the anchored,
  family-filtered neighbours of a and b.
- **Evaluable** pair: n_a > 0 and n_b > 0.
- **Channels:**
  - **primary `own`** = experiments ≥ 150 OR coexpression ≥ 400.
  - Secondary: `exp150`, `exp400`, `coexp400`.
  - Circular reference: `db400_CIRCULAR`.
  - A species' network is used only if it has ≥ 100 nodes.
- **Abundance score:** −|q_a − q_b|, where q is the within-species percentile rank of PaxDb ppm.
- **Shuffled controls**, seeds 11/23/47:
  - Network: node labels are permuted within the species among the channel's nodes. This keeps the
    degree sequence and the coverage.
  - Abundance: values are permuted among each species' covered genes.

### 1.4 Step (a): standalone gate (OrthoBench dev, 34 RefOGs, all 66 species pairs)
- **Truth, candidates and hard queries** come from S5/S9.
  - Truth: `work/S5/refog_pair_labels.parquet`.
  - Candidates: within-RefOG and deep paralog pairs.
  - Sequence score: S7 canonical DIAMOND bits.
  - Hard: the best paralog scores ≥ 0.9 × the best ortholog. This is recomputed on the evaluable
    candidates, as in S9.
  - Decision rule: the S5 `query_decisions` (comb@r).
- **Reported:**
  - Coverage: the fraction of dev candidate pairs, ortholog pairs and globally-hard pairs that are
    evaluable.
  - Pooled ROC-AUC on evaluable pairs, with a RefOG-bootstrap 95% CI and the 3 shuffled AUCs.
  - On evaluable hard queries: the top-1 decisions (correct, wrong, abstain), the precision when the
    layer decides, and comb@0.9 / comb@0.75 − seq with CI. Each is compared with the shuffled seeds,
    including the paired real − shuffled difference.
- **Gate (a) passes for the primary channel/score only if all hold:**
  1. The pooled AUC CI lower bound is > the maximum shuffled AUC.
  2. On evaluable hard queries, precision when deciding is > 0.5 and > every shuffled seed's precision.
  3. comb@0.9 − seq ≥ 0 (point estimate), and real − shuffled comb@0.9 > 0 for all 3 seeds.
- Same gate for abundance (H2). If abundance passes, it is reported as a candidate only. Its
  integration is not run this round.
- **Coverage caveat (pre-stated).** Hard cases are dominated by vertebrate–teleost ohnologs, and
  zebrafish, Tetraodon, chicken, dog, opossum and Ciona have almost no own-species interaction data.
  If the evaluable hard pairs are < 10% of all dev hard pairs, any integration gain is bounded
  accordingly. That will be stated, but it does not by itself stop step (b).

### 1.5 Step (b): integration (only if gate (a) passes for `own`)
- **Base (`off`):** the frozen `accurate` preset with gene order, re-merged from the cached stages of
  `work/ob_final/main_accurate_go`: constrained synteny, veto and R3. `off` must reproduce the frozen
  dev pair counts exactly (flag default off → byte-identical).
- **Arms** (network adjacency at each node = the union of the child species' networks mapped to units,
  as synteny builds its unit adjacency):
  - **NET_sum (primary):** support = synteny support + network interolog support. It enters the
    unchanged constrained R1 bonus (multiplicative, λ = 0.3, cap 3), the R2/R3 veto and R3 partner
    choice.
  - **NET_add (secondary):** synteny unchanged, plus a separate multiplicative network bonus
    (λ_net = 0.3, cap 3) in R1 only. Veto support stays synteny-only.
  - Each arm is run with the network shuffled at seeds 11/23/47.
- **Primary dev metric:** OrthoBench dev **macro pair F1** (`scripts/orthobench_dev_pairs_eval.py`).
- **Keep rule (protocol §6b)** for NET_sum:
  - Δ(real − off) has a paired RefOG bootstrap 95% CI (B = 2000, seed 11) with lower bound > 0.
  - real − shuffled > the seed spread.
  - Guards: micro pair PPV no worse than off by > 0.01; dev OG-F (official benchmark) no worse by
    > 1.0; `verify_orthoxml` valid.
- **TreeFam-A dev:** only if the OrthoBench arm passes and time allows. This would need STRING mapping
  for the QfO proteomes and about 1 h of accurate merge per arm. Otherwise it is listed as not run.
- No parameter search: λ, cap, thresholds and channels are fixed above. Any later change is logged
  in §8 as a deviation.

### 1.6 CF-MS
hu.MAP 3.0 (human) and Wan et al. 2015 (metazoan CF-MS) are added only if quick to obtain. They would
enter as an extra channel `cfms` in the human network (step (a) only), reported as secondary. STRING
`experiments` already includes much of BioGRID/IntAct, including AP-MS such as BioPlex.


## 2. Data (all under `data/extra/l2_network/`, gitignored; md5 lists copied to `results/l2_network/data_md5/`)

| source | version | files | URL | retrieved | licence | md5 (example) |
|---|---|---|---|---|---|---|
| STRING | v12.0 | `{taxid}.protein.links.full`, `.protein.sequences`, `.protein.info`, `.protein.aliases` for the 12 OrthoBench taxa (6239, 7227, 7719, 7955, 9031, 9598, 9606, 9615, 10090, 10116, 13616, 99883) and 4 more QfO taxa (3702, 4932, 284812, 8364) | https://stringdb-downloads.org/download/ | 2026-10-27 18:02 (OB) / 18:13 (QfO) CET | CC BY 4.0 | 9606 links.full `35e0baa7c65d921a3b501e496e5a11f0` |
| PaxDb | 5.0 | whole-organism integrated for 9606, 10090, 10116, 9031, 7955, 7227, 6239; chimp `9598-Chimp_iBAQ_Khan_2013`, dog `9615-GPM_2012_09` (no integrated set) | https://pax-db.org/downloads/5.0/datasets/ | 2026-10-27 18:02 CET | CC BY 4.0 (PaxDb terms) | 9606 `0131f9d6e02cd854dce66047104144a4` |
| hu.MAP 3.0 | 2024-09-29 complex pairs (+ 2023-05 all PPIs, downloaded, unused) | `humap3_complex_pairs_20240929.pairsWprob` | https://humap3.proteincomplexes.org/download | 2026-10-27 18:08 CET | not stated on the download page: **local use only, flagged** | `f213ad02ee410237462c86a784bfa52d` |

Wan et al. 2015 and McWhite et al. 2020 CF-MS (per-species co-elution networks) were **not obtained**. Without
a quick bulk download they fell outside the "only if quick" rule. In any case the binding constraint turned
out to be the partner species' networks, not human (§4.3).

**ID mapping and coverage** (`scripts/L2_network_prep.py`, `results/l2_network/mapping_coverage.json`).
Mapping is by Ensembl protein ID, then sequence md5, then gene. Mapped fraction of OrthoBench proteins:
- 1.00 for Ciona, fly, opossum, Tetraodon and chimp; 0.999 for worm (md5) and dog; 0.994 for chicken;
  0.955 for mouse; 0.839 for human; 0.838 for zebrafish.
- Human and zebrafish are lower because OrthoBench keeps alt-haplotype genes that STRING does not.

The own-species network after mapping (`own` = experiments ≥ 150 OR coexpression ≥ 400) is dense only for
human, mouse, rat, fly and worm. In nodes:

| species | experiments ≥ 150 | experiments ≥ 400 | coexpression ≥ 400 |
|---|---|---|---|
| human | 470k edges | 13,079 | 6,382 |
| mouse | 62k edges | 3,400 | 5,933 |
| fly | 45k edges | 4,216 | 5,794 |
| worm | 57k edges | 2,054 | 6,860 |
| zebrafish | – | 178 | 238 |
| dog | – | 118 | 30 |
| chicken | – | 133 | 344 |
| opossum | – | 4 | 0 |
| Ciona | – | 2 | 0 |
| Tetraodon | – | 0 | 0 |

The circular `database` channel, by contrast, is dense everywhere: 55k to 219k edges per species, including
Ciona, opossum and Tetraodon. That is what orthology-transferred pathway membership looks like.

QfO species for TreeFam-A dev (mapping through UniProt_AC aliases; nodes in the `own` network):

| species | nodes |
|---|---|
| human | 17,767 |
| mouse | 12,396 |
| Arabidopsis | 12,556 |
| worm | 11,884 |
| fly | 10,830 |
| S. cerevisiae | 5,820 |
| S. pombe | 4,244 |
| chicken | 680 |
| zebrafish | 375 |
| Xenopus | 5 |
| Ciona | 22 |

PaxDb mapped to OrthoBench proteins, with ppm > 0:

| species | proteins |
|---|---|
| mouse | 19,351 |
| human | 19,007 |
| fly | 12,656 |
| worm | 12,622 (via alias fallback) |
| zebrafish | 11,887 |
| chicken | 8,387 |
| dog | 4,303 |
| rat | 14,467 |
| chimp | 2,827 |

## 3. Methods

**Scripts** (branch `l2-network`; no `ancestra/` file changed, so frozen presets are untouched; `pytest -q`:
152 passed, 5 skipped):
- `scripts/L2_network_prep.py`: STRING/PaxDb/hu.MAP → OrthoBench IDs.
- `scripts/L2_network_standalone.py`: step (a) on OrthoBench dev. `L2_RUN=humap` gives the hu.MAP run.
- `scripts/L2_network_treefam.py`: step (a) on TreeFam-A dev.
- `scripts/L2_network_error_coverage.py`, `scripts/L2_network_treefam_bound.py`, `scripts/L2_network_rescue_precision.py`: diagnostic oracle bounds and rescue precision; provenance in `results/l2_network/manifest.json`.

**Score kernel.**
- The interolog count is the synteny `context_support` kernel with the PPI neighbourhood in place of the
  gene-order window.
- Anchors are 1:1 predicted orthologs of the frozen run: OrthoBench 8k–17k per species pair (e.g. 15,439
  for human–mouse); QfO from `cand2_accurate`.
- Anchors are family-excluded with the frozen S0 `families_final.npy`.

**Truth, candidates and decision rule.** Exactly S5/S9: 26,049 OrthoBench dev candidate pairs over the 66
species pairs, of which 4,485 are orthologs. The TreeFam-A dev candidate pairs (22,266, of which 13,717
are S relations) come from the reference relations of the 359 dev families. Their sequence score is a
DIAMOND `--sensitive` all-vs-all over the dev-family proteins.

## 4. Results

### 4.1 OrthoBench dev, standalone (`results/l2_network/standalone.json`, `standalone_humap.json`)

Coverage = evaluable (both genes have ≥ 1 anchored, family-filtered network neighbour). AUC on evaluable
pairs, with a 95% RefOG-bootstrap CI and three shuffled seeds. Hard queries: hard flag recomputed on the
evaluable candidates; "precision" = correct / (correct + wrong) when the layer decides.

| channel | evaluable cand. / orth. / hard pairs | AUC all [CI] | shuffled | AUC hard [CI] | hard queries n | decides: correct / wrong (precision) | seq top-1 | comb@0.9 − seq [CI] |
|---|---|---|---|---|---|---|---|---|
| **own (primary)** | 3.4% / 4.0% / 1.0% | **0.630 [0.532, 0.699]** | 0.50 / 0.49 / 0.50 | 0.501 [0.344, 0.687] | 17 | 4 / 5 (**0.44**) | 0.706 | **−0.176 [−0.237, 0]** |
| exp150 | 2.4% / 3.3% / 0.8% | 0.626 [0.547, 0.689] | 0.49–0.50 | 0.490 | 14 | 6 / 2 (0.75) | 0.857 | −0.071 [−0.083, 0] |
| exp400 | 0.6% / 0.8% / 0.4% | 0.665 [0.526, 0.804] | 0.49–0.50 | 0.614 | 6 | 2 / 1 | 0.833 | 0.000 |
| coexp400 | 0.6% / 0.7% / 0.2% | 0.624 [0.337, 0.732] | 0.50–0.53 | 0.667 | 0 | – | – | – |
| own + hu.MAP 3.0 (human) | 3.4% / 4.1% / 1.0% | 0.626 [0.530, 0.696] | 0.49–0.50 | 0.500 | 17 | 4 / 5 (0.44) | 0.706 | −0.176 [−0.237, 0] |
| hu.MAP only (human) | 2.7% / 3.3% / 0.9% | 0.586 [0.509, 0.656] | 0.50 | 0.565 | 15 | 2 / 4 (0.33) | 0.733 | −0.133 [−0.162, 0] |
| db400_CIRCULAR (reference only) | 13.7% / 23.2% / 14.5% | 0.673 [0.519, 0.815] | 0.49–0.50 | 0.523 | 237 | 47 / 26 (0.64) | 0.620 | +0.004 [−0.025, +0.071] |
| **abundance (PaxDb)** | 19.2% / 22.8% / 14.5% | **0.676 [0.604, 0.732]** | 0.50 / 0.50 / 0.51 | 0.689 [0.591, 0.753] | 122 | 39 / 77 (**0.34**) | 0.557 | **−0.131 [−0.246, −0.046]** |

**Controls.**
- Shuffled networks almost never decide a hard query (0 to 1 decision per seed), so their comb@0.9 − seq is
  0.00 and the "precision" entry is undefined, not 0.
- Shuffled abundance decides every query. Its precision is 0.23 to 0.29, which is the base rate of an
  ortholog among a hard query's candidates. Real − shuffled comb@0.9 for abundance is +0.08 to +0.12, with
  CIs that include 0. Abundance is therefore real signal, but it is below the sequence baseline.

**Per species pair** (own, evaluable pairs only; ≥ 30 pairs):

| species pair | AUC (seq AUC) | count > 0: orthologs / paralogs |
|---|---|---|
| human–mouse | 0.724 (0.870) | 58% / 18% |
| human–rat | 0.688 (0.842) | – |
| mouse–rat | 0.688 (0.781) | – |
| human–fly | 0.563 (0.750) | – |
| human–worm | 0.495 (0.714) | – |
| mouse–worm | 0.456 (0.825) | – |

Shuffled is 0.50 to 0.53 in every pair. The signal lives in mammal–mammal pairs, where sequence is already
right on almost every query (seq top-1 0.97 on all evaluable queries). It is absent in the
vertebrate–invertebrate pairs, which are the frozen run's weakest class (pair F1 0.744), and there is no
data for the teleost ohnolog pairs that make up most hard cases.

**Gate (a), primary `own`:**

| criterion | outcome | verdict |
|---|---|---|
| (1) AUC CI lower bound 0.532 > max shuffled 0.503 | holds | pass |
| (2) hard-query precision 0.44 > 0.5 | does not hold | fail |
| (3) comb@0.9 − seq = −0.176 ≥ 0 | does not hold | fail |

→ **FAIL**. Abundance meets (1) and fails (2) and (3) → **FAIL**.

### 4.2 TreeFam-A dev, standalone (extension, deviation 1; `results/l2_network/standalone_treefam_dev.json`)

| channel | evaluable cand. / orth. / hard pairs | AUC all [CI] | shuffled | AUC hard [CI] | hard queries n | correct / wrong (precision) | shuffled precision (decisions) | seq top-1 | comb@0.9 − seq [CI] |
|---|---|---|---|---|---|---|---|---|---|
| **own (primary)** | 25.8% / 29.2% / 18.7% | **0.622 [0.559, 0.663]** | 0.50 / 0.50 / 0.50 | 0.569 [0.501, 0.625] | 161 | 27 / 25 (**0.52**) | 0.27 (15), 0.67 (9), 0.33 (3) | 0.286 | **−0.025 [−0.075, 0]** |
| exp150 | 23.6% / 26.3% / 17.2% | 0.617 [0.547, 0.657] | 0.50 | 0.578 | 142 | 20 / 20 (0.50) | – | 0.282 | −0.042 [−0.107, −0.007] |
| exp400 | 6.3% / 7.6% / 4.0% | 0.639 [0.573, 0.692] | 0.50 | 0.547 | 35 | 6 / 10 (0.38) | – | 0.286 | −0.057 [−0.148, 0] |
| coexp400 | 3.3% / 4.8% / 0.6% | 0.629 [0.541, 0.681] | 0.50–0.51 | 0.631 | 2 | 1 / 0 | – | – | – |
| db400_CIRCULAR | 45.6% / 44.4% / 50.4% | 0.504 [0.456, 0.595] | 0.49 | 0.477 | 447 | 71 / 160 (0.31) | – | 0.304 | +0.009 [−0.022, +0.062] |

**Per species pair** (own, evaluable, ≥ 50 pairs; 23 pairs): AUC 0.44 to 0.77, and shuffled 0.46 to 0.52
throughout.

| species pair | AUC |
|---|---|
| yeast–Arabidopsis | 0.77 |
| yeast–human | 0.73 |
| human–Arabidopsis | 0.73 |
| worm–yeast | 0.72 |
| fly–yeast | 0.72 |
| human–chicken | 0.71 |
| mouse–Arabidopsis | 0.70 |
| mouse–worm | 0.69 |
| human–mouse | 0.68 |
| mouse–yeast | 0.68 |
| fly–Arabidopsis | 0.66 |
| human–fly | 0.65 |
| fly–worm | 0.55 |
| mouse–fly | 0.53 |
| mouse–S. pombe | 0.44 |
| fly–S. pombe | 0.49 |
| worm–S. pombe | 0.50 |
| S. pombe–human | 0.44 |

The deep yeast–animal pairs where Bandyopadhyay et al. (2006) found network signal do show it in the pooled
AUC, but not in the within-query decision.

**Gate (a), TreeFam-A:**

| criterion | outcome | verdict |
|---|---|---|
| (1) | holds | pass |
| (2) precision 0.52 > 0.5 but not > shuffled seed 23 (0.67 on 9 decisions) | partly holds | fail |
| (3) comb@0.9 − seq = −0.025 | does not hold | fail |

→ **FAIL**.

### 4.3 Why step (b) was not run, and what it could have gained (diagnostic, post hoc)

**OrthoBench dev** (`error_coverage.json`). The frozen accurate run makes 953 FN and 774 FP among the 19,146
dev truth pairs. Only 21 FN and 39 FP are network-evaluable.
- 7 of those FNs have support > 0, and 7 of those FPs have support 0.
- Fixing all 14 would move micro pair F1 from 0.8035 to 0.8051, i.e. **+0.0016**. That is an oracle
  ceiling, far below the ±0.01–0.03 needed to clear the CI rule.
- The circular `database` channel, if it were allowed, would reach +0.024. This shows that the bound is
  set by coverage, not by the kernel.

**TreeFam-A dev** (`treefam_bound.json`, `treefam_rescue_precision.json`). Of 5,450 dev FNs of the frozen
QfO accurate run, 1,333 are evaluable and 431 have support > 0. The oracle ceiling is TPR 0.603 → 0.634 and
PPV 0.940 → 0.951. It assumes the network fixes only errors and never breaks a correct pair, which it
does: hard-case precision is 0.52.
- **Used as a recall rescue,** adding unpredicted pairs with support ≥ k, the precision is 0.48 (k = 1),
  0.52 (k = 2), 0.53 (k = 3) and 0.54 (k = 5).
- Shuffled networks almost never reach k ≥ 2.
- Each rescue would therefore add about one false pair per true pair, against a merge that works at PPV
  0.94. It would fail the PPV guard immediately.
- Cause: paralogs share interaction partners (complex paralogs, shared pathways). This is the same failure
  mode as 3Di ("paralogs share folds") and expression (S9).

### 4.4 Cost
- **Downloads:** 1.9 GB in about 3 min.
- **Mapping and edge extraction:** 313 s CPU, 12 species (`mapping_coverage.json`).
- **Standalone runs:** OrthoBench 195 s, hu.MAP 100 s, TreeFam-A 581 s (including DIAMOND over 3,581
  proteins), diagnostics under 30 s.
- **Total:** ≈ 0.4 CPU-h at 1 to 4 threads. 0 GPU-h. Peak RAM < 8 GB.
- **No merge runs.**

## 5. Verdict and recommendation
- **Network neighbourhood conservation (interologs): do not integrate.** The signal is real (AUC 0.62 to
  0.63 on two benchmarks, shuffled 0.50). But it does not choose correctly among near-equal paralogs, and on
  OrthoBench it covers 1% of hard pairs. On QfO, coverage would be adequate for model eukaryotes, but a
  rescue at precision 0.5 is unusable.
- **Abundance: do not integrate.** It repeats S9 (expression), with a hard-query precision of 0.34.
- **No frozen config** is produced, because no arm passed. Nothing in `ancestra/` changed.
- **Annotation-only use is possible.** Per-HOG interolog conservation could be reported in the output
  (as S9 does for expression). It must not feed merge decisions.

## 6. Literature (verified on PubMed, 2026-10-27)
- Walhout AJ et al. (2000) Protein interaction mapping in C. elegans using proteins involved in vulval
  development. *Science* 287:116–122. doi:10.1126/science.287.5450.116. Coined "interologs".
- Matthews LR et al. (2001) Identification of potential interaction networks using sequence-based searches
  for conserved protein-protein interactions or "interologs". *Genome Res* 11:2120–2126.
  doi:10.1101/gr.205301.
- Bandyopadhyay S, Sharan R, Ideker T (2006) Systematic identification of functional orthologs based on
  protein network comparison. *Genome Res* 16:428–435. doi:10.1101/gr.4526006.
- Singh R, Xu J, Berger B (2008) Global alignment of multiple protein interaction networks with application
  to functional orthology detection (IsoRank). *PNAS* 105:12763–12768. doi:10.1073/pnas.0806627105.
- Szklarczyk D et al. (2023) The STRING database in 2023. *Nucleic Acids Res* 51:D638–D646.
  doi:10.1093/nar/gkac1000.
- Huang Q et al. (2023) PaxDb 5.0. *Mol Cell Proteomics* 22:100640. doi:10.1016/j.mcpro.2023.100640.
- Fischer SN et al. (2025) hu.MAP3.0: atlas of human protein complexes by integration of >25,000 proteomic
  experiments. *Mol Syst Biol* 21:911–943. doi:10.1038/s44320-025-00121-5.
- Wan C et al. (2015) Panorama of ancient metazoan macromolecular complexes. *Nature* 525:339–344.
  doi:10.1038/nature14877. Not used: data not obtained.
- McWhite CD et al. (2020) A pan-plant protein complex map reveals deep conservation and novel assemblies.
  *Cell* 181:460–474. doi:10.1016/j.cell.2020.02.049. Not used: data not obtained, and there are no plants
  in OrthoBench.

## 7. Failed or not pursued (for S1 Text)
- **Interolog support** (STRING v12 own-species experiments and coexpression; anchors = the method's own
  1:1 orthologs) separates orthologs from paralogs, with AUC 0.63 on OrthoBench dev and 0.62 on TreeFam-A
  dev against 0.50 shuffled. Among near-equal paralogs it is right only 44 to 52% of the time, and it
  lowers top-1 accuracy (−0.18 and −0.03). Not integrated.
- **Coverage.** On OrthoBench, only 1% of dev hard pairs have own-species network data: teleosts, chicken,
  dog, opossum and Ciona lack interactomes. The oracle ceiling for the frozen run is +0.002 pair F1.
- **As a recall rescue on TreeFam-A dev,** network-supported unpredicted pairs are orthologs only 48 to 54%
  of the time, against PPV 0.94 for the merge.
- **PaxDb abundance:** AUC 0.68, but right only 34% of the time when it decides a hard query
  (comb@0.9 −0.13 [−0.25, −0.05]). Not integrated.
- **Not pursued:**
  - CF-MS beyond hu.MAP 3.0. hu.MAP adds no coverage, because the partner species lack networks.
    Wan 2015 and McWhite 2020 were not obtained.
  - IsoRank-style global network alignment, made pointless by the local-signal result.
  - STRING `database`, `*_transferred`, textmining and homology channels: orthology-circular by
    construction.

## 8. Deviations from the pre-registration
1. **TreeFam-A dev was used for step (a)**, with the identical score, channels, controls and gate. The
   pre-registration mentioned TreeFam-A only for step (b). It was added because OrthoBench lacks yeast and
   has too few evaluable deep pairs to test the classic network-ortholog setting. It changes no decision:
   both benchmarks fail.
2. **hu.MAP 3.0** was run as the secondary channel described in §1.6 (human only), and Wan 2015 was not.
3. **Post-hoc diagnostics** (§4.3: oracle bounds, rescue precision) were not pre-registered. They are
   explanatory only and were not used to decide anything.
4. **Wording of gate criterion (2):** when a shuffled seed makes no decision, its precision is undefined.
   This was treated as "not above", which favours the real layer.
