# S7: isoforms (canonical vs best-isoform-pair scoring)

Status: 2026-10-26, T+2h50 to T+3h12 (`scripts/clock.sh`), inside the 1.5 h box. Standalone sub-benchmark on OrthoBench DEV only. No core file was edited.

All numbers come from `results/S7_isoforms/isoform_benchmark.json`. Checksums are in `results/S7_isoforms/manifest.json`. The data note is DATA.md §12.

**Verdict: KILL as a scoring layer (null result).**
- Best-isoform-pair scoring changes almost nothing. On hard queries the gain is +0.013 [−0.018, +0.044]; on the S5 hard set it is +0.022 [0.000, +0.070].
- It beats both decoy controls, but its gain over canonical-only is within noise and below any useful size.
- It would cost 2 to 5× more search, since human has 111k Ensembl proteins against 23k genes.

## 1. Data
**Ensembl 99 `pep.all`.** The 12 files were fetched; Ensembl checksums match.
- Every OrthoBench input protein maps to its Ensembl 99 gene, and its sequence equals one isoform of that gene.
- **It is the longest isoform in 100% of genes in all 12 species.** So OrthoBench "canonical" means "longest".

**Genes with more than one distinct isoform, genome-wide:**

| species | genes with >1 isoform |
|---|---|
| human | 16,806 / 23,485 |
| mouse | 13,481 / 22,866 |
| chimp | 11,854 / 23,534 |
| zebrafish | 10,878 / 30,313 |
| dog | 9,741 / 20,257 |
| Ciona | 558 / 16,671 |

Source: `work/S7/prep_stats.json`.

**DEV-RefOG gene set.** These are all genes in the S5 labels: 957 genes, 391 of them (41%) with at least one distinct alternative isoform. That gives 1,993 proteins after collapsing identical sequences within a gene; the most isoforms for one gene is 53.

**Search.** A single all-vs-all `diamond blastp --sensitive -k 0 --max-hsps 1 -e 1e-3 -p 6` took 13.6 s wall. Canonical and isoform scores come from the same run, so the comparison is fair.

## 2. Protocol
**Truth.** S5 labels (`work/S5/refog_pair_labels.parquet`): 26,049 cross-species pairs in 66 species pairs, of which 4,485 are orthologs.

**Scores** (`ancestra/isoforms.py`):
- `canonical`: max bits(canonical a, canonical b), over both directions.
- `best_iso`: max bits over all isoform pairs of the two genes.

**Queries** (S5 definition): a gene with at least one ortholog and at least one paralog candidate in the other species.
- Top-1 = 1 if the max ortholog score is greater than the max paralog score, 0.5 on a tie, and 0 otherwise.
- *Hard*: best paralog ≥ 0.9 × best ortholog, by canonical score.
- *S5 hard*: the exact 157 S5 hard queries, defined by whole-proteome bits on the 6 S5 species pairs.

**Controls** (seeds 11, 23, 47):
1. **Same-family decoy** (as specified). Each gene's k non-canonical isoforms are replaced by k non-canonical isoforms drawn from *other genes of the same homology component*, same species first.
2. **Cross-family decoy** (added). The k decoys are drawn from genes of *other* components of the same species. This is the "max over more sequences" null.

**CIs:** paired bootstrap over RefOGs, B=2000.

## 3. Results
| query set | n (RefOGs) | canonical | best_iso | **best_iso − canonical [95% CI]** | changed: fixed / broken | same-family decoy − canonical | best_iso − same-family decoy | cross-family decoy − canonical | best_iso − cross-family decoy |
|---|---|---|---|---|---|---|---|---|---|
| all queries | 4,680 (29) | 0.886 | 0.888 | +0.002 [−0.004, +0.012] | 28 / 14 | −0.23 [−0.31, −0.15] | +0.23 [+0.15, +0.31] | −0.03 [−0.05, −0.01] | +0.03 [+0.01, +0.06] |
| hard | 1,009 (26) | 0.470 | 0.483 | **+0.013 [−0.018, +0.044]** | 28 / 11 | −0.10 to −0.11 | +0.11 to +0.12 (CI lower ≥ +0.04) | −0.015 to −0.021 | +0.028 to +0.034 (CI lower −0.004 to +0.002) |
| hard, vertebrate pairs | 903 (24) | 0.482 | 0.501 | +0.019 [−0.012, +0.053] | 28 / 7 | −0.11 to −0.13 | +0.13 to +0.14 | −0.017 to −0.023 | +0.036 to +0.042 (CI lower ≥ +0.006) |
| S5 hard (6 pairs) | 157 (17) | 0.487 | 0.510 | +0.022 [0.000, +0.070] | 4 / 0 | −0.12 to −0.18 | +0.15 to +0.20 | −0.02 to −0.035 | +0.04 to +0.06 |

Control columns give the range over the 3 seeds.

**AUC (ortholog vs paralog).**

| score | pooled AUC (all pairs) | within-query mean AUC (4,680 queries) |
|---|---|---|
| canonical | 0.759 [0.576, 0.853] | 0.928 |
| best_iso | 0.760 [0.579, 0.854] | 0.931 |
| same-family decoy | 0.70 to 0.71 | 0.80 to 0.81 |
| cross-family decoy | 0.753 to 0.755 | 0.919 to 0.920 |

On vertebrate pairs only, the pooled AUC is 0.808 for canonical and 0.809 for best_iso.

**Orthologous-isoform evidence: how often the best isoform pair is not the canonical pair.**

| | orthologs (n=4,260 scored) | paralogs (n=14,633) |
|---|---|---|
| either gene has more than 1 isoform | 66.0% | 58.3% |
| best pair ≠ canonical pair | **21.2%** | 17.8% |
| … among pairs where either gene has more than 1 isoform | 32.1% | 30.5% |
| gain over 5% of bits | 6.1% | 3.4% |
| median relative gain when it differs | +1.7% | +1.4% |
| canonical pair missed, but an isoform pair hits | 0.2% | 0.3% |

**Reading.**
- For about one ortholog pair in five, a non-canonical isoform pair aligns slightly better. The gain is usually tiny (median +1.7% bits), and paralogs gain almost as often.
- Isoform choice therefore adds almost no ortholog-vs-paralog information beyond the longest isoform.
- On the hard cases, the 39 changed decisions split 28 fixed / 11 broken. That is a slight positive tendency, but the CI includes 0.
- The same-family decoy control (as specified) is *adversarial*, not a null. Isoforms of the query's own paralogs carry real paralog signal, so it lowers accuracy by 0.10 to 0.23. Beating it is therefore a weak test.
- The cross-family decoy is the fair null for "max over more sequences". It also lowers accuracy (−0.015 to −0.03), through spurious domain-level hits. Best-isoform scoring beats it by +0.03 to +0.04, with CIs mostly above 0.
- So real isoforms are "not harmful, and better than random extra sequences", but not better than canonical in any significant way.

## 4. Keep / kill (DESIGN §7 keep rule)
The keep rule needs Δ(real − shuffled) > 0 *and* Δ(real − off) above the seed spread on dev.
- **Δ(real − shuffled) > 0:** met for both controls.
- **Δ(real − off):** +0.013 on hard and +0.002 on all queries, with CIs [−0.018, +0.044] and [−0.004, +0.012]. This is not a demonstrated gain. The equivalence bound on all queries is ≤ +0.012.

**Recommendation: KILL as a merge/scoring layer.** Keep canonical (longest) scoring. Reasons:
- The cost is 2 to 5× more proteins searched.
- The cross-family control shows that max-over-isoforms adds noise hits at genome scale. The S7 test was restricted to DEV-RefOG genes, which flatters it.
- QfO `_additional.fasta` would need the same machinery.

Optional, at low priority: export the best-matching isoform pair per ortholog pair as annotation (orthologous isoform evidence). It is cheap if only computed for final ortholog pairs.

## 5. Files
- **Code:** `ancestra/isoforms.py`, `scripts/S7_isoform_prep.py` (data, mapping, DIAMOND), `scripts/S7_isoform_benchmark.py` (it imports S5 helpers) and `scripts/S7_S9_manifest.py`.
- **Tests:** `tests/test_expression_isoforms.py` (pep.all parsing incl. FlyBase and WormBase IDs, canonical/best-isoform pair scoring, decoy sets, set scoring).
- **Work:** `work/S7/{isoforms.parquet, isoforms.fa, hits.tsv, prep_stats.json, pair_scores.parquet, query_rows.parquet}`.
- **Rerun:** `python scripts/S7_isoform_prep.py && python scripts/S7_isoform_benchmark.py` (about 1 + 3.5 min).
