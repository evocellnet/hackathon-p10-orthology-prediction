# Verification of the manuscript (body.tex, S1–S3 Text, S1 Table)

Independent check, 2026-10-28 01:15–02:30 CET, by an agent that did not write the manuscript. Scope:
`manuscript/body.tex`, `manuscript/S1_Text.md` (§1–9 and header), `S2_Text.md`, `S3_Text.md`, `S1_Table.csv` (and
its README), the bibliography, cross-references, and the "Note" and "Agent setup" sections against the record.
Not checked: anything in `\pend{...}`, the pending scaling cells, Table `tab:layers2` and S1 Text §10 (verified in
`docs/VERIFICATION_L2.md`). Numbers were recomputed from raw JSON/TSV where these exist. No manuscript file was
edited. `manuscript/ancestra_manuscript.tex` is identical to body.tex from `\begin{document}` on, so every
body.tex fix must also go into it (or it must be regenerated).

Status codes: **OK**; **MISMATCH** (the text disagrees with the evidence); **UNSUPPORTED** (no evidence found);
**OVERSTATED** (the number is right but the wording claims more than the evidence shows).

## Summary

The headline numbers hold. That includes Tables 1, 2, 3 and 5 (tab:ob, tab:qfo, tab:layers, tab:trees), the
measured cells of tab:scaling, the TreeFam-A and SwissTree bootstrap results, the gene-loss percentages, the
guide-tree results and all CPU-hour figures. Every value in S1 Table equals `results/qfo_final/summary.tsv`. S2
and S3 Text are verbatim copies of their sources. All 58 `\cite` keys resolve.

The problems are in reporting and wording. In order of importance:

1. **Selective reporting of the QfO primary endpoint (most important).** DESIGN.md §5 names Hieranoid 2,
   InParanoid 5 and SonicParanoid2-sens-g as *primary comparators*, together with OrthoFinder and FastOMA. The
   body's bullet list (l.132–138) and S1 §9.3 leave all three out, and Ancestra fails against all three. The
   list shows OMA pairs and OMA HOGs instead, which were not pre-registered comparators (they were added at user
   request). The abstract, author summary and Discussion inherit this framing.
2. **Wrong count of re-scored QfO submissions** (l.106; S1_Table_README). 10 submissions were re-scored, not 20.
3. **Factual errors about prior work.** FastOMA uses MAFFT, not FAMSA (l.57). SYNERGY does not search ancestral
   proteomes, and it already derived ancestral-gene distances without new searches, which weakens the novelty
   sentences (l.36, l.38, l.64).
4. **Agent-setup record.** There were seven human interventions, not six; the missing one is the request for the
   second round of layers. The deadline was made soft, not reset. S2 Text §2 omits the same instruction. The
   Opus/Sonnet roles are partly wrong.
5. **Broken cross-reference.** `Fig~\ref{fig:scaling}` (l.231) has no target; it should point to
   `tab:scaling`, which is otherwise never referenced.
6. **GULO** was not one of the 2,523 apparent losses (l.174). It is absent in both human and chimp, so it fails
   the 5-of-6 rule; it was confirmed by a separate clade-level check.
7. **Overstated wording.** "most accurate method" and "as good as or better than methods from the OMA family";
   "6–7 points" (actual 5.7 and 6.5); "highest G-STD2 discordance" (only within Table 2); "structure made them
   worse" (the CI upper bound is 0.0); "close to linearly" (exponent 1.16 [1.10, 1.21] excludes 1); "Before any
   results ... split by homology component" (the split was replaced after three S0 configurations had been
   scored).
8. **FastOMA diagnosis** (l.103). OMAmer leaves only 4.1% of proteins unplaced, so "mainly from proteins that
   OMAmer did not place" is not what the evidence shows.
9. **Missing configuration statement.** Headline Ancestra results use the gene-order layer; the scaling runs do
   not. The manuscript says neither.
10. **Supplement errors.** S1 §8.3 gives wall hours as CPU-hours; S1 §9.1 says five straddling groups (three
    straddled); S1 §5 has a wrong scope and selective wording; S1 §6.2 contradicts the body; the S1 subtitle is
    out of date.
11. **Two bib entries** need cleaning (lefort2015fastme title, yates2019ensembl year, volume and pages).
12. **Four lower-priority items from VERIFICATION_L2 were not applied** to the GPU bullets and the round-2
    paragraph (#14, #15, #16, #19 there).

The full list of required fixes, with exact old → new text, is in section 2.

---

## 1. Claim table

### 1.1 body.tex: title, Note, Abstract, Author summary, Introduction

| location | claim | evidence (file / key) | status | correct value / replacement |
|---|---|---|---|---|
| l.18 Note | AI PI planned, delegated, ran analyses, wrote the manuscript; humans gave initial task + a few later instructions; no human has verified | README human header; S2 Text; transcript (fork check) | OK | — |
| l.21 Abstract | OB test (36 RefOGs) F 76.2 vs 69.7 OF, 70.4 SP2, 60.4 OMA HOGs, 33.1 FastOMA; 3.9 vs 8.2 CPU-h | `results/orthobench_final/summary.tsv` (76.1851; 69.7027; 70.4460; cpu_h 3.909, 8.225); summary.md OMA section (oma_browser 60.4, fastoma_v2 33.1) | OK | — |
| l.21 | gains over OMA and FastOMA significant; vs OF Δ CI [−6.5, +18.1] fails NI | summary.tsv delta rows (OF 6.4824 [−6.4895, 18.0548]); summary.md (oma_browser +15.7 [+1.8, +29.2], fastoma_v2 +43.1 [+14.0, +55.3]) | OK | — |
| l.21 | QfO: "matched FastOMA and OMA pairs on TreeFam-A recall but had clearly lower recall than OrthoFinder (−0.16)" | `results/verification/qfo_primary_ni.json` | OVERSTATED (incomplete) | It omits the pre-registered primary comparators where Ancestra is lower: SonicParanoid2-sens-g −0.087 [−0.132, −0.038], InParanoid5 −0.042 [−0.081, −0.000] and Hieranoid 2 −0.036 [−0.069, +0.004]. Replace with: "matched FastOMA and OMA pairs on TreeFam-A recall but had lower recall than OrthoFinder ($-0.16$), SonicParanoid2 and InParanoid" |
| l.21 | worst species-tree discordance among the main comparators | `results/qfo_final/summary.md` G_STD2_Eukaryota RF (0.257/0.266 vs ≤0.249 for Table 2 rows) | OK | (G-STD2 Eukaryota only.) |
| l.21 | CPU ∝ proteins^1.16 [1.10, 1.21] | `results/scaling/summary.md` slope table | OK | — |
| l.21 | five layers + eight more in round 2 | tab:layers (5 rows), tab:layers2 (8 rows) | OK | Two of the eight are second attempts at gene order and structure. |
| l.21 | 53% of apparent mammalian losses not real | `results/S8_geneloss/summary.json` all: (609+199+481+28+31)/2523 = 0.534 | OK | — |
| l.24 Author summary | "It was as good as or better than methods from the OMA family" | qfo_primary_ni.json: TF-A TPR vs OMA HOGs inconclusive (+0.002 [−0.058, +0.067]); SwissTree TPR vs FastOMA/OMA HOGs inconclusive; PPV NI not shown vs FastOMA; EC 0.899 vs FastOMA 0.928 / OMA pairs 0.952; G-STD2 worse than all OMA rows | OVERSTATED | "It scored higher than OMA-based methods on one benchmark and had similar recall on the other. It was cheaper than OrthoFinder, but it missed more distant orthologs than OrthoFinder and several pairwise methods." |
| l.24 | adding gene order, structure, expression or isoforms did not improve the final result | tab:layers | OK | Expression and isoforms were never run inside the method; "did not improve" is acceptable. |
| l.29–34 Intro | definitions, method classes, citations | rows_cite (§1.6) | OK | — |
| l.31 | pairwise methods' search cost grows with the square of the number of proteomes | LITERATURE.md | OK | — |
| l.36 | "SYNERGY ... and Hieranoid ... Both search ancestral proteomes against each other" | Wapinski et al. 2007: one all-vs-all of extant genes, then NJ-style distance updates; no ancestral search | MISMATCH | See fix F3 |
| l.38 | "Ancestral genes are never searched again: their similarity scores are aggregated from ... member proteins" (presented as the new idea) | SYNERGY did this too | OVERSTATED | Add "as in SYNERGY~\cite{wapinski2007synergy}" (fix F3) |
| l.40 | layers "have each been proposed as orthology evidence" with Bgee cited | Bgee is an expression database, not a proposal | OK (minor) | "have each been proposed or used as orthology evidence" |
| l.42 | pre-registered protocol and 48-hour budget | DESIGN.md; RESEARCH_TREE.md | OK | — |

### 1.2 body.tex: Overview of the method

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| l.49 | DIAMOND clusters at 30% identity; all proteins vs representatives; CC of hits with ≥30% coverage; Leiden split > 7,500 | `docs/S0_families.md` §5 config (`approx_id 30`, `s2_cov 0.3 shorter`, `s2_gmax 7500`) | OK | — |
| l.49 | one clustering pass ≤ 59%; all vs reps 79% | S0_families.md table (deep_c20 0.588, mm_s7.5_c50 0.588; selected 0.794) | OK | — |
| l.51 | MMseqs2 re-aligns missed pairs, "6–10% of pairs on simulated data" | `docs/core_v0.md` l.293 (DIAMOND reports 90–94%); elsewhere 5–13% on the toy (l.60, C-CORE-2) | OK | — |
| l.51 | score = bits / larger self-score | `ancestra/merge.py` `norm_max` | OK | — |
| l.55 | mutual best tolerance 0.95 | `ancestra/merge.py` MergeConfig.tol = 0.95 | OK | — |
| l.57 | "aligns ... with FAMSA and builds a tree with FastTree, as FastOMA and OrthoFinder 3 also do" | FastOMA paper: MAFFT + FastTree; OF3: FAMSA + FastTree | MISMATCH | Fix F4 |
| l.57 | tree mode only at divergent nodes with an outgroup; fast preset uses NJ; species overlap labels | merge.py (`tree_min_d 0.5`, `tree_require_outgroup`); pipeline.py PRESETS | OK | — |
| l.58 | deep nodes = median seed −ln S ≥ 1.5; species-disjoint merge | pipeline.py PRESETS `deep_min_d 1.5`, `deep_so_absorb` | OK | — |
| l.60 | average linkage, validity check incl. pyHam | merge.py `link avg`; `results/orthobench_final/summary.md` §4 (verify_orthoxml PASS) | OK | — |
| l.62 | synteny re-ranks partners in step 1, only within the gene-order-free root HOG | pipeline.py comment ("constrained to the synteny-off root HOGs; plus the R2/R3 veto and R3 partner choice") | OK (incomplete) | Also vetoes in-paralog/repair joins. Optional: "...re-ranks candidate partners in step 1 and can veto in-paralog joins, but only within..." |
| l.64 | "What is new ... without searching ancestral genes again" | see l.36/38 | OVERSTATED | Fix F3 |
| Fig 1 caption | stages and optional layers | — | OK | — |

### 1.3 body.tex: OrthoBench (Table 1 = tab:ob)

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| l.76 | 70 RefOGs, dev 34 / test 36, split by homology component | `splits/orthobench_split.json` | OK | — |
| tab:ob accurate | 76.2 / 69.5 / 84.2 / 3.9 CPU-h | summary.tsv main accurate_go (76.1851, 69.5434, 84.2294, 3.909) | OK | — |
| tab:ob fast | 74.8 / 68.7 / 82.1 / 2.9; ΔF +1.4 [−3.6, +7.1] | summary.tsv fast_go (74.7896, 68.6799, 82.0925, 2.884); ablation_delta fast_go test −1.3955 [−7.0616, 3.5817] (sign flipped) | OK | — |
| tab:ob OF | 69.7 / 62.1 / 79.5 / 8.2; +6.5 [−6.5, +18.1] | summary.tsv | OK | — |
| tab:ob SP2 | 70.4 / 70.8 / 70.1 / 7.9; +5.7 [−4.3, +14.9] | summary.tsv (7.854; 5.7391 [−4.2922, 14.8623]) | OK | — |
| tab:ob OMA browser | 60.4 / 88.9 / 45.8; +15.7 [+1.8, +29.2] | summary.md OMA table | OK | — |
| tab:ob ProteinOrtho | 47.5 / 97.7 / 31.4 / 1.6; +28.7 [+13.1, +41.9] | summary.tsv | OK | — |
| tab:ob FastOMA patched | 33.1 / 23.9 / 53.8 / 19.9; +43.1 [+14.0, +55.3] | summary.md fastoma_v2 | OK | The 19.9 CPU-h excludes about 0.75 CPU-h of OMAmer placement reused from v1 (summary.md "Cost"); optional footnote. |
| tab:ob caption | 2,000 resamples; single runs; OMA has all 12 species | summary.md header; NOTES_extra.md | OK | Optional: add "Ancestra used 10 threads and the competitors 8." |
| tab:ob / Abstract | the Ancestra rows include the gene-order layer | test_log.tsv "accurate_go (MAIN: preset accurate + gene order)" | UNSUPPORTED in text | Not stated anywhere. Add to the tab:ob caption: "Ancestra used the gene-order layer (Ensembl 99 coordinates); without it, test F was 76.1." The same applies to tab:qfo (cand2 runs had gene order). |
| l.101 | "Ancestra had the highest F" | summary.tsv | OK | Point estimate. |
| l.101 | "significantly better than OMA HOGs, FastOMA and ProteinOrtho" | CI lower bounds +1.8, +14.0, +13.1 | OK | — |
| l.101 | "the point estimates favoured Ancestra by 6–7 points" | +6.5 (OF), +5.7 (SP2) | MISMATCH | "by 6.5 and 5.7 points" |
| l.101 | CIs include losses beyond the 3-point margin | lower bounds −6.5, −4.3 | OK | — |
| l.101 | all 70: +7.8 [+0.5, +15.1] | summary.tsv delta all 7.7764 [0.5278, 15.0692] | OK | — |
| l.103 | FastOMA 0.5.1 bug dropped 5,497 genes | `bench_repo/results/orthobench_competitors/OMA_RESUME_NOTE.md` l.22 | OK | — |
| l.103 | coverage 79.0% → 80.6%; test F 33.5 → 33.1 | SUMMARY.md l.116 | OK | — |
| l.103 | "The low score comes mainly from proteins that OMAmer did not place and from fused groups at the root" | OMA_RESUME_NOTE.md l.21–22: OMAmer leaves only 4.10% unplaced; the 21% absent genes = 9.1% never in a FastOMA family + 9.7% root-level singletons + 2.2% bug; precision 23.9 | MISMATCH (partly) | "The low score comes from low precision (23.9: reference orthogroups fused into larger root HOGs) and from the 20% of proteins absent from FastOMA's groups (9.1% never placed in a family, 9.7% left as root-level singletons; OMAmer itself left only 4.1% unplaced), consistent with..." |
| l.103 | consistent with Emms et al.'s FastOMA fusion rate | OF3 paper: FastOMA fused 10% of RefOGs vs 2.85% | OK | — |
| l.103 | OMA browser HOGs P 88.9, R 45.8 | summary.md | OK | — |

### 1.4 body.tex: QfO (Table 2 = tab:qfo)

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| l.106 | "we re-scored 20 published submissions and reproduced all 60 published values that we checked" | `results/qfo_peek1/competitor_reproduction.tsv`: 10 methods × {SwissTrees, TreeFam-A, VGNC} × {TPR, PPV} = 60 rows, all abs_diff 0 | MISMATCH | Fix F2 |
| l.106 | 78 species, 988,778 proteins | `results/qfo_run/verify_ids.json` n_proteins 988778 | OK | — |
| l.106, l.366 | "S1 Table has all 21 methods" | S1_Table.csv: 22 rows = 20 published + 2 Ancestra presets | OK (wording) | "all 20 published QfO 2022 methods and both Ancestra presets" |
| tab:qfo all cells | TPR/PPV, EC (n, thousands), G-STD2 Euk RF (n) | `results/qfo_final/summary.md`; rounding checked (EC n 180,745→181, 196,729→197; RF .2656→.266, .2385→.239, .2315→.232) | OK | — |
| tab:qfo row label | "SonicParanoid2 (sens.)" | values are SonicParanoid2-sens-g; S1 Table also has SonicParanoid2-sens (0.744/0.864, 0.682/0.913) | MISMATCH (ambiguous) | "SonicParanoid2 (sens-g)" |
| tab:qfo row label | "FastOMA" | published submission FastOMA_v0.1.4_nov; the OrthoBench row is 0.5.1 | OK (minor) | Optional: "FastOMA (v0.1.4)" |
| l.132 | TreeFam-A recall on held-out families pre-registered as primary | DESIGN.md l.64 | OK | — |
| l.134–137 | vs FastOMA +0.043 [−0.005, +0.093] NI; OMA pairs +0.018 [−0.028, +0.069] NI; OMA HOGs +0.002 [−0.058, +0.067] inconclusive; OF −0.160 [−0.210, −0.106] inferior | qfo_primary_ni.json (0.0433 [−0.0045, 0.0933]; 0.0178 [−0.0284, 0.0687]; 0.0022 [−0.0577, 0.0666]; −0.1595 [−0.21, −0.1057]); n_units 749 | OK (numbers) | — |
| l.132–138 | the list of comparisons | DESIGN.md l.62: primary comparators are Hieranoid 2, OrthoFinder, FastOMA, InParanoid 5, SonicParanoid2-sens-g. qfo_primary_ni.json: Hieranoid 2 −0.0358 [−0.0689, +0.0036]; InParanoid5 −0.0422 [−0.081, −0.0001]; SonicParanoid2-sens-g −0.0865 [−0.132, −0.0378] | OVERSTATED (selective) | Three of the five pre-registered primary comparators are missing, and Ancestra fails NI against all three. OMA pairs and OMA HOGs were not pre-registered. Fix F1. |
| l.139 | SwissTree recall vs OF −0.103 [−0.177, −0.029] | qfo_primary_ni.json SwissTrees (−0.1033 [−0.1772, −0.0289]; n_units 18) | OK | — |
| l.139 | "Ancestra had the highest G-STD2 discordance" | summary.md: within Table 2 yes; Ensembl Compara 0.303, pbDOOD 0.350/0.353, SonicParanoid2-sens 0.283, PANTHER all 0.272 and MetaPhOrs 0.263 are higher | OVERSTATED | "Ancestra had the highest G-STD2 Eukaryota discordance of the methods in Table 2" |
| l.139 | lower EC than OF | 0.8985 vs 0.9311 (SEM ~0.001) | OK | — |
| l.141 | "sits with the OMA methods: precise, with moderate recall" | TF-A PPV 0.921 vs OMA 0.920–0.923 | OK | — |
| l.141 | 78% of missed pairs split at the pair's LCA | `results/qfo_dev_diag/diag_summary.json` (ii 2,239 + iii 2,411)/5,984 = 0.777 | OK | — |
| l.141 | "These are deep nodes such as Opisthokonta or Eukaryota, where there is often no outgroup" | diag_summary.json FN_level_by_category (ii + iii): Bilateria 34.1%, Eukaryota 20.6%, Opisthokonta 12.7%, Chordata 7.4% | OVERSTATED (imprecise) | "These are mostly deep nodes (Bilateria 34%, Eukaryota 21% and Opisthokonta 13% of these misses), where an outgroup is often missing or distant." |
| l.141 | rules that recovered recall fused OB orthogroups (S1 §4.7) | S1 §4.7; docs/core_v0.md §7.11–7.12 | OK | — |

### 1.5 body.tex: extra data layers (Table 3 = tab:layers), trees (tab:trees), runtime, GPU

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| tab:layers gene order | AUC 0.92 human–mouse, 0.49 human–fly, shuffled 0.50 | `results/S5_synteny/tables.md` (w20: 0.9223, 0.4878; shuffled ~0.50) | OK | (These are w=20 values; the 70/157 decisions are w=10.) |
| tab:layers gene order | OB ΔF −0.1 [−0.4, 0.0] removed; QfO TF-A −0.012 [−0.028, +0.001] | summary.tsv accurate_nogo −0.1252 [−0.4244, 0.0159]; `results/verification/qfo_ablation_familyboot.txt` a_nogo −0.0124 [−0.0277, +0.0012] | OK | The QfO delta is for the `fast` preset (a_nogo vs cand2_fast); add "(fast)". |
| tab:layers 3Di | completeness 0.794 → 0.853 (shuffled 0.794); OB F −5.0 [−13.1, 0.0]; P 69.5 → 60.1 | `docs/S6_structure.md` l.13; summary.tsv accurate_go_rescue −5.0485 [−13.1003, 0.0000], P 60.07 | OK | — |
| tab:layers expression | AUC 0.667; hard −0.019 [−0.20, +0.15] | `docs/S9_expression.md` l.34, l.40 | OK | — |
| tab:layers isoforms | AUC 0.760 vs 0.759; hard +0.013 [−0.018, +0.044] | `docs/S7_isoforms.md` l.56, l.66–67 | OK | — |
| tab:layers genome | 92% recovery; 0/512 false hits | S8 summary.json controls pooled sens 0.922, neg 0/512 | OK | (Masked-locus controls had 5/600 false intact calls.) |
| l.165 | 70 of 157 hard cases at precision 0.84 | S5 tables.md syn_w10: 59 + 11 = 70, precision 0.8429 | OK | — |
| l.165 | dev gain 1.0 (accurate) to 4.4 (fast) | summary.tsv dev: 79.69 vs 78.73; 75.71 vs 71.34 | OK | — |
| l.167 | 3Di for 94.5% of QfO proteins | `docs/structure_feasibility.md` l.126 (928,383 / 982,789) | OK (minor) | Of the 988,778 run genes it is 93.6% (GPU_v2.md l.168: 6.4% without). |
| l.167 | completeness above the 0.82 sequence ceiling; AUC 0.694 vs 0.751 deep | S6_structure.md l.13, l.20; S0_families.md l.101 | OK | — |
| l.169 | expression wrong more often than sequence on the hard cases; isoforms 2–5× more search | S9 l.40–41; S7 l.10 | OK | 2–5× is a projected cost ("would cost"), not a measured one. |
| l.171–176 | 2,523 losses; 32% / 19% / 2% / 25% / 21%; 53% not losses | S8 summary.json: (609+199)/2523 = 32.0%; 481 = 19.1%; 28+31 = 2.3%; 351+287 = 25.3%; 257+63+217 = 21.3% | OK | — |
| l.174 | "including the known pseudogenes GULO and UOX in primates" (among the 2,523 losses) | `docs/S8_geneloss.md` l.178: "GULO is not among the candidates: it is absent in both human and chimp, so it fails the 5-of-6 rule"; UOX is a human candidate | MISMATCH | "...including the known UOX pseudogene in human; GULO, which is missing in both human and chimp, was confirmed as a pseudogene by a separate clade-level check;" |
| l.178 | round-2 arms and keep rules committed before results; shuffled controls for each | VERIFICATION_L2 #16 | OVERSTATED (L2 #16 not applied) | Apply L2 #16: add "apart from logged post hoc arms" and "shuffled controls with the exceptions noted in S1 Text §10". |
| l.202 | round-2 paragraph | VERIFICATION_L2 | (skipped, verified there) | — |
| tab:trees | conflict 0.000 / 0.120 / 0.173 / 0.987 | `docs/S1_speciestree.md` §3 table (published column) | OK | — |
| tab:trees | OB ΔF ASTRAL/FastME +0.0 [−0.3, +0.4]; NNI −3.3 to +0.3; random −16.4 [−27.1, −5.6] | summary.tsv ablation_delta (0.0291 [−0.3173, 0.4294]; NNI −3.2754 / −0.0020 / +0.2667; random −16.3891 [−27.0928, −5.5664]) | OK | — |
| tab:trees | QfO ΔTPR/ΔEC +0.028/−0.125, +0.022/−0.134, −0.243/−0.060; TPR CIs include 0 | `results/qfo_final/summary.md` ablation (0.028/−0.1248; 0.022/−0.1343; −0.243/−0.0602); familyboot c_astral [−0.0090, +0.0774], d_fastme [−0.0009, +0.0446] | OK | — |
| l.228 | random tree −16 points | as above | OK | — |
| l.228 | data trees put *C. elegans* "outside Bilateria" | S1_speciestree.md l.73: "*C. elegans* falls outside the other bilaterians" | MISMATCH (minor) | "(for example *C. elegans* placed outside the other bilaterians)" |
| l.228 | conflicted with 25–30% of bacterial splits | S1_speciestree.md l.79 | OK | — |
| l.228 | small non-significant QfO recall gain, large EC cost; two policies not run | as above; DESIGN.md l.39 | OK | — |
| l.231 | "Fig~\ref{fig:scaling}" | no `\label{fig:scaling}` exists in body.tex or ancestra_manuscript.tex; the scaling data is `tab:scaling` (never referenced) | MISMATCH (broken ref) | "(8, 16, 32, 64 and 78 species; Table~\ref{tab:scaling})" |
| tab:scaling measured cells | 0.78/2.16/6.02/16.3/18.5; 0.47/0.99/3.20/10.6/11.6; OF 1.22/3.54/14.4; SP2 0.77/2.34/11.1; FastOMA 1.83; slopes 1.16 [1.10, 1.21], 1.21 [1.02, 1.41], 1.42, 1.53 | `results/scaling/summary.md` | OK | — |
| tab:scaling | Ancestra rows | summary.md: "accurate preset, seq-only" | UNSUPPORTED in text | Gene order was off in the scaling runs, unlike the headline runs. Add "(without gene order)" to the Ancestra row labels or the footnote. |
| l.255 | 18.5 CPU-h, 29 GB; 8.6 S2, 7.0 merge (~61,000 jobs), 3.0 S0; fast merge 2.5 CPU-min | scaling summary per-stage anc_78 (30,928 / 25,040 / 10,707 CPU-s; 28.9 GB); anc_78_fast merge 150 s; `docs/qfo_run.md` l.274 (61,216 jobs) | OK | — |
| l.255, l.272 | "grew close to linearly" / "a cost that grows close to linearly" | slope 1.16, 95% CI [1.10, 1.21], excludes 1 | OVERSTATED (mild) | "grew slightly faster than linearly (exponent 1.16)" |
| l.260 | S2 = 60–75% of wall time | `docs/GPU_v2.md` l.25 | OK | — |
| l.262 | MMseqs2-GPU ~9 ms per query; DIAMOND 2.2.8 no GPU | GPU_v2.md l.18, l.45, l.68 | OK | Optional: "In our measurements, its GPU prefilter needed about 9 ms per query" (it is not a figure from the paper). |
| l.263 | 600–640 GCUPS; exact vs CPU reference; deterministic; S2 1.2–1.3×; e2e 1.18× (1,027 vs 1,207 s); OG-F 80.2 vs 79.7, +0.5 [−3.5, +3.3]; TF-A 0.614/0.929 vs 0.615/0.922 | GPU_v2.md l.19, l.133–140, l.305, l.390 | OK | "CPU reference" = a pure-Python Gotoh implementation (683 pairs). |
| l.264 | all-GPU bullet | VERIFICATION_L2 rows 102–127 | OK except L2 #15, #19 not applied | "with both differences within the CIs" → "(the CIs of both differences include 0; this was the best of seven arms)"; "about 15 times fewer CPU-hours" → "an estimated 15 times fewer CPU-hours (plus 0.2–0.3 GPU-hours)" |
| l.265 | "About 93% of each FastTree call" | VERIFICATION_L2 #14 (20 five-sequence jobs only) | OVERSTATED (L2 #14 not applied) | "About 93% of a FastTree call (in a sample of 20 small jobs)" |
| l.265 | batching saved 1.6% of merge time | L2 #10 (applied) | OK | — |
| l.266 | ProstT5 6,100–7,300 res/s, 40–60×, 99.9% agreement, 48,903 proteins in < 2 h | GPU_v2.md l.22, l.150–172 (1 h 56 min) | OK | — |
| l.267 | ESM-2+FAISS ~10 min; up to 85% (k = 25) vs 82%; 42% outside family | GPU_v2.md l.23, l.218–227 (592 s + 4 s; 29/34 = 0.853) | OK | — |
| l.269 | "modest, accuracy-neutral speed-up" | GPU_v2.md (b) | OK | — |

### 1.6 body.tex: Discussion and Limitations

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| l.272 | "On OrthoBench it was the most accurate method we tested" | Δ vs OF and SP2 CIs include 0 | OVERSTATED | "On OrthoBench it had the highest F of the methods we tested (not significantly above OrthoFinder or SonicParanoid2)" |
| l.272 | "significantly better than OMA-based resources and ProteinOrtho, and cheaper than OrthoFinder" | as above; CPU 3.9 vs 8.2 (indicative) | OK | — |
| l.272 | "On QfO it matched the OMA methods but not OrthoFinder" | qfo_primary_ni.json (OMA HOGs inconclusive; fails vs SonicParanoid2-sens-g and InParanoid5) | OVERSTATED | "On QfO its TreeFam-A recall was non-inferior to FastOMA and OMA pairs, but below OrthoFinder, SonicParanoid2 and InParanoid." |
| l.274 | OrthoFinder roots trees with duplication information \cite{emms2017stride} | OF2 paper uses STRIDE | OK (minor) | Add `emms2019orthofinder` to the cite. |
| l.276 | "Neither improved the final groups on held-out data, and structure made them worse" | rescue ΔF −5.05 [−13.10, 0.00]; all-70 −2.69 [−8.38, +0.26]; QfO −0.004 | OVERSTATED (mild) | "...and structure lowered their precision (69.5 to 60.1; F −5.0 [−13.1, 0.0])" |
| l.276 | genome sequence exposed our own fission errors | S8 A_annotated + A_otherHOG = 32% | OK | — |
| l.278 | OB test intervals about ±12 F | OF CI half-width 12.3; SP2 9.6; OMA 13.7 | OK | — |
| l.278 | deviations: heavy dev reuse, TreeFam-A split after first QfO evaluation, OMA saw the species | DISCLOSURES #1, #2; NOTES_extra.md | OK | — |
| l.280 | about 30 approaches tried and dropped | S1 Text §1–8, 10 (~35 subsections) | OK (approx.) | — |
| l.286 | verification found a false disclosure, an omitted primary endpoint and too-narrow CIs | `docs/VERIFICATION.md` §5.1, §5.2, §5.4 | OK | Note that a second omission (3 primary comparators) is still in the text (F1). |
| l.288 | child-order sensitivity, 18 OGs differ | summary.md §5 Guide tree | OK | — |
| l.288 | 1% jitter −3 to −7, two intervals include 0, TF-A barely moved | VERIFICATION_L2 #11 (applied) | OK | — |
| l.288 | "Only 22 of 81 quantitative claims were independently recomputed" | CLAIMS.md: 81 rows; 9 VERIFIED + 13 VERIFIED-WITH-CAVEAT | OK (count) / OVERSTATED (scope) | "Of the 81 claims in our claim ledger, which covers results up to the freeze, 22 were independently recomputed; the others were checked against their evidence files. Post-freeze GPU and second-round results were checked separately." |

### 1.7 body.tex: Materials and methods

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| l.293 | PI did not write method code | transcript: no PI edits to `ancestra/*.py` (fork check) | OK | — |
| l.293 | Opus for design, algorithm, verification; Sonnet for data acquisition, benchmark runs, bookkeeping | transcript: data acquisition and literature were done by Opus agents; Sonnet did the environment, competitor/scaling runs, bookkeeping and reference checks | MISMATCH (minor) | "Claude Opus 5.5 for literature, data acquisition, design, algorithm and verification work, and Claude Sonnet 5.5 for the software environment, competitor and scaling runs, bookkeeping and reference checks." |
| l.297–301 | research tree, lessons, design doc (two rounds, eight reviewers), claim ledger, deviations list | RESEARCH_TREE.md, LESSONS.md, DESIGN.md, `docs/reviews/` (r1: 5, r2: 3), CLAIMS.md, DISCLOSURES.md | OK | — |
| l.303 | key results recomputed by an independent subagent; clock checks | VERIFICATION.md; RESEARCH_TREE.md milestone headers | OK | — |
| l.305 | "Humans intervened six times ... reset the deadline ..." | transcript; RESEARCH_TREE.md M6; `docs/L2_protocol.md` l.3 (user asked for more extra-data work); START.txt (deadline kept, made soft) | MISMATCH | Fix F5 |
| l.305 | 24-core shared node for 19.5 h, then 48 cores / 500 GB / A100 | README.md; RESEARCH_TREE.md M4 (pause 10:10 = T+19h30), M5 | OK | — |
| l.308–313 | "Before any results, the design document fixed: Splits. OrthoBench RefOGs were split by homology component..." | S1 §9.1, DISCLOSURES #13: the first split was random; it was replaced after review, by which time three S0 configurations had been scored on RefOGs that later fell into the test set | OVERSTATED | Add to the Splits item: "(This split replaced a first, random split that a reviewer found leaky; three early S0 configurations had been scored by then, S1 Text §9.1.)" |
| l.310 | QfO at most three evaluations before freeze | DESIGN.md §4; test_log.tsv (2 peeks before freeze) | OK | — |
| l.311–312 | NI rule −0.03; paired bootstrap 2,000, seed 11 | DESIGN.md l.65; summary.md | OK | — |
| l.315 | TreeFam-A dev split (30%) after first QfO evaluation; "18 other deviations" | DISCLOSURES.md items 1–19; `splits/treefamA_split.json` (359/1,108 = 32%) | OK | — |
| l.318 | DIAMOND 2.2.8, MMseqs2 18; deterministic partition; freeze-v1 with hashes | `git tag` freeze-v1; `results/freeze.json` (ancestra_sha256) | OK | — |
| l.318 | presets: accurate FastTree, fast NJ + species-overlap structure | pipeline.py PRESETS (fast so_struct True) | OK | — |
| l.322 | OMA browser All.May2026, cut at Bilateria; local runs of OF 3.2, SP2, PO, FastOMA 0.5.1 | SUMMARY.md l.77 | OK | — |
| l.323–324 | QfO official code, local; scaling nested domain-stratified subsets, `/usr/bin/time` | qfo_run.md; scaling summary.md | OK | — |
| l.329–343 | data-layer sources (Ensembl 99, AFDB/Foldcomp, Bgee, STRING v12 own channels, hu.MAP 3.0, PaxDb 5.0, Pfam-A 38.2 `--cut_ga`, ESM-2/ProtT5, contrastive head on families with no dev/test genes, introns, ProstT5; five agents; seeds 11/23/47) | `docs/L2_network.md` l.38, l.134–135; `docs/L2_domains.md`; `docs/L2_plm.md` l.143; RESEARCH_TREE.md M6 | OK | Seeds: exceptions listed in L2 #16. |

### 1.8 Cross-references and figures

| item | status | note |
|---|---|---|
| `fig:method`, `tab:ob`, `tab:qfo`, `tab:layers`, `tab:layers2`, `tab:trees` | OK | defined and referenced |
| `Fig~\ref{fig:scaling}` (l.231) | MISMATCH | undefined label → `Table~\ref{tab:scaling}` |
| `tab:scaling` | OK (defined) | never referenced (fixed by the line above) |
| "S1 Text, section 7" (l.62), "section 4.7" (l.141) | OK | §7.3 and §4.7 exist and match |
| `docs/figures/fig2_benchmarks.png` | note | exists but is not used in the manuscript |

### 1.9 S1 Text (header and §1–9)

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| subtitle (l.3) | "Supplement to: Ancestra — orthology inference by level-synchronous ancestral-genome reconstruction inside de novo protein families" | body.tex l.10 | MISMATCH | "Supplement to: Ancestra: orthology inference by level-wise reconstruction of ancestral genes inside de novo protein families, and a test of extra data layers" |
| header | dev = 34 RefOGs; TreeFam-A 30%; > 240 dev runs; OG-F ±2–4 | splits/*.json; `results/core_v0_orthobench/runs.tsv` (242 rows) | OK | — |
| §1.1, §1.2 | review numbers (2.0e11 vs 4.9e11 pairs etc.); precedents | `docs/reviews/r1_*.md` | OK | — |
| §2.1–2.3 | 0.794 / 0.722 test; ceiling 0.82; ≤ 0.59 single pass; variants | S0_families.md §3–5; dev_ceiling.json | OK | §2.3: the Leiden-20k row shows "—" but its completeness is 0.765. |
| §2.4 | ARI 0.78; 79.2 → 75.6 | core_v0.md §7.5–7.6 | OK | — |
| §2.4 | "All merge settings tuned on the old partition were re-checked on three partitions." | core_v0.md §7.6: only gene-order arms, shuffled controls and locus repair on three partitions | OVERSTATED | "The key merge arms (gene order, its shuffled controls and locus repair) were re-checked on three partitions, and the others on the deterministic one." |
| §2.5, §2.6 | giant families; 3Di rescue dev and test numbers | S0_giant.md; S6_structure.md; S0_rescue.md; summary.tsv | OK | — |
| §3, §4.1–4.9 | scoring and merge-arm numbers | runs.tsv rows; core_v0.md §3–7 | OK | §4.7 heading: "(DEV; two pre-fixed grids)" (the deep-tolerance rows come from a later grid). |
| §5 text | "projected cost exceeded the bar by 1%" | `docs/S4_trees.md` §3.3: 1.01 CPU-h by direct measurement, 6.7 CPU-h by the pre-registered size fit | OVERSTATED | "because projected cost exceeded the 1 CPU-h bar (1.01 CPU-h by direct measurement, 6.7 CPU-h by the pre-registered size fit)" |
| §5 table | "rejected: 3.3 CPU-h on the subset" | S4_trees.md §3.3: projected over all 15,456 OrthoBench clusters | MISMATCH | "rejected: 3.3 CPU-h projected for OrthoBench (bar 1 CPU-h)" |
| §6.1 | conflicts, ΔF, ΔTPR, EC, random | S1_speciestree.md; summary.tsv; familyboot | OK | — |
| §6.2 | "The design planned three guide-tree policies (strict, delayed merge, per-family) ... Neither was done." | DESIGN.md l.39; body l.228 says "two other policies" | MISMATCH | "Besides the strict universal guide tree that Ancestra uses, the design planned two other policies (delayed merge, per-family guide trees) and a check with an independent simulator (Zombi or ALF). None of these was done." |
| §6.3, §7.1–7.5 | child order; synteny arms; test non-replication | DISCLOSURES #11; S5_synteny.md; core_v0.md §7.2–7.6; summary.tsv | OK | — |
| §8.1, §8.2 | isoforms, expression | S7, S9 docs | OK | — |
| §8.3 | "so QfO would need 490–980 CPU-hours" | `docs/structure_feasibility.md` l.16–35: 492 h / 983 h are wall hours on a 24-core node | MISMATCH | "so QfO would need 490–980 hours (20–41 days) on a 24-core node" |
| §9.1 | "five groups of paralogous RefOGs that straddled dev and test" | old split (`git show 6711002:splits/orthobench_split.json`): only {016,…,058}, {038,052}, {021,068} straddled | MISMATCH | "A reviewer agent found five groups of paralogous RefOGs, three of which straddled dev and test." |
| §9.2 | other deviations | DISCLOSURES | OK | — |
| §9.3 | table of unmet pre-registered endpoints | DESIGN.md l.62–65; qfo_primary_ni.json | OVERSTATED (incomplete) | Add rows: "QfO TreeFam-A TPR vs Hieranoid 2 / InParanoid5 / SonicParanoid2-sens-g \| −0.036 [−0.069, +0.004] / −0.042 [−0.081, −0.000] / −0.087 [−0.132, −0.038]: not shown / fails / fails" and "QfO TreeFam-A PPV \| non-inferiority not shown vs any primary comparator (shown only vs OMA HOGs and OMA pairs)" |
| §9.4 | FastOMA bug, OMA browser rule change, scaling symlink failure | OMA_RESUME_NOTE.md; DISCLOSURES #19; scaling summary.md | OK | — |
| S1 vs body | numbers in both places | — | OK | No numeric conflict except §6.2. |

### 1.10 S2 Text, S3 Text, S1 Table

| location | claim | evidence | status | correct value / replacement |
|---|---|---|---|---|
| S2 §1 | verbatim PROMPT.md | diff: identical | OK | — |
| S2 §3 | verbatim PROMPT_MANUSCRIPT.md | identical | OK | — |
| S2 §2 | list of later human instructions (6 items) | transcript 2026-10-27 16:57Z; RESEARCH_TREE.md M6; L2_protocol.md l.3 | MISMATCH (omission) | Add: "7. While the scaling runs finished, try further data layers (domain architecture, proteomics, co-fractionation, protein language models), rethink structure, expression and synteny, and revisit dropped GPU ideas (second round of layers, S1 Text section 10)." |
| S3 Part A | verbatim docs/DISCLOSURES.md | diff: only the H1 heading differs | OK | — |
| S3 Part B | verbatim docs/VERIFICATION.md | diff: only heading levels | OK | Optional editor note under B §5.1: "(Part A item 1 has since been corrected.)" |
| S1_Table.csv | all 22 rows × 37 columns | `results/qfo_final/summary.tsv`: 0 differences | OK | — |
| S1_Table_README.txt | "published OpenEBench QfO 2022 submissions, re-scored locally (local scoring reproduced 60/60 published values checked)" | competitor rows come from `bench_repo/reference_results/openebench_QfO_2022.tsv`; 10 methods re-scored | MISMATCH | "Other rows: published OpenEBench QfO 2022 values. As a check of the local setup, 10 of these submissions were re-scored locally on SwissTree, TreeFam-A and VGNC and reproduced all 60 published TPR/PPV values." |
| S1_Table_README.txt | "±" = 1.96 × SEM for GO/EC/FAS/STD | VERIFICATION.md §5.6 checked GO/EC only | UNSUPPORTED (partly) | "1.96 × SEM for GO and EC (checked); FAS and STD use the same field" or check them. |

### 1.11 Citations (spot-check of 25 load-bearing keys; all 58 keys resolve)

| location | key / claim | status | correction |
|---|---|---|---|
| l.57 | majidian2025fastoma: FastOMA uses FAMSA | MISMATCH | FastOMA uses MAFFT + FastTree (F4) |
| l.36 | wapinski2007synergy: SYNERGY searches ancestral proteomes | MISMATCH | F3 |
| l.103 | emms2026orthofinder3: FastOMA fusion rate on OrthoBench | OK | optional "(10% of RefOGs fused vs 2.85%)" |
| l.274 | emms2017stride: OF roots gene trees with duplication information | OK (minor) | add emms2019orthofinder |
| l.29–33, 49, 55–57, 60, 76, 106, 171, 262–267, 311 | altenhoff2016qfo, altenhoff2013gethogs, sarton2025hogreview, cosentino2024sonicparanoid2, klemm2023proteinortho, altenhoff2024oma, rossier2021omamer, schreiber2013hieranoid, kaduk2017hieranoid2, buchfink2026deepclust, remm2001inparanoid, huertacepas2007overlap, schmitt2011orthoxml, train2019pyham, trachana2011orthobench, emms2020orthobench, altenhoff2024qfo, li2023miniprot, kallenborn2025mmseqs2gpu, heinzinger2024prostt5, douze2026faiss, piaggio2012consortnoninferiority | OK | metadata matches PubMed/Crossref; each supports its sentence |
| bib | lefort2015fastme title ends "...Program: Table 1." | MISMATCH | drop ": Table 1." |
| bib | yates2019ensembl: year 2019, no volume, pages gkz966 | MISMATCH (minor) | year 2020, volume 48, number D1, pages D682--D688 |

---

## 2. Required fixes (most important first; exact old → new text)

**F1. QfO primary comparators (body l.132–141, Abstract l.21, Author summary l.24, Discussion l.272; S1 §9.3).**
In body.tex, replace the bullet list and the sentence after it:
- Old: "A paired bootstrap over the 749 held-out TreeFam-A families gave these differences in recall for \texttt{accurate}:" + four bullets.
- New: keep the intro sentence and the four bullets, and add three bullets for the remaining pre-registered primary comparators:
  `\item vs Hieranoid 2: $-0.036$ \CI{$-0.069$}{$+0.004$}, not shown non-inferior;`
  `\item vs InParanoid 5: $-0.042$ \CI{$-0.081$}{$-0.000$}, lower, not non-inferior;`
  `\item vs SonicParanoid2 (sens-g): $-0.087$ \CI{$-0.132$}{$-0.038$}, inferior.`
  After the list, add: "OrthoFinder, FastOMA, Hieranoid 2, InParanoid 5 and SonicParanoid2 were the pre-registered primary comparators; OMA pairs and OMA HOGs were added after a user asked for an OMA resource in every comparison. On TreeFam-A precision, non-inferiority was shown only against OMA pairs and OMA HOGs."
- Abstract, old: "Ancestra matched FastOMA and OMA pairs on TreeFam-A recall but had clearly lower recall than OrthoFinder ($-0.16$)" → new: "Ancestra matched FastOMA and OMA pairs on TreeFam-A recall but had lower recall than OrthoFinder ($-0.16$), SonicParanoid2 and InParanoid"
- Discussion, old: "On QfO it matched the OMA methods but not OrthoFinder." → new: "On QfO its recall was non-inferior to FastOMA and OMA pairs, but below that of OrthoFinder, SonicParanoid2 and InParanoid."
- Author summary, old: "It was as good as or better than methods from the OMA family and cheaper than OrthoFinder, but it missed more distant orthologs than OrthoFinder." → new: "It scored higher than OMA-based methods on one benchmark and had similar recall on the other. It was cheaper than OrthoFinder, but it missed more distant orthologs than OrthoFinder and several pairwise methods."
- S1 §9.3: add the two rows given in §1.9.

**F2. Re-scored submissions (body l.106; S1_Table_README.txt).**
- Old: "To validate our local setup, we re-scored 20 published submissions and reproduced all 60 published values that we checked."
- New: "To validate our local setup, we re-scored 10 of the 20 published submissions on SwissTree, TreeFam-A and VGNC and reproduced all 60 published TPR/PPV values exactly."
- README: see §1.10.

**F3. SYNERGY and novelty (body l.36, l.38, l.64).**
- l.36 old: "SYNERGY did this for fungi with sequence and gene order~\cite{wapinski2007synergy}, and Hieranoid did it with an ancestral consensus proteome at each node~\cite{schreiber2013hieranoid,kaduk2017hieranoid2}. Both search ancestral proteomes against each other, which becomes expensive at deep nodes."
- l.36 new: "SYNERGY did this for fungi with sequence and gene order, from one all-against-all comparison of extant genes whose distances it updates as genes are merged~\cite{wapinski2007synergy}. Hieranoid builds an ancestral consensus proteome at each node and searches it against that of the sister clade~\cite{schreiber2013hieranoid,kaduk2017hieranoid2}. SYNERGY's all-against-all search grows quadratically, and Hieranoid's ancestral searches become expensive at deep nodes."
- l.38: after "aggregated from the scores of their member proteins" add ", as in SYNERGY~\cite{wapinski2007synergy}".
- l.64 old: "...without a reference database and without searching ancestral genes again." → new: "...without a reference database, without an all-against-all search and without searching ancestral genes again."

**F4. FastOMA alignment tool (body l.57).**
- Old: "...builds a tree with FastTree~\cite{price2010fasttree2}, as FastOMA and OrthoFinder~3 also do~\cite{majidian2025fastoma,emms2026orthofinder3}."
- New: "...builds a tree with FastTree~\cite{price2010fasttree2}, as OrthoFinder~3 also does~\cite{emms2026orthofinder3}; FastOMA uses MAFFT and FastTree~\cite{majidian2025fastoma}."

**F5. Human interventions (body l.305; S2 Text §2).**
- Old: "Humans intervened six times after the initial prompt (S2 Text): they suggested Foldcomp structure databases, required at least one OMA resource in every comparison, paused the project to move it to a node with a GPU, reset the deadline, asked for a header in the repository README, and requested this manuscript."
- New: "Humans intervened seven times after the initial prompt (S2 Text): they suggested Foldcomp structure databases, required at least one OMA resource in every comparison, paused the project to move it to a node with a GPU, kept the 48-hour deadline but made it soft, asked for a header in the repository README, asked for a second round of data layers, and requested this manuscript."
- S2 Text §2: add item 7 (see §1.10).

**F6. Broken reference (body l.231).** "Fig~\ref{fig:scaling}" → "Table~\ref{tab:scaling}".

**F7. GULO (body l.174).**
- Old: "25\% were confirmed losses (pseudogenes or deletions), including the known pseudogenes \textit{GULO} and \textit{UOX} in primates;"
- New: "25\% were confirmed losses (pseudogenes or deletions), including the known \textit{UOX} pseudogene in human (\textit{GULO}, which is missing in both human and chimpanzee, was confirmed as a pseudogene by a separate clade-level search);"

**F8. FastOMA diagnosis (body l.103).**
- Old: "The low score comes mainly from proteins that OMAmer did not place and from fused groups at the root, consistent with..."
- New: "The low score comes from low precision (23.9, from reference orthogroups fused into larger HOGs) and from the 20\% of proteins missing from FastOMA's groups (9.1\% never placed in a family and 9.7\% left as root-level singletons; OMAmer itself left only 4.1\% unplaced), consistent with..."

**F9. Overstated wording (body).**
- l.101: "favoured Ancestra by 6--7 points" → "favoured Ancestra by 6.5 and 5.7 points".
- l.139: "Ancestra had the highest G-STD2 discordance" → "Ancestra had the highest G-STD2 Eukaryota discordance of the methods in Table~\ref{tab:qfo}".
- l.141: "These are deep nodes such as Opisthokonta or Eukaryota, where there is often no outgroup" → "These are mostly deep nodes (Bilateria 34\%, Eukaryota 21\% and Opisthokonta 13\% of these misses), where an outgroup is often missing or distant".
- l.255: "grew close to linearly" → "grew slightly faster than linearly (exponent 1.16)"; l.272: "at a cost that grows close to linearly" → "at a cost that grows only slightly faster than linearly".
- l.272: "On OrthoBench it was the most accurate method we tested" → "On OrthoBench it had the highest F of the methods we tested (not significantly above OrthoFinder or SonicParanoid2)".
- l.276: "and structure made them worse" → "and structure lowered their precision".
- l.228: "\textit{C. elegans} outside Bilateria" → "\textit{C. elegans} placed outside the other bilaterians".
- l.310: append to the Splits item: "This split replaced a first, random split that a reviewer found leaky; three early S0 configurations had been scored by then (S1 Text, section 9.1)."
- tab:qfo: "SonicParanoid2 (sens.)" → "SonicParanoid2 (sens-g)".

**F10. State the configuration (tab:ob, tab:qfo, tab:scaling captions).**
- tab:ob caption, add: "Ancestra rows use the gene-order layer; without it, test F was 76.1 (\texttt{accurate}) and 74.8 (\texttt{fast})."
- tab:qfo caption, add: "Ancestra rows use the gene-order layer."
- tab:scaling footnote, add: "Ancestra runs are sequence only (no gene order)."
- tab:layers gene-order cell: "QfO TreeFam-A $-0.012$ ..." → "QfO TreeFam-A ($\texttt{fast}$) $-0.012$ ...".

**F11. Methods details (body l.293, l.288).**
- l.293 old: "Claude Opus 5.5 for design, algorithm and verification work, and Claude Sonnet 5.5 for data acquisition, benchmark runs and bookkeeping." → new: "Claude Opus 5.5 for literature, data acquisition, design, algorithm and verification work, and Claude Sonnet 5.5 for the software environment, competitor and scaling runs, bookkeeping and reference checks."
- l.288 old: "Only 22 of 81 quantitative claims were independently recomputed; the others were checked against their evidence files." → new: "Of the 81 claims in our claim ledger, which covers results up to the freeze, 22 were independently recomputed; the others were checked against their evidence files. Post-freeze GPU and second-round results were checked separately."

**F12. S1 Text.** The subtitle, §2.4, §5 (text and table), §6.2, §8.3, §9.1, §9.3 and the §4.7 heading, as in §1.9.

**F13. Bibliography.** `lefort2015fastme` title and `yates2019ensembl` metadata, as in §1.11. Optional: add `emms2019orthofinder` at l.274, and "proposed or used" at l.40.

**F14. Residual VERIFICATION_L2 items (low priority).**
- l.264: "with both differences within the CIs" → "(the CIs of both differences include 0; this was the best of seven arms)"; "It needed about 15 times fewer CPU-hours" → "It needed an estimated 15 times fewer CPU-hours (plus 0.2--0.3 GPU-hours)".
- l.265: "About 93\% of each FastTree call" → "About 93\% of a FastTree call (in a sample of 20 small jobs)".
- l.178: apply L2 #16 (post hoc arms and seed exceptions).

**F15. Housekeeping.** Regenerate `ancestra_manuscript.tex`/PDF and the S1/S2 PDFs after the fixes. In S1 Table's caption (l.366), "all 21 methods" → "all 20 published methods and both Ancestra presets".
