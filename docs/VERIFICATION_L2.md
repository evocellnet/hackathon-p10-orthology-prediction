# Independent verification of the L2 (second round) text

Verifier: independent agent, 2026-10-28. I did not produce these results. Scope:
REPORT.md §3.9 and §3.8 (from "The L2-gpu follow-up" to the end), manuscript/body.tex (paragraph
"Second round of layers", Table `tab:layers2`, GPU bullets "All families on the GPU" and "Local gene trees", the
1% jitter sentence in Limitations), manuscript/S1_Text.md §10.

Method: every number was compared with the JSON files in `results/l2_*/`. Where a JSON holds raw values the
number was recomputed (percentages, ratios, differences). `docs/L2_*.md` and `results/l2_plm/tables.md` were used
only as a pointer, and spot-checked against the JSON they are generated from (they agree). Nothing was edited
except this file.

Status codes: OK = matches the evidence; MISMATCH = number or label is wrong; UNSUPPORTED = no JSON or the data
do not support the claim as worded; OVERSTATED = number is right but the wording claims more than the data show.

Abbreviations: OB = OrthoBench dev (34 RefOGs); TF-A = TreeFam-A dev (359 families); "A" = CPU table without
cross-family rescue (reference of the GPU diagnostics); "frozen" = the frozen CPU path with rescue.

## 1. Mismatches and overstatements (each with its fix)

| # | Location | Problem | Evidence | Exact fix |
|---|---|---|---|---|
| 1 | REPORT §3.9 table header; body.tex `tab:layers2` caption ("AUC is for ortholog vs paralog on hard cases ... unless stated") | The network AUCs (0.630, 0.622), the abundance AUC (0.676) and the pLM "+0.23 to +0.33 above shuffled" are all-candidate / all-triplet values, not hard-case values. Nothing is "stated". Hard-case values are much weaker. | `l2_network/standalone.json` evaluation.own.auc.own_count: all 0.6303 [0.532, 0.6987], **hard 0.5013 [0.344, 0.6869]**. `standalone_treefam_dev.json`: all 0.6221, **hard 0.5687 [0.5009, 0.6253]**. `standalone.json` evaluation.abund: all 0.6763, hard 0.6887 [0.591, 0.7529]. `l2_plm/discrimination_final.json` T_seq/hard and hard_deep: Δ(real − shuffled) is −0.048 to +0.139, most CIs include 0 (e.g. e35_mean +0.002 [−0.068, +0.051]; pt5_mean +0.042 [−0.033, +0.091]). | Mark these cells "(all pairs)": network "AUC 0.630 [0.532, 0.699] on all evaluable OrthoBench pairs (hard pairs: 0.501 [0.344, 0.687]); 0.622 on TreeFam-A (hard: 0.569 [0.501, 0.625])"; abundance "AUC 0.676 [0.604, 0.732] (hard 0.689 [0.591, 0.753])"; pLM "above shuffled by +0.23 to +0.33 on all triplets (not on hard triplets)". Same change in the body.tex cells and caption. |
| 2 | REPORT §3.9 row "Contrastive head"; body.tex row "Contrastive model"; also S1 §10.4 | "+0.040 [+0.010, +0.056] over raw embeddings on easy cases": +0.040 is the all-triplet contrast (n = 17,750, hard triplets included). The easy (shallow) contrast is larger. | `l2_plm/discrimination_final.json` groups.T_seq/all.contrasts."head_e3b_s11:e3b_mean" = 0.0397 [0.0104, 0.0561]; T_seq/shallow = 0.0505 [0.0256, 0.074]; hard −0.0070 [−0.0649, 0.0556]. | "+0.040 [+0.010, +0.056] over raw embeddings on all triplets (+0.051 [+0.026, +0.074] on shallow/easy ones)". In S1 §10.4 write "improved all triplets (+0.040 ...) and easy ones (+0.051 ...)" instead of "easy cases (+0.040 ...)". Note it is seed 11 only. |
| 3 | REPORT §3.9 row pLM: "Where sequence fails on hard cases they are below chance (0.24–0.41)"; body.tex "below chance where sequence fails"; S1 §10.4 "scored 0.24–0.41, below chance" | The range leaves out ProtT5 alignment, which is at chance, not below. | `l2_plm/complement_final.json` hard (n_seq_wrong 1,944): e35_mean 0.293, e650_mean 0.351, e3b_mean 0.366, pt5_mean 0.239, e650_aln 0.411, e3b_aln 0.369, **pt5_aln 0.525** (its shuffled control 0.529). | "0.24–0.41 for six of the seven arms (ProtT5 alignment: 0.53, at chance)". In body.tex: "below chance for six of seven arms". |
| 4 | REPORT §3.9 "Main lesson" ("Every new layer carries real signal: each beats its shuffled control"); body.tex paragraph after the table ("Every layer beat its shuffled control") | Not true for the decision that matters in several layers. See rows 5a–5e below. | See rows 5a–5e. | Replace with: "Most layers carry real signal when all pairs are pooled (AUC 0.58–0.68 against 0.50 for shuffled controls), but on the hard cases, where sequence is ambiguous, several do not beat their control: networks (hard AUC 0.50), domain architecture within an orthogroup (0.54, and below a within-family shuffle), tissue specificity (45% against 36–39% shuffled, overlapping CIs), and the structure rescue (half of its gain is reproduced by shuffled 3Di)." |
| 5a | (basis of #4) Network | Hard AUC 0.501 [0.344, 0.687] against shuffled 0.4875–0.50. TF-A decisions: real precision 0.519 (27/52) against shuffled seed 23 0.667 (9 decisions). | `standalone.json` own.q.hard.own_count: shuf_precision_when_decides [0.0, 0.0, 0.0] with shuf_decides [0, 1, 0]; `standalone_treefam_dev.json`: shuf_precision_when_decides [0.2667, 0.6667, 0.3333]. | See #4. |
| 5b | (basis of #4) Domains | Within-orthogroup AUC is below the within-family shuffle, and TF-A is at 0.5. | `l2_domains/standalone.json` orthobench_dev.within_only.real aucJ_all 0.5416 vs shufFam 0.5677 / 0.579 / 0.5839; hard: 0.4882 [0.446, 0.5823] vs shuf_max_hard 0.4997, `gate_hard` false; treefamA_dev aucJ_all 0.4928. The 0.979 between-orthogroup value does beat shuffled (0.4976–0.5012). | State "0.54 within one orthogroup (hard pairs 0.49; no better than a within-family shuffle)". |
| 5c | (basis of #4) Tau | Real precision 0.4476 [0.3798, 0.554]; shuffled 0.394 [0.303, 0.548] (seed 11), 0.367, 0.360. CIs overlap, so "above shuffled" is not shown. | `l2_genome/expr_tau.json` real.hard.precision(_ci), shuf11.hard.precision(_ci). | Do not claim tau beats its control. |
| 5d | (basis of #4) Structure rescue | R4-full minus shuffled-3Di: Δmacro +0.0021 [0.0000, +0.0064]; the shuffled arms keep +0.0021 of the +0.0042. | `l2_gpu/t3_ob_R4full_vs_shuf11.json` dmacro_ci [0.0, 0.0064]; `t3_ob_R4.json` R4shuf11 dmacro 0.0021. | See #4. |
| 5e | (basis of #4) pLM on hard cases | See #1 (hard Δ vs shuffled mostly includes 0). | `discrimination_final.json` T_seq/hard, minus_shuf. | See #4. |
| 6 | REPORT §3.9 "Main lesson" ("Synteny ... repairs real fission errors"); body.tex "it fixed real errors"; body.tex row "removed 29–40% of known fission cases"; S1 §10.5 "fission cases found in S8" | The 29–40% is a drop in the count of S8-style fission *candidates* (a gene sits in a target-only HOG) on the OrthoBench accurate run (706 → 426–499). It is not a count of confirmed fission errors, and a wrong join also removes a candidate. The 706 is a recomputed OrthoBench figure, not the 808 errors of S8. The direct evidence that joins are correct is tiny: 7 labelled executed joins on OrthoBench (F_c 1 of 2, F_u 4 of 5), 24 on TF-A (79% / 88% speciation joins). Real support AUC over all labelled decisions on TF-A is 0.46–0.47. One RefOG (RefOG007, 5 groups → 1) is repaired fully at OG level in every real arm. | `l2_genome/keep_ob.json` off.s8 706; F_c.s8 499; F_u.s8 446; R_c.s8 476; R_u.s8 426; shuffled 706–709. `gsr_standalone.json` F_c executed_labelled 2 / executed_correct 1; F_u 5 / 4. `gsr_standalone_qfo.json` fast_R_u_log executed_labelled 24, executed_precision 0.7917, auc_s1 0.4688; fast_R_c_log 24, 0.875, 0.459. | Say: "reduced the number of S8-style fission candidates by 29–40% (706 → 426–499 on OrthoBench dev; shuffled gene order: none). Only 7 executed joins touch a labelled OrthoBench pair (5 correct), so whether the joins are correct is largely untested." Replace "repairs real fission errors" / "fixed real errors" by "removes candidate fission errors". In S1 replace "the fission cases found in S8" by "the S8-style fission candidates (706 on the OrthoBench accurate run)". In body.tex replace "known fission cases" by "fission candidates". |
| 7 | REPORT §3.9 row "Intron positions ... fails"; body.tex "tie-breaking lowered accuracy by 0.07"; S1 §10.5 "it lowered hard-case accuracy from 0.473 to 0.403" | The CI of the change includes 0 (−0.0705 [−0.1488, +0.0065]); "lowered" is a point estimate. The coded pre-registered stop rule returned PASS (hard-AUC CI lower bound 0.5137 > shuffled max 0.4986); integration was skipped under a "conservative reading" (docs/L2_genome.md §5.3). | `l2_genome/introns_standalone.json` hard_queries.real gain −0.0705, gain_ci [−0.1488, 0.0065]; stop_rule {pass_auc true, pass_gain false, verdict "PASS"}. | Write "tie-breaking changed accuracy by −0.071 [−0.149, +0.007]" and in the REPORT verdict cell "not integrated (stopped by a conservative reading of the stop rule; the coded rule returned PASS)". |
| 8 | S1 §10.3 "A bonus arm and a veto arm had no effect" | The bonus arm lowered OrthoBench OG-F by 3.0 and pair F1 by 0.0085, and failed the OG-F guard. The veto arm fired on 1,923 of about 596k QfO seeds and changed nothing on OrthoBench; QfO ΔF1 +0.0013 [0.0000, +0.0044]. | `l2_domains/integration.json` ob.bonus dOG_F −3.01, og_guard_pass false, dmacro −0.0085; ob.veto dOG_F 0.0; qfo.veto dF1_ci [0.0, 0.0044], stats.arch_veto 1923. | "A bonus arm (re-ranking) lowered OrthoBench OG-F by 3.0 and a veto arm almost never fired (1,923 of about 596k QfO seeds)." |
| 9 | S1 §10.6 "We ran eight arms, a random-sparsification control and a score-jitter floor" | There are 7 arms (K39, K10, R207, F0528, MK39, SP3000, SP5000) plus the C-rand control and the jitter floor. "Eight" counts the control. docs/L2_gpu.md repeats "8". | `l2_gpu/task1_summary.json` rows with kind = "arm": 7; kind "control": 1; kind "noise": 3 (J11/J23/J47). | "We ran seven arms, a random-sparsification control and a score-jitter floor (three seeds)." |
| 10 | REPORT §3.8 table row "GPU for the local trees" ("saves 1.5% of merge time"); body.tex "saved only 1.5% of merge time"; S1 §10.6 "saved 1.5% of merge time" | Merge wall went 1,898.0 → 1,867.0 s = −1.63%. The −1.5% is the FAMSA + FastTree share (1,846 → 1,818 s). It is one run each on a shared, loaded node. The 1,000-job micro-benchmark gives −2.8% (359.74 → 349.67 s). | `l2_gpu/t2_ob_batch_vs_ref.json` stages.merge.wall_s 1897.999 vs 1866.968; `t2_bench_ob1000.json` fasttree_single_s 359.74, fasttree_batch50_s 349.67. | "saved 1.6% of merge wall time (one run each; 2.8% in a 1,000-job micro-benchmark)". Conclusion "not worth enabling" is unaffected. |
| 11 | REPORT §3.8 (bullet 2), body.tex Limitations and GPU bullet, S1 §10.6: "A 1% random jitter ... costs −3.2 to −7.3 OG-F"; "sharp optimum"; "robustness limitation of the method"; "The jitter floor explains why" | Point estimates are right (−3.22, −7.30, −4.98 versus A; mean −5.17). But two of the three seed CIs include 0 (−7.30 [−15.42, +0.14]; −4.98 [−12.08, +0.85]). It is one benchmark (34 RefOGs; OG-F dominated by four RefOGs). The frozen parameters were tuned on dev, which alone predicts regression under any perturbation. On TF-A dev the 1% jitter moves TPR/PPV by −0.009 [−0.035, +0.013] / +0.001 only (seed 11). A 3-join shuffled-GSR perturbation also cost −3.65 OG-F (genome doc), so part of this is OG-F chaos through RefOG021. "Explains why" is a post hoc reading (D2 control added after seeing the arms). | `l2_gpu/t1_ob_jitter_vs_cpuA.json` J11 dOG_F −3.22 [−8.37, −0.04], J23 −7.3 [−15.42, 0.14]; `t1_ob_jitter47_vs_cpuA.json` J47 −4.98 [−12.08, 0.85]; macro −0.0106 [−0.0197, −0.0032], −0.0136, −0.0098 [−0.022, 0.0001]. `t1_qfo_vs_cpuA.json` J11 dTPR −0.0088 [−0.0352, 0.0128], dPPV 0.001. `l2_genome/keep_ob.json` F_c.dOG_F_shuf [−3.64, 0, 0]. | Limitations sentence: "a 1% random jitter of the DIAMOND scores lowered OrthoBench development OG-F by 3.2 to 7.3 points (mean 5.2; two of three seed-level intervals include 0) and macro pair F1 by 0.010 to 0.014, whereas TreeFam-A development TPR/PPV moved by −0.009/+0.001 (one seed). The frozen parameters were tuned on this benchmark, so the merge may sit on a dev-specific optimum." Replace "explains why" by "shows how large a perturbation floor is". |
| 12 | body.tex GPU bullet: "The merge averages scores over the pairs that are present, so extra or rescaled edges change which units look closest, and ambiguous clusters are then formed differently." | Stated as the general mechanism, but it rests on one family trace (RefOG021). The max-link control, which should remove the dilution, did not shrink the CPU/GPU gap: avg-link GPU − CPU −4.38 [−11.0, −0.2]; max-link GPU − CPU +4.86 [−7.42, +15.11]. Swapping only GPU self-scores costs −17 OG-F (score values, not presence, matter). | `l2_gpu/t1_ob_diag_maxlink.json` results[1].dOG_F 4.86; `t1_ob_diag_avg_ref_cpuA.json` results[1].dOG_F −4.38; `t1_ob_diag_X3_vs_cpuA.json` −16.97; `trace_fam77_*.json`. | "In a trace of the largest orthogroup (RefOG021), the merge's average over present pairs let a unit leave its cluster; a max-link control did not settle whether this is the general mechanism." |
| 13 | REPORT §3.9 domain row "dominated by a sequence-only setting"; S1 §10.3 "An existing sequence-only arm dominated it" | The comparison is across bases. The gate arm is on the QfO `fast` base; the dominating arm (C2: accurate + so_struct + deep absorption, .6548/.9136) is on the `accurate` base. The doc says so. | `docs/L2_domains.md` §4 finding 1; arm values in `l2_domains/integration.json` qfo.gate TPR_official 0.6534, PPV_official 0.8974. | Add "(different base: accurate versus fast)" or drop "dominated". |
| 14 | REPORT §3.8 table "≈93% of each FastTree call is ML optimisation"; body.tex "About 93% of each FastTree call"; S1 "93% of each call" | No JSON holds it. The doc quotes 20 five-sequence jobs (0.058 s with ML, 0.004 s with -noml; 1 − 0.004/0.058 = 93.1%). The 1,000-job benchmark averages 0.36 s per call, so those 20 jobs are not representative of the mean call. | `l2_gpu/t2_bench_ob1000.json` has no ML/no-ML split. | "about 93% in a sample of 20 small jobs". |
| 15 | REPORT §3.8 bullet 3, body.tex, S1: "0.59 CPU-h" for R207 | Not a measured run. It is the measured unsparsified all-GPU 0.55 CPU-h plus the post hoc sparsification time (134–142 s single thread), about 0.59; the wall time is marked "≈ 2,220 s". It also excludes 0.18–0.33 GPU-h. | `docs/L2_gpu.md` §2.3; `docs/GPU_v2.md` line 406 (0.55). | "about 0.59 CPU-h (estimate: measured 0.55 plus sparsification time), plus 0.2–0.3 GPU-h". |
| 16 | REPORT §3.9 intro, body.tex ("Each layer's arms and keep rule were committed before its integration results, and each was compared with shuffled controls"), S1 §10 intro ("shuffled controls used seeds 11/23/47") | Overgeneral. Post hoc arms were added and logged as deviations: pLM D2 heads, GPU D1 arms (MK39, SP3000, SP5000) and the D2 jitter control, domains post hoc ungated arms, network TF-A standalone extension, genome idea 4 (committed after ideas 1–3 were seen). Shuffled controls: genome F arms on TF-A not run; GPU structure rescue on QfO seed 11 only; QfO jitter seed 11 only; the sparsification arms have C-rand and jitter, not shuffled-score controls. | `docs/L2_plm.md` §3 D2; `docs/L2_gpu.md` §0.4 D1–D3; `docs/L2_genome.md` §1.7; `docs/L2_network.md` §8; `l2_genome/keep_qfo.json` fast_F_c d_real_minus_shuf null; `docs/L2_gpu.md` D3. | Add: "apart from the deviations logged in each doc (post hoc arms: pLM D2, GPU D1/D2, domain ungated arms, network TreeFam-A standalone)" and "shuffled controls used seeds 11/23/47 except where noted (QfO structure rescue and QfO jitter: seed 11 only; genome F arms on TreeFam-A: none)". |
| 17 | REPORT §3.9 row network: "right when it decides 44% / 52%" (and S1 §10.1, S1 §10.2 "34%") | Numbers are right, but the OrthoBench 44% is 4 of 9 decisions (17 hard queries); no CI is given. TF-A 27/52; abundance 39/116. | `standalone.json` own.q.hard.own_count syn_correct 4, syn_wrong 5; `standalone_treefam_dev.json` 27/25; abund 39/77. | Add the counts: "44% (4 of 9 decisions)". |
| 18 | REPORT §3.8 bullet 1 vs table row "All families on the GPU" | Mixed baselines in one paragraph. −9.1 (X1) and −17 (X3) are against A (CPU without cross-family rescue); −5.3 [−12.0, −0.7] is against the frozen CPU path with rescue. Against A the all-GPU graph is −4.38 [−10.99, −0.23], so X1 (−9.1) is worse than the full GPU graph. Jitter −3.2 to −7.3 is also against A (against frozen: −4.2 to −8.2). | `l2_gpu/task1_summary.json` rows ob_vs_cpuA; `t1_ob_diag_X_vs_cpuA.json` results[1].dOG_F −9.1 [−14.32, −4.04]. | Add "(versus the CPU path without rescue)" after −9.1, −17 and the jitter range. |
| 19 | REPORT §3.8 bullet 3, body.tex: "with both differences inside the CIs" | Ambiguous wording: it means the CIs of the differences include 0. R207 ΔTPR +0.0215 [−0.0029, +0.0503], ΔPPV +0.0075 [−0.012, +0.0303]. R207 was the best of 7 arms on this set, so the point estimates are selected maxima. | `l2_gpu/task1_summary.json` R207 qfo. | "(the CIs of both differences include 0; R207 was the best of seven arms on this set)". |
| 20 | body.tex Discussion: networks "failed for the same reason: the evidence ... is also shared by close paralogs" | For networks the main failure is coverage (3.4% of OrthoBench candidates, 1.0% of hard pairs; oracle ceiling +0.0016 micro pair F1). Paralog sharing is shown only on TF-A (precision 0.52, rescue precision 0.48–0.54). | `l2_network/standalone.json` coverage.own frac_hard_pairs_evaluable 0.01; `error_coverage.json` own.oracle_upper_bound delta 0.0016. | Add "or, for networks, almost no coverage". |

## 2. Claim-by-claim check

### 2.1 REPORT.md §3.9

| Location | Claim | Evidence (file: key) | Status | Correct value / note |
|---|---|---|---|---|
| §3.9 | Five agents, shared protocol, dev only, same keep rule (CI lower bound > 0; PPV −0.01; OG-F −1.0) | `docs/L2_protocol.md` rules 4, 6 | OK | See #16 for "committed before" and seed wording. |
| §3.9 | "No arm passed" | `l2_domains/integration.json` verdict.*.keep_before_orthoxml_check false for all four; `l2_genome/keep_ob.json` keep_rule_ob false; `keep_qfo.json` keep_rule_tf false; `l2_gpu/task1_summary.json` E1/E2 fail for every all-GPU arm | OK | True for integration. For intron the coded stop rule returned PASS (#7). Domain `gate` passes the primary CI and fails two guards. |
| Network | AUC 0.630 [0.532, 0.699], 3.4% covered | `l2_network/standalone.json` evaluation.own: 0.6303 [0.532, 0.6987]; coverage.own.frac_candidates_evaluable 0.0339 | OK | All-pairs AUC (#1). |
| Network | AUC 0.622 on TreeFam-A | `standalone_treefam_dev.json` evaluation.own.auc.own_count.all 0.6221 [0.559, 0.6633] | OK | All-pairs (#1). |
| Network | right when it decides 44% / 52% | OB 4/9 = 0.4444; TF-A 27/52 = 0.5192 (`syn_precision_when_decides`) | OK | Add counts (#17). |
| Network | "As recall rescue the added pairs were orthologs only 48–54% of the time" | `treefam_rescue_precision.json` own_count>=1 0.479, >=2 0.518, >=3 0.534, >=5 0.538 | OK | Shuffled rescues are tiny (0–63 pairs), so no control for the precision itself. |
| Network | cost 0.4 CPU-h | `standalone*.json` wall_s 194.7 + 100.3 + 580.8 (+313 s mapping per doc) = about 0.33 h wall | OK (approx.) | Doc says ≈ 0.4; no CPU-time key in the JSON. |
| Abundance | AUC 0.676 [0.604, 0.732]; right 34% | `standalone.json` evaluation.abund: 0.6763 [0.6039, 0.7319]; hard syn_precision 0.3362 | OK | All-pairs AUC (#1). Real minus shuffled comb@0.9 CIs include 0 (+0.082 [−0.030, +0.194]) as the doc says. |
| Domains | hard AUC 0.637 [0.547, 0.738] vs sequence 0.743 | `l2_domains/standalone.json` orthobench_dev.all_pairs.real aucJ_hard 0.6368, aucJ_hard_ci [0.547, 0.7384], aucSeq_hard 0.743 | OK | Within-family shuffle reaches 0.556–0.581 (#5b). |
| Domains | between orthogroups 0.979 | family_fusion.aucJ_sameRefOG 0.9787 [0.8392, 0.9949], aucSeq 0.7083, shuf 0.4976–0.5012 | OK | |
| Domains | within one orthogroup 0.54 | within_only.real aucJ_all 0.5416 (hard 0.4882) | OK | Hard value is 0.49 (#5b). |
| Domains | re-ranking +0.001 [−0.008, +0.011] | all_pairs.real top1_delta 0.001, top1_delta_ci [−0.0082, 0.0111] | OK | |
| Domains | best arm TF-A F1 +0.018 [+0.003, +0.034], TPR +0.038, PPV −0.025 | `integration.json` qfo.gate dF1 0.0182 [0.0032, 0.0336]; dTPR 0.0383; dPPV −0.025 [−0.0458, −0.0062] | OK | |
| Domains | OrthoBench OG-F −2.9 | ob.gate dOG_F −2.86 (79.69 → 76.82) | OK | Rounds to −2.9. |
| Domains | cost 48 CPU-h, annotation alone 31 | doc §6: 4.9 + 26.4 = 31.3; 31.3 + 16.5 + 0.02 = 47.8 | OK | No JSON holds CPU-h. |
| pLM | above shuffled by +0.23 to +0.33 | `discrimination_final.json` T_seq/all minus_shuf over the 7 pretrained arms: +0.228 (pt5_aln) to +0.328 (e3b_aln) | OK | On all triplets (#1). |
| pLM | deep 0.61–0.72 vs sequence 0.757 | T_seq/deep acc: 0.608 (e35_mean) to 0.720 (e650_aln); seq 0.7567 | OK | "below sequence everywhere" is loose: on hard triplets e650_aln (+0.001) and pt5_aln (+0.003) are equal to sequence, hard∩deep +0.022 / +0.035 (CIs include 0). The doc itself says "reach sequence level on hard". |
| pLM | below chance 0.24–0.41 where sequence fails | `complement_final.json` hard | MISMATCH | See #3. |
| pLM | best tie-breaker +0.022 [−0.001, +0.046] | T_seq/hard_deep arms.comb_pt5_aln.minus_seq 0.0224 [−0.0007, 0.0462] | OK | The CI includes 0, correctly shown as not passing. S1 "one of about 40 comparisons": doc §4.1. |
| Head | easy cases +0.040 [+0.010, +0.056] | T_seq/all contrast 0.0397 [0.0104, 0.0561] | MISMATCH | See #2. |
| Head | hard −0.007 [−0.065, +0.056]; deep −0.077 | T_seq/hard contrast −0.0070 [−0.0649, 0.0556]; head_e3b_s11 minus seq, deep −0.077 [−0.122, −0.031] (`tables.md`) | OK | |
| pLM cost | 3 GPU-h | `docs/L2_plm.md` §5 | OK | Wall on a shared device. |
| kNN rescue | guard passed 90% of links; 229k of 251k genes chained | `l2_plm/knn_guard.json` real: n_accepted 337,831 / n_mutual_cross 376,238 = 89.8%; largest_family 229,289 | OK | 251,378 genes (doc). The shuffled controls do the same (largest 238,716–238,749). |
| Genome | removes 29–40% of S8 fission candidates (shuffled none) | `keep_ob.json` s8: 706 → 499/446/476/426 = −29.3% to −39.7%; shuffled 706–709 | OK (number) | Wording: #6. |
| Genome | OrthoBench macro pair F1 +0.0059 [−0.0015, +0.0159]; TF-A TPR +0.0042 [−0.0049, +0.0134]; PPV within 0.002 | `keep_ob.json` R_u dmacro_ci; `keep_qfo.json` fast_R_u dTPR 0.0042 [−0.0049, 0.0134]; dPPV −0.002 | OK | These are the best arm (R_u). PPV deltas over arms: −0.0020 to +0.0021. R_u OG-F is +3.2 (not mentioned). One shuffled seed alone moves TF-A TPR by +0.0097 (`keep_qfo.json` fast_R_u shuffled_TPR 0.6248). |
| Genome | fewer than 1% of its decisions touch a labelled pair | `gsr_standalone.json` coverage 0.0078 (F_c), 0.0055 (F_u); TF-A R_u 235/66,229 | OK | |
| Genome | cost ~4 CPU-h | doc §7 table: 2.7 + 1.0 + 0.03 + 0.22 + < 0.05 | OK | |
| Intron | AUC 0.578 [0.514, 0.661] vs sequence 0.654 | `introns_standalone.json` auc.hard real 0.5777, CI [0.5137, 0.6614], seq_auc 0.6536 | OK | |
| Intron | accuracy −0.071 [−0.149, +0.007] | hard_queries.real gain −0.0705, gain_ci [−0.1488, 0.0065] | OK (number) | Wording: #7. |
| Tau | right when it decides 45%; −0.076 [−0.131, −0.010] | `expr_tau.json` real.hard precision 0.4476 (128/286); gain −0.0755 [−0.1308, −0.0096] | OK (numbers) | #5c. |
| Structure | OrthoBench macro +0.0042 [−0.0005, +0.011], identical to AFDB-only; TF-A +0.002 TPR | `l2_gpu/t3_ob_R4.json` R4full = R4afdb dmacro 0.0042 [−0.0005, 0.011]; `t3_qfo_R4.json` R4full dTPR 0.0021 | OK | |
| Structure | ProstT5 5 GPU-h | doc §4: 5 h 01 min wall on a shared A100 | OK (approx.) | Wall time on a shared device. No JSON. |
| All-GPU sparsified | TF-A 0.637/0.930 vs 0.615/0.922; 0.59 vs 9.05 CPU-h; OG-F −8.4 | `task1_summary.json` R207 qfo 0.6366/0.93; frozen 0.6151/0.9224; ob dOG_F −8.40 | OK | #15, #19. 9.05/0.59 = 15.3. |
| Main lesson | all layers beat shuffled | | OVERSTATED | #4. |
| Second lesson | jitter costs 3–7 dev OG-F | | OK (number) | #11. |

### 2.2 REPORT.md §3.8 (from "The L2-gpu follow-up")

| Claim | Evidence | Status | Correct value / note |
|---|---|---|---|
| GPU scores on exactly DIAMOND's pair set cost −9.1 OG-F | `t1_ob_diag_X_vs_cpuA.json` results[1].dOG_F −9.1 [−14.32, −4.04] | OK | Versus A (#18). Versus frozen −10.04. |
| Swapping in only the GPU self-scores costs −17 | `t1_ob_diag_X3_vs_cpuA.json` dOG_F −16.97 [−30.06, −4.61] | OK | Versus A (#18). |
| 1% jitter costs −3.2 to −7.3 (seeds 11/23/47) | `t1_ob_jitter_vs_cpuA.json`, `t1_ob_jitter47_vs_cpuA.json` | OK (numbers) | #11. |
| "as much as the full GPU graph" | all-GPU −4.38 vs A (−5.32 vs frozen); jitter mean −5.17 | OK | |
| Rank-ratio keep if bits ≥ 0.207 × query's best hit in that species | `density_match_ob.json` alpha 0.20736 | OK | Density-matched label-free on OrthoBench, not tuned on accuracy. |
| QfO 0.637/0.930 vs 0.615/0.922 | `task1_summary.json`; `t1_qfo_vs_cpuA.json` R207 0.6366/0.93 | OK | |
| 0.59 vs 9.05 CPU-h | | OK (estimate) | #15. |
| −8.4 OrthoBench dev OG-F | `task1_summary.json` R207 dOG_F −8.4 [−14.5, −2.4] | OK | Versus frozen. |
| Table: "All families on the GPU ... OG-F −5.3 [−12.0, −0.7]; QfO PPV −0.016" | `task1_summary.json` all-GPU dOG_F −5.32, ci [−11.95, −0.66]; dPPV −0.0163 | OK | −11.95 rounds to −12.0 (GPU_v2.md) or −11.9 (L2_gpu.md); both are fine. |
| Table: "QfO 16× fewer CPU-h" | `docs/GPU_v2.md` line 406–407: 0.55 vs 9.05 CPU-h = 16.5× | OK | |
| Table: "≈93% ... not start-up" | | UNSUPPORTED (small sample) | #14. |
| Table: "batching saves 1.5% of merge time" | | MISMATCH | #10 (1.6% of merge wall). |
| Table: "4 of 1.39 M pairs differ; dev metrics identical" | `t2_ob_batch_vs_ref.json` dOG_F 0.0, dmacro 0.0, macro 0.8026 both; pairs 1,394,142 vs 1,394,138 (doc) | OK | OrthoBench only; no QfO run of the batch arm. |
| "sharp optimum ... robustness limitation of the method" | | OVERSTATED | #11. |

### 2.3 manuscript/body.tex

| Location | Claim | Evidence | Status | Note |
|---|---|---|---|---|
| Paragraph | "No arm passed, so none was run on test data" | see §2.1 | OK | |
| `tab:layers2` | Network 0.63 / 0.62, 3.4%, 48–54% | see §2.1 | OK | Caption mislabels (#1). |
| `tab:layers2` | Abundance 0.68; 34% | see §2.1 | OK | #1. |
| `tab:layers2` | Domains 0.64 (seq 0.74); 0.98 / 0.54; +0.018 [+0.003, +0.034]; PPV −0.025; OG-F −2.9 | see §2.1 | OK | |
| `tab:layers2` | pLM "above shuffled, below sequence; below chance where sequence fails" | | OVERSTATED | #1, #3. |
| `tab:layers2` | Contrastive +0.040 "on easy cases"; −0.007 [−0.065, +0.056] | | MISMATCH | #2. |
| `tab:layers2` | Intron AUC 0.58 (seq 0.65); "lowered accuracy by 0.07" | | OK (numbers) / OVERSTATED (wording) | #7. |
| `tab:layers2` | Gene order: "29–40% of known fission cases (shuffled: none)"; macro +0.006 [−0.002, +0.016]; TPR +0.004 [−0.005, +0.013] | | OK (numbers) / OVERSTATED (wording) | #6. |
| `tab:layers2` | Structure rescue macro +0.004 [−0.001, +0.011] | | OK | |
| After table | "Every layer beat its shuffled control" | | OVERSTATED | #4. |
| After table | "it fixed real errors; fewer than 1% of its decisions touched a labelled pair" | | OVERSTATED / OK | #6. The 1% figure is OK. |
| GPU bullet | 16-fold; 0.55 vs 9.05; OG-F −5.3 [−12.0, −0.7]; PPV −0.016; 9.1 OG-F | | OK | #18. |
| GPU bullet | mechanism sentence; "A control explains why no fix met the equivalence rule" | | OVERSTATED | #11, #12. |
| GPU bullet | jitter 3.2–7.3 (three seeds) | | OK (number) | #11. |
| GPU bullet | R207 0.637/0.930 vs 0.615/0.922; 15 times fewer CPU-h; lost 8.4 OG-F | | OK | #15, #19. |
| Local trees bullet | ~93%; 1.5% | | UNSUPPORTED / MISMATCH | #14, #10. |
| Limitations | "1% jitter ... lowered development OG-F by 3–7 points"; "sharp optimum" | | OK (number) / OVERSTATED (inference) | #11. |

### 2.4 manuscript/S1_Text.md §10

| Location | Claim | Evidence | Status | Note |
|---|---|---|---|---|
| §10 intro | shared protocol, controls, no test data | | OK / loose | #16. |
| §10.1 | AUC 0.630 [0.532, 0.699] and 0.622 vs 0.50 shuffled; 3.4%; right 44% / 52%; rescue 48–54% | see §2.1 | OK | #1, #17. |
| §10.1 | "Most hard cases lie in species with almost no interaction data (fish, chicken, dog, opossum, Ciona)" | `standalone.json` coverage.own.frac_hard_pairs_evaluable 0.01 (1.0%) | OK | |
| §10.2 | AUC 0.676 [0.604, 0.732]; right 34%; lowered accuracy −0.131 [−0.246, −0.046] | `standalone.json` abund q.hard comb0.9_minus_seq −0.1311 [−0.2462, −0.0463] | OK | The CI excludes 0 (stated correctly). Real − shuffled CIs include 0. |
| §10.3 | annotated 58k OrthoBench and 446k TF-A proteins | `docs/L2_domains.md` §1: 58,130; 3,581 + 442,858 = 446,439 | OK | Doc-level only. |
| §10.3 | hard AUC 0.637 [0.547, 0.738] vs 0.743; 0.979; 0.54; +0.001 | see §2.1 | OK | |
| §10.3 | best arm: tolerance 0.8 for identical architectures; F1 +0.018 [+0.003, +0.034], TPR +0.038, PPV −0.025, OG-F −2.9; beat its shuffled control | `integration.json` qfo.gate real_minus_shuf 0.0194 vs seed_spread 0.01 | OK | "beat shuffled" holds on F1 (0.0194 > 0.0100). |
| §10.3 | "dominated by an existing sequence-only arm" | | OVERSTATED | #13. |
| §10.3 | "bonus and veto no effect" | | MISMATCH | #8. |
| §10.3 | 31 CPU-h annotation; 60 CPU-h for all QfO | doc §6 (0.22 CPU-s/protein × 988,778 = 60 CPU-h) | OK | |
| §10.4 | beat shuffled by 0.23–0.33; deep 0.61–0.72 vs 0.757; 0.24–0.41; +0.022 [−0.001, +0.046] | see §2.1 | OK / MISMATCH (#3) | 0.24–0.41 omits ProtT5 alignment. |
| §10.4 | head: "easy cases (+0.040 ...)"; hard −0.007 | | MISMATCH | #2. |
| §10.4 | "A post-hoc variant ... did not beat its random-winner control" | `discrimination_final.json` D2_e3b_s11:D2_e3b_ctrl: all +0.013 [−0.015, +0.043]; hard +0.010 [−0.050, +0.095]; shallow +0.027 [+0.016, +0.051] | OK (mostly) | On shallow triplets all three D2 seeds exceed the control with CIs above 0; not on all, deep or hard. |
| §10.4 | 229k of 251k genes chained | `knn_guard.json` | OK | |
| §10.5 | removed 29–40% of the fission cases found in S8 | | OVERSTATED | #6. |
| §10.5 | macro +0.0059 [−0.0015, +0.0159]; TPR +0.0042 [−0.0049, +0.0134]; PPV within 0.002; fewer than 1% | `keep_ob.json`, `keep_qfo.json` | OK | |
| §10.5 | intron AUC 0.578 [0.514, 0.661] vs 0.50 and 0.654; accuracy 0.473 → 0.403 | `introns_standalone.json` | OK (numbers) | #7. |
| §10.5 | tau right 45%; −0.076 [−0.131, −0.010] | `expr_tau.json` | OK | #5c. |
| §10.6 | "eight arms" | | MISMATCH | #9. |
| §10.6 | jitter 3.2–7.3; R207 0.637/0.930 at 0.59 vs 0.615/0.922 at 9.05; −8.4 | | OK | #11, #15. |
| §10.6 | FastTree -n 1.5%; 4 of 1.39 M pairs; 93% | | MISMATCH / UNSUPPORTED | #10, #14. |
| §10.6 | pyfamsa changed 8 of 998 alignments | `t2_bench_ob1000.json` inproc_identical_alignments 990, inproc_rejected 2 (1,000 jobs) | OK | |
| §10.6 | ProstT5 filled 63k QfO proteins in 5 h; OrthoBench macro +0.0042 [−0.0005, +0.011] | doc §4 (63,237; 5 h 01 min); `t3_ob_R4.json` | OK | |
| §10.6 | Foldseek GPU 47.6 vs 29.0 s for 3,000 queries; different hit set | `t3_foldseek_gpu_bench_qfo.json` gpu_search_s 47.64, cpu_search_s 28.98, hits 4,743 vs 8,602, jaccard 0.3587 | OK | |

## 3. Summary

All headline numbers in the four texts trace to the JSON files; I found no transcription error larger than a
rounding step (the all-GPU CI lower bound is −11.95 in the JSON, printed as −12.0 in the REPORT and −11.9 in
docs/L2_gpu.md; both are acceptable).

Four items are wrong as written and need a one-line fix:
1. Table label: network, abundance and pLM numbers are all-pair values, but the table header and caption say
   "hard cases" (#1).
2. The head's "+0.040 on easy cases" is the all-triplet contrast; the easy-case value is +0.051 (#2).
3. "Below chance (0.24–0.41)" leaves out ProtT5 alignment, which is at chance (0.525) (#3).
4. "Eight arms" in S1 §10.6 should be seven arms plus a control (#9), and "1.5% of merge time" should be 1.6% of
   merge wall time (#10).

The main wording problems, which overstate the evidence:
- "Every new layer ... beats its shuffled control" is not supported on hard cases for networks, domain
  architecture within an orthogroup, tissue specificity and the structure rescue (#4, #5).
- "Synteny repairs real fission errors" rests on a drop in S8-style candidate counts and on 7 labelled executed
  joins on OrthoBench (#6).
- The "sharp optimum" and "explains why" statements about the 1% jitter come from one benchmark, two of three seed
  CIs include 0, TF-A dev barely moves, and the frozen parameters were tuned on that benchmark (#11, #12).
- Intron tie-breaking "lowered accuracy" (CI includes 0) and "fails" (the coded stop rule returned PASS) (#7).
- S1 says the bonus and veto arms had "no effect"; the bonus arm cost 3.0 OG-F (#8).

Not checked (no JSON or outside scope): CPU-h and GPU-h cost entries (doc-level only), the first-round GPU rows of
the §3.8 table, ESM-2/FAISS row, and literature references.
