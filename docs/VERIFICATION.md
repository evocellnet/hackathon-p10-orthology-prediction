# Independent verification (post-freeze, freeze-v1)

Verifier: independent agent (produced none of the results), 2026-10-27, T+17h29 to T+17h45 (scripts/clock.sh).
Scope: OrthoBench final test, reproducibility, toy end-to-end, 20-HOG spot check, and one recomputed number for
each of S0, S1, S4, S5 and S8. All verifier code and outputs are in `results/verification/`. Nothing in `ancestra/`,
`CLAIMS.md` or other agents' results was modified. Threads ≤ 6 throughout.

## 1. Claim table

| claim (CLAIMS id) | reported | recomputed | match? | how / notes |
|---|---|---|---|---|
| OB **test** F, Ancestra accurate (final) | 76.2 (P 69.5, R 84.2, 10/36 exact) | 76.19 (P 69.54, R 84.23, 10 exact) | yes | `ob_rescore.py`: official `benchmark.py` readers + `calculate_benchmarks_pairwise` called on the 36 test RefOG lists; does not import `bench_repo/scripts/orthobench_eval.py`. Input: raw `work/ob_final/main_accurate_go/orthobench_ogs.txt` |
| OB dev / all-70 F, accurate | 79.7 / 77.9 | 79.69 / 77.90 | yes | Official CLI on the raw file, all 70: 77.9 F, 72.5 P, 84.1 R, 21 exact (`official_cli_all70.txt`) |
| OB test F, fast | 74.8 (P 68.7, R 82.1, 9/36) | 74.79 (68.68 / 82.09, 9) | yes | dev 75.71, all 75.24 (CLI 75.2) |
| OB test F, OrthoFinder 3.2 | 69.7 | 69.70 | yes | from raw `Orthogroups.og`; CLI on the raw OrthoFinder `Orthogroups.tsv`: all-70 70.1 / 19 exact = bench_repo |
| OB test F, SonicParanoid2 / ProteinOrtho | 70.4 / 47.5 | 70.45 / 47.49 | yes | raw `.og` files |
| ΔF test accurate − OrthoFinder | +6.5 [−6.5, +18.1], SE 6.4, NI not shown | +6.48 [−6.49, +18.05], SE 6.38 | yes | own paired bootstrap (B=2000, seed 11). Robustness: seeds 23/47/1234 with B=10,000 give [−5.8…−6.0, +18.6]; leave-one-RefOG-out Δ ranges +3.2…+9.5; per-RefOG 15 wins / 12 ties / 9 losses. **Verdict "not shown non-inferior" is robust** |
| ΔF test accurate − SonicParanoid2 | +5.7 [−4.3, +14.9] | +5.74 [−4.29, +14.86] | yes | other seeds: lower bound −3.8 to −4.0 (still below −3) |
| ΔF test accurate − ProteinOrtho | +28.7 [+13.1, +41.9] | +28.70 [+13.12, +41.94] | yes | |
| ΔF test fast − OF / SP2 / PO | +5.1 [−7.5,+16.6] / +4.3 [−6.0,+14.1] / +27.3 [+12.5,+40.1] | +5.09 [−7.49,+16.61] / +4.34 [−6.01,+14.07] / +27.30 [+12.53,+40.12] | yes | |
| Test-split IDs | 36 test RefOGs | per_refog.tsv test labels == `splits/orthobench_split.json` test list; dev ∩ test = ∅ | yes | `sha256sum -c splits/SHA256SUMS`: all 3 OK |
| No OB test score before freeze except S0 completeness | claimed | none found | yes (one caveat, §2.2) | `git grep` of the freeze-v1 tree and `git log -p` of the history for all 36 test RefOG IDs; scan of untracked pre-freeze files in results/ work/ logs/. Only `results/S0_families/2sA_dc30_cc_sh30_RA.TEST.json` (logged). `orthobench_dev_eval.py` refuses non-dev subsets. Final per-config scores were written 08:05, after the 06:40 test_log entries |
| freeze.json hashes == working tree | 24 files | 24/24 identical; `git diff freeze-v1 HEAD -- ancestra/` empty | yes | see §2.3 for a commit-id inconsistency |
| pytest (not slow) | 156 tests (tree) | **156 passed, 1 deselected**, 86 s | yes | `pytest_not_slow.log` |
| Toy e2e (C-CORE-1, pair F1 0.94–0.98) | 0.94–0.98 | seed 23, 8 sp × 300 genes: accurate 0.961 (P 0.993, R 0.931), fast 0.963; with gene order 0.961 / 0.963 | yes | `toy_e2e.json`. hogs.orthoxml **byte-identical** across 2 runs per preset (with and without gene order) and across 1/4/6 threads |
| OB output validity (C-DOC-2 analogue) | PASS | `verify_orthoxml.py` on main_accurate_go: all 8 checks PASS (0 TaxRange≠LCA, 0 nesting violations, QfO flattening identical on 308,802 sampled pairs, pyham full load 15,246 top-level HOGs) | yes | own check: top-level groups (13,669 orthologGroups + 1,632 paralogGroups) re-derived from the OrthoXML == the 15,301 multi-gene OGs in `orthobench_ogs.txt` exactly; 20-HOG subset loads in pyham **without** the compat rewrite (20/20) |
| C-S0-2 S0 dev completeness | 0.794 (27/34) | 0.7941 (27/34) | yes | `s0_completeness.py`, own code on `work/S0/2sA_dc30_cc_sh30_RA/families.tsv` (sha256 = manifest). Keep-uncertain: 0.765 |
| C-S0-4 S0 test completeness | 0.722 | 0.7222 (26/36) | yes | |
| C-S5-2 synteny AUC human/mouse w=20 | 0.92 | 0.9223 | yes | `s5_auc.py`, Mann–Whitney AUC from `work/S5/pair_scores_*.parquet`; also human/chicken 0.928, zebrafish/Tetraodon 0.817, human/zebrafish 0.696, human/fly 0.488, shuffled (seed 11) 0.48–0.51; all RefOGs in dev. Labels themselves not re-derived |
| C-S8-2 hidden-gene recovery | 92.2% [89.7, 94.1] | 553/600 = 92.17%, Wilson [89.7, 94.1]; 0/512 random windows intact | yes | `s8_controls.json`; per species 95/96/98/90/89/85 as doc |
| S4 label accuracy NJ vs FastTree (c3di, matched) | 0.928 vs 0.955, Δ +0.027 [+0.004, +0.060] | 0.9277 vs 0.9551, Δ +0.0274 [+0.0038, +0.0591] | yes | `s4_acc.py`, own pooling + bootstrap over 32 sets from `per_set.tsv` (full scope: +0.034 [+0.019, +0.062]); all sets dev. Pair labels were not re-derived from the Newick files |
| C-S1-1 NCBI vs published QfO tree | conflict 0.000, nRF 0.128 | conflict 0.000 (0/51 splits), nRF 0.1277 | yes | `s1_conflict.py`, own bipartition code (ete3 only for parsing) |

Not rerun (cost or scope): full OrthoBench/QfO pipeline runs (>30 min), S4 tree inference, S5 label building,
S8 genome searches, QfO final scores (not in this task).

## 2. Discrepancies and concerns

1. **No numerical discrepancy found.** Every recomputed number matches to the reported precision.
2. **Undisclosed S0 split-v1 exposure (minor; disclosure gap).** Before the homology-aware v2 split, three S0 configs
   (`dc_a30_c50`, `deep_c20`, `deep_c50`, commit 6711002) were scored per-RefOG on the v1 "dev" half. 16 of those
   35 RefOGs are now in v2 **test** (RefOG009, 014, 016, 017, 018, 022, 027, 031, 032, 036, 037, 047, 049, 056, 064, 069).
   Their per-RefOG completeness was in git. `docs/S0_families.md` §1 records the v1→v2 switch and the re-scoring,
   but `docs/DISCLOSURES.md` does not say that test RefOGs were seen. The impact is probably small: it affects only
   S0 stage-1 choices, and v2 was built by a fixed-seed mechanical rule. It should still be disclosed.
3. **Bootstrap reproduction is not an independent Monte Carlo draw.** My seed-11 CIs match the reported ones to
   about 1e-6 because both implementations use `default_rng(11).integers(0, n, (B, n))`. Other seeds and B=10,000
   move the bounds by ≤0.7 pp, and no verdict changes.
4. **Commit-id inconsistency.** `results/freeze.json` gives `git_commit` 6c593a8 (the parent of the tag). The tag
   `freeze-v1` is 71117bf, and test_log entries cite f7bcddc. `ancestra/` is identical in all of them and in HEAD, so
   the frozen code is unaffected, but freeze.json should name the tag commit. freeze.json hashes only `ancestra/*.py`.
   `ancestra/schemas/`, `ancestra_eval/` and the harness scripts (`results/orthobench_final/scripts/run_cached.py`)
   are not hashed (git shows no ancestra/ change).
5. **Runtime comparison is indicative only.** Ancestra ran with 10 threads and the competitors with 8. DESIGN §5 asks
   for 24 threads, and the node was shared. I did not check whether gene-order table building (GFF parsing) is
   included in Ancestra's 3.91 CPU-h (DESIGN §5: "all preprocessing counts"). "3.9 vs 8.2 CPU-h" should carry these
   caveats.
6. **Giant mixed root HOGs.** The OB accurate output has 6 root HOGs with ≥500 genes (5,840 genes in total; largest
   1,910 at Amniota, likely olfactory receptors). Some have implausible copy profiles, e.g. HOG0010073 (Bilateria,
   1,193 genes: opossum 532, dog 376, chicken 189, human 3). None contains a RefOG gene, so the OrthoBench score is
   unaffected, but pair-level precision in such families (QfO VGNC over-merging) is a known weakness. 20 random root
   HOGs (seed 11) look plausible: 11/20 have at most 2 copies per species, sizes 2–176, and one is an exact
   RefOG069 group (8/8). The rest are one larger multi-copy vertebrate family (176 genes, up to 31 copies) and
   in-paralog expansions (fly 17 copies, chimp 9).
7. **Wording (minor over-claim).** RESEARCH_TREE 8b2 (QfO peek 1) says "Precision comparable" and "ties FastOMA" on
   point estimates without a CI. Under DESIGN §5 this should read "inconclusive".
8. **S4 vs presets.** S4's pre-registered verdict was "keep NJ in localtree" (cost rule). The frozen `accurate` preset
   uses FastTree, adopted later on dev OG-F. That is legitimate, but REPORT should say the S4 cost verdict was
   superseded and why.
9. **Test-set reuse for ablations.** 11 configs were scored on OB test. All were logged before scoring, and the primary
   endpoint (accurate) was fixed at the freeze, so this is acceptable. The ablation CIs are unadjusted (no Holm
   correction, DESIGN §5 "Layers"). Report them as exploratory.

## 3. Files
`results/verification/`: `ob_rescore.py/.json` (task 1), `ob_bootstrap_sensitivity.json`,
`official_cli_all70.txt`, `pytest_not_slow.log`, `toy_e2e.json`, `verify_orthoxml_accurate_go.{txt,json}`,
`spotcheck_hogs.{py,txt,json}`, `spotcheck_20hogs.orthoxml`, `s0_completeness.{py,json}`, `s5_auc.{py,json}`,
`s8_controls.json`, `s4_acc.{py,json}`, `s1_conflict.{py,json}`.
Suggested CLAIMS status: the rows above marked "yes" can move to `verified` (verifier, 2026-10-27). That covers
C-S0-2, C-S0-4, C-S1-1, C-S5-2 (human/mouse and the other pairs listed), C-S8-2 (recovery and negatives),
C-CORE-1 (toy range) and the OrthoBench final rows. The PI decides.

---

# Pass 2: final QfO results (results/qfo_final/, 2026-10-27 T+19h15 to T+19h22)

Same independence rules as pass 1. Code and outputs: `results/verification/qfo_*`.

## 4. Claim table (QfO)

| claim | reported | recomputed | match? | how / notes |
|---|---|---|---|---|
| Provenance cand2_accurate / cand2_fast | freeze-v1 code | manifest `driver.code_sha256`: 24/24 == `freeze.json` for both; md5 of `hogs.orthoxml` == manifest outputs; the final QfO run read exactly these files (`work/qfo_final/<run>.log` input path) | yes | Manifest commits are d61da8b (accurate) and 7dec881 (fast). test_log says 394da46 / d61da8b. Same code, but the ids are inconsistent (see §5.5) |
| SwissTree TPR/PPV accurate | .721/.926 | .7215/.9255 (±.037/.034) | yes | `qfo_refphylo.py`: own reimplementation of RefPhyloTest (truth parsed from `ReconciledTrees_*.drw`, own HOG flattening of `hogs.orthoxml`, protIds mapped with `mapping.json.gz`). Per-pair TP/FP/FN/TN counts are identical to the official raw file (3198/35/865/6318) |
| SwissTree fast | .731/.923 | .7307/.9230 | yes | counts identical to the official raw file |
| TreeFam-A accurate / fast | .583/.921, .608/.916 | .5831/.9213, .6080/.9165 | yes | counts identical (accurate 25,615 TP / 2,187 FP / 18,310 FN) |
| GO / EC accurate | .4748 (n 171,884) / .8985 (n 180,745) | .47484 / .89846, same n | yes | mean over the official raw per-pair files. 300/300 randomly sampled scored pairs per test are orthologs in my flattening of `hogs.orthoxml` |
| GO / EC fast | .4705 (181,857) / .8961 (196,729) | .47049 / .89606, same n | yes | 300/300 pair-provenance check |
| Published competitor rows (20 tools incl. FastOMA, OMA HOGs / Pairs / Groups) | from OpenEBench | 638/638 values in `summary.tsv` == `openebench_QfO_2022.tsv` (all benchmarks incl. STD/G-STD2 and n); no tool missing | yes | `qfo_competitors.py`. OMA rows are present in the main, STD and ablation tables |
| Ablations used freeze-v1 code | yes | 7/7 ablation manifests: 24/24 code hashes == freeze.json, git f7bcddc, not dirty; tree md5 in manifest == the arm's tree file; hogs md5 == manifest | yes | |
| Cache reuse only tree-independent stages | S0 + S2 only | `genes.tsv`, `families.npy`, `scores.npz`, `families_final.npy` byte-identical (the same inode) in det_nogo, cand2_* and all 7 ablations. merge/output are per-run (nlink 1). Driver `transplant()` re-keys only families/scores, checks that the key formula reproduces the source keys, and aborts on a code-key mismatch. The gene table is tree-independent (`load_proteomes` sorts by file name), and scoring takes no tree | yes | |
| Ablation c_astral SwissTree/TreeFam-A | .782/.929, .636/.914 | .7821/.9288, .6360/.9138; counts identical to the official raw file | yes | |
| All QfO final evaluations logged before scoring | 24 logged | g_accurate_nogo logged 09:48:13, its assess files written 09:53–09:54; FAS rerun logged as a rerun | yes | |

## 5. Discrepancies and concerns (QfO)

1. **DISCLOSURES #1 is factually wrong.** It says test-only TreeFam-A families were "never scored separately before
   the freeze". `results/qfo_peek2/treefamA_dev_test.json` (commit 6c593a8, 06:05) contains TreeFam-A **test**
   subset values for both cand2 presets (accurate .580/.913, fast .605/.914). RESEARCH_TREE 8.0 quotes them. The freeze
   followed one minute later. The impact is small: the presets were declared before peek 2, and peek 2's pooled
   TreeFam-A already includes those families. The disclosure still needs correcting.
2. **The pre-registered primary endpoints are not the ones 8.2 reports.** DESIGN §5 names TreeFam-A **test**, SwissTree,
   **G-STD2-Eukaryota (RF with n)** and EC (with n). 8.2 reports pooled TreeFam-A (it includes the 30% dev families used
   for tuning) and STD-Eukaryota. It does not report G-STD2-Eukaryota at all. On G-STD2-Eukaryota Ancestra is worst
   among the main comparators: RF .257 (n 3,379) accurate and .266 (4,753) fast, against OrthoFinder .249 (6,629),
   FastOMA .239 (3,092), OMA HOGs .237 (2,981) and Hieranoid .222 (3,643). On EC accurate .899 (n 181k) is below
   OrthoFinder .931 (n 204k) on both axes.
3. **Non-inferiority under the pre-registered rule** (my paired bootstrap from the official raw files, B=2000, seed 11;
   over the 749 TreeFam-A test families, or over the 18 SwissTree cases; `qfo_primary_ni.json`). Accurate on
   TreeFam-A test: TPR .580 vs OrthoFinder .739, Δ −0.160 [−0.210, −0.106], so **it fails**. Vs FastOMA +0.043
   [−0.005, +0.093] (NI on TPR). Vs OMA HOGs +0.002 [−0.058, +0.067] (inconclusive). Vs OMA Pairs +0.018
   [−0.028, +0.069] (NI). PPV is not shown NI vs any comparator except OMA HOGs and OMA Pairs. SwissTree TPR vs
   OrthoFinder is −0.103 [−0.177, −0.029] (fails), and inconclusive vs the rest. 8.2 says "not on the Pareto front;
   TreeFam-A ≈ FastOMA/OMA; behind OrthoFinder on recall", which is consistent. But the formal verdict (fails vs
   OrthoFinder on both recall endpoints; NI vs FastOMA / OMA Pairs on TreeFam-A TPR only) should be stated.
4. **The ablation CIs are too narrow.** The official ± (TreeFam-A ±0.006) is a pair-level binomial CI that ignores the
   clustering of pairs within families. A paired bootstrap over the 1,108 families gives: ASTRAL-Pro tree +0.028
   [−0.009, +0.077], FastME +0.022 [−0.001, +0.045], no gene order −0.012 [−0.028, +0.001], 5 NNIs +0.003
   [−0.020, +0.026] (`qfo_ablation_familyboot.txt`). So RESEARCH_TREE 8.2's "gene order helps QfO recall a little" and
   any reading of the ASTRAL/FastME recall gains as significant are **not supported**. The random-tree control
   (−0.243) is far outside any CI.
5. **Over-statements in 8.2.** "STD comparable" has no CI. Accurate STD-Eukaryota .113±.006 is close to OrthoFinder
   (.118) but clearly worse than OMA HOGs .086, OMA Pairs .085, Hieranoid .074 and InParanoid .079. Under DESIGN §5 it
   should read "inconclusive vs OrthoFinder, worse than OMA / Hieranoid / InParanoid". "VGNC precision weak (vertebrate
   over-merging)": PPV .982 equals OMA Pairs (.980) and OMA HOGs (.979), and the over-merging cause is a hypothesis.
   The explanation for the EC drop with ASTRAL/FastME trees ("Bacteria+Archaea merged") is not tested here.
6. **Mislabelled error bars.** `results/qfo_final/summary.md` calls the GO/EC "±" values SEM. They are 1.96×SEM
   (GO accurate: raw SEM 0.000567, reported 0.00111).
7. **Commit ids** (as in §2.4). test_log peek-2 entries name 394da46 / d61da8b, but the manifests say d61da8b / 7dec881.
   The ablation entry for g_accurate_nogo names git=b310546, but its manifest says f7bcddc. Code hashes match
   freeze-v1 in every case, so this is a bookkeeping issue only.

Not rerun: the official `run_qfo_local.sh` (not needed, because the per-pair counts reproduce the official raw output
exactly), STD/G-STD2 and FAS (not recomputed; competitor values checked against OpenEBench only).
