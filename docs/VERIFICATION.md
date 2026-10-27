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
