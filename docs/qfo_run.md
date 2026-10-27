# QfO 2022 run (Ancestra core, current defaults) — run engineering

Status, 2026-10-26 (`scripts/clock.sh`):
- **Dry run** (old, non-deterministic S0): started at T+3h08m; S0 and S2 ran on the code of 17:49. The merge was rerun
  with the code of 19:56, in both arms. It is kept as timing and memory data only, in `work/qfo2022/full_{nogo,go}`.
- **Deterministic run** (S0 from commit acb1f43+; `work/qfo2022_det/det_{nogo,go}`): started at 20:13 CET.

**No QfO benchmark has been run on any output.** Peeks are PI-authorised only and logged in `results/test_log.tsv`
(DESIGN §4).

## 1. How to run

```bash
bash scripts/run_qfo_full.sh full                       # 78 species, 12 threads, both merge arms (default)
bash scripts/run_qfo_full.sh full --workdir work/qfo2022_det --name det   # the deterministic run (§6)
bash scripts/run_qfo_full.sh full --gene-order on|off   # only one arm
bash scripts/run_qfo_full.sh sub10                      # seeded 10-species subset (cost extrapolation)
bash scripts/run_qfo_full.sh full --stop-after scores   # stop cleanly after a stage
bash scripts/run_qfo_full.sh full --config cfg.json --name myarm   # {'scoring':{…},'merge':{…}} overrides
```

**Arms.**
- `<WORKDIR>/<PREFIX>_go` (arm A) is the current default merge: `--gene-order data/inputs/qfo2022/gene_order`, i.e.
  constrained synteny (`ancestra.pipeline.SYNTENY_DEFAULT`).
- `<WORKDIR>/<PREFIX>_nogo` (arm B) is the same merge without gene order.
- With `--gene-order both` (the default), B runs first and computes S0 and S2. A then hard-links B's
  input, families and scores stages (`cp -al`, driver `--reuse-stages-from`), so S0 and scoring run only once.

**Resumable, and code-aware.**
- The pipeline caches each stage under `<run>/stages/<stage>/` with a `DONE` hash of the stage config and the upstream
  keys.
- The driver (`scripts/qfo_pipeline_driver.py`) adds `stages/code_keys.json`. For each stage it holds a sha256 over:
  - the `ancestra/` files that stage runs: input = io, newick; families = families, io; scores = scoring;
    merge = merge, merge_numba, synteny, genome; output = orthoxml, merge;
  - the effective families config (`asdict(FamilyConfig(**DEFAULT_FAMILY_CONFIG))`);
  - the gene-order dir (merge stage only);
  - the upstream code key.
- On a mismatch, or a missing record, the driver deletes that stage and everything downstream before the pipeline
  runs. A fix in `ancestra/families.py` therefore reruns S0, scoring, merge and output; a fix in `merge.py` reruns only
  merge and output.
- `pipeline.py` is not hashed. Its stage structure is covered by the pipeline's own keys, and hashing it would
  invalidate S0 on every orchestration edit.
- `manifest.json` gets a `driver` block with:
  - the git commit and `ancestra/` dirty flag;
  - the sha256 of each code file;
  - the stage code keys;
  - the effective families config and its hash;
  - the gene-order dir;
  - the invalidated and adopted stages.
- `--adopt-unrecorded` accepts cached stages that predate code keys. It is for dry runs only and is logged; it was used
  once, for the dry-run merges in §5.
- A killed stage leaves only a `<stage>.tmp<pid>` dir, which the driver deletes.
- A lock file stops two drivers from running on one dir.
- **Do not edit `scripts/run_qfo_full.sh` while it runs.** Bash reads scripts incrementally; on the dry run I stopped
  the wrapper twice for this reason (the pipeline kept running) and ran the verification by hand.

Steps (logged to `<WORKDIR>/logs/`, gitignored):

| step | script | output |
|---|---|---|
| 1. inputs | `scripts/prepare_qfo_inputs.py` | `data/inputs/qfo2022/{proteomes,gene_order}/` (78 symlinks each), `guide_tree_ncbi.nwk`, `subsets/sub10_s11/`, `inputs.json` |
| 2. input check | `scripts/verify_ids.py --gene-order …` (aborts the run on failure) | `results/qfo_run/verify_ids.json` |
| 3. pipeline | `scripts/qfo_pipeline_driver.py` → `ancestra.pipeline.run_pipeline` under `/usr/bin/time -v` | `<run>/{hogs.orthoxml,pairs.tsv,orthobench_ogs.txt,manifest.json,stage_accounting.json}` |
| 4. output check | `scripts/verify_orthoxml.py` | `results/qfo_run/<run>_verify_orthoxml.json` |
| 5. summary | `scripts/summarize_qfo_run.py` | `results/qfo_run/<run>_summary.{json,md}` |

`stage_accounting.json` records, per stage:
- wall time;
- user+sys CPU of the process and of its finished children (getrusage deltas);
- peak RSS: own VmHWM from `stage.json`, plus the maximum child RSS so far;
- node load at the start and end of the stage.

## 2. Inputs (DESIGN §4)
- **Proteomes:** the 78 canonical `UPxxxxxxxxx_<taxid>.fasta` files of `qfo_release-2022_02_with_updated_UP000000437`
  (`data/qfo/proteomes_2022_02_updZEBRAFISH/{Archaea,Bacteria,Eukaryota}`), symlinked. Each file's md5 is checked against
  the release `md5sums.txt`. `_additional` (isoforms) and `_DNA` files are never linked. Total: 988,778 proteins.
- **Species names:** the file stem (`UP000005640_9606`). This is exactly the leaf naming of the S1 QfO trees
  (`results/S1_speciestree/qfo/*.nwk`); `prepare_qfo_inputs.py` asserts that the leaf set equals the stem set.
- **Guide tree:** `results/S1_speciestree/qfo/ncbi.nwk` (S1 default: NCBI taxdump, k-way polytomies), copied to
  `data/inputs/qfo2022/guide_tree_ncbi.nwk` (md5 in `inputs.json`). After loading it has 52 internal nodes (named NCBI
  clades, root `cellular organisms` = Bacteria/Archaea/Eukaryota trichotomy).
- **`verify_ids.py` result (PASS):** 78 species; canonical files with matching md5; 988,778 bare accessions (as
  `ancestra.io` reads them), no duplicates; 0 of 394,954 isoform entries of the `_additional` files present; **coverage
  of the benchmark mapping 100%** (988,778 / 988,778); 0 input IDs unknown to the mapping; species map 1:1 to the 78
  benchmark species codes; tree leaves equal the species. Negative control: adding the human `_additional.fasta` to the
  human file fails `canonical_files`, `bare_accessions` (21,958 `-n` isoform IDs), `no_isoforms` and `coverage` (7.6%
  unknown IDs).
- **Gene order (S5):** `data/inputs/qfo2022/gene_order/<species>.parquet` links the 78 tables of
  `data/genomes/qfo2022/` (`ancestra/genomes_qfo.py`, docs/S5_synteny.md §1.2); their md5s are in `inputs.json`.
  `verify_ids.py --gene-order` checks that:
  - every species has a table;
  - every `protein_id` is an input protein of that species and appears once;
  - (contig, order) positions are unique.

  It passes. Coverage is 0.931 overall; the lowest species are *Parastagonospora* (UP000001055) 0.0001, amphioxus
  0.078, zebrafish 0.644, *Leptospira* (UP000002521) 0.699 and rat 0.829. Those species get little or no synteny
  evidence.
- `mapping.json.gz` comes from the QfO reference data. Only the verification scripts read it (`verify_ids.py`, and
  `verify_orthoxml.py` for the protId → protein-number mapping of the QfO extractor). The method never reads it, which
  `tests/test_no_leak.py` checks.

## 3. Verification scripts
- **`scripts/verify_orthoxml.py RUN_DIR`** runs these checks (DESIGN §6); each failure is reported, not raised:
  - XSD validity;
  - geneRef integrity (unique gene ids, protIds and group ids; no dangling geneRef);
  - every gene referenced at most once, and the declared genes equal the input gene table exactly, each in its own species;
  - the TaxRange of every orthologGroup equals the LCA of its member species in the guide tree, and its `taxonId` agrees;
    nested orthologGroups lie strictly below their parent; a paralogGroup's level contains its members and lies strictly
    below the enclosing orthologGroup;
  - the number of ortholog pairs implied by the HOGs equals the number of `pairs.tsv` lines;
  - **flattening vs QfO:** on a seeded sample (3,000 random top-level orthologGroups plus the 5 largest, seed 11), the
    QfO benchmark's own `map_relations.py` (`PairwiseOrthologRelationExtractor` + `parse_orthoxml`, run on the file,
    including its id-mapping checks) gives the same pairs as `ancestra.orthoxml.hogs_to_pairs` and as `pairs.tsv`;
  - pyham loads a `pyham_compat` copy (the full file, up to 5M genes; QfO with 741k grouped genes loads in 27–29 s
    and 3.5 GB), and its number of top-level HOGs is as expected;
  - the official Open_Orthobench reader (`read_orthogroups`) accepts `orthobench_ogs.txt`: no duplicates, no unknown
    genes, and OGs equal the multi-gene top-level groups. On OrthoBench inputs it also runs `check_orthogroups`.
- **Controls:**
  - It passes on the OrthoBench run `work/core_v0_ob/run_v2` (251k genes; the official `check_orthogroups` reports
    100% of genes) and on a toy run.
  - Corrupted toy copies fail as expected:
    - a wrong TaxRange fails `taxrange_lca`;
    - a duplicated geneRef fails `genes_once`, `pairs_count`, `flattening_vs_qfo` and `orthobench_export`;
    - a dangling geneRef fails `xsd`, `generef_integrity` and `read_orthoxml`.
  - A first version mis-inverted the QfO mapping, which has several xrefs per protein, and failed on the subset. It was
    fixed by inverting through the file's protIds.
- **`tests/test_no_leak.py`** has four tests:
  - a positive control for the scanner;
  - a **static scan** of `ancestra/` (identifiers, string literals and imports; comments and docstrings are listed but
    allowed) for `reference_data, refog, swisstree, treefam, lineage_tree, species_tree_samples, participants, vgnc,
    bench_repo`;
  - a report of comment and docstring mentions;
  - a **runtime audit**: the pipeline runs with default S0 on a toy simulation under `sys.addaudithook`, and no opened
    path or subprocess command line may match the list.

  **Current status: passes.** The RefOG evaluator (`read_refogs`, `evaluate_refogs`, the `evaluate` sub-command and the
  RefOGs path) was moved out of `ancestra/families.py` into the separate evaluation package `ancestra_eval/`
  (`ancestra_eval/families_eval.py`; CLI `python -m ancestra_eval.families_eval evaluate ...`). The static scan covers
  `ancestra/` only. Repro: `pytest tests/test_no_leak.py::test_static_scan`.

## 4. Cost extrapolation (10-species subset, before the full run)
Subset `sub10_s11` (seed 11, domain-stratified: 6 Eukaryota, 3 Bacteria, 1 Archaea; 129,119 proteins = 13.1% of N),
12 threads, node load 27–41 on 24 CPUs (wall times are upper bounds). `results/qfo_run/extrapolation.json`:

| stage | subset wall | CPU | full, linear in N | full, power law (slope 1.57*) | full, quadratic bound |
|---|---|---|---|---|---|
| families (S0) | 241 s | 1,589 s | 0.5 h | 1.6 h | 3.9 h |
| scores (S2) | 303 s | 2,078 s | 0.6 h | 2.0 h | 4.9 h |
| merge | 4.2 s | 4 s | – | – | 0.07 h |
| output | 2.2 s | 2 s | – | – | 0.04 h |

\*slope fitted between the subset and the OrthoBench run (251k proteins, scores 860 s at 12 threads). Even the quadratic
bound is below the 6 h limit, so the full run was launched without stopping.

## 5. Dry run (old non-deterministic S0): timing, memory, verification
These are the `work/qfo2022/full_{nogo,go}` runs, at 12 threads on a shared 24-CPU node with load 28–65; wall times are
upper bounds. Sources: `results/qfo_run/full_{nogo,go}_summary.{json,md}` and `full_{nogo,go}_verify_orthoxml.json`.

| stage | wall | CPU (self+children) | peak RSS (own / child) | notes |
|---|---|---|---|---|
| input | 0.6 s | 0.6 s | 0.6 GB | 988,778 genes, 78 species, 52 internal nodes |
| families (S0) | 2,137 s (36 min) | 11,565 s (3.2 CPU-h) | 6.2 / 6.9 GB | stage 1 (DIAMOND cluster) 795 s; stage 2 (N×R search) 848 s; CC + Leiden 479 s |
| scores (S2) | 5,404 s (1.5 h) | 31,805 s (8.8 CPU-h) | 29.9 GB | 9 DIAMOND blocks 3,886 s; 414M raw hits; MMseqs2 rescue of 21.9M missed pairs 153 s; cross-family rescue: 371 joins, 188k genes re-scored, 1,098 s |
| merge, no gene order | 131 s | 130 s | 7.1 GB | single-threaded (numba) |
| merge, gene order | 348 s | 346 s | 8.2 GB | constrained synteny |
| output | 34–36 s | 34–35 s | 5.3 GB | XML + XSD + pairs + OGs |

The whole S0→output chain takes about 2.1 h at 12 threads. Extrapolations (§4) were linear 1.1 h, power law 3.6 h,
quadratic 8.8 h; the measured scoring time (1.5 h) lies between the linear and power-law estimates.

**Verification: both arms PASS all checks.** XSD; geneRef integrity; 988,778 genes declared = input, each referenced at
most once; TaxRange = LCA for all ~379k orthologGroups (0 nesting or paralog-level violations); 15.3M pairs = the
`pairs.tsv` lines; QfO `map_relations.py` = ours = `pairs.tsv` on 3,005 groups (88k genes, 2.08M pairs, including a
7,029-gene HOG); pyham loads the full file (35k top-level HOGs, 29 s); the OrthoBench reader accepts the export. Run
times: `verify_orthoxml` about 2.5 min and 7 GB; `verify_ids` about 25 s.

Dry-run content (old S0, for orientation only): 163,092 S0 families (132,825 singletons = 13.4% of genes, largest
7,464 genes, Σg²/N² = 0.0023); 36.3k / 36.5k root HOGs (nogo / go) covering 741k genes (75%). The largest root HOG has
7,029 / 6,528 genes from only 15 species at Bilateria. It fills most of a Leiden-capped S0 family, which is a precision
risk to inspect.

## 6. Deterministic run (`work/qfo2022_det`): arms A (gene order) and B (no gene order)
- **Code.** Commit 34909cd plus uncommitted edits in `ancestra/merge.py` and `ancestra/localtree.py`
  (`ancestra_dirty=True`). The exact sha256 of every `ancestra/*.py` file is in `manifest.json → driver.code_sha256`.
  - A and B ran identical code. S0 and S2 were computed once (in B) and hard-linked into A.
  - The first merge pair (21:45 / 21:50) straddled a `merge.py` edit by the core agent. The driver had keyed the merge
    stage on too few files (`localtree.py` was missing), so it did not catch that. The merge and output stages are now
    keyed on all `ancestra/*.py` files, and both merges were rerun at 22:03–22:17 with identical code.
- **Families config hash:** `dea7e48e5e4da45b`.
- **Load and verification.** 12 threads; node load 22–45 on 24 CPUs. `verify_ids` PASS (including gene order).
  `verify_orthoxml` PASS for both arms, on every check; pyham loads the full files.

| stage | wall | CPU | peak RSS | load start→end |
|---|---|---|---|---|
| S0 families | 1,534 s | 10,786 s (3.0 CPU-h) | 6.4 GB (child 6.8) | 30.9→21.7 |
| S2 scores | 3,886 s | 31,979 s (8.9 CPU-h) | 29.5 GB | 21.7→41.2 |
| merge B (no gene order) | 144 s | 143 s | 8.4 GB | 37.6→40.8 |
| merge A (gene order) | 389 s | 370 s | 9.4 GB | 23.3→39.6 |
| output (A / B) | 36 s / 36 s | 34 s / 36 s | 5.6 / 5.5 GB | – |

End to end, B takes 1.56 h wall and A adds 7 min.

- **S0:** 163,187 families (30,300 multi-gene; 132,887 singletons = 13.4% of genes; largest 7,454; Σg²/N² = 0.0023).
- **S2:** 9 DIAMOND blocks (2,570 s, 407M raw hits) and 74.2M scored pairs. The cross-family rescue made 358 joins,
  re-scored 176k genes (938 s) and left 162,829 families.

| | A (gene order) | B (no gene order) |
|---|---|---|
| ortholog pairs (QfO flattening) | 15,104,847 | 15,591,893 |
| orthologGroups / paralogGroups | 379,792 / 114,791 | 377,609 / 115,639 |
| root HOGs (≥2 genes) / genes in them | 36,628 / 740,445 (74.9%) | 36,399 / 740,516 |
| root HOG size median / mean / p90 / p99 / max | 6 / 20.2 / 48 / 167 / 6,154 | 6 / 20.3 / 48 / 176 / 6,183 |
| root HOGs of size 2–5 / 6–20 / 21–78 / 79–200 / 201–1000 / >1000 | 17,652 / 10,509 / 7,033 / 1,155 / 249 / 30 | 17,494 / 10,500 / 6,964 / 1,152 / 262 / 27 |

**Arm A, top 20 root HOGs** (genes/species @ level): 6154/57 @cellular organisms, 6026/14 @Chordata, 4668/49,
3824/52, 3148/35, 3138/51, 3043/51, 2863/28 (all @cellular organisms), 2759/33 @Eukaryota, 2361/49, 2296/29 (@cellular
organisms), 2131/26 @Eukaryota, 1930/44 @Eukaryota, 1824/27 @Eukaryota, 1777/49, 1536/28 (@cellular organisms),
1488/28 @Eukaryota, 1433/36 @cellular organisms, 1309/39 @cellular organisms, 1282/36 @Eukaryota.
- 30 root HOGs exceed 1,000 genes. They come from Leiden-capped S0 families (about 7.5k genes).
- The 6,026-gene / 14-species Chordata HOG (about 430 genes per species) is a likely over-merge, probably a
  zinc-finger or olfactory-receptor-like superfamily, and a precision risk.

HOGs per level: `results/qfo_run/det_{go,nogo}_summary.md`. Size distribution and top 20: `results/qfo_run/det_root_hog_size_distribution.json`.

## 7. QfO peek 1 (PI-authorised; fast tests only)
- **Logging and inputs.** Logged in `results/test_log.tsv` before scoring: one line per arm, with config hashes
  6452b3441294 (A) and 92e5bf6dcf4b (B) and git 855db34. The code sha256s equal the manifest's.
- **Scoring.** `bench_repo/scripts/run_qfo_local.sh` scored `work/qfo2022_det/det_{go,nogo}/hogs.orthoxml` on GO, EC,
  VGNC, SwissTrees, TreeFam-A and SP. That took about 5 min per arm; outputs are in `work/qfo_peek1/<arm>/`.
- **Competitors.** Their public prediction files were re-scored locally on SwissTrees, TreeFam-A, VGNC and SP
  (`scripts/qfo_peek1_rescore_competitors.sh`) to get evaluated-pair counts (TP+FP from the raw files). This
  reproduces the published TPR/PPV exactly (60/60 values; `results/qfo_peek1/competitor_reproduction.tsv`).
- **Table.** `scripts/make_qfo_peek1_table.py` writes `results/qfo_peek1/summary.{tsv,md}`.
- **Not run:** STD and G_STD2 tests (kept for the final run).

## 8. TreeFam-A dev split and dev-only diagnosis (after peek 1)
- **Split** (`scripts/make_treefamA_split.py` → `splits/treefamA_split.json`, sha256 in `splits/SHA256SUMS`).
  - The QfO 2022 reference stores TreeFam-A as **one case** with no family IDs.
  - Families are therefore the connected components of the reference relation graph: 1,108 families, 11,182
    proteins, 1,013 of them complete cliques.
  - Family ID = `TFA:<smallest UniProt accession>`; dev = `int(sha256(id),16) % 10 < 3`.
  - Result: **359 dev families** (32.4%; 13,993 S / 11,420 D relations) and 749 test families.
- **Subset scoring** (`scripts/qfo_treefam_subset.py RUN_DIR…`) runs on the official raw RefPhyloTest pairs with
  the official formula.
  - Note: RefPhyloTest counts every reference pair as 1/2, on top of a pseudocount of 1.
  - Pooled dev+test equals the official TPR/PPV to 1e-8 for arms A and B and for OrthoFinder.
  - It refuses test-only values without `--allow-test`.
  - Dev results: A TPR 0.572 / PPV 0.928; B 0.575 / 0.931; OrthoFinder 0.732 / 0.920
    (`results/qfo_dev_diag/treefamA_dev_scores.json`).
- **Diagnosis** (`scripts/qfo_dev_diag.py`; dev pairs only, test pairs are skipped at parse time): see
  `results/qfo_dev_diag/diag_summary.json` and `diag_pairs.tsv`.

## 9. Frozen candidate cand1 (presets with so_struct + deep_so_absorb), 2026-10-27
- **Runs.** Merge-only reruns on the cached det S0/S2 (driver `--preset fast|accurate --reuse-stages-from det_nogo`,
  gene order on): `work/qfo2022_det/cand1_{fast,accurate}`. `verify_ids` PASS.
  - fast: merge 374 s / 371 CPU-s / 9.5 GB.
  - accurate: merge 4,073 s / 7.2 CPU-h / 9.5 GB, 61,216 FastTree jobs.
- **verify_orthoxml FAILS (structure).**
  - About 8.9–9.0k orthologGroups are nested directly in an orthologGroup of the **same** TaxRange.
  - About 300 paralogGroups have a level that is not strictly below the enclosing orthologGroup, or that does not
    contain their members.
  - pyham crashes ("NoneType has no attribute children").
  - Pair-level checks pass: pairs count, flattening vs QfO and OrthoBench export.
  - `det_go_accurate` (presets before so_struct) had 0 violations.
  - Minimal repro: `results/qfo_run/repro_so_struct_nesting/HOG0000410_cand1_fast.orthoxml` (4 genes:
    OG@Methanobacteriota{OG@Methanobacteriota, OG@Methanobacteriota}).
- **Peek 2 not run** pending the PI's decision (verification gate, DESIGN §6). Dev-only TreeFam-A (not a peek):
  cand1_fast TPR 0.616 / PPV 0.920; cand1_accurate 0.657 / 0.914 (`results/qfo_dev_diag/treefamA_dev_scores_cand1.json`).

## 10. Frozen candidate cand2 and QfO peek 2, 2026-10-27
- **Reuse refused.** The core agent's `work/qfo_recall` merges have no manifest or code hashes, so they could not be
  reused under the PI's rule. Both arms were rerun merge-only into `work/qfo2022_det/cand2_{accurate,fast}`.
  - Their code sha256 equals HEAD = b504c94 for every `ancestra/*.py` file (`ancestra_dirty=False`).
  - accurate: merge 3,003 s / 7.0 CPU-h / 9.6 GB, 61,216 FastTree jobs.
  - fast: merge 375 s / 372 s CPU / 9.5 GB.
- **Verification.** `verify_ids` PASS. `verify_orthoxml` PASS for both arms: 0 nesting and 0 paralog-level violations,
  pyham loads.
- **Peek 2.** `scripts/qfo_peek2.sh` gates on that verification and logs both arms in `results/test_log.tsv` before
  any scoring. It then runs the fast tests and writes `results/qfo_peek2/summary.{tsv,md}` (with the peek-1 rows and
  the published competitors) and `results/qfo_peek2/treefamA_dev_test.json`.

## 11. Final QfO evaluation (freeze-v1), 2026-10-27
**Code.** The working-tree sha256 of all 24 files in `results/freeze.json` equals both cand2 manifests. No `ancestra/`
file was changed during the evaluation. Every evaluation is logged in `results/test_log.tsv` with peek=final.

**Full suite** (`work/qfo_final/cand2_{accurate,fast}`; `bench_repo/scripts/run_qfo_local.sh`, 13 tests, `-j 4` per arm):
- GO, EC, VGNC, SP, SwissTrees and TreeFam-A.
- STD_{Eukaryota,Fungi,Bacteria} and G_STD2_{Luca,Eukaryota,Fungi,Vertebrata}. These took about 2.2 h per arm;
  STD_Eukaryota and STD_Fungi took 2.0–2.3 h each.

**FAS** (`scripts/qfo_final_fas.sh`). The local runner does not wire FAS, so this runs the main.nf
`fas_benchmark` step: upstream `fas_benchmark.py`, unmodified.
- Setup:
  - greedyFAS 1.18.14 (pip wheel) in a separate venv, `conda_envs/fas_venv`, layered on `conda_envs/data`.
  - `fas.setup -t conda_envs/fas_tools --noAnno` writes `pathconfig.txt` without annotation tools; the precomputed
    annotations in the reference data are used.
  - `TMPDIR`/`SQLITE_TMPDIR` point to the run dir, because the sqlite sort overflows the 64 MB overlay root.
- Runs:
  - The first attempt hit "disk full" (the overlay).
  - The second computed 0 of the missing pairs (no `pathconfig.txt`). It is kept as
    `assess_FAS_firstrun_precomputed_only.json` and not reported.
  - The third run is the official procedure: precomputed sample plus 9,000 computed missing pairs. It took 42 min at
    `--cpus 2`.
- FAS is stochastic: upstream uses an unseeded shuffle.

**Ablations** (`scripts/qfo_final_ablations.sh`; fast preset; merge-only on the cached det S0/S2; each arm passes
verify_orthoxml before scoring).
- **Guide-tree arms (c–f)** use the driver's `--transplant-from`. The pipeline keys its input stage on the tree md5
  and chains every later key to it. S0/S2 do not read the tree, so their stages are hard-linked and re-keyed. The key
  formula is first checked to reproduce the source run's DONE keys exactly.
- **Arm (b)** uses `--family-rescue both --family-rescue-guard --structure-3di data/structure/qfo_3di_best.fa`:
  rescue stage 39 min and merge 13.5 min at 3 threads.
- **Arm (g)** is the accurate preset without gene order (FastTree, 3 workers).

**Summary** (`scripts/make_qfo_final_summary.py` → `results/qfo_final/summary.{tsv,md}`, `ablation.tsv`,
`pareto.png` + `pareto_data.tsv`).
- It covers all 20 published 2022 participants. The OMA resources (FastOMA, OMA HOGs, OMA Pairs, OMA Groups) are
  always shown, including in the ablation context and the plot (user request).
- Pareto colours are the validated dataviz palette (blue / orange / aqua, all-pairs check). Aqua is below 3:1 on the
  surface, so it gets direct callout labels, a square marker and the table (relief rule).
