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

## 6. Deterministic run (`work/qfo2022_det`)
(filled in when the run completes)
