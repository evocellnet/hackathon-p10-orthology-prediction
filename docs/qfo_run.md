# QfO 2022 run (Ancestra core, current defaults) — run engineering

Status: first full run, 2026-10-26 (started T+3h08m, `scripts/clock.sh`). **No QfO benchmark has been run on the output**:
peeks are PI-authorised only and logged in `results/test_log.tsv` (DESIGN §4).

## 1. How to run

```bash
bash scripts/run_qfo_full.sh full            # 78 species, 12 threads (THREADS=… or --threads N)
bash scripts/run_qfo_full.sh sub10           # seeded 10-species subset (cost extrapolation)
bash scripts/run_qfo_full.sh full --stop-after scores   # stop cleanly after a stage
bash scripts/run_qfo_full.sh full --config cfg.json --name myarm  # {'scoring':{…},'merge':{…}} overrides
```

The script is **resumable**: the pipeline caches each stage under `work/qfo2022/<name>/stages/<stage>/` with a `DONE`
file holding a hash of the stage config and the upstream keys. A rerun skips every finished stage, and changing only the
merge config reruns only merge and output. A killed stage leaves only a `<stage>.tmp<pid>` dir, which the driver deletes.

Steps (all logged to `work/qfo2022/logs/<name>.<timestamp>.log`, gitignored):

| step | script | output |
|---|---|---|
| 1. inputs | `scripts/prepare_qfo_inputs.py` | `data/inputs/qfo2022/proteomes/` (78 symlinks), `guide_tree_ncbi.nwk`, `subsets/sub10_s11/`, `inputs.json` |
| 2. input check | `scripts/verify_ids.py` (aborts the run on failure) | `results/qfo_run/verify_ids.json` |
| 3. pipeline | `scripts/qfo_pipeline_driver.py` → `ancestra.pipeline.run_pipeline` under `/usr/bin/time -v` | `work/qfo2022/<name>/{hogs.orthoxml,pairs.tsv,orthobench_ogs.txt,manifest.json,stage_accounting.json}` |
| 4. output check | `scripts/verify_orthoxml.py` | `results/qfo_run/<name>_verify_orthoxml.json` |
| 5. summary | `scripts/summarize_qfo_run.py` | `results/qfo_run/<name>_summary.{json,md}` |

`qfo_pipeline_driver.py` is a thin wrapper and does not change the pipeline. It records per stage the wall time, the
user+sys CPU of the process and of its finished children (getrusage deltas), the peak RSS (own VmHWM from `stage.json`,
plus the maximum child RSS so far), and the node load at stage start and end. These go to `stage_accounting.json`.

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
  - pyham loads a `pyham_compat` copy (the full file up to 300k genes, else the sample), and its number of top-level
    HOGs is as expected;
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

## 5. Full run results
(filled in from `results/qfo_run/full_summary.json` / `full_verify_orthoxml.json`, see below)
