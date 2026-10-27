# OrthoBench competitor baselines (12 species, 251,378 proteins, 70 RefOGs)

Scores: official Open_Orthobench `benchmark.py` (pairwise, RefOG-size normalised; percent). dev/test = sum of per-RefOG TP/FP/FN over the
34/36 RefOGs of `splits/orthobench_split.json` (dev+test reproduces the all-70 numbers; verified identical to benchmark.py). Runs: 8 threads,
`/usr/bin/time -v`, SHARED node (24 cores, load avg 20-55 during runs) -> wall times are inflated vs an idle node; CPU h = user+sys.

| method | F all | F dev | F test | precision | recall | exact RefOGs | wall (h) | CPU (h) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder | 70.1 | 70.5 | 69.7 | 61.4 | 81.7 | 19/70 | 1.29 | 8.2 | 5.6 |
| orthofinder_mcl | 58.7 | 54.2 | 64.4 | 45.7 | 82.0 | 18/70 | 1.29 | 8.2 | 5.6 |
| sonicparanoid_default | 47.4 | 34.2 | 70.4 | 37.3 | 64.7 | 18/70 | 1.31 | 7.9 | 8.1 |
| proteinortho | 45.1 | 42.6 | 47.5 | 97.0 | 29.4 | 6/70 | 0.23 | 1.6 | 0.8 |
| fastoma_luca | 31.7 | 29.6 | 33.5 | 24.0 | 46.7 | 19/70 | 3.50 | 12.7 | 13.6 |
| fastoma_luca_ogs | 31.7 | 27.9 | 35.3 | 95.1 | 19.0 | 4/70 | 3.50 | 12.7 | 13.6 |

Precision/recall on dev / test:

| method | P dev | R dev | P test | R test |
|---|---|---|---|---|
| orthofinder | 60.8 | 84.0 | 62.1 | 79.5 |
| orthofinder_mcl | 39.9 | 84.6 | 54.1 | 79.5 |
| sonicparanoid_default | 24.0 | 59.3 | 70.8 | 70.1 |
| proteinortho | 96.1 | 27.4 | 97.7 | 31.4 |
| fastoma_luca | 23.6 | 39.8 | 24.3 | 53.7 |
| fastoma_luca_ogs | 95.2 | 16.3 | 95.0 | 21.7 |

## Paired bootstrap vs OrthoFinder (dF = F(method) - F(OrthoFinder), percentage points; RefOGs resampled, B=2000, seed 11)

| method | subset | dF | 95% CI | SE |
|---|---|---|---|---|
| orthofinder_mcl | all | -11.4 | [-18.8, -3.4] | 4.0 |
| orthofinder_mcl | dev | -16.3 | [-28.5, -2.8] | 6.6 |
| orthofinder_mcl | test | -5.3 | [-11.4, +0.0] | 3.0 |
| sonicparanoid_default | all | -22.8 | [-39.7, +1.6] | 12.4 |
| sonicparanoid_default | dev | -36.4 | [-49.1, -1.7] | 15.4 |
| sonicparanoid_default | test | +0.7 | [-8.9, +11.6] | 5.2 |
| proteinortho | all | -25.0 | [-33.3, -16.1] | 4.5 |
| proteinortho | dev | -27.9 | [-37.2, -16.9] | 5.3 |
| proteinortho | test | -22.2 | [-34.8, -7.1] | 7.1 |
| fastoma_luca | all | -38.4 | [-52.5, -16.3] | 9.2 |
| fastoma_luca | dev | -40.9 | [-57.6, -15.2] | 11.6 |
| fastoma_luca | test | -36.2 | [-56.3, -4.2] | 14.0 |
| fastoma_luca_ogs | all | -38.5 | [-47.6, -28.3] | 4.9 |
| fastoma_luca_ogs | dev | -42.7 | [-53.6, -31.1] | 5.7 |
| fastoma_luca_ogs | test | -34.4 | [-48.6, -18.2] | 7.8 |

## Notes

- **orthofinder**: OrthoFinder 3.2 (conda, bioconda) ; diamond 2.2.8; famsa 2.4.1; FastTree 2.2.0; mcl 22.282 | orthofinder -f <12 proteomes> -t 8 -a 8 -mog (default: DIAMOND, MSA/famsa+FastTree gene trees, MCL I=1.2); node load avg ~20-55 on 24 cores
- **orthofinder_mcl**: same run as orthofinder (MCL_Orthogroups.tsv before gene-tree splitting) | orthofinder -f <12 proteomes> -t 8 -a 8 -mog (default: DIAMOND, MSA/famsa+FastTree gene trees, MCL I=1.2); node load avg ~20-55 on 24 cores
- **sonicparanoid_default**: SonicParanoid 2.0.9 (bioconda), bundled diamond; mode default | sonicparanoid -m default -t 8 -I 1.5 (diamond very-sensitive, domain-arch step on). Includes one-time Pfam profile DB indexing at first run after install (few minutes). Node load avg ~20-55 on 24 cores
- **proteinortho**: Proteinortho 6.3.6 (bioconda) diamond 2.2.8 | proteinortho -p=diamond -cpus=8 (defaults: e=1e-5, identity 25, coverage 50, conn 0.1, purity; PoFF off). Input proteomes had internal '*' stop symbols stripped (sed, required by ProteinOrtho)
- **fastoma_luca**: FastOMA 0.5.1, omamer 2.1.0, LUCA.h5 (zenodo 20814376, 2026-06-23); species tree hand-built from NCBI taxonomy; variant=tophogs | FastOMA 0.5.1 nextflow pipeline (local, no container), OMAmer 2.1.0 + LUCA.h5 (OMAmer DB, OMA-derived: LEAKAGE - OMA contains these species), nextflow executor cpus=8. wall includes omamer search (12 x 1 cpu), nextflow overhead + report notebook; peak RSS = largest single process
- **fastoma_luca_ogs**: FastOMA 0.5.1, omamer 2.1.0, LUCA.h5; variant=OrthologousGroups | FastOMA 0.5.1 nextflow pipeline (local, no container), OMAmer 2.1.0 + LUCA.h5 (OMAmer DB, OMA-derived: LEAKAGE - OMA contains these species); same run as fastoma_luca; peak RSS = largest single process

## Caveats and provenance

- **Leakage (FastOMA):** the OMAmer LUCA database (`LUCA.h5`, Zenodo 20814376, 2026-06-23) is derived from OMA, which contains these Ensembl species and proteins; FastOMA's placement of genes into root HOGs is therefore not independent of the reference. Treat FastOMA as an informed upper-bound-style comparator, not a clean one. FastOMA assigns only 198,615 of 251,378 proteins (OMAmer leaves the rest unplaced); unplaced genes count as missing (FN) in the benchmark. `fastoma_luca` = top-level HOGs of `FastOMA_HOGs.orthoxml` (identical to `RootHOGs.tsv`); `fastoma_luca_ogs` = FastOMA `OrthologousGroups.tsv` (strict 1-per-species-style ortholog groups, deliberately fine-grained: very high precision, low recall) and is not the OrthoBench-appropriate granularity. The species tree given to FastOMA was hand-built from NCBI taxonomy (root = Bilateria).
- **Coverage:** OrthoFinder/SonicParanoid/ProteinOrtho/FastOMA outputs contain only genes placed in a group of size >=2; singletons are simply absent (equivalent to singletons for the pairwise score).
- **Timing:** runs were executed one at a time with 8 threads on a shared 24-core node (load average 20-55 during all runs, see `<method>/load.log`, uptime every 60 s), so wall times are 1.5-3x what an idle node would give; CPU = user+sys from `/usr/bin/time -v`. `peak RSS` = `Maximum resident set size` of the largest process in the tree (for FastOMA/nextflow this is the largest single task, not the sum). FastOMA wall includes nextflow orchestration and the final report notebook; LUCA.h5 download (9.9 GB) and a failed first attempt (species-tree typo, 3 s) are excluded.
- **OrthoFinder 3.2** default = MSA gene trees (famsa + FastTree), DIAMOND, MCL I=1.2, then HOG-derived `Orthogroups.tsv` (`orthofinder`). `orthofinder_mcl` scores the pre-tree MCL orthogroups from the same run (`-mog`).
- **ProteinOrtho** needed proteome files with internal `*` stop symbols removed (sed on sequence lines only; headers untouched, same protein count).
- Scoring: `scripts/score_competitor.sh` -> `scripts/to_orthobench_og.py` / `orthoxml_to_orthobench_og.py` -> `run_orthobench.sh` (official benchmark.py) and `orthobench_eval.py` (per-RefOG, dev/test). Raw outputs: `data/orthobench/runs/<method>/` (not committed).
- **SonicParanoid2 `-m sensitive` NOT completed:** MMseqs2 s=7.5 all-vs-all needed roughly 5-6 h on the shared node (48/144 alignments after ~2 h, in two attempts); both attempts were stopped (first to prioritise FastOMA, second on PI request because the node is needed for scaling runs). No score is reported for it. Partial outputs: `data/orthobench/runs/sonicparanoid_sensitive_aborted_*`.
- Only OrthoBench `Input/*.fa` were read; a failed ProteinOrtho first attempt had written three `*.fa.len` files into `Input/`, removed afterwards (all 12 FASTA md5 verified against `data/orthobench/md5_inputs.txt`).
