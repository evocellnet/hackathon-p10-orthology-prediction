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
| fastoma_v2 | 29.9 | 26.7 | 33.1 | 21.6 | 48.6 | 19/70 | 2.99 | 19.9 | 12.2 |
| fastoma_v2_luca_released | 29.1 | 25.1 | 33.1 | 21.1 | 47.1 | 19/70 | 2.99 | 19.9 | 12.2 |
| oma_browser | 55.8 | 50.5 | 60.4 | 92.5 | 39.9 | 6/70 | 0.26 | 0.3 | 1.5 |
| oma_browser_root | 7.3 | 7.3 | 7.2 | 3.9 | 47.1 | 6/70 | 0.26 | 0.3 | 1.5 |

Precision/recall on dev / test:

| method | P dev | R dev | P test | R test |
|---|---|---|---|---|
| orthofinder | 60.8 | 84.0 | 62.1 | 79.5 |
| orthofinder_mcl | 39.9 | 84.6 | 54.1 | 79.5 |
| sonicparanoid_default | 24.0 | 59.3 | 70.8 | 70.1 |
| proteinortho | 96.1 | 27.4 | 97.7 | 31.4 |
| fastoma_luca | 23.6 | 39.8 | 24.3 | 53.7 |
| fastoma_luca_ogs | 95.2 | 16.3 | 95.0 | 21.7 |
| fastoma_v2 | 19.3 | 43.4 | 23.9 | 53.8 |
| fastoma_v2_luca_released | 18.2 | 40.4 | 23.9 | 53.8 |
| oma_browser | 97.8 | 34.1 | 88.9 | 45.8 |
| oma_browser_root | 4.0 | 44.8 | 3.9 | 49.5 |

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
| fastoma_v2 | all | -40.2 | [-54.6, -15.3] | 10.1 |
| fastoma_v2 | dev | -43.8 | [-61.3, -13.2] | 14.3 |
| fastoma_v2 | test | -36.6 | [-56.2, -4.0] | 13.9 |
| fastoma_v2_luca_released | all | -41.0 | [-55.1, -16.1] | 10.0 |
| fastoma_v2_luca_released | dev | -45.4 | [-62.4, -14.6] | 13.9 |
| fastoma_v2_luca_released | test | -36.6 | [-56.2, -4.0] | 13.9 |
| oma_browser | all | -14.3 | [-22.7, -5.0] | 4.5 |
| oma_browser | dev | -20.0 | [-30.0, -9.1] | 5.3 |
| oma_browser | test | -9.3 | [-23.2, +5.5] | 7.1 |
| oma_browser_root | all | -62.9 | [-71.4, -45.7] | 6.3 |
| oma_browser_root | dev | -63.2 | [-74.1, -42.6] | 8.3 |
| oma_browser_root | test | -62.5 | [-74.2, -11.3] | 14.9 |

## Notes

- **orthofinder**: OrthoFinder 3.2 (conda, bioconda) ; diamond 2.2.8; famsa 2.4.1; FastTree 2.2.0; mcl 22.282 | orthofinder -f <12 proteomes> -t 8 -a 8 -mog (default: DIAMOND, MSA/famsa+FastTree gene trees, MCL I=1.2); node load avg ~20-55 on 24 cores
- **orthofinder_mcl**: same run as orthofinder (MCL_Orthogroups.tsv before gene-tree splitting) | orthofinder -f <12 proteomes> -t 8 -a 8 -mog (default: DIAMOND, MSA/famsa+FastTree gene trees, MCL I=1.2); node load avg ~20-55 on 24 cores
- **sonicparanoid_default**: SonicParanoid 2.0.9 (bioconda), bundled diamond; mode default | sonicparanoid -m default -t 8 -I 1.5 (diamond very-sensitive, domain-arch step on). Includes one-time Pfam profile DB indexing at first run after install (few minutes). Node load avg ~20-55 on 24 cores
- **proteinortho**: Proteinortho 6.3.6 (bioconda) diamond 2.2.8 | proteinortho -p=diamond -cpus=8 (defaults: e=1e-5, identity 25, coverage 50, conn 0.1, purity; PoFF off). Input proteomes had internal '*' stop symbols stripped (sed, required by ProteinOrtho)
- **fastoma_luca**: FastOMA 0.5.1, omamer 2.1.0, LUCA.h5 (zenodo 20814376, 2026-06-23); species tree hand-built from NCBI taxonomy; variant=tophogs | FastOMA 0.5.1 nextflow pipeline (local, no container), OMAmer 2.1.0 + LUCA.h5 (OMAmer DB, OMA-derived: LEAKAGE - OMA contains these species), nextflow executor cpus=8. wall includes omamer search (12 x 1 cpu), nextflow overhead + report notebook; peak RSS = largest single process
- **fastoma_luca_ogs**: FastOMA 0.5.1, omamer 2.1.0, LUCA.h5; variant=OrthologousGroups | FastOMA 0.5.1 nextflow pipeline (local, no container), OMAmer 2.1.0 + LUCA.h5 (OMAmer DB, OMA-derived: LEAKAGE - OMA contains these species); same run as fastoma_luca; peak RSS = largest single process
- **fastoma_v2**:  | FastOMA 0.5.1 + LUCA.h5 + gene-loss patch; 3 segments (pause/resume; segment 3 re-ran all hog_rest because nextflow cache of /tmp workdirs was not reused: 76 tasks run, 14 cached); omamer step reused v1 hogmaps (add ~0.75 CPU-h, 13.6 GB RSS from v1 trace for a full run); local NVMe work dir; LEAKAGE
- **fastoma_v2_luca_released**:  | FastOMA 0.5.1 + LUCA.h5 + gene-loss patch; 3 segments (pause/resume; segment 3 re-ran all hog_rest because nextflow cache of /tmp workdirs was not reused: 76 tasks run, 14 cached); omamer step reused v1 hogmaps (add ~0.75 CPU-h, 13.6 GB RSS from v1 trace for a full run); local NVMe work dir; LEAKAGE
- **oma_browser**: OMA browser All.May2026 (Zenodo 10.5281/zenodo.22869443, published 2026-10-21; HOGs 'infered May 2025'), oma-hogs.orthoXML.gz md5 cf41ae36b554f6a60270297f6ab95d66, CC BY 4.0; Bilateria-level HOGs restricted to OrthoBench genes; ID mapping oma-ensembl.txt.gz (md5 86a09a24...) + exact/near-identical sequence (oma-seqs.fa.gz md5 417fa34d...); unmapped = singletons | PRECOMPUTED OMA browser HOGs (All.May2026, inferred by OMA on its own infrastructure; inference cost not measurable/not included). Runtime = local conversion only (stream oma-hogs.orthoXML.gz, ID/sequence mapping incl. DIAMOND near-identical step on 2 threads, OG export); Bilateria-level HOGs. LEAKAGE: OMA contains all 12 OrthoBench species.
- **oma_browser_root**: OMA browser All.May2026 (Zenodo 10.5281/zenodo.22869443, published 2026-10-21; HOGs 'infered May 2025'), oma-hogs.orthoXML.gz md5 cf41ae36b554f6a60270297f6ab95d66, CC BY 4.0; root-level HOGs restricted to OrthoBench genes; ID mapping oma-ensembl.txt.gz (md5 86a09a24...) + exact/near-identical sequence (oma-seqs.fa.gz md5 417fa34d...); unmapped = singletons | PRECOMPUTED OMA browser HOGs (All.May2026, inferred by OMA on its own infrastructure; inference cost not measurable/not included). Runtime = local conversion only (stream oma-hogs.orthoXML.gz, ID/sequence mapping incl. DIAMOND near-identical step on 2 threads, OG export); root-level HOGs. LEAKAGE: OMA contains all 12 OrthoBench species.

## Caveats and provenance

- **Leakage (FastOMA):** the OMAmer LUCA database (`LUCA.h5`, Zenodo 20814376, 2026-06-23) is derived from OMA, which contains these Ensembl species and proteins; FastOMA's placement of genes into root HOGs is therefore not independent of the reference. Treat FastOMA as an informed upper-bound-style comparator, not a clean one. FastOMA assigns only 198,615 of 251,378 proteins (OMAmer leaves the rest unplaced); unplaced genes count as missing (FN) in the benchmark. `fastoma_luca` = top-level HOGs of `FastOMA_HOGs.orthoxml` (identical to `RootHOGs.tsv`); `fastoma_luca_ogs` = FastOMA `OrthologousGroups.tsv` (strict 1-per-species-style ortholog groups, deliberately fine-grained: very high precision, low recall) and is not the OrthoBench-appropriate granularity. The species tree given to FastOMA was hand-built from NCBI taxonomy (root = Bilateria).
- **Coverage:** OrthoFinder/SonicParanoid/ProteinOrtho/FastOMA outputs contain only genes placed in a group of size >=2; singletons are simply absent (equivalent to singletons for the pairwise score).
- **Timing:** runs were executed one at a time with 8 threads on a shared 24-core node (load average 20-55 during all runs, see `<method>/load.log`, uptime every 60 s), so wall times are 1.5-3x what an idle node would give; CPU = user+sys from `/usr/bin/time -v`. `peak RSS` = `Maximum resident set size` of the largest process in the tree (for FastOMA/nextflow this is the largest single task, not the sum). FastOMA wall includes nextflow orchestration and the final report notebook; LUCA.h5 download (9.9 GB) and a failed first attempt (species-tree typo, 3 s) are excluded.
- **OrthoFinder 3.2** default = MSA gene trees (famsa + FastTree), DIAMOND, MCL I=1.2, then HOG-derived `Orthogroups.tsv` (`orthofinder`). `orthofinder_mcl` scores the pre-tree MCL orthogroups from the same run (`-mog`).
- **ProteinOrtho** needed proteome files with internal `*` stop symbols removed (sed on sequence lines only; headers untouched, same protein count).
- Scoring: `scripts/score_competitor.sh` -> `scripts/to_orthobench_og.py` / `orthoxml_to_orthobench_og.py` -> `run_orthobench.sh` (official benchmark.py) and `orthobench_eval.py` (per-RefOG, dev/test). Raw outputs: `data/orthobench/runs/<method>/` (not committed).
- **SonicParanoid2 `-m sensitive` NOT completed:** MMseqs2 s=7.5 all-vs-all needed roughly 5-6 h on the shared node (48/144 alignments after ~2 h, in two attempts); both attempts were stopped (first to prioritise FastOMA, second on PI request because the node is needed for scaling runs). No score is reported for it. Partial outputs: `data/orthobench/runs/sonicparanoid_sensitive_aborted_*`.
- Only OrthoBench `Input/*.fa` were read; a failed ProteinOrtho first attempt had written three `*.fa.len` files into `Input/`, removed afterwards (all 12 FASTA md5 verified against `data/orthobench/md5_inputs.txt`).
- **oma_browser / oma_browser_root (added at user request, post-freeze, 2026-10-27):** precomputed OMA browser HOGs, release All.May2026
  (Zenodo 10.5281/zenodo.22869443, CC BY 4.0), restricted to OrthoBench genes; `oma_browser` = HOGs at the Bilateria level (primary,
  pre-declared in `oma_browser/PREREGISTRATION.txt`), `oma_browser_root` = root-level HOGs (deep LUCA/Eukaryota families -> huge groups,
  precision 3.9). Mapping Ensembl 99 -> OMA: 95.1% of proteins mapped, 89.2% to an entry in a HOG (RefOG genes 93.2%; chicken only 60%
  because OMA uses RefSeq GRCg7b), see `oma_browser/mapping_stats.txt`; unmapped = singletons. LEAKAGE: OMA contains all 12 species.
  Runtime = local conversion only (OMA inference itself not measurable). Scripts: `oma_hogs_paths.py`, `oma_browser_map_orthobench.py`,
  `oma_browser_map_nearident.py`, `oma_browser_to_orthobench_og.py`, `oma_comparators_vs_ancestra.py`.
- **FastOMA v1 diagnosis (2026-10-27):** inputs (12 proteomes, species tree, names) and export level were correct (top-level OrthoXML HOGs ==
  RootHOGs.tsv == Bilateria-level HOGs). OMAmer leaves only 4.1% unplaced (Metazoa.h5: 4.2%, so a more specific DB does not help). The 21%
  absent genes = 9.1% never in a FastOMA family + 9.7% root-level singletons + 2.2% (5,497 genes) silently LOST by a FastOMA 0.5.1 bug:
  `_hog_class.split_hog()` rebuilds a split HOG only from the subHOGs that contain a representative and drops the others (also on the
  `dev` branch). The loss concentrates in big families: RefOG005/006/008/013 almost entirely missing. The v2 run (LUCA + patch) is complete: see the section below.

## OMA comparators (added at user request, post-freeze, 2026-10-27)

ΔF = F(Ancestra accurate) − F(comparator), paired bootstrap over RefOGs (B=2000, seed 11), NI margin 3 pp. Rows: `fastoma_v2` = FastOMA 0.5.1 + LUCA.h5 + gene-loss patch (pre-declared primary); `fastoma_v2_luca_released` = same run minus rescued HOGs (= unpatched 0.5.1); `oma_browser` = precomputed OMA browser HOGs (All.May2026) cut at Bilateria (primary), `oma_browser_root` = root level; `fastoma_luca` = v1 (unpatched) run.

| comparator | test F | P | R | exact | dev F | all-70 F | ΔF test (Ancestra − comparator) [95% CI] | SE | NI verdict (Ancestra vs comparator) | ΔF all-70 [95% CI] | CPU-h | wall-h | peak RSS GB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fastoma_v2 | 33.1 | 23.9 | 53.8 | 10/36 | 26.7 | 29.9 | +43.1 [+14.0, +55.3] | 11.1 | comparable (non-inferior; CI > 0) | +48.0 [+24.3, +60.6] | 19.86 | 2.99 | 12.2 |
| fastoma_v2_luca_released | 33.1 | 23.9 | 53.8 | 10/36 | 25.1 | 29.1 | +43.1 [+14.0, +55.3] | 11.1 | comparable (non-inferior; CI > 0) | +48.8 [+25.3, +61.2] | 19.86 | 2.99 | 12.2 |
| oma_browser | 60.4 | 88.9 | 45.8 | 3/36 | 50.5 | 55.8 | +15.7 [+1.8, +29.2] | 7.2 | comparable (non-inferior; CI > 0) | +22.1 [+13.2, +30.7] | 0.31 | 0.26 | 1.5 |
| oma_browser_root | 7.2 | 3.9 | 49.5 | 2/36 | 7.3 | 7.3 | +69.0 [+22.3, +75.3] | 12.7 | comparable (non-inferior; CI > 0) | +70.6 [+54.6, +77.2] | 0.31 | 0.26 | 1.5 |
| fastoma_luca | 33.5 | 24.3 | 53.7 | 10/36 | 29.6 | 31.7 | +42.7 [+14.1, +55.2] | 11.2 | comparable (non-inferior; CI > 0) | +46.2 [+26.3, +57.7] | 12.74 | 3.50 | 13.6 |

- **LEAKAGE:** OMA (browser HOGs, and the OMAmer LUCA DB used by FastOMA) contains all 12 OrthoBench species/proteins, so these are informed upper-bound-style comparators, not independent ones. Ancestra did not use OMA.
- **Patch verification:** genes in predicted groups rose to 202,691/251,378 (80.6%) vs 198,615 (79.0%) in v1 (unpatched re-run 199,302); 3,795 subHOGs (630 of size>=2, 7,464 genes) were rescued. Big families partly recovered (RefOG005 TP 0.03 -> 3.8, RefOG006 0.06 -> 4.3, RefOG008 0.08 -> 4.3; RefOG013 unchanged). Test F is essentially unchanged (33.1 vs 33.5 v1): the remaining FastOMA gap is not the gene-dropping bug (most absent genes are OMAmer-unplaced / root singletons).
- **Cost (fastoma_v2, 3 segments, 8 threads, shared node):** 19.9 CPU-h, 2.99 h wall, peak RSS 12.2 GB (+~0.75 CPU-h omamer step reused from v1). oma_browser: conversion only, 0.31 CPU-h, 0.26 h, 1.5 GB (OMA inference not measurable). Both are not comparable to Ancestra's 3.9 CPU-h.
- Reading: Ancestra's point estimate is above every OMA row on test; vs oma_browser ΔF +15.7 [+1.8, +29.2] (CI above 0), with the caveat that OMA saw these species.
