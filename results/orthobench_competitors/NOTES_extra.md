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
  `dev` branch). The loss concentrates in big families: RefOG005/006/008/013 almost entirely missing. v2 run (LUCA + patch) in progress,
  see `OMA_RESUME_NOTE.md`.
