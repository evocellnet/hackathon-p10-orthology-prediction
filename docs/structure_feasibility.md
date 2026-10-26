# 3Di structural-alphabet feasibility (QfO 2022 ~1M proteins, OrthoBench ~250k)

Date: 2026-10-26. Work dir: `data/structure/`. Node: 24 CPUs (AMD EPYC 7402, AVX2 only, no GPU). **The node was heavily shared during all timings (load average 40-90 on 48 logical CPUs; my jobs got 10-16 effective cores)**, so wall-clock numbers are pessimistic; CPU-seconds per residue are the robust quantity.

## Dataset sizes (computed from the files)
| set | proteins | residues | mean len |
|---|---|---|---|
| QfO 2022_02 (78 proteomes, canonical `UP*_taxid.fasta`) | 982,789 | 442.5 M | 450 |
| OrthoBench (12 Ensembl-99 proteomes, `BENCHMARKS/Input`) | 251,378 | 127.8 M | 508 |

## 1. ProstT5 on CPU (Foldseek built-in) -- works, but far too slow for bulk use
- `foldseek databases ProstT5 prostt5_weights tmp` downloads a 2.0 GB `prostt5-f16.gguf` (fp16 weights, 2.4 GB on disk) in 39 s. The weights are already in `data/structure/prostt5_weights/`. `createdb --prostt5-model` works. foldseek 10.941cd33.
- There are no CPU-speed options: only `--threads` (`--gpu` is CUDA-only; the weights are already fp16 GGUF; there is no batch or quantisation flag). So the HF-transformers fallback was not needed. Its fp32 weights are 11 GB, and Zen2 has no bf16 or fp16 arithmetic, so no gain is expected.
- Random human proteins (seed 1), `/usr/bin/time -v`:

| run | proteins | residues | wall | CPU (user+sys) | prot/s | res/s | core-s/res |
|---|---|---|---|---|---|---|---|
| `--threads 24` | 50 | 25,831 | 202 s | 2,672 s (1322% CPU) | 0.25 | 128 | 0.103 |
| `--threads 24` | 67 | 32,442 | 201 s | 3,239 s (1612%) | 0.33 | 161 | 0.100 |
| `--threads 12` | 50 | 25,831 | 214 s | 2,240 s (1046%) | 0.23 | 121 | 0.087 |

- The planned 500- and 2000-protein runs were not completed. The 500-protein run (278k residues) was killed by my 30-minute background timeout after about 30 min, which is consistent with the rate above (about 36 min expected). The 2000-protein run was skipped; the time-box does not allow it and the extrapolation is linear.
- 12 and 24 threads give the same wall time under the current load. More threads do not help on a contended node.
- Cost is about 0.09-0.10 core-seconds per residue, which is about 10 residues per core-second. That matches the FLOP estimate: the ProtT5-XL encoder is about 1.5 B parameters, so about 3 GFLOP per residue.
- Realistic rate on 24 uncontended cores: about 250 res/s (about 0.5 protein/s at 500 aa). The measured contended rate is about 125-160 res/s.

| job | residues | hours at 250 res/s (ideal) | hours at 125 res/s (measured) |
|---|---|---|---|
| OrthoBench 251k proteins | 127.8 M | 142 h (5.9 d) | 284 h (11.8 d) |
| QfO 983k proteins | 442.5 M | 492 h (20.5 d) | 983 h (41 d) |
| generic 250k proteins @500 aa | 125 M | 139 h | 278 h |
| generic 1M proteins @500 aa | 500 M | 556 h (23 d) | 1,111 h (46 d) |
| QfO gap only (95,119 proteins without an AFDB accession) | 92.1 M | 102 h | 205 h |

**ProstT5-on-CPU for the whole of either benchmark is not feasible. Use it only for gaps or small subsets.**

## 2. AlphaFold DB coverage and cost of getting 3Di cheaply
**EBI FTP** (`ftp.ebi.ac.uk/pub/databases/alphafold/latest/`, AFDB v6 = UniProt 2025_03, about 241 M structures):
- 46 per-proteome tars: 16 "proteome" (49.9 GB) + 30 "global_health" (34.7 GB), about 84.6 GB total. Examples: human 4.8 GB, mouse 3.5 GB, yeast 1.0 GB, E. coli 0.46 GB, soybean 7.1 GB.
- Plus Swiss-Prot CIF 37 G / PDB 27 G (550k entries), `msa_depths.csv` 3.5 G, `accession_ids.csv` 8.1 G, `sequences.fasta` 110 G (top level).
- **Only 19 of the 78 QfO proteomes have a tar** (282k of 983k canonical proteins, 44.8 GB). The other 59 proteomes (about 700k proteins) have none.
- Download speed from ftp.ebi.ac.uk (1 GB partial fetch, aborted): **3.6 MB/s single stream** (1 GB in 4 min 58 s). Parallel range requests scale: 16 streams give about 38 MB/s; 24 streams gave 117 MB/s (8.7 GB `accession_ids.csv` in 74 s).
- The 44.8 GB for the 19 QfO tars would take 6-7 min with 24 streams, or about 3.5 h single-stream. There is no advantage over route (a) below for those proteomes.

**Foldseek pre-built DBs** (`foldseek databases`, Cloudflare R2, v6 2025_03). Single stream: 37.5 MB/s. aria2 with 5 connections (as foldseek does): 112 MiB/s.

| DB | tar.gz size | est. download at 112 MB/s | contents |
|---|---|---|---|
| `Alphafold/Proteome` | **1.8 GB** (2.8 GB unpacked) | **69 s measured** | 590,435 entries, 46 proteomes, AA + 3Di (`_ss`, 259 MB) + CA + taxonomy |
| `Alphafold/Swiss-Prot` | 1.6 GB | 15 s | Swiss-Prot only |
| `Alphafold/UniProt50` | 123 GB | 18 min | cluster reps only, not suitable |
| `Alphafold/UniProt` | **491 GB** | about 75 min (3.6 h single stream) | all ~241 M entries. First tar member `afdb` (AA) = 79.8 GB, so `afdb_ss` ~80 GB, plus index/lookup/_h. |

- `Alphafold/UniProt` is gzip, so it is not seekable. It must be streamed once, e.g. `aria2c -x16 ... | tar -xz --wildcards 'afdb_ss*' 'afdb.lookup' 'afdb.index' 'afdb_h*'`. Total time is about 1.5-2 h, with 0.1 TB of disk kept and 0.4 TB transient. **Not run** (the brief was to measure only).
- Per-accession HTTPS (`https://alphafold.ebi.ac.uk/files/AF-<acc>-F1-model_v6.pdb`, **v6 is current; v4 URLs give 403**) needs a non-curl User-Agent; the default curl UA gets 403.
  - Rate: 16 parallel connections gave about 60 files/s (117 files in 2 s, 320 KB per PDB on average). 64 parallel connections gave HTTP 429 for 90% of requests.
  - 623k files (the QfO accessions not in the Proteome DB) would take about 2.9 h and about 200 GB, plus `foldseek createdb` on the PDBs (about 250 structures/s on 4 threads). About 2.5% of requests returned 404.
  - The API `https://alphafold.ebi.ac.uk/api/prediction/<acc>` works (0.1 s per call), but is not needed.
- Coverage of QfO 2022 accessions in AFDB v6 (`accession_ids.csv` lookup of the 982,789 QfO accessions):
  - **887,670 (90.3%)** have an AFDB entry.
  - 95,119 proteins (92.1 M residues) do not. They are concentrated in UP000000539 (chicken, 66% missing), UP000000437 (zebrafish, 54%) and UP000009136 / UP000006906 (about 19%).
  - These are probably accessions since removed or merged in UniProt. An exact-sequence match against the full AFDB AA file may recover some of them. This is untested.
- For accessions that are in AFDB, QfO and AFDB sequences are identical in **99.4%** of the cases (262,660 of 264,149 in the Proteome DB; the rest are mostly isoform or version drift). The identity test must be applied.
- Coverage in the 1.8 GB Proteome DB: **264,149 QfO proteins (27%)** across the 19 proteomes, of which 262,660 have an identical sequence. Extracted to `data/structure/qfo_3di_afdbproteome.fa` (`>UPID|acc`, 3Di string, F1 fragment only; proteins >2700 aa are fragmented in AFDB, so F1 covers only part of them).
- OrthoBench uses Ensembl-99 IDs, not UniProt. Exact-sequence match against the Proteome DB gives 97,647 of 251,378 (39%): worm 18,621/20,191, fly 13,024/13,947, human 19,283/23,485, mouse 18,264/22,866, rat 11,518/22,250, zebrafish 13,064/30,313, chimp 3,277/23,534, and about 0 for dog, chicken, opossum, tetraodon and ciona (no proteome in the DB).
  - Extracted to `data/structure/orthobench_3di_exactmatch.fa`.
  - The remaining 8 species and the unmatched isoforms need either an exact-sequence lookup against the full AFDB AA DB (e.g. the 79.8 GB `afdb` file, hashed with `foldseek convert2fasta` + md5) or ProstT5.
  - The fraction recoverable via the full DB is unknown; UniProt reference proteomes of dog, chicken, etc. are Ensembl/NCBI-derived, so a high fraction is plausible but unverified.

## 3. Sanity check: predicted (ProstT5) vs true 3Di (AFDB v6 structures)
- 117 human proteins (random, from the 200-protein sample) were fetched as v6 PDBs (`afdb50_pdb/`).
- True 3Di comes from `foldseek createdb` on the PDBs. It is **identical to the Alphafold/Proteome DB `_ss` for 117/117**, so the DB route is equivalent to computing 3Di from the structures. All 117 sequences equal the AFDB sequences.
- ProstT5 predictions (`pred50_ss.fasta`, `pred100_ss.fasta` (the 67-protein batch; file name is historical); same length as the protein) were compared per residue (scripts `cmp.py`, `cmp2.py`):
  - 117 proteins: mean per-protein identity **0.691** (median 0.697, min 0.35, max 0.93). Residue-weighted identity **0.688**. The first 50 proteins alone gave 0.700 / 0.696.
  - Baseline (always predict the most common state `D`): 0.356.
  - By pLDDT bin:

| pLDDT | residues | identity |
|---|---|---|
| <50 | 17,340 | 0.69 |
| 50-70 | 5,614 | 0.57 |
| 70-90 | 10,936 | 0.67 |
| >=90 | 24,383 | 0.72 |

  - Per-state recall varies from 0.86 (K) and 0.84 (V) to 0.35 (N) and 0.38 (H).
- Predicted 3Di is a noisy proxy (about 70% identity even in confident regions). Real 3Di should be preferred wherever an AFDB entry exists.
- ProstT5 predicted 3Di also lacks Cα coordinates (no TM-score/LDDT).

## Recommended route
1. **Right now (done, 1.8 GB / 69 s):** use `Alphafold/Proteome` 3Di for the 19 QfO proteomes (262,660 proteins) and the 6 OrthoBench species by exact-sequence match (97,647 proteins). The extracted FASTAs are in `data/structure/`.
2. **For the remaining ~625k QfO proteins with an AFDB accession (and the other OrthoBench species):** stream `Alphafold/UniProt` (491 GB, about 75-120 min at 112 MB/s, about 80 GB of `_ss` kept) and pull 3Di by accession through `afdb.lookup`. This is cheaper and more robust than 623k per-file HTTPS downloads (about 3 h, rate-limited at 429, 200 GB). It is not yet downloaded; it needs disk of about 0.5 TB transient plus a decision from the PI.
   - For OrthoBench, hash AFDB AA sequences and match Ensembl proteins by exact sequence.
3. **ProstT5 only for the gap:** about 95k QfO proteins (92 M residues, 4-8.5 days at 24 cores) plus unmatched OrthoBench proteins. Subsample it, or restrict it to proteins that matter for the analysis (e.g. OrthoBench RefOG members and their orthogroup neighbours; those are about 1-2k proteins, so hours). Re-matching these proteins to a current AFDB entry by sequence identity might shrink the gap.
4. Do not plan on ProstT5 for all 1M proteins (20-46 days) or all OrthoBench (6-12 days) without a GPU.

## Files in `data/structure/`
`prostt5_weights/` (2.4 GB), `afdb_proteome_db/` (foldseek Alphafold/Proteome v6 + converted FASTAs), `qfo_fasta/` (the 78 canonical QfO FASTAs, extracted), `qfo_accessions.json`, `qfo_afdb_cov.json`, `qfo_missing_afdb.fasta` (95,119 proteins), `qfo_3di_afdbproteome.fa`, `orthobench_3di_exactmatch.fa`, `afdb_download_metadata.json`, `afdb_latest_listing.txt`, `sanity50.fa` / `sanity_next.fa` / `sanity150.fa`, `pred*_ss.fasta` (ProstT5 3Di), `true3di_ss.fasta`, `time_*.txt` (`/usr/bin/time -v` output), `cmp.py`, `cmp2.py`, `cov.py`, `cov2.py`, `extract_3di_from_proteome_db.py`.

---

# 4. Foldcomp prebuilt databases, and the final 3Di tables (added 2026-10-26)

Work dir: `data/structure/foldcomp/`. Reproduce everything with `scripts/build_3di.sh` (helpers: `scripts/fc3di.py`, `scripts/get_foldseek_afdb_uniprot.sh`, `data/structure/foldcomp/download_v6.sh`; about 1 h when the big files are on disk, plus about 1.7 h and 0.5 h of downloads). **This supersedes "Recommended route" steps 1-3 above: AFDB-wide 3Di is now on disk, and QfO is 94.5 % and OrthoBench 80.5 % covered with exact sequence matches.**

## 4.1 Inventory of the prebuilt Foldcomp DBs
Server: `https://foldcomp.steineggerlab.workers.dev/<name>` (a Cloudflare worker that 302-redirects to `https://steineggerlab.s3.amazonaws.com/foldcomp/<name>`; the HTML listing is incomplete, e.g. it lacks v6 and the organisms). A DB = `<name>` (data) + `.index` + `.lookup` + `.dbtype`. Entry counts are index line counts (computed). Foldcomp skipped structures with discontinuous residues or other issues, so counts are slightly below AFDB's.

| DB (`foldcomp.setup(name)`) | data size | entries | content / version | upload |
|---|---|---|---|---|
| **afdb_uniprot_v6** | **1.26 TB** (+6.7 GB index, 10.8 GB lookup) | **241,070,489** | AFDB v6 (same 241,070,489 entries as Foldseek's v6 = UniProt 2025_03), all of UniProt, fragments F1..Fn; ids `AF-<acc>-F1-model_v6` | 2026-08-26 |
| afdb_uniprot_v4 | 1.10 TB (+5.6 GB index) | 214,683,829 | AFDB v4 (UniProt 2021_04) | 2025-02-27 |
| afdb_swissprot_v4 | 3.05 GB | 542,378 | Swiss-Prot, v4 | 2025-02-26 |
| afdb_rep_v4 | 10.2 GB | 2,278,845 | AFDB cluster representatives (v4) | 2025-02-26 |
| afdb_rep_dark_v4 | 2.3 GB | 705,931 | dark-cluster representatives (v4) | 2025-02-26 |
| a_thaliana / c_albicans / c_elegans / d_discoideum / d_melanogaster / d_rerio / e_coli / g_max / **h_sapiens** / m_jannaschii / m_musculus / o_sativa / r_norvegicus / s_cerevisiae / s_pombe / z_mays | 177 / 47 / 126 / 104 / 106 / 202 / 21.8 / 347 / 234 / 8.2 / 173 / 215 / 166 / 46.6 / 37.7 / 245 MB | 27,434 / 5,974 / 19,694 / 12,622 / 13,458 / 24,664 / 4,363 / 55,799 / 23,391 / 1,773 / 21,615 / 43,649 / 21,270 / 6,039 / 5,128 / 39,299 | the 16 AFDB model-organism proteomes, **v4 models** (`AF-..-model_v4.pdb` names) | 2025-08-14 |
| highquality_clust30 | 120 GB (+0.9 GB index) | 37,192,348 | ESM Atlas high-quality, 30 % clusters (metagenomic, MGYP ids) | 2025-03-17 |
| esmatlas_v2023_02 / esmatlas | 137 GB / 1.82 TB | not counted (index 1.2 GB / 15.9 GB) | ESM Atlas metagenomic (MGYP), not relevant for QfO/OrthoBench | 2023-03 / 2026-10 |

Also on the server: `afdb_swissprot*` legacy files and `e_coli_folddisco.tar.gz`. Only the 16 organisms above exist (no dog, chicken, etc.); there is no v6 Swiss-Prot or organism DB.

**Best coverage:** only `afdb_uniprot_v6` is useful. The organism / Swiss-Prot / rep DBs are v4 or a tiny subset: human v4 has 20,431 of 20,598 QfO accessions, but its 3Di agrees with the v6 3Di at only 89.8 % per residue (0.3 % of proteins identical), because the v4 models differ from the v6 ones. The 16 organism DBs cover only 6 of the 12 OrthoBench species (human, mouse, rat, zebrafish, worm, fly), so they were not used.

## 4.2 Coverage before / after (only proteins whose AFDB-v6 sequence is IDENTICAL to ours)
| set | previous (Foldseek `Alphafold/Proteome`, exact acc) | Foldcomp v6 only | **best = Foldseek `Alphafold/UniProt` v6, exact sequence** |
|---|---|---|---|
| QfO 2022_02, 982,789 proteins | 262,660 (26.7 %) | 923,018 (93.9 %) | **928,383 (94.5 %)** |
| OrthoBench, 251,378 proteins | 97,647 (38.8 %) | 190,770 (75.9 %) | **202,475 (80.5 %)** |

- QfO: 916,429 matched through their own UniProt accession, 11,954 through another AFDB entry with the same sequence (merged/renamed accessions). 74 of 78 proteomes have >= 90 %. Lowest: chicken UP000000539 37 % (18,112 proteins), cow UP000009136 72 %, zebrafish UP000000437 80 %, Chlamydomonas 85 %. **Gap: 54,406 proteins, 75.6 M residues.**
- OrthoBench per species (before -> after, of n): worm 18,621 -> 18,920 /20,191; fly 13,024 -> 13,300 /13,947; human 19,283 -> 22,123 /23,485; mouse 18,264 -> 21,451 /22,866; rat 11,518 -> 16,028 /22,250; zebrafish 13,064 -> 16,177 /30,313; chimp 3,277 -> 22,137 /23,534; dog 359 -> 12,016 /20,257; chicken 84 -> 5,398 /16,878; opossum 126 -> 20,125 /21,384; tetraodon 26 -> 18,408 /19,602; ciona 1 -> 16,392 /16,671. **Gap: 48,903 proteins, 38.9 M residues**, mostly non-canonical isoforms (AFDB models only the UniProt canonical sequence), plus dog/zebrafish/chicken/rat Ensembl-only proteins.
- Ensembl -> UniProt mapping for OrthoBench is by exact sequence, as in the first probe. The best table matches against the sequences of the whole AFDB-v6 AA DB (241 M entries), so no UniProt id mapping is needed. When several AFDB entries have the same sequence, UniProt accessions that carry the same sequence in that species' UniProt proteome are preferred (`source.tsv` column `match_kind`: `uniprot_candidate` / `seq_other_acc`).
- The Foldcomp route needs the candidate UniProt accessions from UniProt REST (`taxonomy_id:<tax>` FASTAs in `data/structure/uniprot_taxon_fasta/`, release 2026_03) plus `sec_ac.txt` (merged accessions); 246,888 candidate accessions for OrthoBench, +7,210 merged QfO accessions.
- Remaining gap options: ProstT5 on CPU for 75.6 M (QfO) residues is 84-168 h on 24 cores (section 1); restrict to proteins that matter, or for OrthoBench use the canonical isoform's 3Di for its isoforms (identical for the shared domains; not done).

## 4.3 Verification
1. **Sequence identity.** The AA of every kept entry equals our proteome sequence (exact match only). Foldcomp v6: 924,375 of 982,789 QfO accessions exist as F1 entries (the previous probe's `accession_ids.csv` lookup gave 887,670; the difference was not investigated), 923,018 of them have an identical sequence (99.85 %), +6,589 via merged accessions.
2. **3Di identity.** (a) `Alphafold/UniProt` 3Di equals `Alphafold/Proteome` 3Di for 262,660/262,660 shared QfO entries, i.e. the true 3Di. (b) **Foldcomp-derived 3Di is NOT identical to the true 3Di.** Foldcomp stores a lossy backbone (torsion quantisation, about 0.1 A). After `foldseek createdb --input-format 5` and comparison with the true 3Di of the same AFDB entry: mean per-residue identity **97.9 %** (median 98.0 %, residue-weighted 97.9 %), only 4.7 % of proteins fully identical, 98.7 % of proteins >= 95 % (QfO, 922,763 entries; OrthoBench, 190,770 entries: 97.97 % mean, 99.1 % >= 95 %). The error is uniform over pLDDT (97.5-98.8 % in every bin of a 979-protein human check), so it is coordinate noise, not low-confidence regions. This is still far better than ProstT5 (69 %). Use the Foldseek DB where possible; the Foldcomp tables (`*_3di_foldcomp.fa`) are a cross-check or fallback, and are a strict subset of the best tables (0 proteins only in Foldcomp).
3. The Foldcomp data file passes `index max(offset+length) == file size` (1,261,454,421,092 B) and 971,288 entries decoded without error.

## 4.4 Speeds, tooling, gotchas
- Download speed from S3 (aria2c -x16): 150-250 MiB/s alone, 110-130 MiB/s while competing with another download; `afdb_uniprot_v6` (1.26 TB) took 1 h 47 min (avg 186 MiB/s). Foldseek `Alphafold/UniProt` (491 GB tar.gz): 39 min at about 200 MiB/s. `foldseek databases` then extracts with a single-threaded `tar xvfz` at 30-90 MB/s (> 3 h), so `scripts/get_foldseek_afdb_uniprot.sh` uses `rapidgzip -P 12` (about 340 MB/s; extraction about 35 min), excluding `afdb_ca*` and `afdb_h*`.
- `foldseek createdb <foldcomp_db> <out> --input-format 5` reads a Foldcomp DB directly (23,391 human entries in 22 s on 8 threads). For a subset use `foldseek createsubdb --id-mode 0 --subdb-mode 1 keys.txt <big_db> <sub>` (soft link, 23 s for 971k keys; `--id-mode 1` fails on the v6 names because the organism DBs carry a `.pdb` suffix and v6 does not) and `createdb` on that: 971,288 entries in 8 min 43 s on 16 threads (random reads from the 1.26 TB file). Foldseek writes AA (`<out>`), 3Di (`<out>_ss`), Calpha and headers; `fc3di.py` reads them directly.
- Python `foldcomp` 1.0.0 (`pip install foldcomp`) was only used for the pLDDT check. S3 range requests for single entries work (about 70 entries/s with 32 threads, `scripts/fc_range_fetch.py`) but are slower than the bulk download for 971k entries.
- The SHA-256 of the small files and the S3 ETags are below; there is no end-to-end checksum for the 1.26 TB file other than the S3 multipart ETag.

## 4.5 DATA.md entries
| item | source URL | version / date | size | checksum | license |
|---|---|---|---|---|---|
| Foldcomp `afdb_uniprot_v6` (+`.index`, `.lookup`, `.dbtype`) | `https://foldcomp.steineggerlab.workers.dev/afdb_uniprot_v6` -> `https://steineggerlab.s3.amazonaws.com/foldcomp/afdb_uniprot_v6` | AFDB v6 (entry count identical to Foldseek v6 = UniProt 2025_03), uploaded 2026-08-26 | 1,261,454,421,092 B; index 6,685,778,582 B; lookup 10,779,949,313 B | S3 ETag data `02ea7fb0a0215a68ea5f77a2226dda0f-9399`, index `a05aaf24d7c731e6c23b13c713151b18-50`, lookup `ce0842989cc198b94e06c71406d1875f-81`; sha256 index `4d4a0053be1b8e450bcd6af6118aac6205cfff19d0276fc1f384cfcca77286d2`, lookup `30f3cfbcf9c7f82bb50932967f6340976bfe3e899814d3d64b6e0ad67a16f037` | AFDB data CC-BY-4.0 (DeepMind / EMBL-EBI; cite Jumper 2021, Varadi 2022/2024); Foldcomp tool GPL-3.0 (Kim et al. 2023, Bioinformatics) |
| Foldseek `Alphafold/UniProt` (`afdb.tar.gz`) | `https://foldseek.steineggerlab.workers.dev/afdb.tar.gz` (+`afdb.version`) | AFDB v6, UniProt 2025_03; file dated 2025-10-08 | 491,470,141,753 B (extracted AA 79.8 GB, 3Di 79.8 GB, index 6.1 GB, lookup 10.8 GB, mapping 4.0 GB, taxonomy) | S3 ETag `5fd3e5f7a8750b5f5b2a7b91ed79bf49-7324` | CC-BY-4.0 (data); Foldseek GPL-3.0 |
| Foldseek `Alphafold/Proteome` | see section 2 | v6 | 1.8 GB | -- | CC-BY-4.0 |
| UniProt taxon FASTAs (12 OrthoBench taxa), `sec_ac.txt` | `https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:<tax>&format=fasta` ; `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt` | UniProt 2026_03 (2 Sep 2026) | 8-48 MB each; 49.5 MB | -- | CC-BY-4.0 |

Outputs (generated, `data/structure/`): `qfo_3di_best.fa` (`>UPID|acc`, 928,383 entries, 387 MB), `orthobench_3di_best.fa` (`>ensembl_id|uniprot_acc`, 202,475 entries, 95 MB), each with `.source.tsv` (source, matched accession, match kind) and `.coverage.tsv` (per proteome / species); `qfo_3di_foldcomp.fa`, `orthobench_3di_foldcomp.fa` (Foldcomp-derived only; approximate 3Di). Intermediates in `foldcomp/work/`.

## 4.6 Download status and how to finish / free disk
- Both downloads are **complete**: `data/structure/foldcomp/v6/afdb_uniprot_v6` (1.26 TB, no `.aria2` file) and `data/structure/afdb_uniprot_fs/ex/` (extracted Foldseek UniProt DB, 175 GB; `ex/.extracted` marks completion; the 491 GB tar.gz was deleted after extraction).
- If a download must be resumed: `data/structure/foldcomp/download_v6.sh` (resumable `aria2c -c`; done when `afdb_uniprot_v6.aria2` is gone and the size is 1261454421092) and `scripts/get_foldseek_afdb_uniprot.sh`.
- `foldcomp/v6/afdb_uniprot_v6` (1.2 TB) is only needed for the Foldcomp fallback / coordinates. It can be deleted without affecting `*_3di_best.fa`; the Foldseek DB `afdb_uniprot_fs/ex/` (AA + 3Di for all 241 M entries) is the one to keep for further lookups.
