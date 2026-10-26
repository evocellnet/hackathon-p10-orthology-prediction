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
