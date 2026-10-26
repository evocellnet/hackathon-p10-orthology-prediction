# ENVIRONMENT

Single env (no conflicts, no aux env needed; `conda_envs/orth_aux` does not exist):
`/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/conda_envs/orth` (6.2 GB, Python 3.11.18, channels conda-forge > bioconda, strict priority)

## Activate
```
source /opt/conda/etc/profile.d/conda.sh && conda activate /scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/conda_envs/orth
export HF_HOME=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/data/models   # HF cache (models already downloaded)
```
Lock files: `envs/orth.yml` (env export), `envs/orth.lock.txt` (explicit, with md5; `mamba create -p X --file envs/orth.lock.txt`), `envs/orth.export.txt` (`mamba list --export`).
Smoke test: `bash scripts/smoke_test_tools.sh` (25 checks, all PASS).

## Versions
Binaries: diamond 2.2.8, mmseqs2 18.8cc5c, foldseek 10.941cd33, mafft 7.526, famsa 2.4.1, FastTree 2.2.0, iqtree 3.1.4, muscle 5.3, hmmer 3.4, trimal 1.5.1, blast 2.17.0, cd-hit 4.8.1, mcl 22.282, orthofinder 3.2, fastoma 0.5.1, omamer 2.1.0, sonicparanoid 2.0.9, proteinortho 6.3.6, generax 2.1.3, aster 1.25 (astral, astral-pro3, astral4, wastral), newick_utils 1.6, snakemake 9.27.1.
Python: numpy 1.26.4, scipy 1.12.0, pandas 2.3.3, numba 0.68.0, scikit-learn 1.9.1, networkx 3.6.1, ete3 3.1.3 (ete4 4.4.0 also present), biopython 1.88, dendropy 5.1.3, pyham 1.1.11, lxml 6.1.3, faiss-cpu 1.14.3, pyarrow 23.0.0, matplotlib 3.11.2, seaborn 0.13.2, tqdm 4.70.1, pytest 9.1.1, python-igraph 1.0.0, leidenalg 0.12.0, pytorch 2.7.1 (cpu), transformers 5.20.2, tokenizers 0.23.3, sentencepiece 0.2.1, fair-esm 2.0.0, tiktoken, protobuf.
(Everything else: see envs/orth.yml.)

## Notes / gotchas
- Nothing skipped, nothing failed. All optional tools installed from bioconda.
- transformers 5.x needs `tiktoken` to load the ProstT5 tokenizer (installed; AutoTokenizer -> T5Tokenizer works, `do_lower_case=False`).
- Models in `data/models` (HF_HOME, 149 MB): facebook/esm2_t6_8M_UR50D, facebook/esm2_t12_35M_UR50D (load + forward OK on CPU). ProstT5: tokenizer only; weights `pytorch_model.bin` are 11.28 GB (fp32) -> NOT downloaded (exceeds 5 GB limit). foldseek itself is available (`foldseek version`); its ProstT5 weights (`foldseek databases ProstT5`) not downloaded yet.
- OMAmer databases (LUCA.h5 etc.) and OMA/FastOMA data are not downloaded.

## Timing calibration (24 threads, self all-vs-all, `scripts/calib_timing.sh`)
Data in `data/calib/`: UP000002311_559292.fasta (S. cerevisiae, 6,066 proteins), UP000005640_9606.fasta (human reference proteome, canonical, 20,652 proteins; UniProt current release). Outputs in `data/calib/out_*/`.
Diamond: `-k 500 -e 1e-3`; mmseqs easy-search: `--max-seqs 500`. Wall seconds (single run; yeast numbers are noisy, sub-10 s):

| mode | yeast 6k | human 20.6k | human hit rows |
|---|---|---|---|
| diamond makedb | 0.0 | 0.1 | |
| diamond --fast | 0.4 | 1.6 | 226,602 |
| diamond --mid-sensitive | 1.9 | 7.0 | 624,354 |
| diamond --sensitive | 3.5 | 12.1 | 756,562 |
| diamond --very-sensitive | 2.4 | 12.0 | 885,417 |
| diamond --ultra-sensitive | 8.7 | 34.4 | 925,098 |
| mmseqs easy-search -s 5.7 | 4.7 | 37.3 | 1,013,792 |
| mmseqs easy-search -s 7.5 | 14.8 | 88.1 | 1,166,448 |

Foldseek: binary present and runs; no structure search timed (no structures/weights downloaded).

## Addendum (S8 gene-loss validation, 2026-10-26)
- Installed into `conda_envs/orth`: `mamba install -p conda_envs/orth -c bioconda -c conda-forge miniprot samtools`
  -> **miniprot 0.18-r281** (bioconda build h577a1d6_0), **samtools 1.24** (h9dcdb79_1). Used by `ancestra/geneloss.py`
  (protein-to-genome alignment; `.fai` random access). tblastn (blast 2.17.0) was available as fallback but not needed.
