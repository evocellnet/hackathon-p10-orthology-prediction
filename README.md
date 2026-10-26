# Ancestra

Ancestra infers **hierarchical orthologous groups (HOGs)** from a set of proteomes and writes them as OrthoXML.
It is a research prototype built in a 48-hour autonomous project (clock in `START.txt`). Method development and all numbers
quoted here are **dev-only** and not yet independently reviewed; see "Status" below and `CLAIMS.md`.

## What Ancestra is

Ancestra splits the problem in three steps. (1) A near-linear **family partition** (DIAMOND clustering, then every gene against the
cluster representatives, then connected components with a Leiden cap) sorts all proteins into homologous families, with no
all-vs-all search and no reference database. (2) Sequence similarity is computed **only inside each family** (DIAMOND in
blocks, plus an MMseqs2 re-alignment of the pairs DIAMOND missed). (3) A **level-synchronous merge** walks a guide species tree
bottom-up. At each node, all families merge their child ancestral genes into the ancestral genes of that node, and each of these
is a HOG. Ambiguous many-to-many clusters at divergent nodes are resolved with a small local neighbour-joining tree and
species-overlap labelling. When gene order is supplied, ancestral gene-order (synteny) context chooses among candidate partners,
and vetoes in-paralog absorptions that the context contradicts. Output is OrthoXML with named taxonomic levels.

**Honest positioning.** Ancestra is not a new kind of method; it is a scalable, synteny-aware descendant of earlier
bottom-up approaches. **SYNERGY** (Wapinski 2007) already merged orthologs bottom-up along a species tree using sequence plus synteny,
and **Hieranoid** did progressive InParanoid merges; both rely on all-vs-all searches. **OrthoFinder** infers orthogroups from
all-vs-all DIAMOND plus gene-tree reconciliation, and **FastOMA** places proteins into a reference HOG database with OMAmer and then
processes HOGs along the species tree. What Ancestra tries to add is (a) a reference-free, near-linear family step feeding a bottom-up
merge over all families at once, (b) ancestral gene-order context inside the merge decision, (c) a hybrid of score rules and local
trees, and (d) optional structure (3Di) rescue of families. Whether this is competitive on the community benchmarks is **not yet
established**. So far it has only been scored on the OrthoBench *dev* split (34 RefOGs, comparator point estimates only, no CIs) and
no QfO benchmark has been run. See `DESIGN.md` section 2 for the novelty discussion and `LITERATURE.md` for prior work.

## Status (at milestone M2, 2026-10-26)

- Dev only: OrthoBench dev OG-F 75.7 with gene order and 71.3 without (deterministic partition), against OrthoFinder 70.5 on the same
  34 RefOGs (`CLAIMS.md` rows C-CORE-12, C-BASE-1). Partition sensitivity is +-2-4 points and no confidence interval is available yet.
- Test-split and QfO results: **not yet available** (placeholders in `REPORT.md`).
- Layers tested with shuffled-layer controls: synteny kept; structure only as an optional S0 rescue; expression and isoforms rejected by
  pre-registered rules (`RESEARCH_TREE.md`, tree item 7).

## Install

Everything runs in one conda environment (Python 3.11; DIAMOND, MMseqs2, Foldseek, FAMSA, FastTree, ASTER, OrthoFinder, FastOMA,
SonicParanoid and others). Details: `ENVIRONMENT.md`. Lock files are in `envs/`.

```bash
# exact rebuild from the explicit lock file (conda-forge + bioconda, linux-64)
mamba create -p conda_envs/orth --file envs/orth.lock.txt
# or from the environment export
mamba env create -p conda_envs/orth -f envs/orth.yml

source /opt/conda/etc/profile.d/conda.sh && conda activate "$PWD/conda_envs/orth"
bash scripts/smoke_test_tools.sh       # 25 tool checks
```

Ancestra itself is plain Python in `ancestra/` (numpy, scipy, lxml, pandas, numba; see `pyproject.toml`). We ran it from the repository
root (`python -m ancestra ...`); `pip install -e .` should also work but was not tested.

## Quick start on the toy simulator

These exact commands were run on 2026-10-26 (8 species, 300 root genes, 2,340 proteins, 4 threads). `simulate` takes about 1 s
and `run` about 14 s.

```bash
# 1. simulate 8 species x 300 genes with exact orthology truth (seeded, byte-reproducible)
python -m ancestra.simulate --out toy --seed 11 --species 8 --genes 300

# 2. run the pipeline (sequence only)
python -m ancestra run --proteomes toy/proteomes --tree toy/species_tree.nwk --out toy_out --threads 4

# 3. score against the simulator truth
python -m ancestra eval --truth toy/truth/hogs.orthoxml --pred toy_out/hogs.orthoxml
```

Observed on this run: `ortholog pairs: P=0.9929 R=0.9344 F1=0.9628 (tp=7569 pred=7623 true=8100)`, with per-level HOG
tables for N0..N6. This is a toy check, not a benchmark result: the toy has no domain shuffling, no gene-model errors and a
WAG-only substitution model (`docs/simulator.md`, "Explicit simplifications").

With synteny. The simulator writes `toy/gff/<species>.gff3`; turn each into a gene-order table, then pass the directory:

```bash
python - <<'EOF'
import glob, os
from ancestra.genomes import build_ensembl_species, write_gene_order
os.makedirs("toy_gene_order", exist_ok=True)
for fa in sorted(glob.glob("toy/proteomes/*.fa")):
    sp = os.path.basename(fa)[:-3]
    df, cov = build_ensembl_species(fa, f"toy/gff/{sp}.gff3")
    write_gene_order(df, f"toy_gene_order/{sp}.parquet")
EOF
python -m ancestra run --proteomes toy/proteomes --tree toy/species_tree.nwk \
    --gene-order toy_gene_order --out toy_out_go --threads 4
python -m ancestra eval --truth toy/truth/hogs.orthoxml --pred toy_out_go/hogs.orthoxml
```
On this run, F1 was 0.9618 with gene order vs 0.9628 without. The toy has mostly tandem duplicates, where synteny abstains, so no gain
is expected here; the gain is measured on OrthoBench (below).

`python -m ancestra run --help` lists all flags: `--proteomes`, `--tree`, `--out` (required), `--families`, `--families-cols`,
`--family-config`, `--threads`, `--backend {hybrid,diamond,mmseqs}`, `--arm`, `--config`, `--gene-order`, `--force`.

## Inputs

| Input | Required | Format |
|---|---|---|
| `--proteomes DIR` | yes | One FASTA per species (`.fa .fasta .faa .fas .pep`, `.gz` allowed). The species name is the file stem, and it must equal a leaf name in the guide tree. `sp\|ACC\|NAME` headers are reduced to `ACC`. Protein IDs must be unique across species. The pipeline expects one protein per gene (e.g. the canonical isoform); it does not choose isoforms itself. |
| `--tree FILE` | yes | Rooted Newick guide species tree. Internal node names become the HOG levels (OrthoXML `TaxRange`). Multifurcations are merged k-way. Leaves without a proteome are pruned; a proteome missing from the tree is an error. |
| `--gene-order DIR` | optional | `<species>.parquet` per species with columns `protein_id, gene_id, contig, start, end, strand, order, contig_n, primary, source` (schema in `ancestra/genomes.py`). Enables the synteny layer (constrained re-ranking plus veto; default w=10, lambda=0.3). Species without a table simply get no synteny evidence. |
| `--families TSV` | optional | Gene-to-family table (protein ID or integer gene ID in one column, family in another; `--families-cols g,f` picks 0-based columns). Skips S0. Unlisted genes become singleton families. |
| `--family-config JSON` | optional | Overrides of `ancestra.families.FamilyConfig` for S0. Default is `2sA_dc30_cc_sh30_RA` (`docs/S0_families.md`). |
| `--config JSON` | optional | `{"scoring": {...}, "merge": {...}}` overrides of the defaults. `--arm NAME` selects a named merge configuration (default: `MergeConfig` defaults). |

**Guide tree.** For real data the default is the NCBI taxonomy tree restricted to your species (it had 0% split conflict with the
published QfO reference; data-derived trees had 12-21%; `docs/S1_speciestree.md`). Download `taxdump.tar.gz` from NCBI (we used the
2026-10-26 dump, md5 in `DATA.md` section 11), then:

```bash
python -m ancestra.speciestree ncbi --proteomes proteomes/*.fa --taxdump data/taxonomy \
    --leaf-naming ensembl --out guide_tree.nwk
```
This was run on the 12 OrthoBench proteomes and produced a binary tree with 12 leaves and 11 named clades. Use
`--leaf-naming qfo` for UniProt `UPxxxxxxxxx_<taxid>.fasta` file names (NCBI taxid taken from the file name); `ensembl` takes
the `Genus_species` part of Ensembl-style file names and looks it up in the taxdump `names.dmp`. Check the result: leaf names must equal your proteome file stems.

**Building a gene-order directory.** You need genomic coordinates for the proteins in your FASTA files (a GFF3 that links protein IDs
to genes, or UniProt/Ensembl `gcoord` XML). The table builders are:
- `ancestra.genomes.build_ensembl_species(fasta, gff3_full, gff3_primary=None)` for Ensembl-style GFF3 (protein to transcript to gene, as in the toy
  snippet above). `scripts/build_gene_order_orthobench.py` is the full multi-species example (alt-haplotype genes are projected onto the
  primary chromosome).
- `ancestra.genomes_qfo` for the QfO 2022 proteomes (2022 `_gcoord.xml`, then 2026 GFF where the sequence MD5 matches). Driver:
  `bash scripts/build_gene_order.sh` (needs the data under `data/`; see `DATA.md`).
- Proteins without coordinates are left out of the table and get no synteny evidence. Coverage is reported by the builders
  (`coverage.tsv`). `scripts/verify_ids.py --gene-order DIR` checks the tables against the proteomes.

**Families.** Optional precomputed families (for example from your own clustering) go in as a TSV. Without it, S0 runs inside the
pipeline.

## Outputs (in `--out`)

| File | What |
|---|---|
| `hogs.orthoxml` | **Primary output.** OrthoXML 0.5, XSD-validated at write time. Nested `orthologGroup`/`paralogGroup` elements with `TaxRange` = the lowest guide-tree level of the members; `<taxonomy>` block included. pyham loads it (with `write_orthoxml(..., pyham_compat=True)` for the rare top-level paralogGroup case). |
| `pairs.tsv` | Secondary export derived from the same HOGs: ortholog pairs `gene_a, gene_b, level` using the QfO `map_relations.py` flattening rule (deepest common group is an orthologGroup). Not a separate prediction. |
| `orthobench_ogs.txt` | Root-level groups in the Open OrthoBench input format (`Orthogroup` line header, one group per line). |
| `manifest.json` | Provenance: tool/package versions, git commit, argv, full config, input md5s, per-stage wall time and peak RSS, per-node merge statistics, output hashes. |
| `stages/<stage>/` | Checkpoints (`input`, `families`, `scores`, `merge`, `output`), each with a `DONE` hash of its config and upstream stages. A rerun skips unchanged stages; changing only the merge config reruns only merge and output. `--force` recomputes everything. |

Verification: `python scripts/verify_orthoxml.py RUN_DIR` (XSD, geneRef integrity, TaxRange equals the LCA of the member species,
flattening equals the QfO extractor, pyham loads, OrthoBench reader accepts). `docs/qfo_run.md` section 3 lists the checks and their
negative controls.

## Pipeline stages and documentation

| Stage | Code | What | Documentation |
|---|---|---|---|
| input | `ancestra/io.py`, `ancestra/newick.py` | gene table with integer IDs, species table, pruned guide tree | `docs/core_v0.md` section 1 |
| S0 family partition | `ancestra/families.py` | DIAMOND cluster, genes vs representatives, connected components, Leiden cap 7,500, reassignment; deterministic | `docs/S0_families.md` |
| S1 guide species tree | `ancestra/speciestree.py` | NCBI taxdump tree (default) and data-derived/perturbed/random arms | `docs/S1_speciestree.md` |
| S2 scoring | `ancestra/scoring.py` | in-family DIAMOND blocks plus MMseqs2 rescue; large-family caps | `docs/core_v0.md` section 2 |
| S3 merge | `ancestra/merge.py`, `merge_numba.py`, `localtree.py` | level-synchronous k-way merge, in-paralog absorption, fission repair, local trees | `docs/core_v0.md` section 3 |
| S5 synteny | `ancestra/genome.py`, `genomes.py`, `genomes_qfo.py`, `synteny.py` | gene-order tables, ancestral adjacency, constrained re-ranking and veto | `docs/S5_synteny.md` |
| S6 structure (optional, off) | `ancestra/structure.py` | 3Di rescue of S0 families; not in the default pipeline | `docs/S6_structure.md`, `docs/structure_feasibility.md` |
| S7 isoforms (rejected) | `ancestra/isoforms.py` | best-isoform-pair scoring; killed by the pre-registered rule | `docs/S7_isoforms.md` |
| S8 gene-loss validation (add-on) | `ancestra/geneloss.py` | miniprot search of the syntenic locus to check apparent losses; not part of `run` | `docs/S8_geneloss.md` |
| S9 expression (rejected) | `ancestra/expression.py` | Bgee profile similarity; killed | `docs/S9_expression.md` |
| output | `ancestra/orthoxml.py` | OrthoXML writer/reader, HOG-to-pair flattening, XSD in `ancestra/schemas/` | `docs/core_v0.md`, `docs/simulator.md` |
| QfO-scale run | `scripts/run_qfo_full.sh`, `scripts/qfo_pipeline_driver.py` | resumable, code-aware driver for the 78-species QfO 2022 set | `docs/qfo_run.md` |

Toy simulator (exact truth, WAG model, DTL, WGD, inversions): `ancestra/simulate.py`, `docs/simulator.md`. Toy scorer:
`ancestra/evaluate_toy.py`.

## Repository map

| Path | Content |
|---|---|
| `ancestra/` | The method (what `python -m ancestra` runs). It reads only declared inputs: proteomes, gene order, guide tree, optional structure and expression files. `tests/test_no_leak.py` enforces that no benchmark reference data is touched. |
| `ancestra_eval/` | Evaluation helpers that read benchmark references (RefOG scorer for S0). Kept apart from `ancestra/` on purpose. |
| `scripts/` | Experiment drivers, summarisers and verification scripts (per sub-task `S0_*` ... `S9_*`, `run_core_v0_toy.py`, `verify_*.py`, `scaling_*.py`, `build_gene_order*.sh`). |
| `tests/` | pytest suite (14 files, 132 tests; includes toy-with-known-answer, determinism and leak tests). |
| `results/` | Small committed result files with manifests, one directory per sub-task (`S0_families` ... `S9_expression`, `core_v0`, `core_v0_orthobench`, `qfo_run`, `scaling`) and `test_log.tsv` (every test-set peek). |
| `docs/` | Sub-task write-ups (`S0` ... `S9`, `core_v0`, `qfo_run`, `simulator`, `structure_feasibility`) and design reviews in `docs/reviews/`. |
| `bench_repo/` | **Separate benchmark repository** (nested git repo, ignored by this repo) with local runners for Open OrthoBench and the QfO 2022 benchmark, controls, published competitor results and competitor baselines. See below. |
| `splits/` | OrthoBench dev/test split by homology component (seed 11) with SHA256 sums. |
| `envs/` | Conda lock files. |
| `PROMPT.md` `DESIGN.md` `LITERATURE.md` `DATA.md` `ENVIRONMENT.md` `LESSONS.md` `RESEARCH_TREE.md` `CLAIMS.md` `REPORT.md` | Project record (see "Where to read next"). |
| `data/`, `work/`, `conda_envs/` | Not committed: raw data (checksums in `DATA.md`), large intermediates, the environment. |

## Tests

```bash
python -m pytest -q -m "not slow"     # 132 passed in 92 s on 2026-10-26 (24-CPU shared node)
python -m pytest -q tests/test_no_leak.py tests/test_s0_determinism.py   # leakage and determinism only
```
Some tests call DIAMOND, MMseqs2, Foldseek or miniprot, so they need the conda environment. One test marked `slow` (simulator timing)
is deselected by default in the command above.

## Benchmarking (separate repository)

`bench_repo/` holds the scoring code needed to reproduce every comparison step by step; the method repository holds no benchmark
reference data. Quick pointers (full instructions in `bench_repo/README.md`):

```bash
cd bench_repo
bash scripts/download_orthobench.sh
bash scripts/run_orthobench.sh my_run/orthobench_ogs.txt my_method          # OrthoBench; my_run = an Ancestra run on the 12 OrthoBench proteomes
bash scripts/run_qfo_local.sh -i hogs.orthoxml.gz -n Ancestra               # QfO 2022, local, no docker
```
OrthoBench input proteomes are the 12 files in `Open_Orthobench/BENCHMARKS/Input`; QfO input is the 78 canonical
`qfo_release-2022_02_with_updated_UP000000437` proteomes (`scripts/prepare_qfo_inputs.py` links them, and
`scripts/verify_ids.py` must pass before a run). Competitor baselines (OrthoFinder, SonicParanoid2) and the 20 published QfO 2022
results are under `bench_repo/results/` and `bench_repo/reference_results/`. Test-set scoring is limited: at most 3 logged QfO peeks
before the config freeze (`DESIGN.md` section 4; `results/test_log.tsv`).

## Where to read next

- `RESEARCH_TREE.md`: status summaries and the full tree of experiments, positive and negative.
- `REPORT.md`: final report (skeleton now; test-set numbers are TODO until the frozen-config runs).
- `CLAIMS.md`: every quantitative claim with its evidence file, dev/test label and review status.
- `DESIGN.md` (v2, pre-registered rules), `LITERATURE.md`, `DATA.md` (sources, versions, checksums, licences), `ENVIRONMENT.md`, `LESSONS.md`.

## For phase-2 projects

Other projects will build on Ancestra's outputs and code. What they can rely on, and what to be careful about:

- **Stable interface:** `hogs.orthoxml` (OrthoXML 0.5, levels = guide-tree node names) and `manifest.json`. Read HOGs with
  `ancestra.orthoxml.read_orthoxml`, cut at a level with `hogs_at_level(ox, level)`, flatten to pairs with `hogs_to_pairs(ox)`, or
  load in pyham. `ancestra/geneloss.py` (`HogIndex`) is an example consumer.
- **Cite numbers from `CLAIMS.md`, not from this README.** Anything labelled dev-only or pending was tuned on 34 OrthoBench RefOGs and is
  not independently reviewed.
- **Do not read the test split** (`splits/orthobench_split.json`, test RefOGs; QfO benchmark references) from method code; use the
  leak test as a template.
- **Known weaknesses to design around:** over-merging in the largest families (precision 63-75% at OG level on dev), deep (Bilateria)
  decisions are score-based without an outgroup, S0 cannot reach members with no detectable sequence similarity (completeness ceiling
  about 0.82 on dev), synteny gives no help at deep levels (human vs fly AUC 0.49), and the gene-loss validator found that about 32% of
  apparent single-species losses are fission errors of our own HOGs.
- **Reproducibility:** fix `--threads` only for speed (S0 is deterministic across thread counts and input order), keep seeds from
  `DESIGN.md` section 5 (11, 23, 47) for stochastic steps, and record `manifest.json` with every result.
- Runtime numbers in the docs come from a shared node (load 20-90 on 24 CPUs) and are upper bounds.

## License and data

Code license: not yet chosen (to be set when the repository is published). Datasets keep their own licences (`DATA.md`); OrthoBench has
no explicit licence and must not be redistributed; no raw data is committed.
