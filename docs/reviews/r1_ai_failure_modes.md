# Review R1: AI-agent failure modes (reviewer: failure-modes specialist)

**Clock check.** `date` = 15:06 and START = 14:40, so we are at T+0h26. DESIGN.md is stamped "T+1h15". Agents should stamp time from `date` vs START.txt with a script (`scripts/clock.sh`), never estimate it.

## Verdict
The idea is coherent, but the design invites five failure modes: (1) dev/test leakage into the QfO tests it calls held-out, (2) silent ID and coordinate mapping errors, (3) an over-claimed novelty and linearity story with no dominance check against Hieranoid 2, (4) a 9-layer scope with no timeboxes or kill criteria, and (5) numbers in REPORT.md that cannot be traced to a file.

## 1. Leakage between development and the "held-out" test (MUST)
- **S3, S5 and S9 use VGNC for development, but VGNC is a QfO 2022 test.** S2 measures recall on SwissTree/TreeFam-A pairs, yet §4 says these families "are not inspected during development".
  - Fix: develop on OrthoBench, the toy simulator, and non-test QfO signals (consensus of public competitor pairs).
  - QfO full evaluations before the freeze are a **budgeted peek**: at most 3, pre-declared in RESEARCH_TREE.md, all reported.
- **VGNC is partly synteny-derived** (HGNC/VGNC ortholog assignment uses conserved location). A synteny gain on VGNC is therefore partly circular. It cannot be the S5 headline evidence.
- **Benchmark truth files are on disk and readable by any agent.** This covers `data/qfo/controls/oracle_orthologs.tsv`, `data/qfo/reference_data/2022.1/*` (including `species_tree_samples_*.nwk`, `lineage_tree.phyloxml`, `ReconciledTrees_*`, `vgnc-orthologs`), and OrthoBench `RefOGs`.
  - Fix: a lint test, `tests/test_no_leak.py`, that fails if any method-side code or config references `reference_data`, `controls`, `RefOG`, `vgnc`, `SwissTree`, `TreeFam`, `species_tree_samples` or `lineage_tree`.
  - The method may read only `data/inputs/` (symlinks to proteomes, GFF and structures).
  - S1's "NCBI taxonomy" must come from the NCBI taxdump via proteome taxids, **not** from QfO's lineage/species-tree files, which define the STD test.
- **OrthoBench is both tuning set and final benchmark**: label its numbers "development-set". Ablation deltas need a seeded paired bootstrap over the 70 RefOGs; no claim on Δ < 2·SE (SwissTree SE ≈0.025).

## 2. ID and coordinate mapping (MUST)
QfO silently drops unmappable pairs, so a mapping bug looks like a "recall problem" that agents will try to tune away.

- **`verify_ids.py`, run on every QfO submission before scoring.** It checks:
  - every `protId` is a bare accession (regex `^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$`), not `tr|X|NAME`;
  - every accession belongs to the **updated-zebrafish** tarball (both the plain 2022_02 and the updated tarball sit in `data/qfo`; delete or rename the non-updated one);
  - ≥99.5% of canonical accessions appear in the output;
  - no isoform accessions from `_additional.fasta` leak in;
  - the fraction of submitted pairs that QfO maps (from its logs/DB counts) is ≥99%.
- **Synteny coordinates.** The 2022 tarball already holds `*_gcoord.xml` for 19 eukaryotes; I checked the tar listing. Use those first. Use 2026 GFF only for the rest, and only when the accession maps **and the sequence MD5 (or at least the length) is identical**. Report per-species coverage and the number of discarded mismatches.
- **OrthoBench/Ensembl 100 GFF.** Coverage looks good (99.8–100%), but:
  - restrict to the primary assembly (`chr_patch_hapl_scaff` adds false neighbours);
  - map protein→gene via the CDS `protein_id`;
  - collapse overlapping/nested genes;
  - log genes on contigs shorter than 3 genes (Tetraodon, Ciona).
- **Toy round-trip test.** Build a 3-species fake genome with known order. `neighbours(gene)` must return the expected IDs after mapping. A shuffled-order negative control must drive the S5 gain to ~0 on OrthoBench.

## 3. OrthoXML and export validity (MUST)
`verify_orthoxml.py` must pass on every output. It checks:
- (a) `xmllint --schema orthoxml-0.4.xsd` passes;
- (b) every `geneRef` resolves, each gene is referenced at most once, there are no empty groups, and species `NCBITaxId` matches the proteome taxid;
- (c) every `TaxRange` equals the LCA of the species it contains, and child levels are strict descendants of parent levels;
- (d) **the pairs from QfO's own `map_relations.py` equal the pairs from our internal pair export**, so that the internal metrics and QfO agree;
- (e) pyHam loads the file;
- (f) the OrthoBench export (root HOGs) passes `benchmark.py` with no rejected genes.

The positive controls (perfect = 1, random ≈ 0) must also run through **our** exporter, not only through hand-made files. Random-pairs GO Schlicker is 0.32, not 0: state baselines.

## 4. Over-claiming and metric gaming (MUST)
- **Pre-register a single primary output.** The method is the OrthoXML HOGs (the prompt requires them). The pair export is a secondary, clearly labelled variant. Choosing the better export per test after the fact is forbidden.
- **Hieranoid 2 is the mandatory baseline in every QfO plot.** Kill/pivot criterion at T+20h: if Ancestra v0 is dominated by Hieranoid 2 on both SwissTree and TreeFam-A (lower TPR at equal or lower PPV, by more than 2 SE), stop layer work and fix the core first.
- **Rescore all competitors locally.** Rescore their B2SHARE files with our local runner and compare like with like. First verify the reproduction: FastOMA and OrthoFinder must be within |Δ| ≤ 0.005 of the OpenEBench numbers, or no comparison is reported.
- **Linearity claim.** Measure, do not assert:
  - log-log slope of wall time and CPU-h vs S over ≥5 points, with total genes reported;
  - per-merge instrumentation (ancestral proteome size, hits, seconds), because "unmatched passes up" can make the root merges grow. Report the largest merge.
  - Competitors run with identical threads and a hard timeout. A timed-out run is reported as "> T", never extrapolated as a measurement.
- **CLAIMS.md ledger.** Each claim → result JSON path → commit → reviewer who reproduced it. REPORT.md tables are generated by `make_report_tables.py` from `results/**/*.json`, with no hand-typed numbers.

## 5. Scope, timeboxes and kill criteria (MUST)
The design has 9 sub-tasks plus a simulator plus competitor runs in about 40 h. Timebox each layer:

| Item | Timebox | Kill / downgrade criterion |
|---|---|---|
| Toy simulator | 2 h; tree + dup/loss + adjacency only, sequences from numpy LG substitutions | If not done in 2 h, use OrthoBench subsets as the only positive control |
| Core v0 (seq only) on OrthoBench | first OrthoBench F by **T+10h** | If F < 0.5 at T+12h, or a run takes > 2 h, simplify (RBH+in-paralog, 1 rep) before any layer |
| QfO local runner | already partly working | Must reproduce FastOMA/OF by T+6h |
| S5 synteny | 5 h | No ΔF > 2 SE on OrthoBench, and shuffle control fails → negative result |
| S6 3Di | 4 h; ProstT5 gguf exists, so time 1k proteins first | Projected full-QfO cost > 6 h wall → restrict to OrthoBench + deep nodes |
| S7 isoforms | 2 h | — |
| S8 DNA loss / S9 Bgee | 3 h each, OrthoBench vertebrates only | Report effect or negative; no QfO-scale run |

Hard dates: **code freeze T+36h, final runs started by T+38h, REPORT.md skeleton from T+24h**.

## 6. Sub-agent outputs and context loss (SHOULD)
- Every agent run writes `results/<id>/manifest.json` (command, git SHA, input/output checksums, seed, wall/CPU time, 5-line summary).
- Reviewer agents must re-run `verify_*.py` and recompute one number independently.
- Long jobs: Snakemake + per-merge checkpoints; `RESUME.md` with the next 3 actions at each milestone.
- Seeds {11,23,47}: name the stochastic steps and report seed spread for them.

## 7. Could
- Spot-check 20 random HOGs (MSA + tree) per milestone; hierarchy check of duplication levels vs reconciled SwissTree trees at final evaluation only.
