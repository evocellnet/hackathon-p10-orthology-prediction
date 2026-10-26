# S1 — guide species trees (DESIGN.md v2, §3 row S1)

Status: topology arms done and benchmarked on their own. The policy arms (strict, delayed, per-family) need the merge code and belong to full-pipeline ablation.
Author: S1 sub-agent, 2026-10-26 16:32–17:25 CET (T+1h52 to T+2h45).

- Code: `ancestra/speciestree.py` (module + small CLI: `ncbi`, `perturb`, `random`, `compare`).
- Driver: `scripts/S1_run.py` (stages `ncbi search armA armB hybrid perturb`).
- Evaluation: `scripts/S1_evaluate.py` (metrics, table, manifest).
- Tests: `tests/test_speciestree.py` (16 tests). They use toy trees with known answers: an additive tree for NJ, hand-worked RF and quartet values, polytomy refinement, rooting, a toy taxdump with a merged taxid and a duplicated name, determinism of the NNI and random generators, and random trees as a negative control.
- Results: `results/S1_speciestree/{ob,qfo}/*.nwk`, `metrics.json` (all numbers below), `summary.md` / `summary.tsv` (generated from metrics.json), `manifest.json` (input checksums, code sha256, commands, seeds, runtimes, output sha256), and `*/runtime.json`.
- Large intermediates (hit tables, gene trees, supermatrix) are in `work/S1/` (gitignored, reproducible).

## 1. Leakage control (DESIGN §4)
Method inputs are limited to three things: proteomes, S0 families, and the NCBI taxdump. The taxdump was downloaded on 2026-10-26, md5 `1c0b66d2…`, matching NCBI's `.md5`; see DATA.md §11.

- `ancestra/speciestree.py` and `scripts/S1_run.py` contain none of the forbidden tokens (reference_data, RefOG, SwissTree, TreeFam, lineage_tree, species_tree_samples). This was checked with grep.
- The published QfO tree is read only by `scripts/S1_evaluate.py`. That file is `data/qfo/reference_data/2022.1/lineage_tree.phyloxml`, the "SwissTree speciestree for QfO 2018 dataset", which covers all 78 species of 2022_02. The evaluation script writes a Newick copy to `results/S1_speciestree/eval_only/`, a path the method never reads.
- The QfO species-tree samples (STD) were not used at all.
- **OrthoBench:** the NCBI tree of the 12 species is fully binary (11 clades). It matches the accepted phylogeny (Ecdysozoa; Olfactores is not testable with one tunicate; Euarchontoglires; Clupeocephala). So "vs NCBI" doubles as "vs published" there, and no separate tree was hand-entered.

## 2. Arms (all rooted, with named internal nodes, leaves = proteome file stems as in `ancestra.io`)
| arm | file | how |
|---|---|---|
| NCBI | `ncbi.nwk` | Induced taxdump subtree (merged taxids followed; OrthoBench names resolved via names.dmp, e.g. *Canis familiaris* → 9615). Unary nodes are collapsed and only nodes with ≥2 children are kept. Names are NCBI scientific names; a duplicate or leaf-clashing name gets `_<taxid>`. `ncbi_nodes.tsv` lists taxid, rank, n_children and n_leaves. QfO: 51 internal nodes for 78 leaves (binary max 77), with a root trichotomy Bacteria/Archaea/Eukaryota. OrthoBench: 11 (binary). |
| A: best-hit NJ | `armA_{bionj,fastme}_{ncbiroot,midroot}.nwk` | Seeded sample (seed 11) of proteins ≥50 aa: 500 per eukaryote and 1,000 per prokaryote (QfO: 53,370 queries = 5.4% of 980,536); 500 per species for OrthoBench (6,000). DIAMOND `--mid-sensitive -k 5 --max-hsps 1 -e 1e-5` against **each proteome separately**, so best and second-best hits are known per species. d_ij = −ln(mean identity of the best hits of i's queries in j), symmetrised; pairs with fewer than 5 hits get 1.5 × the maximum. Trees are built with FastME 2: `-m B` (BIONJ), and BIONJ plus balanced-ME NNI+SPR (`-s -n`). |
| B: ASTRAL-Pro | `armB_astralpro_ncbiroot.nwk` | 400 near-single-copy families, then FAMSA, FastTree `-lg -gamma`, and ASTRAL-Pro3 (ASTER 1.25, seed 11). **OrthoBench** uses S0 `2sA_dc30_cc_sh30_RA` families with ≥50% occupancy, ≤2 copies per species and ≤1.3 mean copies; 2,295 candidates, 400 sampled with seed 11; variant `_sc1` = strictly single-copy, 1,177 candidates. **QfO** has no S0 run yet: the full S0 would take about 1–2 h at 12 threads on 989k proteins, given that OrthoBench's 251k took 263–421 s wall and 0.55 CPU-h. So families were built from the arm-A hits instead. Each seed query forms a family with the best hit in every species where that hit has ≥50% query coverage and bits(2nd) < 0.9·bits(1st), i.e. a clear single best copy. Seeds with ≥39 species were kept and filtered greedily for <20% member overlap: 5,977 candidates, 400 kept, mean occupancy 57/78. |
| B: concatenation | `armB_concat_fasttree*_ncbiroot.nwk` | The same 400 alignments (first copy per species; columns with >50% gaps removed), concatenated and run through FastTreeMP `-lg -gamma`. QfO supermatrix: 179,660 sites. |
| hybrid | `hybrid_ncbi_refined_by_{armB_astralpro,armA_fastme}.nwk` | The NCBI tree with its polytomies resolved by the data-tree clades that are compatible with every NCBI clade (`refine_polytomies`). The root polytomy is left unresolved because the data trees' root is itself chosen by NCBI agreement. No NCBI clade is ever removed. QfO: +13 clades (ASTRAL) or +11 (FastME) on top of 51. OrthoBench: +0. |
| perturbed | `perturbed_ncbi_nni{k}_s{11,23,47}.nwk` | k random NNIs on the NCBI tree, with k ≈ 5/10/25% of internal nodes: QfO k = 3, 5, 13; OrthoBench k = 1, 3. `nni_perturb` swaps a random child of u with a random sibling of u. |
| random | `random_s{11,23,47}.nwk` | Random sequential joining (binary), seeded. This is the negative control. |

**Rooting.** Data trees are rooted on the edge that maximises the number of NCBI clades recovered, with ties broken by balance (`root_by_reference`). This uses only the taxonomy, which is an allowed input; `midroot` variants use midpoint rooting.

- QfO: the NCBI root is a trichotomy, so all three domain-level root positions tie. Balance then picks Prokaryotes|Eukaryota. This is an arbitrary convention, not an inference. Midpoint rooting puts the root on the Archaea branch instead (split 7|71).
- OrthoBench: NCBI-rooting gives Ecdysozoa|Chordata. Midpoint rooting puts the root on the *C. elegans* branch.

**Internal names.** A data-tree node whose leaf set equals an NCBI clade takes that clade's name. Any other node is named `<smallest containing NCBI clade>~<k>`; these are synthetic and are dissolved on output (DESIGN §8). All 41 trees load with `ancestra.io.load_guide_tree` and `TreeIndex`, and every node name is unique (this was checked).

## 3. Results
Metrics are unrooted, computed on common leaves (`compare_trees`):

- **nRF** = (|Q\R|+|R\Q|)/(|Q|+|R|) over non-trivial splits.
- **recall** = |Q∩R|/|R|, the share of resolved reference splits that are recovered.
- **conflict** = the share of query splits that are *incompatible* with some reference split. A split that only refines a reference polytomy is not counted, which makes this the "resolved splits only" measure.
- **qd_res** = quartet distance over quartets resolved in both trees (all C(78,4) = 1.42 M quartets).

The data trees are binary (75 splits); NCBI has 51 and the published tree 43. Plain nRF therefore also counts legitimate refinements of polytomies: even a perfect binary tree would score nRF > 0.

**QfO 2022_02 (78 species).** Columns are vs NCBI and vs the published SwissTree/QfO tree (evaluation only). Per-domain cells are vs the published tree. All three domains are monophyletic in every data tree.

| tree | NCBI: nRF / recall / conflict / qd_res | published: nRF / recall / conflict / qd_res | per domain vs published, nRF (conflict): Bact / Arch / Euk |
|---|---|---|---|
| NCBI itself | – | 0.128 / 0.953 / **0.000** / 0.000 | – |
| A BIONJ | 0.444 / 0.686 / 0.400 / 0.065 | 0.458 / 0.744 / 0.213 / 0.031 | 0.76 (0.30) / 1.00 (0.25) / 0.34 (0.18) |
| A FastME (BIONJ+SPR) | 0.397 / 0.745 / 0.347 / 0.050 | 0.390 / 0.837 / 0.173 / 0.020 | 0.76 (0.30) / 0.60 (0.00) / 0.24 (0.09) |
| B ASTRAL-Pro | **0.333 / 0.824 / 0.267 / 0.045** | **0.373 / 0.860 / 0.120 / 0.013** | 0.76 (0.25) / 0.60 (0.00) / 0.24 (0.09) |
| B concat FastTree | 0.349 / 0.804 / 0.267 / 0.045 | 0.390 / 0.837 / 0.133 / 0.012 | 0.76 (0.25) / 0.60 (0.00) / 0.24 (0.09) |
| hybrid NCBI+ASTRAL | 0.113 / 1.000 / 0.000 / 0.000 | 0.234 / 0.953 / 0.016 / 0.000 | 0.65 (0.08) / 0.60 (0.00) / 0.12 (0.00) |
| hybrid NCBI+FastME | 0.097 / 1.000 / 0.000 / 0.000 | 0.219 / 0.953 / 0.000 / 0.000 | 0.60 (0.00) / 0.60 (0.00) / 0.12 (0.00) |
| NNI k=3 (mean of 3 seeds) | 0.052 / 0.948 / 0.052 / 0.001 | 0.170 / 0.907 / 0.039 / 0.001 | 0.54 (0.00) / 0.50 (0.00) / 0.09 (0.06) |
| NNI k=5 | 0.092 / 0.908 / 0.092 / 0.024 | 0.213 / 0.860 / 0.078 / 0.025 | 0.59 (0.04) / 0.50 (0.00) / 0.11 (0.08) |
| NNI k=13 | 0.229 / 0.771 / 0.229 / 0.036 | 0.333 / 0.729 / 0.196 / 0.037 | 0.64 (0.12) / 0.50 (0.00) / 0.25 (0.21) |
| random (3 seeds) | 0.995 / 0.007 / 0.987 / 0.662 | 1.000 / 0.000 / 0.987 / 0.660 | 1.00 (0.77) / 1.00 (0.83) / 1.00 (0.99) |

Per-domain nRF vs the published tree is high for Bacteria (0.76) because that tree is almost unresolved for Bacteria, so most of the difference is refinement. The conflict column is the meaningful one.

**OrthoBench (12 species), vs NCBI (= published).** All six data trees have the same topology: A BIONJ/FastME with either rooting, B ASTRAL-Pro, B concat, and the `_sc1` strict single-copy variants of B. Each scores nRF 0.111, recall 0.889, conflict 0.111 (1 of 9 splits) and qd_res 0.057. Their only error is the clade (*dog*, (*human*, *chimp*)), which excludes the rodents. ASTRAL local PP for that branch is 1.0 for the near-single-copy set and 0.999 for strict single-copy. The hybrids equal NCBI (nRF 0). Controls: NNI k=1 gives nRF 0.111, k=3 gives 0.296, and random gives 0.815 (mean of 3 seeds).

**What the discordances are** (`metrics.json → conflicting_splits`, QfO vs published):

- **Eukaryotes.** Every data arm (A FastME, B ASTRAL, B concat) has the same classic long-branch artefacts:
  1. Laurasiatheria (dog, cow) grouped with primates instead of Euarchontoglires; rodents are long-branched. This is the same signal as in OrthoBench.
  2. *C. elegans* falls outside the other bilaterians, so Ecdysozoa is broken: arthropods group with the leech *Helobdella* (a lophotrochozoan).
  3. The cnidarian *Nematostella* is drawn into Bilateria, next to the chordates.

  In addition, *Ciona* sits outside (amphioxus + vertebrates) instead of in Olfactores. Neither reference resolves Chordata, so that error is not scored, but it is visible in the trees.

  NCBI encodes the accepted resolution of 1–3.
- **Bacteria.** 25–30% of the data-tree splits conflict with the published tree, and 60% (ASTRAL, concat) to 70% (FastME) conflict with NCBI. Deep bacterial phylum order is poorly determined both by 400 genes and by the taxonomy (NCBI's 2024 Pseudomonadati/Bacillati kingdoms; the SwissTree 2018 tree is a near-polytomy). The one hybrid conflict vs published (Helicobacter grouped with Pseudomonadota) reflects old vs new taxonomy: the 2018 tree keeps delta/epsilon proteobacteria inside Proteobacteria. Current NCBI places *Helicobacter* in Campylobacterota, so the hybrid has 0 conflicts vs NCBI.
- **Archaea.** ASTRAL, concatenation and FastME show no conflict vs the published tree; BIONJ has 25%.

**Robustness.** On OrthoBench, arm B gives the identical topology with near-single-copy and with strict single-copy families. On QfO, ASTRAL and concatenation differ by RF 8 within Eukaryota (4 of 45 splits) but share the artefacts above. Arm A (FastME) makes the same eukaryote errors at a fraction of the cost: it conflicts with the published tree on 0.173 of its splits, against 0.120 for ASTRAL.

## 4. Runtime
Wall and CPU times were measured at ≤12 threads on a node shared with other agents (load average 21–54 on 24 CPUs), so wall times are upper bounds. Source: `*/runtime.json` and `armB_info.json`.

| step | OrthoBench (251k proteins) wall / CPU | QfO (981k proteins) wall / CPU |
|---|---|---|
| NCBI tree (parse taxdump) | 6 s / 6 s | 4 s / 4 s |
| A: sample + per-species DIAMOND | 33 s / 280 s | 757 s / 3,007 s (0.84 CPU-h) |
| A: distances + BIONJ/FastME | 1 s | 9 s |
| B: family selection + 400 × FAMSA+FastTree | 25 s / 163 s | 234 s / 1,685 s |
| B: ASTRAL-Pro3 | 0.5 s | 78 s / 24 s |
| B: concatenation FastTreeMP | 98 s / 217 s | 926 s / 1,869 s |
| hybrid, perturbed, random | < 0.1 s | < 0.1 s |
| evaluation (incl. 1.42 M quartets per pair) | – | ~45 s total |

**Scaling.**

- **Arm A.** Cost is (sample per species × S) queries × N targets. The sample is 5% of an all-vs-all, but this is still quadratic in the number of species S. For thousands of genomes the per-species targets would also need subsampling, or a k-mer sketch, or reuse of the S0 representative hits.
- **Arm B.** Cost is linear in the number of families. With S0 families available it adds only the gene trees: about 0.05 CPU-h on OrthoBench and about 0.5 CPU-h on QfO.
- **NCBI arm.** Free, apart from a 6 s taxdump parse.

## 5. Recommendation
1. **Default guide tree: NCBI taxdump tree (`ncbi.nwk`), polytomies merged k-way (DESIGN §8).**
   - It has zero conflict with the published QfO tree (recall 0.953), against 12–21% for any data-derived tree.
   - Its nodes are real taxonomy clades, so HOG TaxRange levels come directly from it.
   - It costs nothing and is reproducible from a dated, md5-checked dump.
   - Every data-derived arm repeats known long-branch artefacts (nematode, cnidarian, rodent, *Ciona*) that would push wrong merge orders into exactly the vertebrate levels where OrthoBench and SwissTree are scored.
2. **Arms to carry into full-pipeline ablation** (DESIGN §3, S1 downstream metric = dev F):
   - (a) **NCBI** (default).
   - (b) **hybrid NCBI + ASTRAL-Pro** (`hybrid_ncbi_refined_by_armB_astralpro.nwk`). This tests whether resolving polytomies helps; it has 0 conflict with NCBI and 1 split vs published (a taxonomy-version difference). It matters only for QfO, which gets 13 new clades:
     - 4 in Bacteria, 1 in Archaea, and 3 fungal pairs.
     - Human + chimp.
     - Amorphea (Opisthokonta + *Dictyostelium*) and Choanozoa (Metazoa + *Monosiga*), both plausible.
     - A eukaryote root at Metamonada, which is debated.
     - The *Ciona* artefact (`Chordata~r4` = amphioxus + vertebrates). It is wrong per the literature but is not detectable against either reference.
   - (c) **B ASTRAL-Pro** (`armB_astralpro_ncbiroot.nwk`). This is the best purely data-derived tree, for the "no taxonomy available" use case and as a realistic error model.
   - (d) **A FastME** (`armA_fastme_ncbiroot.nwk`). This is the cheapest data-derived tree; keep it only if (c) shows that data trees matter.
   - (e) **NNI-perturbed NCBI**, k ≈ 10% (QfO k=5, OrthoBench k=1), seeds 11/23/47, for sensitivity.
   - (f) **random**, seeds 11/23/47, as the negative control. If F does not drop under (f), the guide tree is not doing its job.

   Topology vs downstream F should be read against the error sizes measured here: NNI k=5 ≈ 9% conflicting splits vs NCBI, and arm B ≈ 27%.
3. **For the policy arms (S1 discordance).** The family-level discordance count is not computed here, because it needs the merge code. The data-tree conflicts suggest where the delayed-merge policy should matter: deep Bacteria (HGT plus a poorly resolved phylum order), and four eukaryote nodes (Boreoeutheria, Chordata, Bilateria, Ecdysozoa).
4. **Scale-out path.** If a species has no usable taxid, place it with arm A distances (insert at the NCBI node of its nearest neighbours) rather than building a whole data tree.

## 6. Limitations
- The published reference exists only for QfO, and it is a 2018 SwissTree species tree with polytomies and older bacterial taxonomy.
- The QfO arm-B families come from sampled best hits, not from S0, because S0 had not been run on QfO.
- The rooting of the data trees in QfO (Prokaryotes|Eukaryota) is a convention.
- Branch lengths of the ASTRAL tree are in coalescent units, and FastME branch lengths are in −ln(identity) units. Neither is used by the pipeline.
- Single runs: no bootstrap over the arm-A sample. ASTRAL local PP is 1.0 on all OrthoBench branches.
