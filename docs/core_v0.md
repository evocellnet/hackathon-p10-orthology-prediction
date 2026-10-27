# Ancestra core v0 → v0.2: S2 scoring + S3 level-synchronous merge (+ local trees, synteny hook)

> **Update v0.2 (T+2h–T+3h30, OrthoBench dev phase; §7):**
> - New merge default: `min_cov` 0 and the hybrid local-tree mode (§3.7) at divergent nodes with an outgroup.
> - Numba engine (§3.6) and relative-divergence tests (§3.6, off by default).
> - Root parsimony merge R4 (off by default).
> - Synteny via the context hook: R1-only, constrained to synteny-off root HOGs, plus a veto (§3.5, §7). It is on in the pipeline when `--gene-order` is given.
> - §1–§6 describe v0. Their toy tables use the v0 default (`v0_default` arm).

Status: v0, toy-verified, 2026-10-26 (built T+0h50 to T+1h55, read from `scripts/clock.sh`).

- **Code:** `ancestra/io.py`, `ancestra/scoring.py`, `ancestra/merge.py`, `ancestra/pipeline.py` and `ancestra/__main__.py`.
- **Tests:** `tests/test_merge.py` has 24 tests. They cover exact hand-made HOGs with duplications, a polytomy, the coverage guard, the cap rule, the family join, and the pipeline end to end with resume.
- **Benchmark:** `scripts/run_core_v0_toy.py` writes `results/core_v0/toy_grid.jsonl` (600 records). `scripts/summarize_core_v0.py` turns it into `results/core_v0/toy_tables.md`.
- **Sweep:** `scripts/sweep_core_v0.py` writes `results/core_v0/sweep.tsv`.
- **Cross-family rescue experiment:** `scripts/exp_family_rescue.py` writes `results/core_v0/exp_family_rescue.tsv`.
- **Error analysis:** `scripts/diagnose_core_v0.py` writes `results/core_v0/diagnose_*.txt`.
- **Example manifests:** `results/core_v0/example_def32_s11_*_manifest.json`.

```
python -m ancestra run --proteomes DIR --tree tree.nwk [--families fam.tsv --families-cols 0,1] \
    --out OUTDIR [--threads 24] [--backend hybrid|diamond|mmseqs] [--arm NAME] [--config cfg.json]
python -m ancestra eval --truth SIM/truth/hogs.orthoxml --pred OUTDIR/hogs.orthoxml
```

Outputs:
- `hogs.orthoxml`: OrthoXML 0.5. It is XSD-validated in the output stage and pyham loads it.
- `pairs.tsv`: the QfO flattening of the same HOGs. This is a labelled secondary export, and a test checks that it equals the flattening of the XML.
- `orthobench_ogs.txt`: root-level groups.
- `manifest.json`: git commit and dirty flag, full config, md5 of every input, and per-stage wall time and peak RSS. Own-process VmHWM is reset per stage; for child processes it records the maximum so far. It also has per-node statistics.
- `stages/{input,families,scores,merge,output}/`: one checkpoint per stage. Each is written to a temp dir and then renamed. A `DONE` file holds the hash of the stage config plus the upstream keys. A rerun skips unchanged stages, and changing only the merge config reruns only merge and output (tested).

## 1. Inputs (`io.py`)
- **Proteomes:** one FASTA per species (`.fa .fasta .faa .fas .pep`, gz allowed), and the species name is the file stem.
  - Species are sorted by name and genes keep their file order.
  - Gene IDs are contiguous `0..N-1`, and every tool file uses only the integer.
  - `sp|ACC|NAME` headers are reduced to `ACC`.
  - Protein IDs must be unique.
- **Guide tree:** rooted Newick, where internal names are the taxonomic levels (they become TaxRange).
  - Multifurcations are kept, and the merge is k-way. Nothing is binarised, so there are no synthetic levels.
  - Unnamed internal nodes become `node<k>`, leaves without a proteome are pruned, and unary nodes are collapsed.
  - A proteome species missing from the tree is an error.
- **Families:** a gene→family TSV, where genes are given by protein ID or integer ID and the columns are selectable. Genes not listed become singleton families.
  - If no TSV is given, the pipeline runs S0 through `ancestra.families`. `DEFAULT_FAMILY_CONFIG` in `pipeline.py` is `2sA_dc30_cc_sh30_RA`, selected by the S0 owner on OrthoBench DEV RefOGs (docs/S0_families.md):
    - DIAMOND cluster at 30% identity / 50% member cover;
    - all genes vs reps (`--sensitive`, e ≤ 1e-5), keeping edges with ≥30% coverage of the shorter sequence;
    - connected components, with a Leiden split of components larger than 7,500 genes and reassignment.
  - If that module is not importable, it falls back to DIAMOND all-vs-all plus CC.

## 2. S2 within-family scoring (`scoring.py`)
Only pairs inside one family are scored. All backends return one `ScoreTable` with:
- unordered pairs `q<t`;
- `bits`, `qcov`, `tcov`, `pident` and `evalue`;
- `self_bits` per gene.

| backend | what it runs | cost |
|---|---|---|
| `diamond` (PI r2 choice) | Whole families are packed, largest first, into blocks of ≤100k proteins; a family larger than a block is its own block. Each block runs `diamond blastp` block-vs-block with `-k 0 --max-hsps 1 -b 4 -c 1 --sensitive`, and hits between different families are dropped. | ~N·B |
| `mmseqs` (exhaustive reference) | One `createdb`. A **custom prefilter DB** (`tsv2db --output-dbtype 7`) lists exactly the within-family pairs q≤t, including self pairs, which give the self-scores. Then `mmseqs align -a --alignment-mode 3` (full SW) and `convertalis`. | Σ_f g_f(g_f+1)/2 alignments, with no seeding loss |
| **`hybrid` (default)** | DIAMOND blocks as above, then `mmseqs align` **only on the in-family pairs that DIAMOND did not report** (families ≤ `cap_all`). | DIAMOND plus the missed fraction (5–13% of pairs on the toy) |
| `foldseek3di` | Hook only (S6). It has the same interface and raises `NotImplementedError`. | – |

Other scoring rules:
- **E-value:** 10 inside families. S0 has already established homology, and the coverage guard in S3 removes local junk. Compared with e=1e-3, it gains 0–0.4 recall points with DIAMOND (§5.3) and 0.5–1.3 points with exhaustive MMseqs2 (2-sim screen, def8_s11 and hard8_s11), at no precision cost.
- **Large families** (r2 algorithms review, M3):
  - g ≤ 5,000: keep all in-family hits.
  - g > 5,000: keep, per query, the top 2 hits per other species plus ≤20 same-species hits; a pair survives if either direction keeps it.
  - Dropped edges are logged per family (`pairs_dropped_per_family`).
  - g > 20,000 should already have been split by S0; such families are counted (`n_families_over_gmax`) and get the same reduction.
  - Tested only by a unit test so far; no toy family exceeds 132 genes.
- **Cross-family rescue** (`family_rescue`, default on; r2 external-PI M1): the DIAMOND pass also sees hits *between* S0 families in the same block.
  - Two families are joined when some hit between them has min(qcov, tcov) ≥ 0.5 and e ≤ 1e-5.
  - Only the joined families are re-scored (`score_with_rescue`), and the cross-family hits are saved as `stages/scores/cross_family_hits.npz`.
  - On QfO, block packing is by size, not similarity, so blocks should be ordered by centroid similarity (r2 algorithms review, item 7) for the rescue to see the relevant hits.

## 3. S3 level-synchronous merge (`merge.py`)

### 3.1 State and similarity graph
- **Ancestral genes.** After node X is processed, the genes of clade(X) are partitioned into the *ancestral genes of X*, which have global integer IDs.
  - Each ancestral gene has a nested structure: a leaf gene, `("O", X, children)` (orthologGroup), or `("P", c, children)` (paralogGroup).
  - Leaf genes are the ancestral genes of their species.
- **Leaf score.** s(a,b) = `bits / max(self_a, self_b)` (`norm_max`, default).
  - Pairs with min(qcov, tcov) < `min_cov` = 0.5 are dropped. This is the bidirectional coverage guard against multi-domain chaining (DESIGN §9).
- **LCA bucketing (r2 algorithms review, item 5).** Each leaf pair is put in the bucket of the species-tree LCA of its two species; same-species pairs go to the leaf.
  - A pair enters the graph exactly once, at the node where it becomes a cross-child pair.
  - All families are processed together, so a node is one vectorised job. Families need no explicit handling, because no edge ever crosses families.
- **Contraction.** After X is decided, the graph H_X over the children's ancestral genes is contracted along the merge map, keeping (sum, count, max) per edge. This gives the graph G_X over X's ancestral genes, which is passed to the parent. Therefore:
  - max-link S(A,B) = max over a∈A, b∈B of s(a,b);
  - avg-link S(A,B) = mean over the **observed** member pairs. A missing pair is unknown, not zero (r2 external PI, M2).
  - Both are exact, and nothing is ever re-searched.

### 3.2 Decision at node X with children C1..Ck (k ≥ 2; any k)
"Units" are the ancestral genes of the children. Cross edges join units of different children; same-child edges join units of the same child and come from G_Ci.

**R1, mutual-best seeds.**
1. For every unit A of child i and every other child j, compute best(A, j) = max over B in Cj of S(A,B).
2. A cross edge (A,B) is *mutual-best* if S ≥ tol·best(A, j) and S ≥ tol·best(B, i), with tol = 0.95 (tol = 1 is strict RBH).
3. Mutual-best edges are unioned greedily in decreasing score.
4. A union that would put two units of the *same* child into one group is allowed only if they are in-paralogs: the best same-child score between them must be ≥ the cross-child score of the edge.

**R2, in-paralog absorption.** Take an unassigned unit u (child c) whose best cross partner v lies in a group G. u joins G if some w in G from child c satisfies S(u,w) ≥ `inpar_factor` · d(w). Here d(w) is w's best score to the members of G in other children. This is InParanoid's rule seen from the main ortholog w (`inpar_rule="ortholog"`, default).
- The step is repeated (≤5 passes, in decreasing S(u,v)) so that chains of in-paralogs are absorbed.
- Interpretation: u and w are closer to each other than w is to its orthologs, so the duplication happened after X, on the branch to c.
- Alternatives exist as arms:
  - `self`: S(u,w) ≥ S(u,v).
  - `either`.
  - `quartet`: a four-point test with the closest outgroup gene o from a sibling clade of X. u is absorbed iff S(u,w)·S(v,o) ≥ S(w,v)·S(u,o). This topology test, ((u,w),(v,o)), uses d = −ln S and falls back to `ortholog` when no outgroup is scored against all three.

**R3, fission repair.** Iterated mutual-best between still-unassigned units and groups that have *no* member of the unit's child.
- The group score is the max over its members.
- A join is accepted only if its score is ≥ `repair_ratio` · (the unit's best cross score); for a singleton partner the test is applied in both directions.
- This guard keeps out-paralogs of an existing group from being attached (hidden paralogy). Default `repair_ratio` = 0.9.

**Output.**
- Every group with ≥2 units becomes one ancestral gene of X, with structure `("O", X, [per child: the unit, or ("P", child_name, units) if the child contributed >1 unit])`. The paralogGroup level is the child, i.e. the node directly below the duplication branch, which follows the convention in `orthoxml.py`.
- Units without a partner **pass up unchanged**. Their HOG at X holds one child's genes, and their group keeps its lower TaxRange; `hogs_at_level` handles this.
- At the root, every multi-gene ancestral gene is one top-level HOG.

### 3.3 Defaults and how they were chosen

| parameter | default | evidence |
|---|---|---|
| link | **avg** | Pre-stated rule (in `run_core_v0_toy.py` before the grid ran): highest mean pair F1 over the 12 truth-family sims, ties within 0.002 broken by precision. avg, max and their variants were all within 0.0045; avg had the best precision among the tied arms. |
| score | norm_max | ≈ raw bits on the toy (±0.001 F1 over the 12 sims). Kept normalised per r2 because of length bias in real data, which the toy lacks. A 2-sim screen gave norm_min −0.001 and **OrthoFinder-style per-species-pair length fit −0.04 to −0.07 F1**, because the fit breaks the comparison of same-child vs cross-child scores that R1 and R2 rely on. |
| tol | 0.95 | Sweep: 0.9 −0.000 F1 but −0.9 P points; 0.98 / 1.0 −0.0015 / −0.003 F1. |
| inpar_rule / factor | ortholog / 1.0 | Sweep: 0.95 / 0.9 / 0.8 give +0.003 / +0.003 / +0.001 F1 at −0.8 / −1.4 / −2.4 P points, which fails the precision floor (pre-stated: P ≥ default − 0.003). |
| repair / ratio | on / 0.9 | Sweep: +0.0003 F1 (negligible on the toy, where the fissions have other causes, §6). |
| min_cov | 0.5 | Toy sweep favours 0.0 (+0.0034 F1). Kept at 0.5 as a documented exception: the toy has no multi-domain proteins, so it cannot measure the failure mode the guard exists for (DESIGN §9). Re-test on OrthoBench dev. |

Named arms are in `merge.ARMS`: `default`, `rbh`, `max`, `max+repair`, `avg`, `avg+repair`, `*/bits`, `*/quartet` and `avg+repair/either`.

### 3.4 Extension hooks
- **Context evidence / synteny (S5): `run_merge(..., context=f, gene_order=go, adjacency_w=5)`.** Specified in §3.5.
- **Relative-divergence tests (added after v0, §3.6):** `MergeConfig(rd_hidden=True, rd_inpar=True)`.
- **Local-tree mode:** `run_merge(..., resolver=g)`, where `g(node, child, lu, lv, score, labels) -> labels` may rewrite the partition. This is the slot for ambiguous clusters (FAMSA + FastTree + species overlap).
- **3Di backend:** `scoring.BACKENDS["foldseek3di"]` has the same `ScoreTable` interface. A second score column can be added to `leaf_scores`.
- **Diagnostics:** `run_merge(..., trace={})` records, per node, each clade gene's ancestral ID before and after the node. `diagnose_core_v0.py` uses it.

### 3.5 Synteny interface (for the S5 agent; data layer done, no scoring implemented)
**Gene order.** `ancestra.genome.load_gene_order_gff3(gff_dir, gene_table) -> GeneOrder` reads Ensembl-style GFF3 (protein → transcript → gene; WormBase `<gene_id>.1` IDs fall back to the gene ID).
- `GeneOrder` holds, per integer gene ID: `contig` (global index into `contig_names = "species:seqid"`), `rank` along the contig among *proteome genes* (gap-tolerant with respect to unannotated loci), `start` and `strand`. Genes without coordinates have `contig = -1`.
- For human, mouse and zebrafish the `chr_patch_hapl_scaff` file is used, because OrthoBench contains alt-haplotype genes. Those genes sit on their own contigs, so they only get neighbours there.
- `GeneOrder.shuffled(seed)` is the pre-registered shuffled-order control: ranks are permuted within each species.
- **OrthoBench (Ensembl 99) coverage:** 100% for 10 species, 99.99% for human and mouse. 3,837 contigs, 1.22M neighbour pairs at w=5, 32 s to load (`work/core_v0_ob/gene_order.npz`).
- **Toy:** `SIM/gff/` gives 100% coverage.

**Neighbours.** `genome.neighbor_pairs(go, w)` returns all same-contig pairs (a, b, rank distance ≤ w).

**Adjacency.** `genome.unit_adjacency(...)` maps those pairs through gene → ancestral-gene maps. Two ancestral genes are adjacent if some member genes are neighbours in at least one species; the weight is the number of such species, plus the minimal distance.

**What the context function receives.** At each node X, `run_merge` calls `f(ctx)` with an `ancestra.merge.NodeContext`:

| field | content |
|---|---|
| `node`, `children` | X and its child names; child index = position in `children` |
| `units`, `child` | local unit i ↔ global ancestral-gene ID; child index of each unit (units = ancestral genes of X's children) |
| `cand_a`, `cand_b`, `cand_score` | every scored cross-child unit pair (local indices) and its sequence score (avg- or max-link of the normalised leaf scores) |
| `anchors` | bool per candidate: pass-1 anchors, i.e. strict mutual best with each unit in exactly one strict mutual-best pair (sequence only) |
| `adjacency` | `(ua, ub, n_species, min_dist)`: within-child adjacencies between local units at window w (built from the members' gene order in that child's species); `None` without `gene_order` |
| `members(i)` | leaf gene IDs of local unit i (for gene-level features such as strand, tandem arrays or contig) |

**What it returns.** A float array of length `len(cand_a)`. It is an additive bonus on the candidate scores, in the same units as `cand_score` (normalised, ≤1); NaN counts as 0.

**When it is applied.** The bonus is applied **before** R1–R3: mutual-best, in-paralog absorption and repair are recomputed on `score + bonus`. That makes the scheme three-pass and non-circular:
1. Anchors come from sequence only.
2. The context is scored against the anchors.
3. The rules run on the adjusted scores.

**Using it for S5.** The DESIGN S5 support for (a,b) is the number of anchor pairs (n_a, n_b) with n_a adjacent to a and n_b adjacent to b. It can be computed from `adjacency` (two adjacency lists, one per child) and `anchors`. Tandem duplicates (members adjacent to each other) should abstain.

**Guarantees** (tested in `test_gene_order_and_context_hook`): a zero bonus reproduces the sequence-only result exactly, adjacency is always within-child, and the context is called once per internal node.

**Cost.** Each call computes adjacency from N·w neighbour pairs, which is O(N·w) per node.

### 3.6 Relative-divergence (RD) tests and the numba engine (added after v0)
**Definition.** r(a,b) = d(a,b) / mean(d(a,o), d(b,o)).
- d = −ln S, where S is the normalised (avg-link) score.
- o is the outgroup gene with the highest S(a,o)+S(b,o). The outgroups are genes of X's sibling clades, with their scores to the units (`out_edges`).
- For true orthologs at X, r ≈ T_X/T_parent. The family's rate multiplies every d, so it cancels.

**`rd_hidden` (veto of hidden paralogy).** At each node, r is computed for all R1 mutual-best seeds that have an outgroup (≥10 needed). Seeds with r > median + `rd_z`·1.4826·MAD are not merged at X.
- An out-paralog pair then meets at the parent.
- There both are typically absorbed as in-paralogs relative to the outgroup gene, which gives the right paralogGroup.

**`rd_inpar` (rate-robust absorption).** In R2, if the InParanoid test fails for (u, w), u is still absorbed when r(u,w) ≤ `rd_inpar_factor` × the node's median seed r. This means the u–w split is younger than X; the rate of a fast duplicate u enters both numerator and denominator.

**Toy sweep** (12 truth-family sims, hybrid+rescue scores, `results/core_v0/sweep_rd.tsv`; same pre-stated rule):

| config | P | R | F1 |
|---|---|---|---|
| default | 0.9895 | 0.9595 | 0.9742 |
| rd_hidden z=3 | 0.9913 | 0.9592 | 0.9749 |
| rd_inpar f=1.0 | 0.9885 | 0.9653 | 0.9767 |
| **rd_both (z=3, f=1)** | **0.9902** | **0.9651** | **0.9774** |
| rd_inpar f=1.2 | 0.9818 | 0.9720 | 0.9768 |

Both precision and recall improve. `rd_both` is arm `rd_both`; OrthoBench dev results are in §7.

**Numba engine.** The per-node rules R1–R3, RD and quartet run in `ancestra/merge_numba.py` (`MergeConfig.engine="numba"`, default).
- It uses CSR adjacencies with binary-search lookups and an array union-find with member lists.
- The pure-Python `decide_node` stays as the reference (`engine="python"`).
- On three toy sims and six arms the two engines give **identical** ortholog-pair sets (`test_numba_engine_matches_python`).
- Ties between exactly equal scores could in principle break differently, so the engines are not guaranteed bit-identical.
- Compilation takes ~10 s once and is cached. On OrthoBench (251k genes, 11 internal nodes) the merge takes 7–13 s.

### 3.7 Hybrid local-tree mode (default since v0.2; `ancestra/localtree.py`)
This is the DESIGN §1(b) hybrid: score rules where they are unambiguous, a small tree where they are not.

**Detection.** It runs at node X after R1–R3.
- Every unit u has best(u), its best edge (cross- or same-child).
- An edge between two *different* rule groups with S ≥ τ·min(best(u), best(v)) (τ = 0.9) is an ambiguity edge.
- Groups connected by ambiguity edges form a cluster. A cluster with ≥2 groups and ≤60 units is re-decided by a tree; larger clusters keep the rule result and are counted as `tree_too_big`.

**Tree.** NJ on unit distances −ln S of the linked (avg) member scores; missing pairs get max + 1.
- One outgroup leaf is added: the sibling-clade gene with the highest summed score to the cluster (`out_edges`).
- NJ is used rather than an MSA + ML tree because it is rate-robust, deterministic and cheap: O(m³) with m ≤ 60, no subprocess per cluster. FAMSA + FastTree is a possible later arm.

**Rooting.** On the outgroup branch. Otherwise at the edge that minimises species-overlap duplications.

**Labelling.** A node whose subtrees share a child of X, and together span ≥2 children, is a duplication at or before X, so it is split. Any other subtree is one ancestral gene of X: a speciation, or duplications confined to one child (in-paralogs).

**Where it runs (pre-stated selection on dev with a toy check, §7).** Only when an outgroup exists (`tree_require_outgroup`) and the node's median mutual-best seed distance −ln S ≥ `tree_min_d` = 0.5, i.e. at divergent nodes.
- On OrthoBench this is Euteleostomi, Chordata and Ecdysozoa.
- At the root (no outgroup) trees collapse precision: dev OG-F 61.
- At shallow nodes they add FPs: dev OG-F 65.7 when used there only.

**Hand cases** (`test_localtree_labelling`): a fast-evolving in-paralog is kept with its sister copy; an out-paralog after loss is kept apart; an old duplication with both copies retained is split correctly without an outgroup.

## 4. Complexity
- **S2:** DIAMOND blocks ~N·B. The rescue cost is proportional to the pairs DIAMOND missed (5–12% of Σ g_f(g_f+1)/2 on the toy). The cross-family rescue re-scores only the joined families (3–14% of genes on the toy, 2–15 s).
- **S3:**
  - E scored leaf pairs are sorted once by LCA, O(E log E).
  - At node X, the work is aggregation of H_X = contracted child graphs ∪ the X bucket, O(|H_X| log |H_X|), plus a decision step O(|H_X| + n_X·k) and contraction.
  - Each leaf pair enters once. Contracted graphs only shrink, because aggregation is massive: on def32_s11, Σ_X |H_X| = 56k edges for 596k leaf pairs.
  - Worst case (nothing merges, many paralogs) is O(E·depth).
  - Memory is O(N + E). Leaf species do no merge work.
- **Measured:** the merge takes 0.5–1 s for 33k genes on 32 species (31 nodes, max 0.16 s per node). The decision step is pure Python (dict and union-find loops over units and cross edges), so QfO deep nodes (1e5–1e6 units, ~1e8 edges) need numba (planned).

## 5. Toy results (in-house simulator; guide tree = true species tree)
**Sims.** Seeds 11, 23 and 47, in four settings:
- `def8`: 8 species × 300 genes, default model.
- `hard8`: dup 0.5, loss 0.5, a WGD on a random branch, neo-functionalisation 2×.
- `def32` and `hard32`: 32 species × 1000 genes (33k proteins) with the same two models.

**Families.** `truth` uses the simulator root families and isolates S3. `S0` uses the current default S0 (`2sA_dc30_cc_sh30_RA`; `s0b` in the JSON). The grid also kept the earlier S0 default (`s0`: reps-vs-reps CC, min-cov 50%) as a comparison.

**Metrics.** Ortholog pairs use QfO flattening against the true pairs. HOGs are scored per level with `evaluate_toy`. All values are mean ± sd over the 3 seeds.

### 5.1 Headline: default configuration (hybrid + family rescue scoring, arm `default`)
| setting | families | P | R | F1 | root HOG exact | all-level HOG exact | mean best Jaccard | deep exact (d≤1) |
|---|---|---|---|---|---|---|---|---|
| def8 | truth | 0.993±0.001 | 0.968±0.008 | **0.980±0.005** | 0.939 | 0.964 | 0.985 | 0.943 |
| def8 | S0 | 0.993±0.001 | 0.936±0.006 | **0.964±0.003** | 0.879 | 0.951 | 0.980 | 0.908 |
| hard8 | truth | 0.987±0.004 | 0.938±0.020 | **0.962±0.012** | 0.826 | 0.894 | 0.955 | 0.828 |
| hard8 | S0 | 0.986±0.004 | 0.902±0.021 | **0.942±0.013** | 0.771 | 0.879 | 0.948 | 0.792 |
| def32 | truth | 0.992±0.002 | 0.978±0.005 | **0.985±0.003** | 0.960 | 0.979 | 0.991 | 0.956 |
| def32 | S0 | 0.992±0.002 | 0.958±0.002 | **0.975±0.002** | 0.920 | 0.975 | 0.990 | 0.933 |
| hard32 | truth | 0.986±0.001 | 0.954±0.002 | **0.970±0.001** | 0.909 | 0.941 | 0.977 | 0.899 |
| hard32 | S0 | 0.986±0.001 | 0.936±0.007 | **0.961±0.003** | 0.876 | 0.938 | 0.976 | 0.880 |

Without the cross-family rescue, F1 with the current S0 is 0.962 / 0.939 / 0.972 / 0.958. With truth families the rescue makes **0 joins**: no spurious joins, which serves as a negative control.

HOG exactness rises from the root to the tips: on hard32 it is 0.909 at depth 0, 0.894 at depth 1, 0.917 at depth 2, 0.946 at depths 3–5 and 0.957 at 6+. The full per-depth table is in `toy_tables.md`.

### 5.2 Merge arms (truth families, DIAMOND-only scoring, mean over the 12 sims)
| arm | P | R | F1 |
|---|---|---|---|
| rbh (strict mutual best only, no in-paralogs) | 0.9945 | 0.6025 | 0.7434 |
| max | 0.9887 | 0.9226 | 0.9543 |
| max + repair(0.5) | 0.9871 | 0.9243 | 0.9544 |
| avg | 0.9895 | 0.9284 | 0.9578 |
| avg + repair(0.5) | 0.9883 | 0.9298 | 0.9580 |
| avg + repair, bits | 0.9881 | 0.9312 | 0.9587 |
| avg + repair, quartet in-paralog rule | 0.9835 | 0.9357 | 0.9588 |
| avg + repair, either rule | 0.9721 | 0.9413 | 0.9563 |

What the arms show:
- In-paralog absorption is essential. RBH-only recall is 0.73 (def8), 0.50 (hard8), 0.70 (def32) and 0.47 (hard32).
- avg-link beats max-link by 0.1–0.7 F1 points, mostly on hard sims (max-link favours large groups).
- Fission repair is nearly neutral.
- The quartet and either rules trade precision for recall; `quartet` is best on hard8 (+0.4 F1) and worst for FP paralog pairs on hard32 (2× default).

With hybrid scoring, the `avg` arm reaches mean F1 0.974 over the 12 sims (sweep).

### 5.3 S2 backends (truth families, avg + repair)
| setting | DIAMOND e10 | DIAMOND e1e-3 | MMseqs2 exhaustive | hybrid |
|---|---|---|---|---|
| def8 F1 | 0.961 | 0.961 | 0.979 | 0.980 |
| hard8 F1 | 0.936 | 0.935 | 0.962 | 0.961 |
| def32 F1 | 0.976 | 0.974 | 0.984 | 0.984 |
| hard32 F1 | 0.960 | 0.958 | 0.970 | 0.970 |

DIAMOND `--sensitive` reports 90–94% of the in-family pairs that MMseqs2 aligns. On the two sims checked, every DIAMOND pair is also in the MMseqs2 set. The missed pairs are the divergent ones the deep merges need, so DIAMOND loses 1–5 recall points. The hybrid recovers all of it.

Wall time (shared node, load 30–80, 3 threads per job): on the toy, exhaustive MMseqs2 is *faster* than DIAMOND, because Σg² is small: def32 71 s vs 92 s, hybrid 115 s. At 24 threads the hybrid on def32_s11 takes 15–22 s.

The PI's QfO measurement found DIAMOND faster than MMseqs2 (126 s vs 457 s at 5e6 pairs), so the hybrid (DIAMOND + rescue of ~7–13% of pairs) is the scalable choice that keeps exhaustive recall.

### 5.4 S0 and the cross-family rescue

| setting | S0 config | truth families exact | fissions | fusions | F1 hybrid | F1 hybrid + rescue |
|---|---|---|---|---|---|---|
| def8 | earlier (`s0`) | 0.840 | 47 | 0 | 0.940 | 0.960 |
| def8 | **current (`s0b`)** | 0.896 | 31 | 0 | 0.962 | 0.964 |
| hard8 | earlier | 0.842 | 45 | 0 | 0.917 | 0.938 |
| hard8 | **current** | 0.895 | 30 | 0 | 0.939 | 0.942 |
| def32 | earlier | 0.838 | 162 | 0 | 0.943 | 0.974 |
| def32 | **current** | 0.932 | 67 | 0 | 0.972 | 0.975 |
| hard32 | earlier | 0.831 | 165 | 0 | 0.924 | 0.960 |
| hard32 | **current** | 0.936 | 63 | 0 | 0.958 | 0.961 |

S0 fissions are permanent unless they are repaired, and with the earlier S0 they were the largest loss (4–8 recall points).

Joining families through the S2 cross-family hits (cov ≥ 0.5, e ≤ 1e-5) roughly halved the earlier S0's fissions (def32: 162 → 77; `exp_family_rescue.tsv`) and added +2.0 to +3.6 F1 points with no fusions. With the current S0 most of that gap is already closed, and the rescue adds +0.2 to +0.3 points. It is kept on because it is cheap (16 joins and 477 re-scored genes on def32_s11) and created no fusions on the toy.

Looser variants (cov 0.3, e 1e-3) add another +0.1–0.2 F1 with the earlier S0 and create 0.3–1.3 fusions per sim. The remaining gap to truth families is 0.3–2 F1 points.

### 5.5 Runtime
`results/core_v0/example_def32_s11_s0fam_manifest.json` records a full CLI run with default S0 on def32_s11 (33,044 proteins, 32 species, 24 threads, node load ~27). Pair F1 is 0.973.

| stage | wall | peak RSS |
|---|---|---|
| input | 0.03 s | 51 MB |
| families (S0 `2sA_dc30_cc_sh30_RA`) | 16.0 s | 85 MB (DIAMOND child: 630 MB) |
| scores (hybrid + rescue) | 16.3 s | 338 MB |
| merge (31 nodes) | 0.5 s | 125 MB |
| output (XML + XSD + pairs + OGs) | 1.8 s | 271 MB |

The end-to-end wall time is 35.1 s with a maximum RSS of 645 MB (`/usr/bin/time -v`). The cost is dominated by S0 and S2.

The grid jobs (3 threads each, 8 in parallel, load 30–80) took 5–7 s (8 species) and 65–115 s (32 species) of scoring per condition, and 0.03–1 s of merge.

## 6. Known failure modes (from `diagnose_core_v0.py`; results/core_v0/diagnose_*.txt)
**Method.** Each missed true ortholog pair is located at its true speciation level L. The script then checks what happened at L to the two "units", i.e. the ancestral genes at L's children. Counts below are distinct unit pairs.

1. **Fast-evolving young duplicates fail the in-paralog test.** This is the `one_unmatched` class: 944 of 1,207 FN unit pairs on hard32_s11, and 312 of 483 on def32_s23.
   - With neo-functionalisation, the copy u evolves faster, so S(u,w) < S(w,v) although u and w are in-paralogs.
   - The `self` rule (S(u,w) ≥ S(u,v)) fixes recall but cannot separate out-paralogs. For an out-paralog, d(u,w) − d(u,v) = d(X,w) − d(X,v), which is a coin flip. That is why `either` costs 2 precision points.
   - The quartet test with an outgroup is the principled fix. With S → d = −ln S it is too crude, because bitscores are not exponential in distance: it gains recall but loses precision.
   - Next: a calibrated distance (Scoredist-style) and the local-tree mode for these clusters.
2. **No detectable similarity (`no_edge`).** On def8_s11, 50 of 73 FN unit pairs have no scored member pair at all. These are saturated fast families at the root and depth 1; family rates are Gamma(2), so some families run at 3× the rate.
   - The hybrid scoring and e=10 fixed the share that came from seeding (DIAMOND) and the e-value.
   - The rest needs profile or structure (3Di) rescue at deep nodes. The coverage guard (0.5) accounts for a further ~0.3 recall points on the toy.
3. **S0 fissions are permanent** unless they are repaired (§5.4). With the current S0 plus the cross-family rescue, the loss against truth families is 0.3–2 F1 points.
4. **Precision errors are hidden paralogy.**
   - 100% of FP pairs are true paralogs; none are cross-family. Most of the FP mass comes from one or two duplication branches per sim: def32_s23 has 53% of its FP pairs on the branch above N15 (length 0.25), and hard32_s11 has 46% on the branch above N2 (0.30).
   - The FPs are predicted as orthologs at the node just below the duplication. This is out-paralogy after a duplication on P→X:
     - complementary loss leaves a clean 1:1 between out-paralogs (case iii of r1_external_pi);
     - or retained copies are mis-paired.
   - Short duplication branches are enriched 4–6× among FPs (13–17% of FP pairs on branches < 0.05, against 2–7% of all paralog pairs). A duplication almost simultaneous with a speciation is ~50/50 for any score rule. Still, most of the FP mass is on ordinary branches.
   - Score rules cannot see case iii. An outgroup does not help as a witness when the duplication happened after the outgroup split.
   - Proposal: a rate-robust **relative-divergence test**. r(a,b) = d(a,b) / mean(d(a,o), d(b,o)) with the closest outgroup gene o, calibrated per node against confident 1:1 pairs. The family rate cancels in the ratio, so out-paralogs show up as high-r outliers. Synteny is the other lever (S5).
5. **HOG level of lineage-specific duplicates without a partner.** Two copies in one child with no partner in a sibling pass up separately and are joined only where they meet a common partner. Their paralogGroup is then placed at that higher level, or at the MRCA if it is unknown, rather than at the true duplication branch. Pair metrics are unaffected; per-level HOG exactness is slightly affected.
6. **Not tested by the toy.** HGT (the sims use none), guide-tree errors (the guide tree is the true tree), multi-domain proteins and fusion/fission genes, fragments and isoforms, and families above 132 genes (the cap rule has a unit test only).

## 7. OrthoBench DEV (v0 → v0.2)
### 7.1 Setup
- **Input:** 12 proteomes (251,378 genes) and the S1 NCBI-taxdump guide tree (`results/S1_speciestree/ob/ncbi.nwk`, clade-named internal nodes).
- **Families:** S0 `2sA_dc30_cc_sh30_RA`, the pipeline default. Iterations reuse the S0 agent's precomputed families TSV, which is the same config.
- **Scoring:** hybrid + capped cross-family rescue (`work/core_v0_ob/run_v2`). Merge iterations reuse these scores (`scripts/ob_merge_iter.py`).
- **Scoring (DEV only; test RefOGs are never scored or read):**
  - **OG level:** `scripts/orthobench_dev_eval.py`. It calls the official `benchmark.py` `calculate_benchmarks_pairwise` on the 34 dev RefOGs of `splits/orthobench_split.json`. It reproduces the controls: positive 100, split-halves 65.1, merged-pairs 47.0.
  - **Pair level:** `scripts/orthobench_dev_pairs_eval.py` (S5 agent; ortholog/paralog labels from the RefOG trees).
- **Log:** every configuration is in `results/core_v0_orthobench/runs.tsv` (shared with the S5 agent's `S5p2` rows).
- **Error analysis:** `scripts/diagnose_orthobench.py` (S0 split vs merge split vs over-merge, per dev RefOG).

### 7.2 Main results (dev, 34 RefOGs; OG = OrthoBench pairwise-on-OGs, pair = ortholog pairs)
| config | OG-F | OG-P | OG-R | exact | pair P | pair R | pair F1 |
|---|---|---|---|---|---|---|---|
| v0 default (cov guard 0.5) | 59.0 | 75.4 | 48.5 | 3 | 0.736 | 0.596 | 0.659 |
| `min_cov` 0 (`seq_ob`) | 68.1 | 71.2 | 65.3 | 11 | 0.756 | 0.705 | 0.730 |
| + local-tree mode (**default v0.2**, `tree_ob`) | 72.5 | 63.2 | 85.0 | 9 | 0.810 | 0.681 | 0.740 |
| + synteny R1-only, unconstrained (w10 λ0.3 core) | 70.2 | 59.7 | 85.2 | 9 | 0.843 | 0.726 | 0.780 |
| + synteny R1-only, **constrained** (two-pass) | 78.5 | 75.0 | 82.3 | 9 | 0.837 | 0.736 | 0.783 |
| + constrained + veto + R3 choice (**pipeline default with `--gene-order`**) | **79.2** | 73.9 | 85.4 | 9 | 0.841 | 0.733 | **0.783** |
| shuffled order, constrained + veto + R3 (seeds 11 / 23 / 47) | 72.5 / 73.3 / 72.5 | | | | | | 0.740 / 0.740 / 0.740 |

Δ(real − shuffled) is +6 to +7 OG-F and +0.043 pair F1. These are point estimates on 34 RefOGs; no CI has been computed for the combined configuration.

### 7.3 What was learned
1. **The coverage guard is the largest single lever (+9 OG-F).**
   - OrthoBench is full of fragments and odd gene models: alt haplotypes, partial genes. A bidirectional ≥50% coverage guard isolates them as singletons.
   - The toy has none of these, so it could not show this. 0.3 gives 66.3; 0 gives 68.1–68.6.
   - The multi-domain fusion risk that the guard was meant for did not show up on dev. Re-test on QfO.
2. **The cross-family rescue needs a size cap.** Uncapped, it rebuilt the 89,721-gene component that S0 had Leiden-split; the top-2 cap then removed most edges in 13 dev RefOGs' family. Capped at 7,500 (= S0 `s2_gmax`), it gives +1.5 OG-F over no rescue (scores_pass1).
3. **Large S0 families** (6 > 5,000 genes, 28M hits dropped by the cap) hold several RefOGs, i.e. deep paralogs. This is where the precision errors live.
4. **Dev error classes (seq_ob):**
   - **Merge splits (20/34):** in vertebrate-only families, the 2R ohnolog subfamilies are not united. The pieces do have weak Ciona partners (not partner-less), and in-family scores sit at the noise floor (RefOG061: median normalised score 0.058 inside vs 0.044 to other family genes).
   - **Over-merges (13/34):** fast Ciona genes bridge a vertebrate RefOG to a deep-paralog vertebrate subfamily at Chordata (RefOG007 +104 genes, RefOG033 +82).
   - **S0 splits:** 7/34.
5. **Score-rule fixes do not transfer from the toy.**
   - RD tests: toy +0.3 F1, dev −1.3 OG-F, and −4.7 combined with trees.
   - Root parsimony merge R4: unrestricted −4 OG-F; mutual-best with t ≥ 0.2 is neutral.
   - max-link: −6 OG-F.
   - All remain as options, off by default.
6. **Local-tree mode** (§3.7) fixes the 2R splits: RefOG021 TP 506 → 7,261, and RefOG019/043/054 repaired.
   - Applied everywhere it over-merges (OG-P 41%). The root and shallow nodes are the bad places.
   - Divergent nodes with an outgroup: +4.4 OG-F; pair P +0.05, pair R −0.02.
7. **Synteny (S5, `ancestra/synteny.py`, applied through the NodeContext hook, R1-only).**
   - Unconstrained, it improves pairs (+0.040) but fuses OGs (−2.3 OG-F).
   - Constrained to the synteny-off root HOGs (two-pass, `merge.run_merge_constrained`), it keeps the pair gain and *raises* OG-F to 78.5. The partner swaps then reorganise HOGs inside a root HOG instead of linking two.
   - The R2/R3 veto (support 0 while a competitor has ≥2) adds +0.8 OG-F. R3 partner choice alone has no effect.
   - Shuffled controls are flat.

8. **Locus-aware fission repair** (S8 proposal; `MergeConfig.locus_repair`, numba engine).
   - Mechanism: a still-alone unit whose genes all come from one species joins the group lacking its child with the highest synteny support (≥ `locus_min_support`), using any in-family-scored partner. Support is unmasked, so the join can cross synteny-off root HOGs.
   - Result (support ≥ 2, vs constrained + veto + R3): pair F1 0.7996 (+0.016), S8 5/6-mammal fission candidates 639 → 350 (−45%), but **OG-F 79.2 → 75.5** (OG-P −4).
   - Shuffled controls (seeds 11, 23): OG-F 72.1–73.1, pair F1 0.740, ~815 candidates. The effect is real.
   - Without the single-species restriction, OG-F falls to 68.1.
   - Off by default (trade-off, not a gain on both metrics); a decision for the PI.

### 7.4 Toy settings resembling OrthoBench (new; `ancestra.simulate`)
**Simulator changes** (backward-compatible): `--wgd-branch` accepts several WGDs (`N@pos,...`), and `--lineage-rates node:factor` sets per-lineage substitution rates.
- `twoR12`: an OrthoBench-shaped 12-species tree (`work/core_v0/twoR_tree.nwk`), two WGDs above Euteleostomi (positions 0.3, 0.7), a teleost WGD, fractionation 0.6, Ciona 2.5× faster; 1000 families.
- `twoRstem12`: the same plus a 0.6 stem branch above Bilateria. Families then hold several Bilateria-level HOGs (deep paralogs), as OrthoBench families do; the OrthoBench analogue is the Bilateria-level HOG partition.

| setting, truth families | arm | ortholog-pair F1 | Bilateria-level HOG pair F1 (P / R) |
|---|---|---|---|
| twoRstem12 | v0_default | 0.881 | 0.851 (0.910 / 0.799) |
| twoRstem12 | seq_ob | 0.887 | 0.863 (0.903 / 0.828) |
| twoRstem12 | **default (tree_ob)** | 0.887 | **0.898** (0.936 / 0.864) |
| twoRstem12 | tree at all nodes | 0.877 | 0.925 (0.896 / 0.956) |
| twoR12 | seq_ob / default | 0.937 / 0.933 | root OG F1 0.923 / 0.938 |

On the standard 12 toy sims, the v0.2 default is within ±0.004 F1 of `seq_ob` and at least v0 in most settings (`results/core_v0/toy_grid.jsonl`, arms `default`, `v0_default`, `seq_ob`).

### 7.5 Runtime on OrthoBench (12 threads, shared node, load 45–60)
**Full end-to-end run with the v0.2 default and synteny** (`python -m ancestra run --proteomes .../Input --tree results/S1_speciestree/ob/ncbi.nwk --gene-order data/genomes/orthobench --threads 12`; manifest in `results/core_v0_orthobench/e2e_default_manifest.json`):

| stage | wall | peak RSS (own) | notes |
|---|---|---|---|
| input | 0.3 s | 0.2 GB | 251,378 genes, 12 species |
| families (S0 `2sA_dc30_cc_sh30_RA`) | 435 s | 0.9 GB (children 4.6 GB) | 22,629 families, max 7,527 |
| scores (hybrid + capped rescue) | 1,186 s | 12.9 GB (children 13.9 GB) | DIAMOND 912 s, MMseqs2 rescue 24 s, family rescue 195 s (85 joins); 15.4M scored pairs |
| merge (two-pass constrained synteny, trees) | 75 s | 2.8 GB | 15,456 tree clusters (79,782 units), 5.4 s in trees |
| output (OrthoXML + XSD + pairs + OGs) | 13 s | 1.7 GB | 15,761 root HOGs |

The total is **28.5 min wall** with **14.2 GB max RSS** (`/usr/bin/time -v`), on a shared node at load 45–60.

**Dev result of this fresh run:** OG-F 75.6 (P 72.8, R 78.6), pair F1 0.636.

**S0 reproducibility caveat.** The iteration runs reused the S0 agent's families for the *same* config (dev OG-F 79.2, pair F1 0.783). Those families differ from the freshly computed ones: ARI 0.78; largest families 7,478 vs 7,527.
- The likely cause is DIAMOND `cluster` thread-count sensitivity (24 vs 12 threads). It is being reported to the S0 owner.
- The dev pair metric is dominated by a few large RefOGs (RefOG021 alone: 7.9k–9.4k FP pairs), so it is volatile under such changes.
- Decisions so far were made on one fixed family partition (run_v2) and are internally consistent. **Absolute numbers carry roughly ±4 OG-F and ±0.15 pair F1 of S0-induced variation.**


### 7.6 S0 determinism (fixed; `ancestra/families.py`, `tests/test_s0_determinism.py`)
**Cause.** The non-reproducibility (ARI 0.78 between two runs of the same config) had two sources.
1. **Input order and ID names, not threads.** DIAMOND `cluster` is greedy, and ties follow input order and IDs. The S0 agent ran on `ob_all.fa` with protein IDs in file order; the pipeline wrote integer IDs in gene-table order. Stage-1 clusters differed.
   - A pure thread check (4 vs 16 threads, 80k proteins; 12 vs 24 on full OrthoBench) gave identical stage-1 output, and repeat runs were identical.
2. **The stage-2 `diamond blastp` search was thread-dependent.** DIAMOND's automatic block and chunk sizing follows `--threads`. 463 of 2.95M hits (0.016%) changed between 12 and 24 threads; this is not a `-k` truncation effect (max 278 hits per query). Connected components amplified it into different families.

**Fix.**
- `run_families` clusters a *canonical* copy of the input:
  - sequences normalised as in `io.load_proteomes` (terminal `*` stripped, internal `*` → X);
  - sorted by (length desc, sequence, original ID);
  - renamed by rank, then mapped back.
- The stage-2 search pins `-b 4 -c 1` (`FamilyConfig.s2_block`, `s2_chunks`).
- Rep edges are sorted before connected components and Leiden; Leiden keeps seed 11.
- The pipeline passes protein names to S0.
- Stage-1 reuse (`s1_reuse`) is refused for pre-fix runs.

**Verification.**
- `test_s0_determinism`: the toy at 2 vs 8 threads, plus shuffled order, renamed IDs and terminal `*`, with and without forced Leiden splits. It gives identical partitions, and identical family IDs for identical input.
- OrthoBench: `families.tsv` is **byte-identical** (md5 5561551b…) for `ob_all.fa` at 24 threads, `ob_all.fa` at 12 threads, and the pipeline-ordered FASTA at 12 threads. Stage-1 clusters and stage-2 hits are also identical.
- The deterministic partition has 22,643 families, the largest 7,269, with 1 Leiden-split component and 940 genes reassigned.

**Partition sensitivity of the dev metrics** (same merge code; three partitions):

| partition | off: OG-F / pair F1 micro / macro | constrained synteny + veto: OG-F / micro / macro | shuffled (3 seeds): OG-F / micro | + locus repair: OG-F / micro / macro |
|---|---|---|---|---|
| S0-agent (old, run_v2) | 72.5 / 0.740 / 0.768 | 79.2 / 0.783 / 0.791 | 72.5–73.3 / 0.740 | 75.5 / 0.800 / 0.800 |
| old pipeline order (run_e2e) | 71.4 / 0.739 / 0.769 | 75.6 / 0.636 / 0.763 | 71.4 / 0.739 | 79.1 / 0.701 / 0.775 |
| **deterministic (run_det)** | 71.3 / 0.674 / 0.752 | **75.7** / 0.722 / 0.772 | 71.3–72.0 / 0.674–0.721 (macro 0.752–0.758) | 74.5 / 0.780 / 0.792 |

- The micro pair F1 is dominated by RefOG021 (125 genes, lineage-specific expansions). On run_e2e, synteny raises its FP pairs from 27 to 1,800 of the 2,392 total; elsewhere it removes FPs (RefOG033 24 → 0, RefOG048 254 → 159).
- Macro pair F1 (mean over dev RefOGs) is far more stable. Even a shuffled control can move micro pair F1 by +0.046 (deterministic partition, seed 23) while macro moves only +0.006.

### 7.7 T+12h gate input (deterministic partition, full default pipeline)
`python -m ancestra run ... --gene-order data/genomes/orthobench --threads 12` (`work/core_v0_ob/run_det`; manifest `results/core_v0_orthobench/det_default_manifest.json`).
- **Dev OG-F 75.7** (P 72.2, R 79.6, 11/34 exact). Pair F1 is 0.722 micro, 0.772 macro.
- Without synteny the dev OG-F is 71.3.
- Runtime: 31.5 min wall, 14.7 GB max RSS (S0 462 s, scores 1,365 s, merge 54 s, output 9 s; 12 threads, shared node).


### 7.9 Local-tree backend: NJ vs FAMSA + FastTree (S4 follow-up; deterministic partition, dev only)
**Option.** `MergeConfig.tree_method="fasttree"`. NJ remains the default, and the checkpoint key is unchanged at default (`key()` omits the new fields when they are at their defaults).
- Each ambiguous cluster becomes one FAMSA (`-t 1`) + `FastTree -lg -nosupport` job. All clusters of a node run as one batch in a pool of `tree_ft_workers` (12).
- The Newick output is converted to an edge list and labelled exactly as the NJ tree is (outgroup rooting, species overlap).
- **Leaves, `tree_ft_mode`:**
  - `medoid`: one sequence per unit, the length-medoid member. This is a cheap, deterministic stand-in for a score medoid.
  - `all`: every member gene when the cluster has ≤60 genes (otherwise medoid). Each unit then goes to the tree group holding most of its genes.
- The outgroup gene's sequence is added as the rooting leaf.

**Selection rule (PI, pre-fixed).** Adopt only if Δ OG-F ≥ +2 or Δ macro pair F1 ≥ +0.01, with no loss in the other metric.

Arms (gene order on, i.e. constrained synteny + veto + R3; `run_det` scores; `results/core_v0_orthobench/runs.tsv`, label prefix `T[det]`):

| tree backend | OG-F | OG-P | OG-R | exact | pair F1 micro | pair F1 macro | extra merge cost (two passes) |
|---|---|---|---|---|---|---|---|
| NJ (default) | 75.71 | 72.2 | 79.6 | 11 | 0.722 | 0.772 | 64 s wall, 61 CPU-s |
| FastTree, medoid leaves | **80.36** | 76.2 | 85.0 | 11 | 0.804 | **0.799** | 1,039 s wall, 7,558 CPU-s (2.1 CPU-h) |
| FastTree, all leaves (≤60 genes) | 79.53 | 74.3 | 85.6 | 12 | 0.794 | 0.781 | 2,160 s wall, 21,909 CPU-s (6.1 CPU-h) |

- FastTree-medoid gives +4.65 OG-F and +0.027 macro pair F1, with OG precision and recall both up. That meets the rule.
- Cost: 9,077 jobs per merge pass (Euteleostomi 5,088, Chordata 2,607, Ecdysozoa 1,382), 480 s wall with 12 workers.
- FastTree-all (gene leaves) is worse than medoid on both metrics at 3× the cost; not recommended.
- The NJ rerun after the refactor reproduces 75.71 / 0.772 exactly.
- **Recommendation:** make `tree_method="fasttree"` with medoid leaves the default for the next frozen config. It is **not** switched in code yet, because the QfO run in progress uses the current default.
- At QfO scale, S4 estimates ~11.6 CPU-h for FastTree on all tree clusters (×2 with constrained synteny's two passes); one option is to keep the off-pass on NJ.

### 7.10 Presets (frozen config) and FastTree determinism
| arm (deterministic partition, gene order on) | OG-F | macro pair F1 | micro pair F1 | merge cost |
|---|---|---|---|---|
| FastTree on both passes | 80.36 | 0.799 | 0.804 | 1,039 s wall, 2.1 CPU-h |
| **on-pass FastTree, off-pass NJ** | 79.69 | 0.803 | 0.804 | 428 s wall, 1.05 CPU-h |

The NJ off pass is within noise by the PI's rule (ΔOG-F −0.67 > −1, Δmacro +0.004 > −0.005) at half the cost, so it is preferred.

**Presets.** `python -m ancestra run ... --preset {accurate,fast,custom}`; `PRESETS` and `apply_preset` are in `ancestra/pipeline.py`.
- **`accurate` (default of the frozen config):** local trees use FAMSA + FastTree -lg on length-medoid leaves. With `--gene-order`, the synteny-off pass of the constrained two-pass uses NJ.
- **`fast`:** NJ everywhere (the v0.2 behaviour); used for scaling claims.
- **`custom`:** the merge config is used as given via `--config`.
- The merge-stage key includes the preset and the off-pass tree method, so caches never mix (`test_pipeline_presets_do_not_share_cache`).
- `tree_ft_workers` follows `--threads` and is excluded from the key, because results do not depend on it.

**Check through the CLI.** A merge-only rerun on `run_det` with `--preset accurate` reproduces OG-F 79.69 (436 s wall). Manifest: `results/core_v0_orthobench/det_accurate_manifest.json`.

**Determinism.** FAMSA runs with `-t 1`, and FastTree is deterministic for a given input (fixed default seed). Pool results are collected in submission order. `test_localtree_fasttree_deterministic_across_workers` runs 2 vs 4 workers plus a repeat and gets identical outputs.

### 7.11 QfO recall: species-overlap labelling (pre-fixed 6-arm grid; dev only)
**Problem** (QfO peek 1, docs/qfo_run.md §8): recall is low (TreeFam-A dev TPR 0.572 vs OrthoFinder 0.732). FN classes: paralog call inside the same root HOG (2,411), different root HOGs (2,239), S0 split (1,286).

**The "coverage guard".** The QfO run used `min_cov` 0, and no other code path filters scored pairs by coverage: hybrid scoring has no coverage filter; S0 stage 2 (`s2_cov`) and the family rescue (`rescue_cov`) act on families, not pairs. So nothing rejects the 743 low-coverage scored FN pairs. A plausible soft cause is `norm_max` (bits / larger self-score), which penalises fragment vs full-length pairs. Arm A5 (`norm_min`) tests this.

**New options** (`MergeConfig`; off by default; omitted from `key()` at their defaults):
- `so_struct` (species-overlap structure): inside a group, same-child units with *disjoint species sets* are combined as a speciation, i.e. an orthologGroup at their species MRCA. Previously they formed a paralogGroup at the child. Overlapping units stay paralogs.
- `so_absorb` (species-overlap absorption): a still-alone unit whose best partner's group contains only species-disjoint units of its own child joins that group, as a missed speciation.
- Test: `test_species_overlap_missed_speciation`.

**Evaluation setup.**
- QfO: merge-only reruns on hard-linked det_go stages (`work/qfo_recall/base`), `fast` preset (NJ).
- TreeFam-A DEV is scored in memory by `scripts/qfo_treefam_dev_fast.py`: the dev reference pairs come from the peek-1 raw RefPhyloTest output, with the official formula. It reproduces the official dev values exactly (TP 8,009, FP 623, FN 5,984, TN 10,797). SwissTree and test families are never read.
- OrthoBench: `run_det` with the `accurate` preset base.
- Rule (PI): TreeFam-A dev TPR up, with PPV drop ≤ 0.01 and OrthoBench dev OG-F ≥ −1.

| arm | TF-A dev TPR | TF-A dev PPV | OrthoBench OG-F | OB macro pair F1 | S8 fission cand. | rule |
|---|---|---|---|---|---|---|
| A0 baseline | 0.5723 | 0.9276 | 79.69 | 0.803 | 706 | – |
| **A1 so_struct** | **0.6088** (+0.037) | 0.9229 (−0.005) | **79.69** (±0) | **0.806** | 512 | **pass** |
| A2 so_struct + so_absorb | 0.6113 | 0.9161 (−0.012) | 79.01 | 0.808 | 238 | fail (PPV) |
| A3 A2 + repair_ratio 0.5 | 0.6125 | 0.9205 (−0.007) | 78.14 (−1.55) | 0.800 | 161 | fail (OB) |
| A4 A2 + tol 0.8 | 0.6479 (+0.076) | 0.9156 (−0.012) | 71.15 (−8.5) | 0.803 | 269 | fail |
| A5 A2 + norm_min | 0.5866 | 0.9245 | 51.05 | 0.820 | 187 | fail (OB) |

**Merge cost on QfO** (`fast`): 6–7 min wall and ~0.1 CPU-h per arm, 9 GB.
- `accurate` on QfO is estimated at ~12 CPU-h, about 1 h wall at 12 workers for the on-pass. The estimate scales from OrthoBench: 177.7k QfO tree clusters vs 15.5k → ~100k FastTree jobs.

**Adopted by the rule:** A1 (`so_struct=True`). It is now in both presets (§7.12).

**Observations.**
- Looser mutual best (A4) recovers the most TreeFam recall but fuses OrthoBench OGs.
- `so_absorb` and the relaxed repair cut the S8 fission candidates by 66–77% but cost PPV or OG-F.

### 7.12 Deep-node rules (final small grid; base = A1)
**Generic definition of "deep".** A node is deep if the median distance −ln S of its mutual-best seeds is ≥ `deep_min_d` = 1.5. No taxon names are used.
- QfO values (det_go): cellular organisms 1.77, Eukaryota 1.81, Opisthokonta 1.68, Sar 1.82, Metamonada 1.68, Alveolata 1.61, Dikarya 1.59, Fungi 1.52 are deep.
- Not deep: Ecdysozoa 1.46, Bilateria 1.36, and every metazoan level below.
- On OrthoBench the maximum is 1.42 (Bilateria), so the rule is inactive there by construction. This was checked: no deep node, and OG-F is identical to A1.

**Options** (`MergeConfig`, omitted from `key()` at their defaults):
- `deep_tol`: the mutual-best tolerance at deep nodes.
- `deep_so_absorb`: species-overlap absorption (§7.11) at deep nodes only.
- Test: `test_deep_node_relaxed_tolerance`.

**Results** (QfO `fast` + so_struct; OrthoBench `accurate` + so_struct; rule vs A1):

| arm | TF-A dev TPR | TF-A dev PPV | OrthoBench OG-F | rule |
|---|---|---|---|---|
| A1 (base) | 0.6088 | 0.9229 | 79.69 | – |
| A6 deep tol 0.8 | 0.6494 (+0.041) | 0.8895 (−0.033) | 79.69 | fail |
| A7 deep tol 0.85 | 0.6347 (+0.026) | 0.8964 (−0.027) | 79.69 | fail |
| **A8 deep so_absorb** | **0.6163** (+0.008) | **0.9198** (−0.003) | 79.69 | **pass → adopted** |

**Adopted (superseded for `accurate` in §7.14).** Both presets set `so_struct=True`, `deep_min_d=1.5` and `deep_so_absorb=True`; after the accurate re-check, `accurate` keeps only the deep absorption. A plain `--preset fast` QfO rerun reproduces A8 (`results/core_v0_qfo/runs.tsv`, label `A8check`).

**Remaining FN, winner A8** (TreeFam-A DEV speciation pairs not predicted; `fn_breakdown_A8.json`; level = LCA of the pair):

| class | total | Bilateria | Eukaryota | Opisthokonta | other levels |
|---|---|---|---|---|---|
| different root HOG | 2,306 | 841 | 512 | 278 | 675 (Chordata 137, Amniota 99, ...) |
| paralog call in same root HOG | 1,776 | 572 | 285 | 169 | 750 |
| S0 split (different family) | 1,286 | 84 | 222 | 746 | 234 |

The largest remaining blocks are Bilateria-level misses (not deep under the rule) and Opisthokonta-level S0 splits.

### 7.13 Fix: OrthoXML nesting with so_struct / deep_so_absorb (docs/qfo_run.md §9)
**Bug.** With `so_struct`, species-disjoint units of one child were wrapped in an orthologGroup at their species MRCA. When the units were themselves orthologGroups at that same node, this gave same-level nesting, as in the repro O@Methanobacteriota{O@Methanobacteriota, O@Methanobacteriota}: two separate ancestral genes of a node, disjoint only through losses.
- In my QfO `fast` output: 8,907 same-level orthologGroup-in-orthologGroup cases and 261 paralogGroups at their enclosing orthologGroup's level. pyham crashed; pairs were unaffected.

**Fix** (`ancestra/merge.py`).
1. **`so_combine` (recursive).** Same-level orthologGroups are flattened into their children. The pieces are bucketed by the child of the MRCA they lie in, and each bucket is combined recursively at its own (lower) MRCA. A paralogGroup at the MRCA level is placed at its members' MRCA when that is lower, otherwise the units stay paralogs (fallback, counted as `so_struct_fallback`). A re-bucketed paralogGroup whose level is not strictly below its new orthologGroup is relabelled to its members' MRCA.
2. **Normalisation in `to_orthoxml` (`_fix_paralog_levels`).** A paralogGroup directly inside a paralogGroup is flattened into it. A paralogGroup whose level does not contain its members, or is not strictly below the enclosing orthologGroup, gets its members' MRCA (latest-possible placement). Ortholog and paralog pairs are unchanged (tested). This caught the remaining 49 cases of P-inside-P with ancestral levels.

**Tests.**
- `test_so_struct_nesting_regression`: the repro topology; exact expected HOG with 0 violations.
- `test_so_struct_rebucketed_paralog_level`.
- `test_fix_paralog_levels_flattens_and_relabels`: pairs unchanged.

**Verification** (`scripts/verify_orthoxml.py`, all checks):

| output | violations | pyham | other checks | dev metric |
|---|---|---|---|---|
| OrthoBench `run_det`, `accurate` | 0 | loads (15,233 HOGs) | pass | OG-F 79.69 |
| QfO `fast` merge-only (`work/qfo_recall/fast_fix3`) | 0 | loads (32,708 HOGs) | pairs count, flattening vs QfO `map_relations` and the OrthoBench export all pass | TreeFam-A dev TPR 0.6151 / PPV 0.9224 |

The QfO dev metric was 0.6163 / 0.9198 before the fix; the fallback to paralogs changes a few pairs.

### 7.14 Adoption re-check on the `accurate` base (QfO merge-only, after the nesting fix)
**Base.** `det_go_accurate`, i.e. `accurate` before the new options. Re-scored with `scripts/qfo_treefam_dev_fast.py`: TreeFam-A dev TPR 0.5854 / PPV 0.9413.
- Rule: TPR up and PPV drop ≤ 0.01. Exception allowed only for PPV drop ≤ 0.015.
- Arms: merge-only on the cached det_go stages, `accurate` (FastTree on-pass, NJ off-pass), gene order on. Each takes 47–58 min wall and ~7.0 CPU-h (61k FastTree jobs).

| arm | TPR | PPV | ΔTPR | ΔPPV | rule |
|---|---|---|---|---|---|
| C1 accurate + so_struct | 0.6441 | 0.9158 | +0.059 | −0.025 | fail |
| C2 accurate + so_struct + deep_so_absorb | 0.6548 | 0.9136 | +0.069 | −0.028 | fail (also outside the 0.015 exception) |
| **C3 accurate + deep_so_absorb** | **0.5908** | **0.9397** | **+0.005** | **−0.002** | **pass** |

**Decision (rule).**
- **`accurate`:** deep-node species-overlap absorption only (`so_struct=False`, `deep_min_d=1.5`, `deep_so_absorb=True`).
- **`fast`:** keeps `so_struct` + deep absorption, as adopted on the `fast` base (§7.11–7.12).
- On OrthoBench both options are inactive or neutral for `accurate`: OG-F 79.69.

**Why `so_struct` costs precision on the accurate base.** FastTree trees already resolve many in-paralog/ortholog decisions. Re-labelling species-disjoint paralog units as speciations then converts more true duplications (with differential loss) into ortholog pairs.

**Verification of the C2 output** (both options, `accurate`): `verify_orthoxml` overall PASS, 0 violations, pyham full load (32,021 HOGs).

### 7.8 Next steps
1. **Precision of tree mode and over-merges in large S0 families** (OG-P 63–75%):
   - check ambiguity clusters against FAMSA + FastTree trees on a sample;
   - add an outgroup-based witness test before tree merges.
2. **Root (Bilateria) decisions** have no outgroup and remain score-based. Options: bring in the deep-node S6 structure evidence, or a second outgroup from S0 families of other clades.
3. **QfO run** with the frozen v0.2 config (peek 1), including the min_cov question on multi-domain proteins.

## 8. Next steps (priority, v0 list; see §7.6 for the current list)
1. Numba for `decide_node` and `_aggregate` before QfO; QfO deep nodes have about 1e8 edges.
2. A rate-robust in-paralog test, i.e. a Scoredist-like calibrated distance plus the quartet with outgroup and the local-tree mode for flagged clusters. This targets failure 1, the largest FN class.
3. Against hidden paralogy (failure 4, all FPs): the relative-divergence test above, using the existing outgroup edges (`out_edges`), plus the OMA-style witness test. Then the synteny `context` hook (S5) with a shuffled-order control.
4. Run on OrthoBench dev (35 RefOGs) for the T+12h gate, and re-test the `min_cov` guard there.
5. S0 and S2 at QfO scale: order families in DIAMOND blocks by centroid similarity, so the cross-family rescue sees the hits that matter (r2 algorithms review, item 7).
