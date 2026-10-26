# Ancestra core v0 (sequence only): S2 scoring + S3 level-synchronous merge

Status: v0, toy-verified, 2026-10-26 (built T+0h50 to T+1h55, read from `scripts/clock.sh`).

- **Code:** `ancestra/io.py`, `ancestra/scoring.py`, `ancestra/merge.py`, `ancestra/pipeline.py` and `ancestra/__main__.py`.
- **Tests:** `tests/test_merge.py` has 16 tests. They cover exact hand-made HOGs with duplications, a polytomy, the coverage guard, the cap rule, the family join, and the pipeline end to end with resume.
- **Benchmark:** `scripts/run_core_v0_toy.py` writes `results/core_v0/toy_grid.jsonl` (504 records). `scripts/summarize_core_v0.py` turns it into `results/core_v0/toy_tables.md`.
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
  - If no TSV is given, the pipeline runs S0 through `ancestra.families` (`DEFAULT_FAMILY_CONFIG` in `pipeline.py`: DIAMOND cluster at 30% identity plus a reps all-vs-all, connected components, min-coverage 50%).
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
- **E-value:** 10 inside families. S0 has already established homology, and the coverage guard in S3 removes local junk. On the toy, e=1e-3 loses 0.3–0.4 recall points, and on 8-species sims 1–2 points (§5.3).
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
| score | norm_max | ≈ bits (+0.0005 F1); norm_min −0.001; **OrthoFinder-style (per species-pair length fit) −0.04 to −0.07 F1**, because it breaks the comparison of same-child vs cross-child scores that R1 and R2 rely on. |
| tol | 0.95 | Sweep: 0.9 −0.000 F1 but −0.9 P points; 0.98 / 1.0 −0.0015 / −0.003 F1. |
| inpar_rule / factor | ortholog / 1.0 | Sweep: 0.95 / 0.9 / 0.8 give +0.003 / +0.003 / +0.001 F1 at −0.8 / −1.4 / −2.4 P points, which fails the precision floor (pre-stated: P ≥ default − 0.003). |
| repair / ratio | on / 0.9 | Sweep: +0.0003 F1 (negligible on the toy, where the fissions have other causes, §6). |
| min_cov | 0.5 | Toy sweep favours 0.0 (+0.0034 F1). Kept at 0.5 as a documented exception: the toy has no multi-domain proteins, so it cannot measure the failure mode the guard exists for (DESIGN §9). Re-test on OrthoBench dev. |

Named arms are in `merge.ARMS`: `default`, `rbh`, `max`, `max+repair`, `avg`, `avg+repair`, `*/bits`, `*/quartet` and `avg+repair/either`.

### 3.4 Extension hooks (not used in v0)
- `run_merge(..., context=f)`: `f(node, a, b, s) -> bonus` is added to every cross-child candidate score before R1. This is the slot for synteny (S5): the per-node unit IDs and the edges are available, and the units' leaf members give their adjacencies.
- `run_merge(..., resolver=g)`: `g(node, child, lu, lv, score, labels) -> labels` may rewrite the partition. This is the slot for the local-tree mode on ambiguous clusters (FAMSA + FastTree + species overlap).
- The outgroup machinery (`out_edges` in `decide_node`) already gives every unit its scores to sibling-clade genes. It is reusable for the OMA-style witness test against hidden paralogy (r2 external PI, S1).
- `scoring.BACKENDS["foldseek3di"]`: a 3Di backend with the same `ScoreTable` interface. A second score column can be added to `leaf_scores`.
- `run_merge(..., trace={})` records, per node, each clade gene's ancestral ID before and after the node. `diagnose_core_v0.py` uses it.

## 4. Complexity
- **S2:** DIAMOND blocks ~N·B. The rescue cost is proportional to the pairs DIAMOND missed (5–13% of Σ g_f(g_f+1)/2 on the toy). The cross-family rescue re-scores only the joined families (≈11–14% of genes on the toy, 2–15 s).
- **S3:**
  - E scored leaf pairs are sorted once by LCA, O(E log E).
  - At node X, the work is aggregation of H_X = contracted child graphs ∪ the X bucket, O(|H_X| log |H_X|), plus a decision step O(|H_X| + n_X·k) and contraction.
  - Each leaf pair enters once. Contracted graphs only shrink, because aggregation is massive: on def32_s11, Σ_X |H_X| = 56k edges for 596k leaf pairs.
  - Worst case (nothing merges, many paralogs) is O(E·depth).
  - Memory is O(N + E). Leaf species do no merge work.
- **Measured:** the merge takes 0.5–1 s for 33k genes on 32 species (31 nodes, max 0.16 s per node). The decision step is pure Python dict code, about 3 µs per edge, so QfO (~1e8 edges) would need numba (planned).

## 5. Toy results (in-house simulator; guide tree = true species tree)
**Sims.** Seeds 11, 23 and 47, in four settings:
- `def8`: 8 species × 300 genes, default model.
- `hard8`: dup 0.5, loss 0.5, a WGD on a random branch, neo-functionalisation 2×.
- `def32` and `hard32`: 32 species × 1000 genes (33k proteins) with the same two models.

**Families.** `truth` uses the simulator root families and isolates S3; `s0` uses the default S0 clustering.

**Metrics.** Ortholog pairs use QfO flattening against the true pairs. HOGs are scored per level with `evaluate_toy`. All values are mean ± sd over the 3 seeds.

### 5.1 Headline: default configuration (hybrid + family rescue scoring, arm `default`)
| setting | families | P | R | F1 | root HOG exact | all-level HOG exact | mean best Jaccard | deep exact (d≤1) |
|---|---|---|---|---|---|---|---|---|
| def8 | truth | 0.993±0.001 | 0.968±0.008 | **0.980±0.005** | 0.939 | 0.964 | 0.985 | 0.943 |
| def8 | S0 | 0.993±0.001 | 0.929±0.016 | **0.960±0.008** | 0.865 | 0.944 | 0.978 | 0.890 |
| hard8 | truth | 0.987±0.004 | 0.938±0.020 | **0.962±0.012** | 0.826 | 0.894 | 0.955 | 0.828 |
| hard8 | S0 | 0.986±0.004 | 0.895±0.026 | **0.938±0.016** | 0.759 | 0.876 | 0.947 | 0.784 |
| def32 | truth | 0.992±0.002 | 0.978±0.005 | **0.985±0.003** | 0.960 | 0.979 | 0.991 | 0.956 |
| def32 | S0 | 0.992±0.002 | 0.956±0.004 | **0.974±0.003** | 0.904 | 0.973 | 0.990 | 0.920 |
| hard32 | truth | 0.986±0.001 | 0.954±0.002 | **0.970±0.001** | 0.909 | 0.941 | 0.977 | 0.899 |
| hard32 | S0 | 0.986±0.001 | 0.935±0.006 | **0.960±0.003** | 0.862 | 0.937 | 0.976 | 0.870 |

Without the cross-family rescue, F1 with S0 families is 0.940 / 0.917 / 0.943 / 0.924. With truth families the rescue makes **0 joins**: no spurious joins, which serves as a negative control.

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
- In-paralog absorption is essential. RBH alone loses 30–60% of the true pairs, and on hard32 its recall is 0.35–0.47.
- avg-link beats max-link by 0.3–0.5 F1 points, mostly on hard sims (max-link favours large groups).
- Fission repair is nearly neutral.
- The quartet and either rules trade precision for recall; `quartet` is best on hard8 (+0.4 F1) and worst for FP paralog pairs on hard32 (2× default).

With hybrid scoring the default reaches mean F1 0.974 (sweep).

### 5.3 S2 backends (truth families, avg + repair)
| setting | DIAMOND e10 | DIAMOND e1e-3 | MMseqs2 exhaustive | hybrid |
|---|---|---|---|---|
| def8 F1 | 0.961 | 0.961 | 0.979 | 0.980 |
| hard8 F1 | 0.936 | 0.935 | 0.962 | 0.961 |
| def32 F1 | 0.976 | 0.974 | 0.984 | 0.984 |
| hard32 F1 | 0.960 | 0.958 | 0.970 | 0.970 |

DIAMOND `--sensitive` reports 90–94% of the in-family pairs that MMseqs2 aligns, and every DIAMOND pair is also in the MMseqs2 set. The missed pairs are the divergent ones the deep merges need, so DIAMOND loses 1–5 recall points. The hybrid recovers all of it.

Wall time (shared node, load 30–80, 3 threads per job): on the toy, exhaustive MMseqs2 is *faster* than DIAMOND, because Σg² is small: def32 71 s vs 92 s, hybrid 115 s. At 24 threads the hybrid on def32_s11 takes 15–22 s.

The PI's QfO measurement found DIAMOND faster than MMseqs2 (126 s vs 457 s at 5e6 pairs), so the hybrid (DIAMOND + rescue of ~7–13% of pairs) is the scalable choice that keeps exhaustive recall.

### 5.4 S0 and the cross-family rescue
Default S0 recovers 83–84% of the truth families exactly, with 0 fusions and 45 (8-species) to 165 (32-species) fissions. That costs 4–8 recall points.

Joining families through the S2 cross-family hits (cov ≥ 0.5, e ≤ 1e-5) roughly halves the fissions (def32: 162 → 77; hard32: 165 → 69), lifts F1 by +2.0 to +3.6 points, and adds no fusions.

Looser variants (cov 0.3, e 1e-3) add another +0.1–0.2 F1 and create 0.3–1.3 fusions per sim (`exp_family_rescue.tsv`).

### 5.5 Runtime
`results/core_v0/example_def32_s11_s0fam_manifest.json` records a full CLI run on def32_s11 (33,044 proteins, 32 species, 24 threads, node load ~30):

| stage | wall | peak RSS |
|---|---|---|
| input | 0.02 s | 51 MB |
| families (S0) | 40.2 s | 70 MB (DIAMOND child: 598 MB) |
| scores (hybrid + rescue) | 22.4 s | 317 MB |
| merge (31 nodes) | 0.5 s | 128 MB |
| output (XML + XSD + pairs + OGs) | 1.7 s | 277 MB |

The end-to-end wall time is 58.6 s with a maximum RSS of 612 MB. The cost is dominated by S0 and S2.

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
3. **S0 fissions are permanent** unless they are repaired. This is the largest loss with S0 families (§5.4). The cross-family rescue halves it.
4. **Precision errors concentrate below very short internal branches.**
   - Examples: def32_s23 N15/N16 (branch length 0.001) carries 1,467 of 2,830 FP pairs; hard32_s11 N2 carries 3,068 of 7,318.
   - A duplication on a near-zero branch just above X is almost simultaneous with the speciation, so any score-based rule is near 50/50 there.
   - FP pairs are 100% true paralogs (never cross-family) and come from a few large unit pairs: on hard32_s11, 777 unit pairs make 7,318 FP pairs. This is where synteny and local trees should help.
   - Hidden paralogy from differential loss (a clean 1:1 between out-paralogs) is not yet tested for. The witness-test hook is ready (§3.4).
5. **HOG level of lineage-specific duplicates without a partner.** Two copies in one child with no partner in a sibling pass up separately and are joined only where they meet a common partner. Their paralogGroup is then placed at that higher level, or at the MRCA if it is unknown, rather than at the true duplication branch. Pair metrics are unaffected; per-level HOG exactness is slightly affected.
6. **Not tested by the toy.** HGT (the sims use none), guide-tree errors (the guide tree is the true tree), multi-domain proteins and fusion/fission genes, fragments and isoforms, and families above 132 genes (the cap rule has a unit test only).

## 7. Next steps (priority)
1. Numba for `decide_node` and `_aggregate` before QfO; QfO deep nodes have about 1e8 edges.
2. A rate-robust in-paralog test, i.e. a Scoredist-like calibrated distance plus the quartet with outgroup and the local-tree mode for flagged clusters. This targets failure 1, the largest FN class.
3. A witness test against hidden paralogy on 1:1 seeds using the existing outgroup edges (failure 4). Then the synteny `context` hook (S5) with shuffled-order control, starting with the precision errors on short branches.
4. Run on OrthoBench dev (35 RefOGs) for the T+12h gate, and re-test the `min_cov` guard there.
5. S0 owner: adopt the cross-family rescue, or order blocks by centroid similarity, so the rescue sees the right hits at QfO scale.
