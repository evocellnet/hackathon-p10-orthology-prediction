# R2 review: algorithms and compute

T+1h02. Tests: `r2_compute_tests/`. Load average 50–100 during the tests, so timings are upper bounds.

## Round-1 items
- **M1: addressed.** The dense-clade scaling set is still missing.
- **M2, M4: moot.** v1 never searches ancestral genes.
- **M3: OK.** §8 says both "k-way" and "binarisation"; pick one.
- **M5: open.** `timings.tsv` still records all 0.0.
- **M6: addressed.**

## Measured on QfO (N = 988,778)
- **DeepClust (defaults):** 21 min, 9.2 GB. 328k clusters, 237k of them singletons. Largest cluster 2,871. Σg² = 5e7.
- **Fission:** only **36.8%** of FastOMA root HOGs lie fully inside one cluster (20.8% of HOGs with ≥10 genes). The gene-weighted largest share is 0.64. The defaults (centroid clustering, 80% member cover) give subfamilies, and fission at S0 cannot be undone.
- **Giant family:** connected components of OrthoFinder pairs include one of **115k genes**, which carries 98% of Σg² = 1.35e10. Single linkage needs a cap.
- **Backend test** (85k proteins, 5.0e6 in-family pairs):

  - `mmseqs align` on a custom prefilter (exhaustive reference): 457 s.
  - DIAMOND, all 85k in one block, `-k 0`, in-family filter: mid-sensitive 36 s (86.6% of pairs), sensitive 126 s (92.9%).

## v1 cost on QfO (24 threads)
| Stage | Estimate |
|---|---|
| DeepClust | 21 min |
| Joining clusters (328k reps × reps, DIAMOND) | ~9 min mid, ~30 min sensitive |
| Scoring: families packed into blocks of B ≈ 100k; cost ≈ N·B, not Σg² | ~6 min mid, ~25 min sensitive |
| Edges | 1–3e8 edges, 1–4 GB |
| Merges | Each leaf edge is consumed once, at the LCA of its two species, so the total is **O(E)**: minutes in numba. No per-node search, so the r1 critical-path problem disappears |

- **Total ≈ 1–1.5 h.** All-vs-all on this node takes ~4.8 h at sensitive and ~1.4 h at mid.
- The real saving at 78 species is **~4× at sensitive and ~1× at mid**; S0 dominates.

## Thousands of proteomes (N ≈ 3e7)
- DeepClust ~10 h (close to linear).
- Block scoring: N·B = 3e12, so ~4 h mid and ~15 h sensitive.
- Joining by reps × reps becomes quadratic; recluster the reps with lower cover instead.
- **The real blocker is Σ_f g_f².** g_f grows with S, so storing every leaf pair is quadratic in S: a 2,000-species family gives ~4e6 edges, a kinase-like family ~1e10. The cap must bound edges and alignments.

## Changes
**MUST**
1. **Two-level S0.**
   - DeepClust with defaults, yielding subfamilies.
   - DIAMOND `--sensitive` on centroids plus singletons.
   - Link two subfamilies only when the hit covers ≥50% of the shorter sequence in both directions (anti-chaining).
   - Take connected components. Split any component larger than G_max = 20k genes with Leiden on the centroid graph, and flag the split parts for root rescue.
   - Tune toward **fusion, not fission**.
   - Metric: gene-weighted RefOG completeness on OrthoBench dev. FastOMA agreement is a labelled diagnostic only.
2. **S2 backend.** DIAMOND on blocks of whole families packed to B ≈ 100k. A family larger than B is its own block.
   - Flags: `-k 0 --max-hsps 1 -b 4 -c 1 --sensitive`, then filter to in-family pairs. Mid-sensitive is an S2 arm.
   - Drop pyopal/parasail: they are not installed.
   - `mmseqs align` on a custom prefilter is used only as the exhaustive recall reference on a sample, and for targeted rescue.
3. **Concrete cap.**
   - g ≤ 5k: keep all in-family hits.
   - 5k < g ≤ G_max: stream-reduce to the top 2 hits per (query, target species), plus up to 20 same-species hits.
   - g > G_max: split as in item 1.
   - Log the edges dropped per family.
4. **Calibration.** Fill `timings.tsv` from these numbers, rerun on an idle node, and log the load average with every timing.

**SHOULD**
5. Implement merges as LCA-bucketed edge aggregation, which is O(E). State this in §1.
6. **Thousands-scale configuration.**
   - Store leaf edges only up to a depth threshold.
   - Above that threshold, search medoids within the family. The cost is Σ_f a_f·b_f in ancestral genes.
   - Add the ~200 Enterobacterales dense set to the scaling runs.
7. Order families within blocks by centroid similarity, so cross-family hits flag S0 fission for free; report the rate.

**COULD**: report S0's share of total CPU time.
