# R1 review: algorithms and compute

Reviewed 2026-10-26 15:10 CET (T+0h30). Evidence: `r1_anc_sizes_fastoma_proxy.{py,tsv}` here.

## 1. The complexity claim fails on QfO

Merge cost is C = Σ_v Σ_{i<j} k²·n_i·n_j. Here n is the ancestral proteome size and k the number of representatives per ancestral gene. With "unmatched passes up", n is not bounded.

As a proxy for n, I counted per node the maximal FastOMA QfO-2022 HOGs plus the clade's singletons, on the QfO lineage tree:

| node | leaf genes | ancestral (proxy) | cross pairs, k=1 |
|---|---|---|---|
| Bilateria | 414k | 258k | 1.3e10 |
| Opisthokonta (3-way) | 559k | 337k | 1.9e10 |
| Eukaryota (8-way) | 866k | 534k | **9.5e10** |
| LUCA (3-way) | 966k | 574k | 2.7e10 |
| total, 44 merges | | | **2.0e11** |

All-vs-all on the same genes is N²/2 = 4.9e11.

- **k=1 saves only ~2.5× over all-vs-all. k=3 costs 3.6× more than all-vs-all.**
- Four deep nodes carry 78% of the cost. The Eukaryota polytomy alone carries 48%.
- Orphans are ~60% of the deep ancestral units: FastOMA leaves 341k of 989k genes as singletons.
- Pre-registered claim (a), linear runtime on 8→78 QfO species, will very likely fail. Adding diverse species grows n roughly linearly, so C grows close to S².

## 2. Wall-time estimates (24 threads)

Throughput from ENVIRONMENT.md, human self-search (4.25e8 pairs): sensitive 3.5e7 pairs/s, mid-sensitive 6e7, fast 2.7e8. These are extrapolated from 20k² and need calibration (M5).

| run | pairs | sensitive | mid | fast |
|---|---|---|---|---|
| QfO all-vs-all (N²) | 9.8e11 | 7.8 h | 4.5 h | 1.0 h |
| leaf self-searches Σ G² | 2.2e10 | 10 min | 6 min | 1 min |
| progressive, k=1, one direction | 2.0e11 | **1.6 h** | 55 min | 12 min |
| progressive, k=1, reciprocal | 4.0e11 | 3.2 h | 1.9 h | 25 min |
| progressive, k=3 | 1.8e12 | 14 h | 8 h | 1.9 h |
| OrthoBench, 12 species (estimate) | ~1e10 | ~5 min | | |

**Parallelism.** Every expensive merge lies on one sequential critical path (…→Bilateria→Opisthokonta→Eukaryota→LUCA). Parallel subtrees therefore only help the ~20% of shallow work and the Python merge logic. Give large DIAMOND runs all 24 threads. Pack small merges at 4–6 threads each.

**Memory.** RAM is fine (<50 GB at -b 4). The real risk is hit volume. At ~37 hits per query, a deep merge yields more than 50M rows. Cap -k and store hits as integer IDs in parquet.

## 3. Algorithmic issues

1. **In-paralog self-scores.** A full ancestral self-search costs at least as much as the cross search. Run leaf self-searches once. At internal nodes, align only candidate pairs (a, a′ that share a best hit b) in a targeted DIAMOND run. That takes seconds.
2. **Reciprocal search.** Search one direction only. Re-check only the hit targets in reverse.
3. **Polytomies.** The tree has degree-3 (7), degree-4 (3), degree-8 (Eukaryota) and degree-16 (Bacteria) nodes. All-pairs merging is quadratic in degree, and multi-way InParanoid is undefined. Binarize with a guide (S1 distances or size order) and mark the synthetic nodes. Dissolve them in OrthoXML so that no unknown TaxRange appears; their duplications move to the real node.
4. **Representatives.** k=3 everywhere is unaffordable. Default to k=1 (medoid). Add representatives only in a rescue pass: unmatched multi-member genes searched against the other side's unmatched genes, capped at 1.2·n.
5. **Orphans.** Use orphans as queries but not as targets at deep nodes. This saves f_A·f_B, about 35%. Run it as an ablation with a deep-recall check.

## 4. Thousands of proteomes, AFDB, GPU

Dense sampling within a clade keeps n near the pan-proteome, so the cost per added species stays bounded. Diverse sampling does not.

At ~2,000 eukaryotes, deep merges reach 1e12–1e13 pairs, which is tens of hours each. Getting there needs:
- an orphan policy;
- a family-restricted candidate search (a linear linclust or deepclust prefilter);
- `diamond --multiprocessing` across nodes;
- binarization of NCBI polytomies with hundreds of children.

Claim a measured cost model rather than "linear".

On GPU: DIAMOND has no reliable GPU path. MMseqs2-GPU and Foldseek-GPU can be swapped in, and ProstT5 is the most GPU-bound step. Put search behind one interface now, and do no GPU work.

## 5. v0 in about 8 h

**Stack.** Python 3.11 with numpy, pyarrow and numba (polars and parasail are not installed). DIAMOND 2.2.8 runs as a subprocess. A ~150-line post-order driver with `ProcessPoolExecutor` and per-job thread budgets; no Snakemake for v0.

**IDs.**
- Global uint32 gene IDs, with FASTA headers holding only that integer.
- An explicit species table (UP ID, taxid, 5-letter code) with a test that all 78 leaves map. My proxy found an empty Murinae child, a real code mismatch.

**Per-node directory `nodes/<id>/`:**
- `anc.parquet` (anc_id, rep_gid, n_members, event, child refs)
- `map_<child>.npy`, an int32 array child→parent. Composing these gives membership at any level.
- `hits.parquet` (qid, sid, bitscore, qcov, scov, pident)
- `reps.faa`
- `stats.json` (sizes, pairs, hits, wall time, peak RSS via `/usr/bin/time -v`)

**DIAMOND flags.** `--outfmt 6 qseqid sseqid bitscore pident qcovhsp scovhsp evalue -e 1e-3 -k 50 -b 4 -c 1 --threads 24`. Fixed -b, -c and thread count keep runs reproducible. Sensitivity per depth is an S2 variant.

**Checkpoints.**
- Each node is written to a temp dir and renamed when complete.
- `DONE` holds a hash of the config plus the children's hashes.
- Queries over 200k are split into 50k chunks, so a deep merge can resume part-way.

**Output.** Stream OrthoXML with `lxml.etree.xmlfile`. Export pairs with numba (cross-product at each speciation node). In CI on the toy data, check XSD and pyham, and check that pairs from the XML equal pairs from the internal tables.

**Hour plan.**
- h0–1: I/O and binarization.
- h1–3: merge kernel.
- h3–4: driver and checkpoints.
- h4–5: OrthoXML and pairs.
- h5–6: OrthoBench end-to-end.
- h6–8: QfO run, k=1, mid-sensitive, ~1–2 h in the background.

## Proposed changes

- **M1 (must).** Replace "linear in S" and claim (a) with the cost model C. Report C, wall time and RSS per node against all-vs-all N². Fit the scaling exponent on nested QfO subsets and on a dense set (~200 Enterobacterales reference proteomes). Reason: §1.
- **M2 (must).** Default k=1, with capped extra representatives in a rescue pass only. Reason: k=3 is worse than all-vs-all.
- **M3 (must).** Binarize polytomies, mark the synthetic nodes, and dissolve them in OrthoXML.
- **M4 (must).** No ancestral self-searches: run leaf self-searches once, align only targeted in-paralog candidates, and search one direction.
- **M5 (must).** Spend 15 min calibrating DIAMOND (100k×100k and 500k×50k, mid and sensitive) before QfO. Reason: current throughput is extrapolated.
- **M6 (must).** Atomic, hash-keyed node checkpoints plus chunked large nodes. Reason: deep merges take hours.
- **S1 (should).** Orphans as queries only at deep nodes, run as an ablation. Reason: saves ~35%.
- **S2 (should).** Sensitivity by node depth.
- **S3 (should).** A search-backend interface (DIAMOND, MMseqs2, Foldseek, later GPU).
- **C1 (could).** Fix `data/calib/timings.tsv`, which records all 0.0.
- **C2 (could).** A linclust family-prefilter variant for the scalability section.
