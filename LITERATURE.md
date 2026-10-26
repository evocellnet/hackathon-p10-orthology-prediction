# LITERATURE.md: Orthology inference, scalable homology, extra data layers, benchmarks

Compiled 2026-10-26 by the literature agent. Search tools: WebSearch/WebFetch, the Crossref API (to check DOIs), the GitHub repos themselves (cloned and read), and the EBI/Zenodo/B2SHARE/OpenEBench APIs (queried directly).
Verification marks:
- **[V]**: DOI, venue and year checked against Crossref or the publisher, or the resource was opened and inspected today.
- **[P]**: preprint only, as far as I could find.
- **[U]**: could not verify. Treat these as uncertain.

Every URL listed under "Key URLs" was opened or HEAD-checked on 2026-10-26.

---

## 0. Executive summary (read this first)

1. **The field leaders now scale by placing sequences against a reference.**
   - **OrthoFinder 3** (Nature Methods 2026) runs a full analysis on about 64 "core" species. It then assigns every other species with DIAMOND against SHOOT-style profiles, and the paper reports near-linear time from 128 to 1024 species (1024 species took 128 h on 16 cores).
   - **FastOMA** (Nature Methods 2025) places proteins into an existing OMA HOG database with OMAmer k-mers, then resolves HOGs bottom-up along the species tree using small subsampled trees. It ran 2,086 genomes in under 24 h on 300 cores.
   - **SonicParanoid2** (Genome Biology 2024) uses a machine-learning model to skip half of each pairwise alignment, but it stays quadratic and produces no HOGs.
   - **eggNOG v7** (NAR 2026) builds Pfam-domain-centric trees and runs a species-overlap algorithm (pbDOOD). It is a database, not a de novo tool.
   - No published method uses 3D structure, synteny or expression *inside* inference, and none reports results on QfO with such layers. That is the open space.
2. **Competitor predictions on the QfO 2022 reference proteomes can be downloaded.** I found 21 participant datasets on EUDAT B2SHARE (OrthoFinder, FastOMA, OMA HOGs/Groups/Pairs, SonicParanoid2, PANTHER 18, Ensembl Compara, Hieranoid 2, InParanoid 5, Domainoid+, OrthoInspector 3.5, PhylomeDB 5, MetaPhOrs, OrthoFFGC, RSD, BBH, and pbDOOD/eggNOG v7 from June 2026).
   - Their per-benchmark metric values are retrievable from the OpenEBench API.
   - I harvested them into the table in §7.4. The script is `scripts/fetch_qfo_public_results.py`.
   - **We do not need to rerun competitors to place ourselves on the QfO plots. We do need to run some of them for runtime comparisons.**
3. **The QfO benchmark can most likely run without Docker.** The pipeline is Nextflow plus two kinds of container: Python, and "Darwin" (a closed-source interpreter that ships as a static binary from biorecipes.com). The repo also has a `singularity` profile.
   - Inside our container there is no singularity, apptainer or docker binary (checked).
   - Practical route: install Darwin plus the Python requirements natively from the Dockerfiles, and run the test scripts directly or with Nextflow and no profile. udocker or apptainer from conda are fallbacks. Details in §7.2.
4. **There is a data-version trap.** The current QfO FTP holds **release 2026_02** (81 species). The benchmark reference data on Zenodo is for **2022** (78 species).
   - We must predict on `QfO_release_2022_02_with_updated_UP000000437.tar.gz` (2.4 GB), which now lives under `previous_releases/`. The URL hard-coded in OpenEBench returns 404.
   - The 2026 release adds per-species **gff3/gbff and `*_gcoord.xml` gene coordinates**, which are useful for synteny. Coordinates for 2022 accessions have to be mapped by UniProt accession.
5. **HOG output is a double-edged sword on QfO.** QfO flattens HOGs into all pairs that the nested orthologGroups imply.
   - On the 2022 benchmarks, OMA_HOGs have much lower precision than OMA_Pairs (GO Schlicker 0.417 vs 0.496, EC 0.829 vs 0.952, FAS 0.411 vs 0.810), with higher recall.
   - Over-merged deep HOG levels are penalised hard. Our HOG construction must be conservative at deep levels, or we must report both a HOG export and a pair export.

---

## 1. State-of-the-art orthology methods and their scaling

The QfO standing numbers below come from the 2022 public results in §7.4, unless stated otherwise.

**OrthoFinder v2**: Emms & Kelly 2019, *Genome Biology* 20:238, doi:10.1186/s13059-019-1832-y **[V via PMC6857279]**. v1: Emms & Kelly 2015, *Genome Biol* 16:157, doi:10.1186/s13059-015-0721-2 **[V]**.
- How it works: all-vs-all DIAMOND, length- and distance-normalised bit scores, MCL orthogroups, a gene tree per orthogroup (DendroBLAST or MSA), STAG species tree, STRIDE rooting, and a duplication-loss-aware species-overlap reconciliation that yields orthologs and per-node "phylogenetic HOGs".
- It is quadratic in the number of species.
- QfO 2022 standing: highest recall among de novo tools (SwissTree TPR 0.825, TreeFam-A 0.737), with mid precision.

**OrthoFinder v3**: Emms, Liu, Belcher, Holmes & Kelly 2026, "OrthoFinder: improved phylogenetic orthology inference with enhanced accuracy and scalability", *Nature Methods*, doi:10.1038/s41592-026-03126-6 **[V]**. Preprint: bioRxiv 10.1101/2025.07.15.664860. Code: github.com/OrthoFinder/OrthoFinder (v3.2).
- (i) Orthogroups are split at "ancient" duplications that predate the root, after reconciling each OG tree. This gives a 7% relative gain in OrthoBench accuracy.
- (ii) `--core`/`--assign` workflow:
  - Users run a normal OF on fewer than 100 (default 64) core species.
  - For each core OG, k-means on the MSA embedding picks representative sequences. These become DIAMOND profiles (the SHOOT approach; Emms & Kelly 2022, *Genome Biol* 23:85, doi:10.1186/s13059-022-02652-8 **[V]**).
  - New genes are assigned to OGs. Unassigned genes are clustered within clades that lack core representation.
  - FAMSA + FastTree `-fastest` build trees, and ASTRAL-Pro builds the species tree.
- Reported results:
  - Near-linear runtime from 128 to 1024 species. It was the only tool to finish 1024 Ensembl species within 7 days.
  - About 4x lower peak RAM than competitors.
  - On OrthoBench it beats Broccoli, OrthoHMM, SonicParanoid2, OrthoMCL, ProteinOrtho, Hieranoid and FastOMA. FastOMA fused 3.5x more RefOGs (10% vs 2.85%).
  - On QfO it sits on the Pareto front for STD Euk/Bac, EC, SwissTree, TreeFam-A and VGNC. STD-Euk recall is 15,721 vs 8,686 for FastOMA, at a similar RF.
- **The OF3 QfO predictions do not appear to be public.** Only the 2023 "OrthoFinder" submission is on B2SHARE.
- No OrthoXML output is documented. Output is Orthogroups.tsv plus per-node HOG tsv.
- *Takeaway:* this is the strongest accuracy and scalability competitor. Its scaling depends on a core reference and profile assignment, which is still DIAMOND search, sub-quadratic only because the core size is fixed.

**OMA / FastOMA**: Majidian, Nevers, Yazdizadeh Kharrazi, Warwick Vesztrocy, Pascarelli, Moi, Glover, Altenhoff & Dessimoz 2025, "Orthology inference at scale with FastOMA", *Nature Methods* 22:269, doi:10.1038/s41592-024-02552-8 **[V]**. Code: github.com/DessimozLab/FastOMA (Nextflow).
- Step 1 (root-HOGs): OMAmer k-mer placement of every protein into the reference OMA HOG database (LUCA.h5 from omabrowser.org). Unplaced proteins are clustered with Linclust.
- Step 2: bottom-up traversal of the *input species tree* (NCBI taxonomy by default). At each node, MAFFT + FastTree runs on subsampled members (default 20 per HOG), and the species-overlap rule labels duplications vs speciations.
- Output: OrthoXML HOGs.
- Runtime: linear in genome count, 2,086 eukaryotes in under 24 h on 300 cores.
- QfO standing (2022, v0.1.4): SwissTree PPV 0.953 (top), TPR 0.699; STD-Euk recall is low.
- Dependencies and caveats:
  - Its accuracy depends on the OMA reference DB, which contains the QfO species (a leakage caveat for benchmarking; my inference, not discussed in the paper as far as I saw).
  - It also depends on the resolution of the input species tree.
- OMA database: Altenhoff et al. 2024 *NAR* 52:D513, doi:10.1093/nar/gkad1020 **[V]**. OMAmer: Rossier, Warwick Vesztrocy, Robinson-Rechavi & Dessimoz 2021, *Bioinformatics* 37:2866, doi:10.1093/bioinformatics/btab219 **[V]**. GETHOGs (HOGs from pairs): Altenhoff et al. 2013, *PLoS ONE* 8:e53786, doi:10.1371/journal.pone.0053786 **[V]**.
- *Takeaway:* "conceptually different from FastOMA" means we must not just place proteins into a pre-built OMA DB and run per-level species overlap.

**SonicParanoid2**: Cosentino, Sriswasdi & Iwasaki 2024, *Genome Biology* 25:195, doi:10.1186/s13059-024-03298-4 **[V]**.
- For each proteome pair, an AdaBoost model predicts which of the two directional alignments is faster. Only that one runs fully, and the reverse is restricted to hits. This cut alignment time by 95% on QfO.
- A Doc2Vec "language model" on Pfam domain architectures adds orthologs (needs at least 70% Pfam coverage). Pair graphs are merged and clustered with MCL.
- Results: 2,000 MAGs in about 1 day; 200 eukaryotes in about 6 h.
- Still O(N²) proteome pairs. Output is flat OGs and pairs, no HOGs.
- QfO: claims top aggregate. In the 2022 table, the sens-g mode is balanced and the sens mode has high recall but lower precision.

**Hieranoid / Hieranoid 2**: Schreiber & Sonnhammer 2013, *J Mol Biol* 425:2072, doi:10.1016/j.jmb.2013.02.018 **[V]**. Kaduk & Sonnhammer 2017, *Bioinformatics* 33:1154, doi:10.1093/bioinformatics/btw774 **[V]**.
- Walks up a guide species tree. At each internal node it runs InParanoid between the two child (pseudo-)proteomes, then builds an ancestral consensus proteome that represents each ortholog group.
- **Number of comparisons = N−1 (linear in species)**, and the output is hierarchical by construction.
- QfO 2022: best STD-Eukaryota RF (0.034), SwissTree PPV 0.950, with moderate recall.
- *Takeaway:* the most underexploited scalable paradigm. A modern "progressive ancestral-profile" approach (profiles, embeddings, 3Di, synteny anchors) would be conceptually distinct from OF3 and FastOMA.

**InParanoid-DIAMOND / InParanoid 5 (QfO label)**: Persson & Sonnhammer 2022, *Bioinformatics* 38:2918, doi:10.1093/bioinformatics/btac194 **[V]**. Domainoid: Persson et al. 2019 *BMC Bioinf* 20:523, doi:10.1186/s12859-019-3137-2 **[V]**.
- Pairwise, so O(N²). High precision. Domainoid+ adds domain-level orthology.

**PANTHER / TreeGrafter**: Thomas et al. 2022 "PANTHER: Making genome-scale phylogenetics accessible to all", *Protein Sci* 31:8, doi:10.1002/pro.4218 **[V]**. TreeGrafter: Tang et al. 2019 *Bioinformatics* 35:518, doi:10.1093/bioinformatics/bty625 **[V]**.
- About 15k curated family HMMs with reference trees. New sequences are grafted onto them.
- Two QfO entries: "all" (high recall, SwissTree TPR 0.811) and "LDO" (least-diverged orthologs; high precision).
- A grafting model is effectively linear per query but needs a curated reference.

**eggNOG-mapper v2**: Cantalapiedra et al. 2021, *MBE* 38:5825, doi:10.1093/molbev/msab293 **[V]**.
- Assigns queries to precomputed eggNOG OGs by DIAMOND/MMseqs2/HMMER. Annotation-oriented.

**eggNOG v7 / pbDOOD**: Hernández-Plaza et al. 2026, *NAR* 54:D402, doi:10.1093/nar/gkaf1249 **[V]**. Code: github.com/AnaHrnndz/pbdood.
- Pfam `cut_ga` pre-clustering, plus MMseqs2 de novo subclustering of huge clusters. MAFFT/FAMSA → FastTree2 → MinVar rooting.
- A noise-tolerant species-overlap algorithm (SO ≥ 0.1, with outlier removal) defines duplication nodes and dated OGs.
- Scale: 59.3M proteins, 12,535 species.
- **pbDOOD was submitted to QfO in June 2026.** It has the highest SwissTree TPR (0.85) at PPV 0.95, and the largest STD sampling counts, but a worse STD RF (G_STD2 Euk 0.35).
- *Takeaway:* domain-centric trees plus species overlap is the current recall leader.

**Broccoli**: Derelle, Philippe & Colbourne 2020, *MBE* 37:3389, doi:10.1093/molbev/msaa159 **[V]**.
- k-mer clustering, then DIAMOND within clusters, then fast NJ trees, then a label-propagation network. Handles chimeric proteins. Quadratic within clusters.

**Proteinortho6**: Klemm, Stadler & Lechner 2023, *Front Bioinf* 3:1322477, doi:10.3389/fbinf.2023.1322477 **[V]**.
- A pseudo-reciprocal best-alignment heuristic (avoids a full all-vs-all).
- Its PoFF extension uses **synteny**: Lechner et al. 2014, *PLoS ONE* 9:e105015, doi:10.1371/journal.pone.0105015 **[V]**.

**PhylomeDB V5**: Fuentes et al. 2022, *NAR* 50:D1062 [pages U], doi:10.1093/nar/gkab966 **[V DOI]**.
- Seed-proteome phylomes with species-overlap orthology. MetaPhOrs merges multiple phylomes. Species overlap originates in Huerta-Cepas et al. 2007, *Genome Biol* 8:R109, doi:10.1186/gb-2007-8-6-r109 **[V]**.

**Ensembl Compara GeneTrees**: Vilella et al. 2009, *Genome Res* 19:327, doi:10.1101/gr.073585.107 **[V]**.
- Clustering, then MSA, then TreeBeST gene/species tree reconciliation, giving orthologs and paralogs.
- Recent releases cluster via an HMM library **[U on exact current details]**.
- Ships GOC and WGA orthology QC scores (§5).
- QfO 2022: lower precision (EC 0.840).

**OrthoMCL**: Li, Stoeckert & Roos 2003, *Genome Res* 13:2178, doi:10.1101/gr.1224503 **[V]**. BLAST + MCL; historical baseline.

**OrthoInspector 3**: Nevers et al. 2019, *NAR* 47:D411, doi:10.1093/nar/gky1068 **[V]**. Pairwise BRH-based; strong on QfO precision.

**OrthoHMM**: Steenwyk et al. 2024, bioRxiv 10.1101/2024.12.07.627370 **[P]**. Profile-HMM all-vs-all with length and phylogeny bias correction. Claims OrthoBench gains, but OF3 reports lower precision.

**TOGA / TOGA2 (genome-alignment-based orthology)**: see §6.

**HOGPROP**: a function/trait propagation algorithm *over* HOGs (QTLSearch). It is not an inference method.
- Source: Warwick Vesztrocy et al. 2020 **[U: exact citation not verified]**. Summarised in Sarton-Lohéac et al. 2025, *J Mol Evol*, doi:10.1007/s00239-025-10277-1 **[V]**.

**PLM/embedding-, structure- or graph-ML-based orthology (2025–2026)**
- Despite targeted searches, **I found no published orthology inference method that uses protein-language-model embeddings, 3Di or graph neural networks as its core and reports QfO/OrthoBench results.**
- Closest work:
  - SonicParanoid2's Doc2Vec on domain architectures.
  - Minotto et al. 2025, "Exploring homology detection via k-means clustering of proteins embedded with a large language model", *Bioinformatics*, doi:10.1093/bioinformatics/btaf472 **[V]**.
  - An MDPI *Proteomes* 2026 paper using ESM-2 distances to characterise (not infer) zebrafish ortholog divergence, doi:10.3390/proteomes14030036 **[V, low relevance]**.
  - Two 2025 *J Mol Evol* QfO perspective papers discuss AI and orthology: Glover et al. 2025 "Perspectives on Orthology During the Quest for Orthologs", doi:10.1007/s00239-025-10290-4 **[V]**; "Quest for Orthologs in the era of Data Deluge and AI", doi:10.1007/s00239-025-10272-6 **[V title; authors/content not read, paywalled]**.
- *Takeaway:* an embedding or structure-based, sub-quadratic HOG method would be genuinely novel. The absence of such methods may also mean they did not beat sequence trees on QfO, so we should not assume a gain.

**Background reviews:** Sarton-Lohéac et al. 2025 *J Mol Evol* (HOG review, see above). It explicitly notes that **no benchmark evaluates hierarchical features** such as duplication timing or ancestral gene content. QfO flattens HOGs to pairs.

---

## 2. Sub-quadratic homology search and clustering

**MMseqs2**: Steinegger & Söding 2017, *Nat Biotechnol* 35:1026, doi:10.1038/nbt.3988 **[U: standard citation from memory, not re-checked]**.
- k-mer prefilter, then ungapped, then gapped alignment. `cluster` is cascaded.

**Linclust**: Steinegger & Söding 2018, *Nat Commun* 9:2542, doi:10.1038/s41467-018-04964-5 **[V]**.
- **Linear time**: each sequence picks m k-mers by hash, sequences sharing a k-mer group under the longest member, and only center–member pairs are aligned.
- Sensitivity is limited to roughly ≥50% identity unless cascaded. Used by FastOMA for unplaced proteins.

**MMseqs2-GPU**: Kallenborn et al. 2025, *Nat Methods*, doi:10.1038/s41592-025-02819-8 **[V]**. Not usable here (no GPU).

**DIAMOND v2**: Buchfink, Reuter & Drost 2021, *Nat Methods* 18:366, doi:10.1038/s41592-021-01101-x **[U: from memory]**.

**DIAMOND DeepClust**: Buchfink, Ashkenazy, Reuter, Kennedy & Drost 2026, "Clustering the protein universe of life using DIAMOND DeepClust", *Nat Methods*, doi:10.1038/s41592-026-03030-z **[V]**. Preprint 10.1101/2023.01.24.525373.
- Cascaded clustering down to about 20% identity. 19 billion proteins became 544M clusters.
- *Takeaway:* best available tool for sensitive, near-linear family-level clustering (`diamond deepclust` / `linclust`). A candidate step for building root-level families.

**Foldseek**: van Kempen et al. 2024, *Nat Biotechnol* 42:243, doi:10.1038/s41587-023-01773-0 **[V]**.
- 3Di alphabet (20 states describing tertiary interactions) plus AA. A k-mer prefilter like MMseqs2. Can also cluster.

**Foldseek cluster / AFDB clusters**: Barrio-Hernandez et al. 2023, *Nature* 622:637, doi:10.1038/s41586-023-06510-w **[V]**. 2.3M structural clusters of AFDB.

**pLM embedding search**
- *DHR (Dense Homolog Retriever)*: Hong et al. 2024, *Nat Biotechnol*, doi:10.1038/s41587-024-02353-6 **[V]**.
  - Dual-encoder retrieval. Claims more than 10% higher sensitivity and 22x faster than PSI-BLAST/DIAMOND at superfamily level.
  - GPU-oriented **[U on CPU speed]**.
- *PLMSearch*: Liu et al. 2024, *Nat Commun* 15:2775, doi:10.1038/s41467-024-46808-5 **[V]**.
  - ESM-1b embeddings plus a structural-similarity predictor. MMseqs2-like speed and about 3x sensitivity for remote homology.
- *TM-Vec / DeepBLAST*: Hamamsy et al. 2024, *Nat Biotechnol*, doi:10.1038/s41587-023-01917-2 **[V]**. Predicts TM-score from sequence embeddings; supports vector indexing.
- *pLM-BLAST*: Kaminski et al. 2023, *Bioinformatics* 39:btad579, doi:10.1093/bioinformatics/btad579 **[V]**. Local alignment on embeddings.
- *ESM-2 positional-embedding alignment*: Johnson, Peshwa & Sun 2024, *eLife* 13:RP91415, doi:10.7554/eLife.91415.3 **[V]**. 1-byte per-residue embeddings give much higher remote-homology sensitivity.
- *EBA*: Pantolini et al. 2024, *Bioinformatics* 40:btad786, doi:10.1093/bioinformatics/btad786 **[V]**.
- ESM-2: Lin et al. 2023, *Science* 379:1123, doi:10.1126/science.ade2574 **[V]**.
- *Caveat for us:*
  - Per-protein embedding of about 1M QfO proteins with ESM-2 or ProtT5 on 24 CPU cores takes hours (rough estimate, **[U]**; needs a timing test). Small models (ESM-2 8M/35M) are feasible.
  - Approximate nearest-neighbour search (FAISS/HNSW) is then O(N log N).
  - Embedding similarity measures *homology / structural similarity*. It does **not** distinguish orthologs from paralogs. Trees or reconciliation are still needed.

**k-mer sketching (MinHash)**
- Mash: Ondov et al. 2016, *Genome Biol* 17:132 (preprint 10.1101/029827) **[V]**.
- Dashing 2: Baker & Langmead 2023, *Genome Res* 33:1218, doi:10.1101/gr.277655.123 **[V]**.
- These are designed for genomes. For proteins, OMAmer's family-specific k-mer tables and Linclust's minimiser grouping are the practical equivalents. **I found no protein-MinHash orthology tool with QfO results.**

**Phylogenetic placement** (alternative to search for adding sequences to families): EPA-ng, Barbera et al. 2019, *Syst Biol* 68:365, doi:10.1093/sysbio/syy054 **[V]**. SHOOT (above).

---

## 3. Gene trees, reconciliation and species trees for orthology and HOGs

**Species-overlap (SO) rule**: Huerta-Cepas et al. 2007 (above).
- A node is a duplication if its two child clades share species. It needs no species tree and tolerates topological noise.
- Used by PhylomeDB, FastOMA (per level), eggNOG v7 (pbDOOD, with an SO threshold and outlier removal) and OF (hybrid).

**Duplication-loss reconciliation**
- Notung: Chen, Durand & Farach-Colton 2000, *J Comput Biol* 7:429, doi:10.1089/106652700750050871 **[V]**. Parsimony DL reconciliation with rearrangement of weak edges.
- ALE: Szöllősi et al. 2013, *Syst Biol* 62:901, doi:10.1093/sysbio/syt054 **[V]**. Probabilistic DTL over tree samples.
- GeneRax: Morel et al. 2020, *MBE* 37:2763, doi:10.1093/molbev/msaa141 **[V]**. Species-tree-aware ML gene trees under DTL.
- SpeciesRax: Morel et al. 2022, *MBE* 39:msab365 **[V]**. Species tree from gene-family trees under DTL.
- AleRax: Morel et al. 2024, *Bioinformatics* 40:btae162, doi:10.1093/bioinformatics/btae162 **[V]**. Joint gene and species tree with DTL.
- *Cost note:* GeneRax/ALE are far too slow for ~10⁵ families × thousands of species. Use them only on benchmark subsets or as a "gold-er" tree-building ablation.

**Species-tree inference from multi-copy gene trees**
- STAG: Emms & Kelly 2018, bioRxiv 10.1101/267914 **[P]**.
- STRIDE (rooting via duplications): Emms & Kelly 2017, *MBE* 34:3267, doi:10.1093/molbev/msx259 **[V]**.
- ASTRAL-Pro: Zhang, Scornavacca, Molloy & Mirarab 2020, *MBE* 37:3292, doi:10.1093/molbev/msaa139 **[V]**.
- ASTRAL-Pro 2: Zhang & Mirarab 2022, *Bioinformatics* 38:4949, doi:10.1093/bioinformatics/btac620 **[V]**. Fast and used by OF3 `--assign`.
- DISCO: Willson et al. 2022, *Syst Biol* 71:610, doi:10.1093/sysbio/syab070 **[V]**. Decomposes multi-copy trees into single-copy trees, then ASTRAL.

**Universal vs per-family species tree**
- OMA, FastOMA, Hieranoid and OF all reconcile against *one* species tree.
- ILS (short internal branches, e.g. great apes, Neoaves) and HGT (bacteria/archaea; 30 of the 78 QfO species are prokaryotes) cause genuine discordance.
- The species-overlap rule is robust to ILS-driven topology differences: it only asks whether species sets overlap, not what the topology is. Reconciliation-based labelling can invent spurious duplication-plus-loss pairs under ILS.
- Options for us:
  - (a) SO-based labelling, which needs no species tree for duplication calls but needs one to name HOG levels.
  - (b) A taxonomy-level HOG definition (HOG at clade X = genes descending from one gene in the MRCA of X). This needs only clade membership, not a fully resolved tree, and is robust to polytomies.
  - (c) Per-family tolerance, e.g. pbDOOD outlier exclusion, or treating SO < threshold as speciation.
- **The QfO STD tests use reference species trees.** A method that hard-codes the same species tree could bias STD only weakly, since STD measures gene-tree vs species-tree agreement of *sampled ortholog sets*. Still, document which tree we used.
- Tree builders at scale:
  - FastTree 2: Price et al. 2010, *PLoS ONE* 5:e9490, doi:10.1371/journal.pone.0009490 **[V]**.
  - FAMSA: Deorowicz et al. 2016, *Sci Rep* 6:33964, doi:10.1038/srep33964 **[V]**. FAMSA2 also exists **[U citation]**.
  - VeryFastTree exists **[U citation]**.

---

## 4. Structure in phylogeny and homology

**Foldseek / 3Di**: see §2.
- The 3Di alphabet and its 20×20 3Di substitution matrix are defined in van Kempen et al. 2024. Combined AA+3Di scoring improves remote-homology sensitivity well beyond sequence-only search.

**FoldTree**: Moi, Bernard, Steinegger, Nevers, Langleib & Dessimoz 2025, "Structural phylogenetics unravels the evolutionary diversification of communication systems in gram-positive bacteria and their viruses", *Nat Struct Mol Biol*, doi:10.1038/s41594-025-01649-8 **[V]**. Code: github.com/DessimozLab/fold_tree.
- Builds NJ trees from Foldseek-derived distances (the best was a local 3Di+AA-based "fident" distance). Trees were scored with a taxonomic congruence score (TCS) on CATH/Pfam families and on "highly divergent" families.
- Structure trees beat sequence ML trees for divergent families and are comparable for closer ones.
- *Takeaway:* useful for deep (root-HOG) splits. For shallow speciation-vs-duplication calls, sequence is likely sufficient.

**ProstT5**: Heinzinger et al. 2024, "Bilingual language model for protein sequence and structure", *NAR Genom Bioinform* 6:lqae150, doi:10.1093/nargab/lqae150 **[V]**.
- ProtT5-XL fine-tuned to translate AA↔3Di. Foldseek integrates a fast ProstT5-encoder + CNN predictor: `foldseek databases ProstT5 weights tmp; foldseek createdb in.fasta db --prostt5-model weights` (README checked).
- The README says it "runs by default on CPU" and is about 400–4000x faster than ColabFold. GPU gives 1–2 orders of magnitude more.
- The HF model card: half precision is GPU-only, so CPU uses fp32.
- **CPU throughput for about 1M proteins is unverified.** It must be timed on our node. The ProtT5-XL encoder is about 1.2B parameters, so expect tens to hundreds of ms per protein per core-group (**[U]**).
- Predicted 3Di lacks Cα, so there is no TM-score/LDDT and no `--alignment-type 1`.

**Q.3Di substitution model**: Garg & Hochberg 2025, "A general substitution matrix for structural phylogenetics", *MBE* 42:msaf124, doi:10.1093/molbev/msaf124 **[V]**.
- The first general 3Di substitution model for ML phylogenetics (IQ-TREE compatible). Used to root the tree of life between archaea and bacteria.
- *Actionable:* this is the model to use for 3Di partitions.

**Tertiary-interaction characters**: Puente-Lelievre, Malik, Douglas, Ascher & Matzke 2023, bioRxiv 10.1101/2023.12.12.571181 **[P; author list from memory: U]**.
- 3Di characters in model-based ML/Bayesian phylogenetics (IQ-TREE, a custom 3Di matrix). Ferritin-like superfamily; resolution beyond the twilight zone.
- Follow-up: Fullmer, Puente-Lelievre & Matzke 2025, bioRxiv 10.1101/2025.06.30.662300 **[P]**. **AA+3Di partitioned alignments resolve better than either alone, but the gain is weak for short or α-helical proteins.**
- Review: Puente-Lelievre et al. 2025 "Protein Structural Phylogenetics", *GBE* 17:evaf139, doi:10.1093/gbe/evaf139 **[V]**.
- A BEAST2 3Di package exists (beast2.org blog, Dec 2025) **[V page exists]**.

**FoldMason**: Gilchrist, Mirdita & Steinegger 2026, "Multiple protein structure alignment at scale with FoldMason", *Science* 391, doi:10.1126/science.ads6733 **[V]**.
- Progressive MSTA using the Foldseek 3Di+AA alphabet. Two orders of magnitude faster than other MSTA tools, handles more than 100k structures, and has an LDDT-based column confidence.
- Produces joint AA+3Di MSAs, the natural input for AA+3Di partitioned trees.
- Needs real structures (Cα) for refinement. **[U]** whether it accepts ProstT5-predicted 3Di only.

**Structure data access**: see §9.

---

## 5. Synteny in orthology

**Ensembl Compara QC scores** (Ensembl docs "Orthology quality-controls" **[V]**):
- **GOC** (gene order conservation): for an ortholog pair, the share of the 2 up- and 2 down-stream neighbours that are also orthologous. Each counts 25%, so the score is 0–100 and tolerates inversions and insertions.
- **WGA coverage**: exon-weighted coverage of the pair by whole-genome alignments.
- Both serve as *post-hoc* "high confidence" filters. They are not used in inference.

**Local-synteny orthology for mammals**: Jun, Mandoiu & Nelson 2009, "Identification of mammalian orthologs using local synteny", *BMC Genomics* 10:630, doi:10.1186/1471-2164-10-630 **[V by search result; authors U]**.
- 93% concordance with sequence orthology. Synteny resolves retrogenes and many-to-many cases.

**Proteinortho PoFF**: Lechner et al. 2014 (above). Synteny-aware graph refinement inside a scalable pipeline.

**OrthoRefine**: Ludwig et al. 2024, *BMC Bioinformatics* 25:163, doi:10.1186/s12859-024-05786-7 **[V]**. Post-processes OrthoFinder HOGs with a synteny window to split paralogs (bacteria).

**SOI (Synteny Orthology Index)**: Zhang R-G et al. 2025, *NAR* (PMC12006799) **[V journal/year; DOI not checked]**.
- An orthology index = fraction of OrthoFinder orthologs inside a syntenic block. It separates orthologous from paralogous synteny after WGDs (plants).

**GENESPACE**: Lovell et al. 2022, *eLife* 11:e78526 (preprint 10.1101/2022.03.09.483468) **[V preprint; eLife DOI from memory]**. Synteny-constrained OrthoFinder orthogroups with pan-genome annotation.

**MCScanX**: Wang et al. 2012, *NAR* 40:e49, doi:10.1093/nar/gkr1293 **[V]**. Collinear-block detection. Needs a BLAST/DIAMOND hit table.

**syntenet**: Almeida-Silva et al. 2023, *Bioinformatics* 39:btac806, doi:10.1093/bioinformatics/btac806 **[V]**. Synteny networks.

**pSONIC**: Conover et al. 2021 (PMC8496325) **[V exists; venue U]**. Ploidy-aware syntenic orthologs.

**Rearrangement-based orthology**: Rubert & Braga 2023, "Efficient gene orthology inference via large-scale rearrangements", *Algorithms Mol Biol* 18:14, doi:10.1186/s13015-023-00238-y **[V]**. Genome-wide family-free / DCJ-based. Small genome counts.

**POInTbrowse**: Siddiqui & Conant 2023, *BMC Bioinf*, doi:10.1186/s12859-023-05298-w **[V]**. Orthology in paleopolyploids.

**OMA browser synteny viewer**: OMA 2024 (above), visualisation only.

**Novabrowse**: 2026 preprint, doi:10.64898/2026.03.27.714371 **[P, not read]**.

**OrthoFinder v3 synteny**: none documented.

- *Takeaway:* synteny is used for QC, for refining existing OGs, or for WGD-rich plants. **No scalable method uses synteny as a first-class signal in HOG inference with QfO results.**
- Gene order decays quickly across QfO's deep taxonomic span. It is informative mainly among vertebrates, among close fungi (Saccharomycetaceae ohnologs) and in operons. Expect impact on the VGNC benchmark and on in-paralog / ohnolog splitting, not on global STD.

---

## 6. Expression, epigenomics, isoforms and gene-loss validation

**Ortholog conjecture and expression**
- Altenhoff, Studer, Robinson-Rechavi & Dessimoz 2012, *PLoS Comput Biol* 8:e1002514, doi:10.1371/journal.pcbi.1002514 **[V]**. Orthologs are weakly but significantly more functionally similar than paralogs.
- Chen & Zhang 2012, *PLoS Comput Biol* 8:e1002784, doi:10.1371/journal.pcbi.1002784 **[V]**. Supported by RNA-seq.
- Kryuchkova-Mostacci & Robinson-Rechavi 2016, *PLoS Comput Biol* 12:e1005274, doi:10.1371/journal.pcbi.1005274 **[V]**. Tissue-specificity diverges slowly between orthologs and quickly between paralogs.
- Dunn, Zapata, Munro, Siebert & Hejnol 2018, *PNAS* 115:E409, doi:10.1073/pnas.1707515115 **[V]**. With a phylogenetic framework, support for the conjecture vanishes on the same data, a statistical warning.
- Song et al. 2024/2025, "Revising the ortholog conjecture in cross-species comparison of scRNA-seq data", bioRxiv 10.1101/2024.06.21.600109 **[P]**.
- *Takeaway:* expression is a weak, noisy tie-breaker for *which paralog copy kept the ancestral role*. It is not an orthology/paralogy discriminator by itself.

**Bgee**: Bastian et al. 2025, "Bgee in 2024: focus on curated single-cell RNA-seq datasets, and query tools", *NAR* 53:D878, doi:10.1093/nar/gkae1118 **[V; page numbers U]**. Bgee suite: Bastian et al. 2021 *NAR* 49:D831, doi:10.1093/nar/gkaa793 **[V]**.
- About 50 animal species; downloads at bgee.org.
- **Circularity risk:** Bgee's comparative expression uses OMA HOGs. Use only per-species expression calls, not Bgee's ortholog-based comparisons.
- Overlap with QfO 2022 is mostly vertebrates plus fly and worm (**[U exact list]**).

**Isoforms**
- QfO 2022 ships `*_additional.fasta` with 394,952 isoforms (42 species). Benchmarks score canonical accessions only, so isoform choice affects which protein represents a gene.
- IsoSel: Philippon et al. 2017, *PLoS ONE* 12:e0174250, doi:10.1371/journal.pone.0174250 **[V]**. Picks the isoform that best fits the family alignment.
- Guillaudeux et al. 2022, *BMC Genomics* 23:236, doi:10.1186/s12864-022-08429-4 **[V]**. Orthologous isoforms between human, mouse and dog; tools from the Ouangraoua lab (github.com/UdeS-CoBIUS/TranscriptOrthology).
- Splicing orthology: Blanquart et al. 2016, *BMC Genomics* (PMC5123294) **[authors U]**.
- Jia et al. 2010, *BMC Genomics* 11(S4):S11, doi:10.1186/1471-2164-11-S4-S11 **[V]**.
- Zambelli et al. 2010, *BMC Genomics* 11:534, doi:10.1186/1471-2164-11-534 **[V]**.
- FastOMA can pick the best-placing isoform with OMAmer (feature described in its README/paper; **[U exact flag]**).
- *Takeaway:* "family-consistent isoform selection" (choose the isoform with the best profile/placement score rather than the longest) is a cheap, measurable layer. An isoform-orthology output (orthologous isoform pairs) is novel but has no standard benchmark.

**Gene-loss validation using genome DNA**
- TOGA: Kirilenko et al. 2023, "Integrating gene annotation with orthology inference at scale", *Science* 380:eabn3107, doi:10.1126/science.abn3107 **[V]**.
  - Projects reference transcripts through whole-genome alignment chains. Machine learning on intronic and intergenic alignment features separates orthologs from paralogs and processed pseudogenes.
  - CESAR re-aligns exons and classifies inactivating mutations (intact / lost / missing).
  - Scale: 488 mammals, 501 birds.
- TOGA2: Malovichko … Hiller 2026, bioRxiv 10.64898/2026.06.30.735536 **[P]**. Adds gene-tree reconciliation, 513x less memory, 2,162 vertebrate assemblies. github.com/hillerlab/TOGA2.
- CESAR: Sharma, Elghafari & Hiller 2016, *NAR* 44:e103, doi:10.1093/nar/gkw210 **[V]**.
- Gene-loss pipeline: Sharma et al. 2018, *Nat Commun* 9:1215, doi:10.1038/s41467-018-03667-1 **[V]**.
- *Lightweight alternative for us:* for a HOG absent from species S, take the syntenic interval between S's orthologs of the HOG's flanking genes. Search the S genome interval (not proteome) with miniprot (Li 2023, *Bioinformatics* 39:btad014 **[V from memory: U]**) or tblastn using the HOG consensus or representative. Classify the result as unannotated-intact / pseudogenised (frameshift or stop) / absent.
- Genomes come from the assemblies listed in `PROTEOME_ASSEMBLY_MAPPING` (QfO 2026 FTP). Note that **QfO `*_DNA.fasta` files are CDS, not genomes.**

---

## 7. Benchmarks: what, where, how to run locally

### 7.1 QfO benchmark service

**Papers**
- Altenhoff et al. 2016, "Standardized benchmarking in the quest for orthologs", *Nat Methods* 13:425, doi:10.1038/nmeth.3830 **[V]**.
- Altenhoff et al. 2020, *NAR* 48:W538, doi:10.1093/nar/gkaa308 **[V]**. Generalised STD, consensus calls.
- Nevers et al. 2022, "The Quest for Orthologs orthology benchmark service in 2022", *NAR* 50:W623 (PMC9252809) **[V title/year; volume U]**.
- Altenhoff et al. 2024, "New developments for the Quest for Orthologs benchmark service", *NAR Genom Bioinform* 6:lqae167, doi:10.1093/nargab/lqae167 **[V]**. 2022 reference set: 78 species (48 Euk, 23 Bac, 7 Arch), 988,778 canonical + 394,952 isoforms. 18 methods. New FAS benchmark.
- Boeckmann et al. 2011, *Brief Bioinform* 12:423 **[V via search]**, for SwissTree.

**Tests in the 2022 event** (from `nextflow.config`, `challenges_ids`, and OpenEBench)

| Category | Tests | Metrics |
|---|---|---|
| Function | GO, EC | avg Schlicker similarity vs NR_ORTHOLOGS |
| Function | FAS | feature architecture similarity |
| Reference gene trees | SwissTrees, TreeFam-A, VGNC | TPR/PPV |
| Species tree discordance, classic | STD_Eukaryota, STD_Fungi, STD_Bacteria | RF_DISTANCE, FRAC_INCORRECT_TREES vs NR_COMPLETED_TREE_SAMPLINGS |
| Species tree discordance, generalised | G_STD2_Luca, G_STD2_Eukaryota, G_STD2_Vertebrata, G_STD2_Fungi | as classic STD |

- **Every test is pairwise.** HOGs are converted to pairs. In `map_relations.py`, every pair of genes from different species that sit in different children of an `orthologGroup` is an ortholog. `paralogGroup` children are not paired. Top-level `paralogGroup`s are handled through their nested orthologGroups.
- Genes **must** carry the UniProt accession in `protId`, and the namespace must be `http://orthoXML.org/2011/`. TSV pair format is also accepted.

### 7.2 Running QfO locally (no Docker)

Repo: https://github.com/qfo/benchmark-webservice (last commit 2025-03-26; cloned and inspected).

**Reference data**: `./fetch_reference_data.py 2022`. 2022 maps to 2022.1, Zenodo **https://zenodo.org/records/10517603** (CC-BY-4.0). Files include:
- ServerIndexed.db (823 MB) and ServerIndexed.db.tree (1.78 GB)
- mapping.json.gz, GOdata, enzymes, ReconciledTrees_SwissTrees/TreeFam-A
- species_tree_samples_*.nwk, vgnc-orthologs, swissprot.txt.gz
- fas_precomputed.json.gz (557 MB), fas_annotations.tgz

Other releases: 2022.0 = zenodo 7919753; 2020.2 = 15087752; 2018 = 15087679; 2011 = 15083829.

**Containers**: the profiles use `qfobenchmark/python:2022.1`, `qfobenchmark/darwin:2022.1` and `qfobenchmark/fas_benchmark:2022.1`. There are `docker` and `singularity` profiles, and a `conf/curnagl.config` (UNIL SLURM, `module load singularity`).

**No-Docker options, in recommended order:**
1. **Native install.**
   - Darwin binary: `https://biorecipes.com/darwin/darwin.linux64` (HTTP 200, 5.2 MB) plus `darwin-lib.tgz`.
   - The Dockerfile then patches `lib/darwinit` so that `datadirname` points to the refset dir, sets `ulimit -s unlimited`, and sets `DARWIN_ORTHOLOG_BENCHMARK_REPO_PATH`.
   - Install `mafft` (conda), the Python requirements (lxml, dendropy, biopython, jsonschema, matplotlib, pandas), and Nextflow from conda (needs Java).
   - Then run `nextflow run main.nf` without a container profile, or call `*.sh` / `*.drw` directly.
   - Expect some path patching. **[U]** whether everything runs cleanly outside the image. Budget about 2–4 h of agent time.
2. `pip install udocker`, then pull the three Docker Hub images and run rootless. **[U]** cluster compatibility.
3. `conda install -c conda-forge apptainer` with the `singularity` profile. Needs unprivileged user namespaces, **[U]**. Our container's `curnagl.config` suggests the host has singularity, but nesting from inside our container is unlikely.
4. Public web submission at https://orthology.benchmarkservice.org (OpenEBench login; results can be kept private) **[U on turnaround]**. Not suited to iteration.

**Reference proteomes (must use these exact ones)**: https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2022_02_with_updated_UP000000437/QfO_release_2022_02_with_updated_UP000000437.tar.gz (2.4 GB).
- Contains, per species: `*.fasta` (canonical), `*_additional.fasta` (isoforms), `*_DNA.fasta` (CDS), `*.gene2acc`, `*.idmapping`, and `*.xml` (gene coordinates).
- The OpenEBench metadata still points to `.../QfO/QfO_release_2022_02_with_updated_UP000000437.tar.gz`, **which now returns 404**.
- The current folder https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/ holds **2026_02** (81 species, UniProt 2026_02 / Ensembl 115). It has per-species `*.gff3.gz` (15 eukaryotes), `*.gbff.gz`, `*_gcoord.xml.gz` and `PROTEOME_ASSEMBLY_MAPPING`. Use it for synteny coordinates and genome assembly IDs, mapped by accession. **No 2026 benchmark reference data exists** (fetch script tops out at 2022.1).

### 7.3 Published competitor predictions (download, do not rerun)

B2SHARE API query: `https://b2share.eudat.eu/api/records?q=%22Quest%20for%20Orthologs%22&size=100` (42 records). File URL pattern, which I tested and which returned HTTP 200: `https://b2share.eudat.eu/records/<id>/files/<filename>?download=1`.

2022 reference-proteome submissions (record id → file):

| Method | Record id | File | Format |
|---|---|---|---|
| OrthoFinder (2023) | `r5jz7-fsr77` (dup `d8nrr-9fs80`) | OrthoFinder.txt | 214 MB |
| FastOMA v0.1.4 | `r2hs8-wcx45` | FastOMA_v0.1.4_qfo2022_17Nov3.orthoxml | OrthoXML, HOGs |
| OMA HOGs | `gmf9z-vn551` | HierarchicalGroups.orthoxml.gz | OrthoXML |
| OMA Groups | `t7a3g-ap796` | MergedOMAGroups.orthoxml.gz | OrthoXML |
| OMA Pairs | `w3r7b-da890` | OMA.2.5.0-VPairs.txt.gz | pairs |
| SonicParanoid2-sens | `e5k5r-vsk83` | qfo22-challenge.sp2-sens.tsv.gz | pairs |
| SonicParanoid2-sens-g | `bvfjw-k3d26` | qfo22-challenge.sp2-sens-g.tsv.gz | pairs |
| Ensembl Compara | `mnbxy-pnx93` | qfo_2022_02.default.orthoxml.xml.gz | OrthoXML |
| PANTHER 18 all | `c9wp1-7dh03` | PANTHER18_all_1 | pairs |
| PANTHER 18 LDO | `rv6b7-bpk11` | PANTHER18_ldo | pairs |
| Hieranoid 2 | `rx86s-23p86` | pairs_hieranoid-diamond.tsv | pairs |
| InParanoid5 | `ca8wf-v9w65` | InParanoid_QFO22.pairs_1.gz | pairs |
| Domainoid+ | `hrmsy-tje90` | mergedResults_domainoid.tsv | pairs |
| OrthoInspector 3.5 | `bvb0m-dkz52` | outputbench2022.tsv | pairs |
| PhylomeDB v5 | `yvdbk-8te46` | orthologs.farthest.30Cl.2022.txt.gz | pairs |
| MetaPhOrs 2.5 | `j0yhx-8g173` | metaphors…2022.tsv | pairs |
| OrthoFFGC | `4hky9-h2996` | ffgc.l2.I1.4.xml | OrthoXML |
| RSD | `kzf58-3z834` | rsd.txt.gz | pairs |
| BBH | `dfwf8-tva07` | bbh.tsv.gz | pairs |
| pbDOOD (eggNOG v7) | `0798b-qh662` | pbdood_pairs.tsv | pairs (2026-06) |
| pbDOOD (midpoint) | `eqm9x-5sr17` | pbdood_pairs_midpoint.tsv | pairs (2026-06) |

The 2020-proteome submissions are also on B2SHARE (OrthoFinder MSA v2.5.2, Sonicparanoid 1.3.5, OMA, PANTHER 16 and others).

- Downloaded predictions give us competitor *pairs*. These let us compute our own extra metrics: overlap with consensus, synteny support of competitor pairs, and Jaccard with our HOGs.
- **No runtime is recorded.** Runtime comparisons need local runs. FastOMA, SonicParanoid2 and OrthoFinder 3 are pip/conda installable. OF3 needed 128 h for 1024 species on 16 cores, so compare on subsets of 16–256 species.

### 7.4 QfO 2022 public results (OpenEBench, harvested 2026-10-26)

Generated with `scripts/fetch_qfo_public_results.py`.

Column key:
- ST = SwissTree; TF = TreeFam-A.
- GO/EC = avg Schlicker. "GO n" = number of orthologs in the GO test.
- G2Euk = G_STD2_Eukaryota RF distance, with n = completed tree samplings. STDEuk = STD_Eukaryota RF, with n samplings.
- Higher TPR/PPV/Schlicker/n is better. Lower RF is better.

| Method | ST TPR | ST PPV | TF TPR | TF PPV | VGNC TPR | VGNC PPV | GO n | GO | EC | G2Euk RF | G2Euk n | STDEuk RF | STDEuk n | FAS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BBH (SW) | .570 | .932 | .483 | .956 | .979 | 1.00 | 88910 | .511 | .951 | .232 | 3646 | .046 | 9912 | .774 |
| Domainoid+ | .759 | .932 | .672 | .923 | .986 | .995 | 169998 | .485 | .916 | .237 | 5162 | .045 | 13464 | .740 |
| Ensembl Compara | .648 | .892 | .637 | .818 | .927 | .987 | 357983 | .444 | .840 | .303 | 2436 | .077 | 13583 | .651 |
| FastOMA 0.1.4 | .699 | .953 | .546 | .919 | .969 | .996 | 121949 | .484 | .928 | .238 | 3092 | .047 | 8686 | .669 |
| Hieranoid 2 | .726 | .950 | .621 | .930 | .953 | .998 | 126270 | .494 | .948 | .222 | 3643 | .034 | 10296 | .773 |
| InParanoid5 | .701 | .944 | .627 | .933 | .958 | .998 | 143167 | .491 | .927 | .227 | 3655 | .036 | 11934 | .765 |
| MetaPhOrs 2.5 | .688 | .921 | .630 | .952 | .828 | .987 | 181562 | .475 | .892 | .263 | 2088 | .074 | 14689 | .723 |
| OMA Groups | .496 | .940 | .383 | .954 | .945 | .999 | 62952 | .517 | .971 | .205 | 1063 | .035 | 7745 | .821 |
| OMA HOGs | .696 | .924 | .580 | .923 | .971 | .979 | 362533 | .417 | .829 | .237 | 2981 | .041 | 7395 | .411 |
| OMA Pairs | .677 | .929 | .565 | .920 | .972 | .980 | 110960 | .496 | .952 | .226 | 2369 | .040 | 10339 | .810 |
| OrthoFFGC | .503 | .921 | .353 | .948 | .985 | .976 | 55571 | .521 | .975 | .199 | 404 | .061 | 5701 | – |
| OrthoFinder (2023) | .825 | .917 | .737 | .918 | .980 | .997 | 167696 | .487 | .931 | .249 | 6629 | .051 | 13707 | .716 |
| OrthoInspector 3.5 | .748 | .931 | .670 | .941 | .979 | .998 | 157365 | .489 | .922 | .235 | 5287 | .039 | 12406 | .743 |
| PANTHER 18 LDO | .553 | .944 | .430 | .957 | .940 | .998 | 77921 | .515 | .961 | .239 | 3679 | .040 | 7515 | .766 |
| PANTHER 18 all | .811 | .943 | .698 | .882 | .954 | .994 | 204051 | .474 | .923 | .272 | 5943 | .055 | 14584 | .712 |
| RSD (SW) | .546 | .919 | .443 | .953 | .969 | .999 | 83538 | .510 | .947 | .242 | 2718 | .061 | 10936 | .754 |
| SonicParanoid2-sens | .744 | .864 | .682 | .913 | .982 | .996 | 208515 | .474 | .864 | .283 | 6181 | .057 | 13350 | .708 |
| SonicParanoid2-sens-g | .717 | .938 | .667 | .936 | .980 | .998 | 177643 | .483 | .906 | .232 | 4878 | .037 | 12831 | .693 |
| pbDOOD (eggNOG v7) | .847 | .952 | .735 | .880 | .973 | .993 | 246469 | .463 | .880 | .353 | 8929 | .072 | 14686 | – |
| pbDOOD (midpoint) | .854 | .950 | .735 | .889 | .973 | .995 | 222212 | .463 | .894 | .350 | 9190 | .072 | 14654 | – |
| PhylomeDB v5 | .731 | .925 | .595 | .931 | .974 | .997 | 119460 | .499 | .943 | .248 | 4551 | .070 | 14481 | .731 |

Reading notes:
- Every test is a precision–recall trade-off. Claims must be Pareto-front statements per test, not a single score.
- OF3's paper reports STD-Euk n = 15,721 at RF 0.06 (scale differs slightly from the 0.051 above; **[U]** which variant).
- The "OrthoFinder" row is the 2023 submission, probably OF 2.5.x.
- The same harvest script also returns G_STD2_Luca/Fungi/Vertebrata, STD_Fungi, STD_Bacteria and FRAC_INCORRECT_TREES.

### 7.5 OrthoBench (Open_Orthobench)

- Papers: Emms & Kelly 2020, "Benchmarking orthogroup inference accuracy: revisiting Orthobench", *GBE* 12:2258, doi:10.1093/gbe/evaa211 **[V]**. Original: Trachana et al. 2011, *BioEssays* **[U DOI]**.
- Repo: https://github.com/davidemms/Open_Orthobench (cloned). The release tarball is at https://github.com/davidemms/Open_Orthobench/releases/download/v1.1/BENCHMARKS.tar.gz.
- Contents:
  - 70 RefOGs.
  - 12 bilaterian Ensembl proteomes in `BENCHMARKS/Input/` (human GRCh38, chimp Pan_tro_3.0, mouse GRCm38, rat Rnor_6.0, dog CanFam3.1, opossum ASM229v1, chicken GRCg6a, zebrafish GRCz11, Tetraodon TETRAODON8, Ciona KH, fly BDGP6.28, worm WBcel235).
  - 251,378 genes in total. Headers are bare IDs (e.g., `ENSP…`, `FBpp…`, `WBGene…`).
- Run: `python3 BENCHMARKS/benchmark.py our_orthogroups.txt` (one OG per line; `-b` for the basic format).
- The script rejects files that have genes in more than one OG or that cover fewer than 50% of genes.
- It reports precision, recall, F-score, fusions, fissions and missing genes/RefOGs (OF3 used 7 measures, including entropy).
- **For HOG methods, export the root-level (Bilateria) HOGs as OGs.** No precomputed competitor results are in the repo. OF3's Supplementary File 1/2 (Nature Methods) contain them **[U accessibility]**.
- GFF for synteny would need the matching Ensembl archive (roughly releases 98–102, given CanFam3.1/Rnor_6.0/GRCg6a; **[U exact release]**).

### 7.6 SwissTree

- Site: https://swisstree.vital-it.ch/cgi-bin/swisst (release 2017.0). FTP: https://swisstree.vital-it.ch/ftp/ (per-family `ST0xx/consensus_tree*.phyloxml|nhx`, `swisstree.nhx`, `swisstree.fst`).
- Citation: Boeckmann et al. 2017, bioRxiv 10.1101/181966 **[P]**.
- **SwissTree is already inside QfO** (`ReconciledTrees_SwissTrees.drw` mapped to 2022 proteomes). The QfO SwissTrees challenge is the practical way to run it. Running the standalone site would need ID mapping to the 2011 reference proteomes.

### 7.7 Other benchmarks

- **LEMOrtho** (https://lemortho.ezlab.org): a container-based continuous benchmark from the OrthoDB group, using the 70 OrthoBench RefOGs with F1 and Variation of Information **[site not opened; U on status]**.
- **No benchmark scores HOG hierarchy correctness** (Sarton-Lohéac et al. 2025). We could add an internal check: duplication-node placement vs reconciled reference trees (SwissTree/TreeFam-A trees are reconciled), as a hierarchy-specific sub-benchmark.

---

## 8. OrthoXML format and tools

**Specification**: Schmitt et al. 2011, "SeqXML and OrthoXML: standards for sequence and orthology information", *Brief Bioinform* 12:485, doi:10.1093/bib/bbr034 **[U: from memory]**. Repo: https://github.com/qfo/orthoxml (cloned).

| Version | Status / changes |
|---|---|
| **v0.3** | Original. Only `orthologGroup` allowed at the top level. |
| **v0.4** | **The current standard.** A top-level `paralogGroup` is allowed (root duplication). Same namespace `http://orthoXML.org/2011/`. Schema at https://orthoxml.org/0.4/orthoxml.xsd (HTTP 200). |
| **v0.5rc1** | Branch 0.5. Adds an optional `<taxonomy>` element (nested `<taxon id name>`), a `taxonId` attribute on `orthologGroup`/`paralogGroup` and on `species`, and makes NCBITaxId optional. https://orthoxml.org/0.5/orthoxml.xsd returns 404; use `orthoxml-0.5.xsd` from the GitHub repo. |

Structure: `<orthoXML>` contains:
- `<species name NCBITaxId>` → `<database name version>` → `<genes>` → `<gene id protId geneId transcriptId>`;
- optional `<scores>` definitions;
- `<groups>` with nested `<orthologGroup id>` / `<paralogGroup>` holding `<geneRef id>`, `<score id value>`, and `<property name value>`.

By OMA/FastOMA convention, the HOG level goes in `<property name="TaxRange" value="Mammalia"/>` (0.3/0.4) or in `taxonId` (0.5).

**For QfO:**
- Put UniProt accessions in `protId`.
- Make sure a top-level paralogGroup still nests orthologGroups.
- Deep over-merged orthologGroups create many false pairs.

**Tools**
- **orthoxml-tools**: Yazdizadeh Kharrazi et al. 2025, *J Mol Evol* 93:800, doi:10.1007/s00239-025-10271-7 **[V via README]**. `pip install orthoxml-tools`, GPL-3. Validates against the declared schema version (0.3/0.4/0.5), gives stats, filters by completeness score, splits, and converts to pairs/OGs/NHX/OrthoFinder CSV.
- **pyHam / iHam**: Train et al. 2019, *Bioinformatics* 35:2504, doi:10.1093/bioinformatics/bty994 **[V]**. HOG traversal, gene gain/loss/duplication per branch. Needs a species tree whose node names match the HOG levels.
- ETE toolkit `etree2orthoxml` (tree to OrthoXML).
- Validation with `xmllint --schema orthoxml-0.4.xsd`.

---

## 9. Data for the extra layers

**3D structure**
- **AlphaFold DB**: current version **v6** (FTP `latest/` 2026-01-12).
  - Per-proteome tars for 46 model and global-health proteomes at https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/ (e.g., `UP000005640_9606_HUMAN_v6.tar`, `UP000000625_83333_ECOLI_v6.tar`, `UP000002311_559292_YEAST_v6.tar`, `UP000000437_7955_DANRE_v6.tar`, ARATH, DROME, CAEEL, MOUSE, RAT, SCHPO, DICDI, …).
  - `accession_ids.csv` lists all covered UniProt accessions, so QfO coverage can be checked offline.
  - `sequences.fasta` and `download_metadata.json` are also on the FTP.
  - Per-entry API: `https://alphafold.ebi.ac.uk/api/prediction/<ACC>` returns `pdbUrl`/`cifUrl`/`bcifUrl` (`AF-<ACC>-F1-model_v6.*`).
  - Bulk for other species: Google Cloud public dataset **[U current bucket name for v6]**.
  - Caveats: canonical isoforms only, and long human proteins are fragmented.
- **Foldseek pre-built DBs**: `foldseek databases Alphafold/Proteome|Alphafold/UniProt50|Alphafold/UniProt|Alphafold/Swiss-Prot|ESMAtlas30|PDB|ProstT5`.
  - Alphafold/Proteome covers the AFDB proteome set. We can extract 3Di for our accessions from Alphafold/UniProt, but it is large (hundreds of GB, **[U size]**).
  - Latest tagged release: 10-941cd33 (2025-01).
- **ProstT5**: HF `Rostlab/ProstT5`, or via Foldseek `--prostt5-model` (CPU default, fp32). **Must be timed** on a 1k-protein sample before committing to 1M proteins.
- **Practical plan**:
  - Use AFDB-derived 3Di where the accession has a model (most QfO canonical proteins; **[U %]**, check against accession_ids.csv).
  - Use ProstT5 3Di for the rest and for OrthoBench Ensembl proteins.
  - Compare real vs predicted 3Di on an overlap as a control.

**Synteny / gene coordinates**
- QfO 2022 tarball `*.xml` (gene coordinates).
- QfO 2026 `*_gcoord.xml.gz` (all species) and `*.gff3.gz` (15 eukaryotes), plus `PROTEOME_ASSEMBLY_MAPPING` (GCA/GCF per proteome).
- Ensembl FTP (https://ftp.ensembl.org/pub/; release 115 matches QfO 2026) and NCBI Datasets for genomes and GFF.
- Prokaryotes: use gbff.

**Expression**
- Bgee downloads (https://www.bgee.org/download), per-species present/absent calls and expression ranks.
- About 50 animals; only vertebrates, fly and worm overlap QfO and OrthoBench.

**Isoforms**: QfO `*_additional.fasta` (2022: 394,952 isoforms) and Ensembl `pep.all`.

**Genomes for loss validation**: assemblies via `PROTEOME_ASSEMBLY_MAPPING` → NCBI/Ensembl. Tools: miniprot / tblastn; TOGA2 for vertebrates only (needs chain alignments, too heavy for 48 h).

---

## 10. Implications for our design

### 10.1 Where a conceptually different, scalable method could compete

1. **Reference-free sub-quadratic inference.**
   - OF3 scales via a fixed core plus DIAMOND profile assignment. FastOMA scales via an *external* curated OMA HOG database. PANTHER and eggNOG are curated databases.
   - A method that builds families de novo in near-linear time (Linclust or DeepClust cascades, 3Di k-mers, or embedding ANN graphs) and then infers HOGs **without any pre-built orthology DB** is conceptually distinct. It also avoids the benchmark-leakage question that applies to OMAmer/PANTHER/eggNOG resources built from the QfO species.
2. **Progressive ancestral reconstruction (Hieranoid-style, N−1 merges).**
   - Hieranoid 2 is already competitive on precision and STD (best STD-Euk RF) with linear species-pair complexity, but it is old, slow per merge, and recall-limited.
   - Modern ingredients could address recall without quadratic cost: MMseqs2/DIAMOND profile-vs-profile search at each internal node, an AA+3Di consensus, synteny-anchored pairing of ancestral loci, and SO-based duplication calls.
   - It also yields HOGs natively, one level per merge.
3. **Structure as a deep-split signal.**
   - Root-HOG delineation (fusion/fission errors, where FastOMA fuses 3.5x more RefOGs) and deep duplications are where sequence trees fail and FoldTree/Q.3Di help.
   - Use AA+3Di partitioned trees, or 3Di-aware distances, only for deep or ambiguous nodes. This keeps cost bounded.
   - Expect little effect for shallow vertebrate splits, and weak effect for α-helical proteins (Fullmer et al. 2025).
4. **Synteny inside inference, not as post-hoc QC.**
   - Use neighbourhood conservation as an extra likelihood term when (a) splitting fused families, (b) choosing between in- and out-paralog placements in ambiguous trees, and (c) validating losses.
   - No published scalable HOG method does this. Measurable on VGNC (vertebrates), SwissTree and OrthoBench (12 bilaterians; needs matching Ensembl GFF).
5. **Isoform-aware representative selection.** Cheap, and changes which canonical protein is compared. Benchmark with and without it.
6. **Explicit loss validation.** An "absent vs lost vs unannotated" label per HOG×species, from genome loci, is a novel output that improves HOG taxonomic-level placement. No standard benchmark exists, so validate on known lost genes (e.g., TOGA calls for mammals, **[U availability]**).

### 10.2 Concrete risks

- **Benchmark mechanics**
  - Predictions must be on the *2022_02_with_updated_UP000000437* proteomes with UniProt `protId`. The 2026 FTP proteomes are a different set and must not be used for QfO.
  - The Darwin-based pipeline outside Docker is untested. Allocate time early and verify against a known public submission (e.g., rerun the OMA Groups orthoxml from the repo `example/` or B2SHARE and reproduce its published numbers within tolerance). That is our positive control.
- **HOG→pair flattening penalises over-merging** (OMA HOGs vs OMA Pairs). Export both the full HOGs and a conservative pair set. Report what each is.
- **There is no single "accuracy".** Results are Pareto fronts per test. Pre-register which tests and which side of the trade-off we target, so that metrics are not cherry-picked.
- **Runtime comparison.** Public predictions have no runtimes. Run FastOMA, SonicParanoid2 and OF3 (and ideally Hieranoid 2) on identical 8/16/32/64/78-species subsets on our 24-core node, and extrapolate scaling. OF3's 1024-species run needed 128 h, so full-scale competitor runs are infeasible within 48 h.
- **Compute for extra layers.**
  - ProstT5 on CPU for 1M proteins may take a long time (unverified).
  - AFDB per-accession downloads of about 1M files are impractical. Prefer proteome tars plus the Foldseek DBs, and fall back to ProstT5 for the uncovered remainder.
  - TOGA, GeneRax and ALE are infeasible at full scale. Restrict them to benchmark subsets.
- **Species-tree dependence.**
  - QfO spans bacteria, archaea and eukaryotes, with HGT, so a universal tree is wrong for many prokaryotic families.
  - Prefer SO-style duplication calls and taxonomy-clade HOG levels.
  - Record whether the NCBI taxonomy or an inferred tree is used. Test both, because FastOMA reports sensitivity to tree resolution.
- **Expression layer.** Weak signal and contested statistics (Dunn et al. 2018). Bgee is circular with OMA. Few species overlap. Treat it as an exploratory ablation with a pre-stated negative-result criterion.
- **Embedding layers detect homology, not orthology.** They can only replace the search or clustering step. Orthology still requires duplication/speciation inference. pLM training data include the benchmark proteins, which is a mild leakage for homology but not for orthology labels.
- **Novelty claim.** OF3 (core + profile assignment), FastOMA (k-mer placement + bottom-up SO) and eggNOG v7 (domain trees + SO) cover most "fast clustering then trees" designs. To be "conceptually different", the core inference principle must differ (for example progressive ancestral profiles, structure+synteny-informed scoring, or embedding-graph speciation inference), not just swap DIAMOND for MMseqs2.

---

## Key URLs (all checked 2026-10-26)

**QfO**
- Benchmark repo: https://github.com/qfo/benchmark-webservice
- Reference data 2022.1: https://zenodo.org/records/10517603
- 2022 proteomes: https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2022_02_with_updated_UP000000437/QfO_release_2022_02_with_updated_UP000000437.tar.gz
- Current (2026_02) QfO FTP with gff/gcoord: https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/
- Competitor predictions: https://b2share.eudat.eu/api/records?q=%22Quest%20for%20Orthologs%22&size=100
- Metrics API: https://openebench.bsc.es/api/scientific/staged/Dataset/<id>
- Darwin binary: https://biorecipes.com/darwin/darwin.linux64

**OrthoBench**
- Repo: https://github.com/davidemms/Open_Orthobench
- Release tarball: https://github.com/davidemms/Open_Orthobench/releases/download/v1.1/BENCHMARKS.tar.gz

**SwissTree**: https://swisstree.vital-it.ch/ftp/

**OrthoXML**
- Schemas: https://github.com/qfo/orthoxml (0.3/0.4/0.5rc1 XSDs)
- orthoxml-tools: https://github.com/DessimozLab/orthoxml-tools

**AlphaFold DB**: https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/ and API https://alphafold.ebi.ac.uk/api/prediction/<ACC>

**Tools**
- Foldseek: https://github.com/steineggerlab/foldseek
- ProstT5: https://huggingface.co/Rostlab/ProstT5
- OrthoFinder: https://github.com/OrthoFinder/OrthoFinder
- FastOMA: https://github.com/DessimozLab/FastOMA
- SonicParanoid2: https://gitlab.com/salvo981/sonicparanoid2
- pbDOOD: https://github.com/AnaHrnndz/pbdood
- FoldTree: https://github.com/DessimozLab/fold_tree
- TOGA2: https://github.com/hillerlab/TOGA2
