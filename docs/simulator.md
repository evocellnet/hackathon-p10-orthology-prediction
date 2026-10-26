# Toy genome-evolution simulator (`ancestra.simulate`)

Purpose: verification data with **exact** orthology truth (DESIGN.md §5). It is used for
positive and negative controls of every Ancestra sub-step, not as a realistic benchmark.
Its results do not replace OrthoBench or QfO.

```
python -m ancestra.simulate --out DIR --seed 11 [--species 8 --genes 300 ...]
python -m ancestra.evaluate_toy --truth DIR/truth/hogs.orthoxml --pred pred.orthoxml [--json r.json]
```
Python API: `simulate(SimParams(...), seed) -> SimResult` and `write_simulation(res, DIR)`.
For the same seed and parameters, the outputs are byte-identical (this is tested).

## Model

| Component | What is simulated | Default |
|---|---|---|
| Species tree | Yule (pure-birth) ultrametric tree with `S` leaves, rescaled so the root-to-tip depth is `tree_height` expected substitutions/site. Or a user Newick (`--tree`, branch lengths required; rescaled only if `--tree-height` is given). Names: internal `N0` (root), `N1`, … in pre-order; leaves `sp0`, `sp1`, …. User names are kept. | S=8, height 1.0 |
| Root genome | `G` genes (one gene = one family) in order on `n_chromosomes`, cut at random points. Each root protein has length U[`len_min`,`len_max`] and residues drawn from the WAG stationary frequencies. Output sequences start with `M` (prepended, never mutated). | G=300, 2 chromosomes, 150–600 aa |
| Gene-content events | A continuous-time process on **all branches at once** (a global clock, so HGT donors and recipients are contemporaneous). Rates are per gene per unit branch length. **Duplication**: tandem (copy placed next to the original), or *dispersed* with probability `dispersed_dup_frac` (random position and strand). **Loss**: a uniformly chosen gene. **Inversion**: a block of Geometric(mean `block_mean`) genes is reversed and its strands flipped. **Translocation**: a block moves to a random position on any chromosome, inverted with p=0.5. **HGT**: a copy is inserted at a random position in a random *contemporaneous* other lineage. | dup 0.2, loss 0.2, inv 0.05, transl 0.02, HGT 0 (off), dispersed 0.2, block 4 |
| WGD (optional) | On the branch above `--wgd-branch` (a node name, or `random`), at relative position `wgd_position`, every gene is duplicated. The copies form new chromosomes `chrNb` in the same order. Fractionation: with probability `wgd_fractionation`, one random copy of each pair is lost at a uniform time later on the same branch. | off; position 0.5, fractionation 0.5 |
| Substitutions | **WAG** exchangeabilities and frequencies (Whelan & Goldman 2001; values copied from the PAML-format `wag.dat` shipped with HMMER/easel). Q is normalised to 1 substitution/site/unit time. Exact `P(t)=exp(Qt)` comes from the eigendecomposition. Among-site rate heterogeneity uses a **discrete gamma** (`gamma_categories` mean-of-category rates, shape `gamma_shape`); each site keeps its category, and inserted sites draw a new one. | shape 1.0, K=4 |
| Rate variation | Per-family multiplier ~ Gamma(`family_rate_shape`, mean 1). Optional neofunctionalisation: a new duplicate evolves `neo_factor`× faster for `neo_duration` branch-length units. | shape 2.0; neo off (1.0), 0.05 |
| Indels | Poisson(`indel_rate`·d·L) events per lineage segment of length d. Each is an insertion or deletion (p=0.5) of Geometric(mean `indel_mean_len`) residues, capped at `indel_max_len`. Inserted residues are drawn from π. A deletion is skipped if the protein would fall below `min_protein_len`. | 0.02, mean 3, max 20, min 50 |

Sequences evolve *lazily*. Along a gene lineage the simulator accumulates d = Σ rate×time. The exact transition is applied whenever the lineage reaches a gene-tree node (speciation, duplication, transfer, leaf). Because the chain is Markov, this is exact for piecewise-constant rates.

## Outputs (`DIR/`)
- `proteomes/<sp>.fa`: protein IDs `<sp>_g0001`, … numbered in genomic order (chromosome by chromosome).
- `gff/<sp>.gff3`: minimal Ensembl-like gene/mRNA/CDS records. Each gene is a single exon, CDS = 3×(len+1) nt (stop codon included), with random intergenic spacing. `protein_id=` matches the FASTA ID, and genomic order and strand match the simulated order and orientation.
- `species_tree.nwk`: rooted, with branch lengths and all nodes named.
- `truth/gene_tree_<fam>.nwk`: one tree per root family, pruned of lost lineages, with unary nodes collapsed. Branch lengths are expected substitutions/site. Internal labels are `S_<node>` (speciation at species node), `D_<node>` (duplication on the branch *above* `<node>`) or `T_<donor>_to_<recipient>`. NHX comments repeat this (`Ev=`, `S=`, `WGD=1`). The root of a transferred subtree carries `XFER=donor>recipient`.
- `truth/hogs.orthoxml`: true nested HOGs in OrthoXML **0.5**, with `<taxonomy>`, `taxonId` and `TaxRange` on orthologGroups. paralogGroups carry `taxonId` = the node below the duplication branch. The file validates against the official XSD. A transferred subtree is cut out of its donor family and becomes its own HOG family (`famXXXX.xN`). Families reduced to one gene are listed only as genes.
- `truth/ortholog_pairs.tsv`: `gene_a gene_b level hog_family`, where `level` is the speciation node of the LCA. Orthologs are pairs whose gene-tree LCA is a speciation node **and** whose paths to the LCA cross no transfer, so xenologs are excluded. This equals `hogs_to_pairs(truth/hogs.orthoxml)` (tested).
- `truth/events.tsv`: each duplication (tandem / dispersed / wgd), loss (including `wgd_fractionation`), inversion, translocation, HGT and WGD. Columns: time, branch (named by its lower node), family, gene-tree node id, details.
- `truth/gene_families.tsv`: protein → root family and HOG family. `params.json`: all parameters, seed, version, runtime and event counts.

## HOG and pair semantics (shared with `ancestra.orthoxml`)
- **Pairs**: these follow the QfO `map_relations.py` flattening. Two genes are orthologs iff their deepest common group is an orthologGroup and they come from different species. They are paralogs iff that group is a paralogGroup. This is tested against the unmodified QfO code.
- **HOGs at level L** (`hogs_at_level`): an orthologGroup whose level is at or below L is one HOG. A paralogGroup whose duplication is strictly below L is one HOG. Any other group is split into its children. A group without a level is placed at the MRCA of its species, the latest-possible placement used by pyham. Top-level groups whose level lies below L are still reported as HOGs at L, restricted to L's clade (this is the "gene family" view). pyham reports only HOGs that existed at L. With that one convention difference, our HOGs equal pyham's at every level (tested).
- **pyham**: pyham 1.1.x crashes on a top-level paralogGroup that directly contains a geneRef. Such a group arises when a family survives only through one duplicated lineage. `write_orthoxml(..., pyham_compat=True)` splits such groups. Ortholog pairs are unchanged.

## Explicit simplifications (do not over-interpret toy results)
1. The species tree is a Yule tree with a global clock. There is no rate variation among lineages except per family and the neofunctionalisation boost.
2. A single empirical model (WAG) is used for all families. There are no site-specific profiles (each site is i.i.d. given its rate category), no structural constraints and no compositional heterogeneity. Proteins are single-domain: no domain shuffling, fusion/fission, gene birth from non-coding DNA, or pseudogenes.
3. Events are uniform over genes. There is no dosage or functional bias, and young duplicates have no elevated loss rate apart from WGD fractionation. Tandem duplications copy one gene, not blocks. There are no chromosome fusions/fissions and no segmental duplications.
4. Indels are uniform along the sequence, with no preference for loops or termini.
5. HGT recipients are uniform among contemporaneous lineages, with no distance bias. A transfer replaces nothing (always additive).
6. There are no isoforms, no genome sequence, and no fragmented or mis-annotated gene models. GFF coordinates are synthetic.
7. Default divergence: root-to-tip 1.0 subst/site (deepest ortholog pairs differ by ~2 subst/site, ~40–50% identity). This is moderate compared with real deep eukaryotic splits. Use `--tree-height` to make it harder or easier.

## Performance (measured 2026-10-26 on the shared 24-CPU node; single thread, wall clock incl. writing)
| Size | Proteins | True pairs | Wall | Peak RSS |
|---|---|---|---|---|
| 300 genes × 8 species (default) | ~2.4k | ~8.6k | 1.9–3.5 s | 80 MB |
| 2000 × 32 | ~64k | ~1.0M | 27–46 s | 0.4 GB |
| 2000 × 32 + HGT 0.02 + WGD + neo 3× | ~66k | ~1.0M | 42 s | 0.4 GB |

`evaluate_toy` takes about 40 s and 0.6 GB on the 2000 × 32 truth (31 levels). Most of the simulation time is spent in `P(t)` and the per-segment sampling (~130k lineage segments).
