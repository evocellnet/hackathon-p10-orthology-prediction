# S8: validating apparent gene losses with genome DNA (demo)

Status: demo, 2026-10-26, T+3h09 to T+3h45 (from `scripts/clock.sh`), against a 3 h time-box. Design: DESIGN.md v2 §3 row S8 and PROMPT.md, paragraph "Genome (DNA) sequence".

- **Code:** `ancestra/geneloss.py` (new). It only reads core outputs; no core file was edited.
- **Scripts:** the driver is `scripts/run_S8_geneloss.py`. Supporting scripts are `scripts/summarize_S8.py`, `scripts/S8_homininae_check.py`, and the GFF parsers `scripts/S8_parse_gff_{genes,cds}.py`.
- **Tests:** `tests/test_geneloss.py` (7 pass). They include an end-to-end miniprot run on a synthetic genome holding an intact gene, a frameshifted gene and a deleted gene.
- **Outputs:** `results/S8_geneloss/`.
  - `final_calls.tsv`: one row per candidate.
  - `{cand,pos,mask,neg,gw,known}_<species>.tsv`: per-run tables.
  - `homininae_check.tsv`.
  - `tables.md` and `summary.json`.
  - `manifest.json`: commands, parameters, and md5 of the inputs and code.
- **Data:** DATA.md §13 "S8 genomes" covers the Ensembl 99 `dna_sm` genomes of the 6 mammals, with md5 and Ensembl `sum` (all match). Tool versions are in the ENVIRONMENT.md addendum (miniprot 0.18-r281, samtools 1.24).
- **Input HOGs:** `work/core_v0_ob/run_v2/hogs.orthoxml`. This is the sequence-only OrthoBench run with S0 families (`results/core_v0_orthobench/runs.tsv`, label `run_v2`). HOGs are cut at **Theria** (= Mammalia in the OrthoBench tree, the 6 mammals) and **Amniota** (outgroup check).
- **Cost:** wall time 4 min (human, mouse, rat) plus 9.5 min (chimp, dog, opossum) at 6 threads. Peak RSS was 33 GB, from the miniprot genome indexes.

## 0. Summary
- **Candidates.** 2,523 Theria HOGs are present in 5 of the 6 mammals and absent in exactly one, with ≤2 copies per present species. By missing species: dog 1,009, opossum 937, rat 301, chimp 164, human 67, mouse 45.
  - An opossum-only absence is kept only if the Amniota HOG has a chicken or fish member. Otherwise it could be a eutherian gain.
- **The locus method works on the controls.** There are 600 positive controls: complete 1:1 HOGs whose target gene is hidden and then searched for at the synteny-predicted locus.
  - The intact ORF on the hidden gene is recovered in **92.2%** [89.7, 94.1] pooled. By species: human 95, chimp 96, mouse 98, rat 90, dog 89, opossum 85.
  - The pre-set ≥90% bar is met pooled and in 5 of 6 species; opossum falls short.
  - False "intact" calls: **0/512** random genomic windows and **5/600** masked loci (the hidden gene replaced by N), giving P(intact | absent) = 0.45%.
- **What the 2,523 apparent losses turn out to be:**
  - **1,348 (53%) are not losses.** The gene is in the genome:
    - 609 are annotated but split off into a target-only HOG (a HOG fission, i.e. an orthology-inference miss);
    - 199 are annotated but assigned to another HOG;
    - 481 are intact unannotated ORFs or gene-model errors;
    - 59 are relocated intact copies.
  - **638 (25%) are confirmed losses:** 351 pseudogenized and 287 deleted with conserved flanks.
  - **537 (21%) are unresolved:**
    - 257 are disrupted alignments over an *annotated* coding model (an assembly indel, or a pseudogene modelled as coding);
    - 63 are partial hits;
    - 217 have no bracketed locus or sit in an assembly gap.
- **Precision of the loss calls.** This is estimated from the controls, using each species' P(loss | present) and treating a disrupted alignment over the hidden gene as a false loss.
  - Pooled precision is **0.80**.
  - It is high in the well-assembled genomes: human 0.97, mouse 0.97, chimp 0.93.
  - It is lower where the assembly or annotation is noisier: rat 0.86, opossum 0.85, dog 0.67.
- **Independent check: GENCODE in human.**
  - **20/20** human pseudogene (`B_pseudo`) calls overlap an Ensembl/GENCODE-annotated pseudogene. Examples: GGTA1P, CYMP (chymosin), TDH, NEPNP, CTF2P, CYP2G2P, SULT1D1P, AOX2P, NUDT16P1, and OR pseudogenes.
  - In mouse, 10/13 do.
  - But 5 of the 6 human "intact unannotated ORF" calls also sit on GENCODE *unitary* pseudogenes, such as SLC22A20P and FMO6P. An intact reading frame is necessary but not sufficient for a functional gene (§5).
- **Known biology.**
  - **GULO** (the vitamin C gene) is pseudogenized in human and chimp. Its locus lies between the flanking orthologs; the mouse and rat proteins align with 3 frameshifts and 1 stop, on GENCODE GULOP.
  - **UOX** (urate oxidase) is pseudogenized in human and chimp, with 3 to 4 in-frame stops. It *is* among our candidates as a human-only loss, because Ensembl annotates a chimp UOX "protein".
  - **AIM2** in dog is called deleted, with conserved flanks and 5 of 5 source species agreeing on the locus. This is consistent with the reported loss of AIM2 in dog.
- **The taxonomic placement of losses changes.** The chimp locus was re-searched with the non-hominid proteins only.
  - **13/20** human "pseudogene" losses (including UOX and GGTA1) are also disrupted in chimp. For 12 of them the disruption sits on the *annotated chimp gene*.
  - The loss is at least Homininae-level, and the chimp model is a pseudogene annotated as coding. The HOG should be absent at Homininae, not at human.
  - The other 7 are intact in chimp, i.e. human-specific losses.

## 1. Method (`ancestra/geneloss.py`)

### 1.1 Candidate selection
1. `HogIndex.from_orthoxml(level="Theria", outgroup_level="Amniota")` uses the core's `hogs_at_level`.
2. `select_apparent_losses` keeps HOGs with a member in exactly 5 of the 6 mammals, at most 2 copies per present species, and an outgroup member when the missing species is opossum.
3. Positive-control HOGs come from `select_complete_single_copy`: 9,829 HOGs with exactly one gene in every mammal.

### 1.2 Expected (syntenic) locus: `syntenic_locus`
**Proposals.** Each source species that has the HOG proposes a locus:
1. Take up to k=5 flanking genes on each side of its gene. Use primary-assembly genes from the S5 gene-order tables; alt-haplotype projections are excluded.
2. Map each flanking gene through its own Theria HOG to a single-copy target member.
3. Take the left/right anchor pair with the smallest kL+kR that lies on one target contig, ≤2 Mb apart.
4. The locus is the anchors' span ±200 kb ("bracketed").
5. If only one side maps, the locus is the anchor ±1 Mb ("one-sided").

**Synteny vote.** The bracketed proposals are clustered (same contig, overlapping), and the cluster supported by the most source species wins.
- This was added after a failure mode in the smoke test. The chimp member of several HOGs was an annotated *retrocopy*, so its flanks pointed to the human retrocopy locus.
- With the vote, all bracketed sources agree for 67% of human candidates and for 87 to 98% of the candidates in the other species.

### 1.3 Alignment
miniprot 0.18 (`--gff --outs=0.5 -p 0.5 -N 20`) aligns all HOG proteins from the present species against the locus DNA. The parsed fields are:
- query coverage;
- `fs` (frameshifts) and `st` (in-frame stops);
- intron count, from the CIGAR;
- CDS blocks, from the GFF.

### 1.4 Hit filters
**Paralog filter.** A hit whose CDS blocks overlap ≥20% with the CDS of an annotated target protein in *another HOG that has members in other species* is a paralog hit and is set aside.
- If the overlapped protein belongs to this HOG or to a target-only HOG, the hit counts as the annotated target gene.

**Retrocopy filter.** An intronless hit in a HOG whose median member has ≥2 coding exons is set aside. Processed pseudogenes and retrocopies are not the lost gene.

### 1.5 Locus classes (`classify_hits`)
- `A_annotated`: intact hit (≥80% of the query aligned, fs=st=0) on the CDS of an annotated target protein from this HOG or a target-only HOG. The gene exists; the HOG inference split it off.
- `A_gap`: intact hit not on any annotated input protein.
  - **Majority rule:** more query proteins must give an intact best hit than a disrupted one. Otherwise the call is `B_pseudo` with `conflict=True`.
  - This catches pseudogenes modelled as coding in one source species. UOX in human is the example: the chimp model aligns intact, but the mouse, rat, dog and opossum proteins all show stops.
  - On the positives, the rule costs 1.4 points of sensitivity (0.922 → 0.908).
- `B_pseudo`: ≥30% of the query aligned, with frameshifts or stops.
- `B_partial`: ≥30% aligned, not disrupted, but not intact.
- `C_deleted`: no hit (after the filters) in a *bracketed* locus with <20% N.
- `D_undetermined`: one-sided or no locus, or the locus is ≥20% N.

### 1.6 Genome-wide search and post-rules (`summarize_S8.py`)
For every B/C/D candidate, the same proteins are searched against the whole genome with a miniprot index. The paralog and retrocopy filters are applied the same way. Final classes:
- `C_relocated` / `D_relocated`: an intact copy exists outside the locus.
- `A_otherHOG`: a C call whose set-aside "paralog" hit is ortholog-like. The hit must have qcov ≥0.8 and identity ≥ the 5th percentile of the identity that recovered positives reach for the same target × query species (thresholds in `summary.json`). The gene is annotated, but the HOG pipeline put it in another multi-species HOG.
- `B_annot`: a B call whose disrupted alignment lies on the CDS of an annotated target protein from this HOG or a target-only HOG. Ensembl models a coding gene there, so this is either an assembly indel or a pseudogene modelled as coding. It is unresolved and *not* counted as a loss.

### 1.7 Controls (100 per species, seed 8)
- **pos:** a complete 1:1 HOG with the target gene hidden. The locus is computed from the other species and searched with their proteins. The control is recovered if the call is `A_annotated` on the hidden gene's CDS.
- **mask:** the same loci with the hidden gene's span replaced by N. Expected outcome: no intact hit.
- **neg:** random windows on the target's gene-bearing contigs, each the size of a candidate locus, searched with that candidate's proteins. Expected outcome: no intact hit.

## 2. Results: classes per species (`tables.md` §1, §8)

| target | n | A_annotated | A_otherHOG | A_gap | B_annot | B_pseudo | B_partial | C_deleted | C/D_relocated | D_undet. | not a loss | confirmed loss |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| human | 67 | 16 | 19 | 6 | 0 | 20 | 0 | 5 | 1 | 0 | 42 (63%) | 25 |
| chimp | 164 | 31 | 19 | 49 | 14 | 29 | 4 | 13 | 4 | 1 | 103 (63%) | 42 |
| mouse | 45 | 14 | 2 | 5 | 1 | 13 | 1 | 2 | 1 | 6 | 22 (49%) | 15 |
| rat | 301 | 63 | 15 | 81 | 16 | 80 | 8 | 24 | 3 | 11 | 162 (54%) | 104 |
| dog | 1009 | 301 | 110 | 234 | 88 | 125 | 28 | 95 | 12 | 16 | 657 (65%) | 220 |
| opossum | 937 | 184 | 34 | 106 | 138 | 84 | 22 | 148 | 38 | 183 | 362 (39%) | 232 |
| **all** | **2523** | **609** | **199** | **481** | **257** | **351** | **63** | **287** | **59** | **217** | **1348 (53%)** | **638 (25%)** |

Notes:
- **Opossum.**
  - 258 of its 937 loci are only one-sided: 1 Mb windows, often next to unplaced scaffolds, so absences there cannot be called as deletions. This is why opossum has 183 `D_undetermined`.
  - It also has 138 `B_annot` (annotated models with frameshifts or stops in the alignment).
- **Dog A_gap.** 178 of the 234 intact "unannotated" ORFs fall *inside the span of an annotated dog protein-coding gene* while overlapping <20% of its CDS. These are mostly CanFam3.1 gene-model errors (merged or partial models, e.g. the AFP/AFM region) rather than truly unannotated genes. In both cases the gene exists.
- **A_otherHOG examples.** In human: CD24, IFITM3, ENSA, SNRPG, ELOB and ATP6V1F (identity 0.94 to 1.0). These are genes with many retrocopies. Here the HOG that lacks human is built around annotated retrocopy models in other species, while the real human gene sits in another HOG.
  - For distant query species (opossum or dog proteins), the identity threshold is low (0.52 to 0.75), so in dog and opossum some `A_otherHOG` calls may be genuine paralogs.

## 3. Controls and precision (`tables.md` §2, §3)

| target | pos recovered | sens [95% CI] | pos → loss (B or C, incl. on the hidden CDS) | mask: intact FP | random windows: intact FP | est. precision of loss calls |
|---|---|---|---|---|---|---|
| human | 95/100 | 0.95 [0.89, 0.98] | 2 | 1/100 | 0/67 | 0.97 |
| chimp | 96/100 | 0.96 [0.90, 0.98] | 3 | 1/100 | 0/100 | 0.93 |
| mouse | 98/100 | 0.98 [0.93, 0.99] | 2 | 0/100 | 0/45 | 0.97 |
| rat | 90/100 | 0.90 [0.83, 0.95] | 8 | 1/100 | 0/100 | 0.86 |
| dog | 89/100 | 0.89 [0.81, 0.94] | 10 | 1/100 | 0/100 | 0.67 |
| opossum | 85/100 | 0.85 [0.77, 0.91] | 9 | 1/100 | 0/100 | 0.85 |
| pooled | 553/600 | **0.922** [0.897, 0.941] | 34 (5.7%) | 5/600 | 0/512 | **0.80** |

How to read the table:
- **Failed positives.** Of the 47 positives that were not recovered:
  - 29 gave a disrupted (23) or partial (6) alignment on the hidden gene: an assembly indel in Rnor_6.0, CanFam3.1 or ASM229v1, or a non-canonical gene model. Human XPO4 is one example, with 1 in-frame stop in all 5 queries.
  - 12 had no bracketed locus.
  - 5 found nothing.
  - 1 was an intact hit on another gene.
  - So in these assemblies about 1 in 10 present-and-annotated genes looks "pseudogenized" at the DNA level. This is the dominant error, and it is why the loss precision is 0.67 for dog.
  - Estimates come from §3 of `tables.md`: expected false losses = (present calls / sensitivity) × P(loss | present).
- **Deletion calls versus pseudogene calls.** The masked controls show that when the gene is absent at a conserved locus, the method almost never invents an intact gene (5/600; mostly a tandem copy in the window).
  - A "disrupted hit" call can still come from a pseudogenized paralog: 11/600 masked loci were called B.
  - So `C_deleted` versus `B_pseudo` is less reliable than "loss" versus "not loss".
- **The "intact ORF ⇒ present" call is precise against chance hits:** P(intact | absent) = 0.45%, giving about 3 expected false present calls among 1,289.
  - It does *not* measure function, as the GENCODE check below shows.

## 4. Independent annotation check (GENCODE/Ensembl biotypes; `tables.md` §4)
For each candidate, the best hit's CDS blocks were intersected with all Ensembl 99 gene features (every biotype).

- **Human.**
  - `B_pseudo`: **20/20** overlap an annotated pseudogene (mostly `unitary_pseudogene` / `transcribed_unitary_pseudogene`).
  - `A_gap`: 5/6 overlap an annotated pseudogene (SLC22A20P, METTL21EP, PRORSD1P, FMO6P, TMEM30CP); 1/6 is protein-coding (SPTSSB).
- **Mouse.**
  - `B_pseudo`: 10/13 overlap an annotated pseudogene (e.g. Zfp862-ps, Dcdc5, Best4-ps).
  - `A_gap`: 2/5.
- **Reading.**
  - The DNA evidence for pseudogenization agrees with the curated reference in every human case.
  - The reverse error is real: an intact, full-length alignment with no frameshifts or stops was found at 5 GENCODE unitary pseudogenes. Their inactivation lies outside what miniprot scores: splice sites, start codon, promoter or transcript level. For example, FMO6P is known to be mis-spliced.
  - Our "annotation gap" calls are therefore an *upper bound* on missing genes. In a well-annotated genome they are more often non-functional loci.
  - For dog, opossum, rat and chimp there is no comparable curated pseudogene set: Ensembl annotates few pseudogenes there.

## 5. Known biology
**Direct tests** use the HOG of the mouse protein, the syntenic locus in the hominid, and the same classifier (`known_*.tsv`):

| gene | target | call | evidence | Ensembl at hit |
|---|---|---|---|---|
| GULO | human | B_pseudo | 4/4 queries disrupted; best (rat): 83% of query, fs=3, st=1 | GULOP, unitary_pseudogene |
| GULO | chimp | B_pseudo | 4/4 disrupted; fs=3, st=1 | (none annotated) |
| UOX | human | B_pseudo (majority rule) | mouse: 100% of query, 3 in-frame stops; only the chimp model aligns intact | UOX, transcribed_unitary_pseudogene |
| UOX | chimp | B_pseudo | 4/4 non-hominid queries: 4 stops, on the annotated chimp "UOX" protein | protein_coding (pseudogene model) |

- GULO is not among the candidates: it is absent in *both* human and chimp, so it fails the 5-of-6 rule.
- UOX is a candidate, because Ensembl gives chimp a UOX protein.

**Homininae check** (`S8_homininae_check.py`, `homininae_check.tsv`). For each of the 20 human `B_pseudo` and 5 human `C_deleted` HOGs, the chimp locus was searched with the mouse, rat, dog and opossum proteins only, ignoring the chimp annotation.
- **13 of the 20 human pseudogenes are also disrupted in chimp**, 12 of them on the annotated chimp gene. Examples: GGTA1 (α-1,3-galactosyltransferase; catarrhine loss), UOX, CYM (chymosin), NEPN, HILS1, MGAT4F, AOX2, SULT1D1, TTC39D and 3 olfactory receptors.
  - These are shared Homininae (or older) losses.
  - The chimp "protein" is a pseudogene modelled as coding, and the HOG is wrongly marked present in chimp.
- 7 are intact in chimp (majority rule, chimp annotation ignored), i.e. candidate human-specific losses: CTF2, TMED11, GUCY1B2, NUDT16, ADAM20-like, CYP2G and TDH.
- The 5 human `C_deleted` HOGs have an intact chimp gene at a *different* locus from the chimp HOG member. They are the retrocopy-model HOGs described in §1.2.

**Other recognisable cases:**
- **Dog AIM2:** `C_deleted`, 5/5 sources agree on the locus, 1.7% N. This is consistent with the reported absence of AIM2 in dog.
- **Dog LEP (leptin):** `A_annotated`. The dog leptin gene is annotated and intact, but the HOG put it in a dog-only group.

## 6. How this feeds HOG presence/absence and taxonomic placement
1. **Presence/absence correction.** Of the apparent single-species absences in our HOGs, **53% are not gene losses.**
   - 609 (24%) are HOG fissions: the gene is annotated but in a target-only HOG. These are direct repair candidates for the merge, because the locus evidence names the target protein (`overlap_gene`).
   - 199 are assignments to another HOG.
   - 481 are annotation or gene-model gaps.
   - Such absences should not be counted as losses. In OrthoXML terms, the HOG's species set at Theria should include the target, either by re-attaching the named protein or by recording "present, unannotated" with the genomic coordinates.
   - Confirmed pseudogenes and deletions (638) are evidence *for* the absence.
2. **Taxonomic placement of losses (and of the HOG).**
   - The Homininae check moves 13 losses from the human branch to the Homininae branch or earlier.
   - The same test applied to all HOGs would also fix gains. For example, when opossum lacks a HOG and `A_*` shows the gene in opossum, the HOG is placed at Theria and not at Eutheria. Here the outgroup rule already assumed that, so all opossum `A_*` calls are confirmations.
   - For deeper levels this needs per-node DNA searches. That is cheap: about 0.1 to 1 s per locus, plus about 1 min per genome index.
3. **Downstream analyses.** These are gene-loss counts per branch, phylogenetic profiles, and the gene sampling in QfO's species-tree-discordance (STD) test.
   - All of them should use `confirmed_loss` (B_pseudo + C_deleted) rather than "missing from the HOG".
   - On this data, the terminal-branch "losses" drop from 2,523 to 638 confirmed plus 537 unresolved. At least 53% are removed as present (65% in dog).
   - The unresolved 21% (B_annot, B_partial, D) are exactly the cases to treat as missing data.
4. **What it does not do.** It does not prove function (§4), and it is weakest in fragmented or indel-prone assemblies (dog, opossum, rat; §3). A production version would:
   - use a second assembly or reads to separate assembly indels from real frameshifts;
   - check splice sites and the start codon;
   - weight the calls by the per-species control error.

## 7. Limitations and deviations
- **Scope.** This is a demo layer and is not wired into the merge. Nothing here feeds back into the HOGs or the benchmarks.
  - It uses one run (`run_v2`), one level (Theria) and absences in one species only. Losses in two species, e.g. GULO in human and chimp, need the clade-level version used in §5.
- **Locus parameters were not tuned.** These are k=5, pad 200 kb (1 Mb one-sided), max span 2 Mb, 0.8/0.3 coverage, and 20% CDS overlap. They were set once.
  - The synteny vote, the retrocopy filter, CDS-level (not gene-span) overlap, and the majority rule were added after the 30-candidate human smoke test.
  - That test was inspected against GENCODE names, so it is not a blind evaluation. The controls were run only after these rules were fixed, and so were the GENCODE check on all calls and the dog, opossum, rat and chimp runs.
- **Thresholds behind two post-rules.** The `A_otherHOG` identity thresholds come from the positive controls of the same run (5th percentile).
- **Truth for losses.** There is no ground-truth set of real losses except the GENCODE pseudogene annotation in human and mouse and a handful of literature cases. Hence "precision" is estimated from controls (§3) plus the GENCODE agreement (§4).
- **miniprot da/do tags.** The non-canonical splice counts were recorded (`splice_nc`) but are not used. They are inconsistent; `da` can be −1.
