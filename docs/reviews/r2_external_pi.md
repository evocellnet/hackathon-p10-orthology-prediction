# Review round 2: external PI (comparative genomics)

Reviewed DESIGN.md v1 at 2026-10-26 15:16 CET (T+0h36, 47h23 left, from `scripts/clock.sh`).

## Verdict
v1 addresses most of my round-1 points. Family partition with a level-synchronous merge is the right fix for the cost problem:
- The 2e11-pair deep merges become Σ_f g_f² of local work.
- Each node parallelises across families, which removes the critical path.
- Synteny can be scored at X, because every family's matches at X are known.

New problems: S0 fissions cannot be repaired, large-family semantics are undefined, and synteny (the distinctness carrier) is scheduled as a disposable layer. The brief's "alternatives to a universal tree" is half covered. The design is achievable if v0 stays sequence-only.

## Must
**M1. S0 is a hard partition, so its fissions are permanent.**
- A fusion at S0 costs compute and the merge fixes it. A fission loses recall for good, and the root rescue cannot join families. Tune S0 for completeness.
- Pre-register acceptance before S2 starts: at least 95% of dev RefOGs within one family. Report the largest family and Σg²/(N²/2) on QfO.
- Add a **cross-family rescue**: search family medoids against each other, and merge only when species sets are complementary and the coverage guard passes.

**M2. Define large families and missing pairs.**
- Rule: if g_f ≤ cap (for example 5k), compute exact all-pairs. If larger, re-cluster at stricter coverage. If the family is still too large, run Hieranoid-style medoid search per node, restricted to the family.
- **A missing pair is unknown, not zero.** Aggregates use observed pairs only.
- For the in-family restriction, pack families into bins of about 50k sequences and run DIAMOND per bin with -k at least the largest family in the bin, or run `mmseqs align` on an explicit pair DB. Do not use a global run with -k 50.

**M3. Synteny is what makes the method distinct, so it is not a 4 h layer.**
- Without it the method is a database-free FastOMA. Prepare the data now, start S5 by T+14h, and give it at least 8 h.
- Define adjacency concretely: A~B at X if their members lie within w genes in some descendant species, weighted by the number of such species. Do not reconstruct a linear ancestral order.
- Pre-register:
  - the target stratum (WGD ohnologs in teleosts and yeast; top-2 sequence scores within X%);
  - a modest effect size;
  - the fallback novelty claim if the result is null.

**M4. Guide-tree alternatives (required by the brief).**
- S1 has only topology arms. Add the discordance arms from round 1:
  - strict universal tree;
  - delayed merge;
  - per-family local labelling of flagged families.
- Report the fraction of families flagged, split into prokaryotes and eukaryotes. Use Zombi-LGT as the truth on dev, and STD-Bacteria and G_STD2_Luca at the freeze.

## Should
**S1. Hidden paralogy.**
- Differential loss produces a clean 1:1 match, so an "ambiguity" trigger never fires on it.
- Outgroup in-family scores are already computed. Apply an OMA-style witness test to 1:1 matches, and send failures to the local-tree path.

**S2. Normalisation.**
- Add an S2 arm with OrthoFinder-style length and bitscore normalisation before aggregation. This addresses rate heterogeneity (r1_prior_work, item 9).

**S3. Dev coverage.**
- The 35 RefOGs cover only bilaterians, while deep merges drive QfO precision.
- Adopt the stats reviewer's 30% TreeFam-A dev split, or record in §10 why it was rejected.
- Use the consensus of competitor pairs as a weak dev signal for bacteria.

**S4. v0 scope.**
- v0 by T+8h: S0, binned in-family DIAMOND, mutual-best plus in-paralogs, fission repair, OrthoXML.
- Witness test and local trees after the T+12h gate.
- Four parallel agents from now: S0+S2; merge kernel on toy data; runners and competitor timings; synteny data plus Zombi (1 h box).

## Could
- In §10, log the rejected round-1 suggestions, and include the per-node cost of cross-family rescue in the scaling claim.

## Competitiveness
Exact in-family scores put paralog resolution on par with OMA and OrthoFinder, ahead of Hieranoid's profiles. QfO precision depends on conservative deep merges, because one false LUCA merge adds n₁·n₂ false pairs. Beyond the items above, no further changes are needed.
