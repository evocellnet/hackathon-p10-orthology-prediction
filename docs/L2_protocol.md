# Layers round 2 (L2): shared protocol

Started 2026-10-27 18:00 CET, after the user asked for more extra-data work while the scaling runs finish
(domain architecture, proteomics, co-fractionation, protein language models, and new uses of structure,
synteny, expression and GPU steps). This is **post-freeze exploration**. The frozen method (`freeze-v1`,
presets `accurate`/`fast`) and every reported result stay unchanged. A layer that passes here is reported as
a candidate for a later version, never folded into the frozen numbers.

## Rules (all L2 agents)

1. **Code.** Work in your own git worktree, `git worktree add work/l2_<name>_wt -b l2-<name> main` (GPU agent:
   branch from `v2-gpu`). New behaviour sits behind flags whose default is off, so the frozen presets give
   byte-identical output. `pytest -q` must pass on your branch. Never edit `bench_repo/` metrics or tests.
2. **Shared files.** Do not edit REPORT.md, RESEARCH_TREE.md, DATA.md, LESSONS.md, CLAIMS.md or `manuscript/` on
   main. Write everything to `docs/L2_<name>.md` and `results/l2_<name>/`. The PI merges.
3. **Data.** Downloads go to `data/extra/l2_<name>/` (gitignored). For each source, record the URL, version,
   retrieval date, licence and md5 in the "Data" section of your doc. Never commit data. Restricted or
   licence-unclear data are used locally only and flagged.
4. **Dev only.** Use the OrthoBench dev split (34 RefOGs; `splits/`), the pair-level dev labels
   (`scripts/S5_refog_labels.py`, `scripts/orthobench_dev_pairs_eval.py`; hard cases as defined in
   docs/S5_synteny.md), and TreeFam-A dev (359 families; `splits/treefamA_split.json`,
   `scripts/qfo_treefam_subset.py`). **No test RefOGs, no TreeFam-A test families, no SwissTree, no other QfO
   test.** Test evaluation is decided by the PI and run by a different agent.
5. **Pre-registration.** Before you look at any dev integration result, write and commit the "Pre-registration"
   section of your doc: hypothesis, arms, primary dev metric, controls and keep rule. Changing it afterwards
   must be logged as a deviation, with the reason.
6. **Two-step test of every layer.**
   - (a) Standalone signal: does the layer separate orthologs from paralogs on dev, above all on hard cases?
     Report AUC with a bootstrap CI, a shuffled-layer control (fixed seeds 11/23/47) and coverage. If the
     standalone signal is not above shuffled, stop there and record it as a negative result.
   - (b) Integration into the pipeline (re-ranking, veto, rescue, repair, new feature). **Keep rule:** on the
     primary dev metric, Δ(real − off) has a paired-bootstrap 95% CI (B=2000, seed 11) whose lower bound is
     above 0; real beats shuffled by more than the seed spread; and no guard metric gets worse (PPV by >0.01,
     OrthoBench dev OG-F by >1.0, OrthoXML validity via `verify_orthoxml`).
7. **Honesty.** Fixed seeds. Report every arm you ran, including failures. Never tune on what you then report
   as held-out. Say clearly when a dev effect sits within the partition or seed noise (±2–4 OG-F on 34 RefOGs;
   prefer macro pair F1 and TreeFam-A dev TPR/PPV).
8. **Literature.** Credit the prior ideas you build on. Every reference is checked (PubMed or Crossref DOI)
   before it goes into your doc.
9. **Resources.** The scaling runs (`work/scale_large`, `work/scaling`, `results/scaling*`) must not be
   disturbed. Respect your thread cap. Run `nice -n 5`. Check `uptime` before big jobs; if the load is above
   44, wait or use fewer threads. Keep RAM under 80 GB. The A100 (40 GB) is shared by two agents, so keep
   GPU memory under 18 GB and do not use exclusive mode.
10. **Time.** Read the time with `scripts/clock.sh`. Final results and docs are due by **2026-10-28 08:30 CET**.
    Update `results/l2_<name>/RESUME_NOTE.md` at least every hour, so the work survives an interruption.
11. **Deliverables.**
    - `docs/L2_<name>.md`, containing: pre-registration; data; methods; results table(s) with CIs and
      controls; cost (CPU-h, GPU-h, preprocessing); verdict; literature credit; and a 3–6 line "failed or
      not pursued" list for the manuscript's S1 Text.
    - `results/l2_<name>/*.json`, each with a manifest (git commit, config, inputs, runtime).
    - Commits on your branch, ending with the attribution trailer below.
    - For any arm that passes the keep rule: the exact frozen config file that a test agent will run once.

Commit trailer:
```
Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01QAWokgYzUPevsvsYAsatDC
```
Do not push.

## Agents

| name | layer | threads | GPU |
|---|---|---|---|
| domains | domain architecture (Pfam/InterPro, TED CATH domains) | 8 | no |
| network | interaction networks, co-fractionation (CF-MS), protein abundance | 4 | no |
| plm | protein language models (embeddings, embedding alignment, a contrastive head trained on self-derived labels) | 4 | yes |
| gpu | GPU ideas dropped in GPU_v2 (sparsified all-GPU S2, GPU structure rescue at full coverage) | 8 | yes |
| genome | synteny-driven fission repair and recall at the LCA, gene structure (intron positions), new uses of expression | 8 | no |

## Known failure modes this round should target (from earlier work)

- **Recall at the merge.** 78% of missed TreeFam-A dev pairs are merge failures at the pair's LCA: lineages
  pass up unmatched and are later absorbed as paralogs (RESEARCH_TREE 8b3). On QfO, recall is our main gap
  against OrthoFinder (TreeFam-A test TPR −0.160).
- **Fission.** S8 found 808 orthology fission errors among 2,523 apparent mammal gene losses (docs/S8_geneloss.md).
- **Family completeness.** S0 is capped near 0.82 dev RefOG completeness by sequence alone. 3Di rescue reaches
  0.853, but without a guard it fuses RefOG010.
- **Earlier layers that failed and why.** 3Di in merge scoring (paralogs share folds: deep AUC 0.694 vs 0.751
  for sequence); expression (right only 40–47% when it decides); synteny (strong on dev, not measurable on
  test, since local trees already fix those cases); isoforms (killed).
