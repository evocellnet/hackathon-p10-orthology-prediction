# L2-domains RESUME NOTE (agent `domains`)
Branch l2-domains, worktree work/l2_domains_wt (all code + docs/L2_domains.md there; results in main-repo
results/l2_domains/ are copied into the worktree at commit time).

## 2026-10-27 18:20 CET
- Pre-registration committed (72a9269), docs/L2_domains.md §0.
- Pfam-A 38.2 downloaded to data/extra/l2_domains (md5 verified, hmmpress done).
- Stage-1 annotation running: work/l2_domains/pfam_stage1 (OB dev families 58,130 + TF-A dev genes 3,581),
  8 single-thread hmmsearch chunks, ETA ~18:50. Log work/l2_domains/pfam_stage1.log; done when runtime.txt has "done".
- Code: ancestra/domains.py, merge.py arch_mode hooks (default off), tests/test_domains.py (5 pass).
  scripts/L2d_standalone.py (step a), scripts/L2d_merge_arm.py (step b).
- Next: run L2d_standalone.py --domtbl work/l2_domains/pfam_stage1/pfam.domtbl --out results/l2_domains;
  gate decides QfO family annotation (446k proteins, ~26 CPU-h) and merge arms.

## 2026-10-27 18:45 CET
- Stage 1 annotation done (4.9 CPU-h). Standalone done: results/l2_domains/standalone.json, doc §2. Gate passes
  (OB hard, via between-RefOG pairs); within-orthogroup and TreeFam-A dev signal null.
- Stage 2 annotation (QfO 442,858 proteins) running, 6 procs: work/l2_domains/pfam_stage2 (ETA ~23:15).
  If interrupted: rerun scripts/L2d_hmmsearch.sh work/l2_domains/fasta/stage2.fa work/l2_domains/pfam_stage2 6.
- OB arms: off running (results/l2_domains/arms/ob_off_none.json), then queue (logs work/l2_domains/logs/ob_queue.log);
  queue script skips finished outputs, so rerunning the same command resumes.
- After stage 2: QfO queue: scripts/L2d_queue.sh qfo 2 'off:none gate:none ...' work/l2_domains/pfam_stage1/pfam.domtbl
  work/l2_domains/pfam_stage2/pfam.domtbl ; then scripts/L2d_analyze.py results/l2_domains/arms.

## 2026-10-27 19:45 CET
- OB off reproduces OG-F 79.69 (pair F1 0.8035), 1894 s merge at 2 FastTree workers. OB gate real: OG-F 76.82,
  pair F1 0.7896 (worse). OB queue continues (veto, bonus, then shuffles), ~31 min per run, ETA ~01:30.
- Stage 2 at ~30% (ETA ~22:10); QfO queues chained (logs/qfo_chain.log) start automatically after it.

## 2026-10-27 20:45 CET
- OB real arms done: gate OG-F 76.82 / pair F1 0.7896; veto 79.69 / 0.8035 (61 seeds vetoed, no dev change);
  bonus 76.67 / 0.7950. Off 79.69 / 0.8035. Shuffles running (gate 11/23/47, veto, bonus), ETA ~01:30.

## 2026-10-27 22:10 CET
- Stage 2 done (26.4 CPU-h, 3.2 h). QfO off reproduces B3 exactly (TP 8607 FP 722 FN 5386).
- QfO real: gate TPR .6534 PPV .8976 (PPV guard fails), veto .6168/.9228 (tiny), bonus .6119/.9246.
- OB queue 1 + queue 2 (bonus:47 bonus:23 veto:47) running; QfO queues A/B running.
- Planned post hoc (labelled): QfO/OB ungated tol 0.8 (is architecture adding anything beyond relaxation?),
  outputs to results/l2_domains/posthoc/.

## 2026-10-27 22:58 CET
- All QfO arms + shuffles done; analysis (scripts/L2d_analyze.py) -> results/l2_domains/integration.json.
  No arm passes: gate primary +0.018 [+0.003,+0.034] but PPV -0.025 and OB OG-F -2.87; gate_deep CI incl 0, PPV -0.032;
  veto +0.0013 [0.0000,+0.0044] (one family; LB not > 0); bonus -0.0017.
- Post hoc QfO ungated tol0.8: .6477/.8869 (gate .6534/.8976 dominates it); ungated deep_tol0.8 .6499/.8754.
- Remaining: OB shuffles (veto 11/23, bonus 11/23) + OB ungated tol0.8 post hoc (results/l2_domains/posthoc/), ETA ~00:15.

## 2026-10-28 00:20 CET — DONE
- All 30 pre-registered arms + 3 post hoc done; docs/L2_domains.md complete (§0–§9). Verdict: no arm passes; no frozen
  config. pytest 157 passed / 5 skipped. Committed on branch l2-domains (not pushed). Nothing left running.
