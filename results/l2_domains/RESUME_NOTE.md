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
