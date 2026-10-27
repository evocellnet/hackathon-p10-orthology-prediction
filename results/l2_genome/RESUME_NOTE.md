# RESUME NOTE (genome, L2)
- 18:15 pre-registration committed (docs/L2_genome.md §1). Next: implement GSR (gsr_* flags in merge.py), intron table.
- 18:15 GSR implemented (ancestra/merge.py gsr_* flags, default off; pytest 156 pass) + opt-in FastTree memo (ANCESTRA_FT_CACHE).
  Running: scripts/L2_gsr_ob.py --arms off F_c F_u + shuffled (log work/l2_genome/logs/ob_F.log). Intron standalone run by a sub-agent (scripts/L2_introns_*.py).
  Next: QfO fast arms (scripts/L2_gsr_qfo.py), compare with scripts/L2_compare.py, standalone scripts/L2_gsr_standalone.py.
- 18:35 OB F arms done (results/l2_genome/ob_arms, keep_ob.json): off reproduces 79.69/0.8026/706. F_c macro -0.0014 [-0.0050,+0.0008], F_u +0.0006 [-0.0046,+0.0061] -> keep FAILS (null); S8 fission cand 706 -> 499/446 vs shuffled 706-708 (positive control real).
  GSR decisions touch dev RefOG genes rarely (label coverage <1%) -> OB dev has little power.
  Introns: STOP (hard re-rank -0.071 [-0.149,+0.007]; AUC hard 0.578 vs shuf 0.499).
  Running: OB R arms (ob_R.log), QfO fast arms (qfo_fast.log). fast_off reproduces .6151/.9224.
- 19:45 ALL DONE. Verdict: no arm passes keep rule (OB macro and TF-A TPR CIs include 0); introns STOP; tau STOP. Doc complete (docs/L2_genome.md §2-§10). No frozen config shipped.
