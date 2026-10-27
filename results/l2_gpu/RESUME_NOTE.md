# L2-gpu resume note
- 18:10 pre-registration committed (docs/L2_gpu.md §0). Task 1 starting: sparsification transforms of the all-GPU tables
  (work/gpuv2/ob_gpusw_acc, work/gpuv2/qfo_allgpu_fast) -> merge-only reruns in work/l2gpu/.
- 18:25 code committed (190d3fe): gpu_sparsify option, batched FastTree option, drivers. Queued (work/l2gpu/):
  batch_ob1.sh (OB arms K39 F0528 R207 K10 Crand + max-link diag) -> batch_qfo1.sh (OB X1 + QfO arms) ->
  batch_t23.sh (task 2 ref/batch/batch_inproc on frozen OB stages; task 3 R4full/R4afdb/R4shuf{11,23,47}).
  QfO ProstT5 gap (63,237 proteins) running on GPU: work/l2gpu/prostt5_qfo/run.sh.
  Eval: scripts/l2gpu_eval.py ob|qfo NAME REF RUNS... (refs: work/ob_final/main_accurate_go, work/gpuv2/qfo_cpu_fast).
- 19:45 Task 1 OB round 1: all arms fail E1 (OG-F −5 to −14 vs frozen; results/l2_gpu/t1_ob_*.json). Mechanism traced
  (results/l2_gpu/trace_fam77_*.json): RefOG021 split arises at Chordata local-tree clustering; GPU graph gives a
  28-gene unit degree 281 (vs ~60) and a non-RefOG best partner -> different ambiguity clusters. X1 (GPU scores on
  DIAMOND pair set) also −9 OG-F: OG-F is hypersensitive to any score perturbation -> D2 jitter noise floor queued.
  Queue: batch_qfo1 (QfO round 1) -> batch_t23 (R2 arms MK39/SP3000/SP5000, jitter J11/23/47 + QfO A/J11, T2, T3).
  ProstT5 QfO gap: ~14% after 85 min (GPU shared) -> ETA ~04:00; decide at 01:00 whether to use it.
