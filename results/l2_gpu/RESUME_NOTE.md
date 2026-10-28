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
- 23:30 Task 1 DONE (doc §2; results/l2_gpu/task1_summary.{json,md}): no arm passes E1-E3; 1% jitter of the CPU table
  costs -3..-7 OG-F / -0.01 macro (OB dev at a tuned optimum); QfO R207 best (TPR +0.021 / PPV +0.008 at 0.59 CPU-h).
  Task 2: micro-bench (results/l2_gpu/t2_bench_ob1000.json): FastTree -n byte-identical on 1000/1000 jobs but only
  ~3% faster; pyfamsa no faster, 990/998 identical alignments. Full OB merge with tree_ft_batch=200 queued
  (batch_late.sh) to confirm hogs byte identity (ref rerun = frozen cc261e73 byte-identical).
  Task 3 OB: R4full = R4afdb on dev (OG-F 79.93, macro +0.0042 [-0.0005, +0.011] vs frozen -> fails keep rule);
  shuffled 3Di +0.0021. QfO R4full / R4afdb / R4shuf11 running (batch_late.sh); foldseek GPU bench running
  (work/l2gpu/t3/fsbench). QfO ProstT5 gap done: 63,237 proteins, 5h01 shared GPU, md5 764665b6.
- 01:00 ALL TASKS DONE. docs/L2_gpu.md complete (summary, pre-reg + deviations D1-D3, results, mechanism, Pareto, cost,
  literature, failed list). Verdicts: task 1 no arm passes (R207 = QfO-best speed mode, documented, no frozen config);
  task 2 FastTree -n equivalent (4/1.39M pairs differ) but -1.5% only; pyfamsa rejected; task 3 R4-full = R4-afdb,
  fails keep rule (OB macro +0.0042 [-0.0005, +0.011]; QfO +0.0021/+0.0024 within CI); Foldseek GPU slower.
  pytest: 166 passed, 1 skipped. Nothing left running.
