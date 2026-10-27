# L2-plm resume note
Updated 2026-10-27 20:24 CET (scripts/clock.sh) (agent `plm`). **Status: DONE, verdict negative.**

- docs/L2_plm.md complete: pre-registration (commit 954d85d, before any dev result), data, methods, deviations
  D1/D2, results, cost (~3 GPU-h), verdict, literature, failed/not-pursued list.
- results/l2_plm/discrimination_final.json (all arms + controls), tables.md, complement_final.json,
  head_*.json (curves + md5 of weights in data/extra/l2_plm/heads/), selflabels*.json, knn_guard.json.
- No arm passed step (a) -> no integration (idea 3) and no frozen config.
- pytest: 157 passed (with conda_envs/orth/bin on PATH for diamond). No changes to ancestra/, bench_repo/, tests/.
- Work data (not committed): work/l2_plm/ (embeddings, alignment caches, triplet dumps), data/extra/l2_plm/.
