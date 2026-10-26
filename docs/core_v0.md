# Ancestra core v0 (sequence only): S2 scoring + S3 level-synchronous merge

Status: v0, toy-verified (2026-10-26, T+1h–T+2h). Code: `ancestra/io.py`, `ancestra/scoring.py`,
`ancestra/merge.py`, `ancestra/pipeline.py` (`python -m ancestra run`). Tests:
`tests/test_merge.py` (15 tests, including exact hand-made HOGs). Benchmark:
`scripts/run_core_v0_toy.py` → `results/core_v0/toy_grid.jsonl`, tables by
`scripts/summarize_core_v0.py` → `results/core_v0/toy_tables.md`. Error analysis:
`scripts/diagnose_core_v0.py`.

```
python -m ancestra run --proteomes DIR --tree tree.nwk [--families fam.tsv --families-cols 0,1] \
    --out OUTDIR --threads 24 [--backend diamond|mmseqs] [--arm max+repair] [--config cfg.json]
```
Outputs: `OUTDIR/hogs.orthoxml` (OrthoXML 0.5, XSD-validated), `OUTDIR/pairs.tsv` (QfO
flattening of the same HOGs, labelled secondary export), `OUTDIR/orthobench_ogs.txt`
(root-level groups), `OUTDIR/manifest.json`, `OUTDIR/stages/<stage>/` (checkpoints).

RESULTS_PLACEHOLDER
