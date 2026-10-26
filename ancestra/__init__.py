"""Ancestra: orthology inference by progressive reconstruction of ancestral proteomes.

Design-independent components available so far:

* :mod:`ancestra.simulate`     -- toy genome-evolution simulator with exact truth
* :mod:`ancestra.orthoxml`     -- OrthoXML writer/reader, QfO-style flattening, OrthoBench export
* :mod:`ancestra.evaluate_toy` -- compare predicted HOGs against simulated truth
* :mod:`ancestra.newick`       -- minimal Newick parser/writer and tree index

Core v0 engine (docs/core_v0.md; ``python -m ancestra run``):

* :mod:`ancestra.io`           -- proteomes -> integer gene IDs, guide tree, families
* :mod:`ancestra.scoring`      -- S2 within-family scores (DIAMOND / MMseqs2 / hybrid)
* :mod:`ancestra.merge`        -- S3 level-synchronous bottom-up merge -> HOGs
* :mod:`ancestra.pipeline`     -- checkpointed end-to-end CLI with manifest
"""
__version__ = "0.1.0"
