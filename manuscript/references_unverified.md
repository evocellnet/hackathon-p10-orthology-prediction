# References not included in references.bib (could not be verified)

Checked 2026-10-27 against Crossref (DOI lookup and `query.bibliographic`) and PubMed.

| Reference | Needed for | Reason it is not in references.bib |
|---|---|---|
| Csardi G, Nepusz T (2006). The igraph software package for complex network research. *InterJournal Complex Systems* 1695. | igraph (python-igraph 1.0.0 is in `ENVIRONMENT.md`) | The paper has no DOI and no Crossref record. The Crossref bibliographic query returned only the CRAN package DOI (10.32614/cran.package.igraph, a different work: the R package) and unrelated items. Not in PubMed. Replaced by the verified igraph 1.0 paper, Antonov et al. 2026 (`antonov2026igraph`), which matches the python-igraph 1.0.0 version used. |
| Douze M et al. (2024). The Faiss library. arXiv:2401.08281 (DOI 10.48550/arXiv.2401.08281) | FAISS | This is a DataCite (arXiv) DOI and returns "Resource not found" at Crossref. It is not in PubMed. Replaced by its published version (IEEE Trans Big Data 12:346, doi 10.1109/tbdata.2025.3618474, `douze2025faiss`), found via Crossref. |

All other requested references were matched in Crossref (two completed from PubMed, see `references_check.tsv`).

Not attempted (outside the requested scope, marked [U] in LITERATURE.md): HOGPROP (Warwick Vesztrocy 2020), FAMSA2, VeryFastTree, Blanquart 2016 splicing orthology, pSONIC (Conover 2021), OrthoFFGC, BBH/RSD, and the 2026 bioRxiv preprints (Novabrowse, TOGA2).
