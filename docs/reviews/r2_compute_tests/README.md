# r2 compute tests (2026-10-26 15:15-15:42 CET, node load avg 50-100 on 24 cpus: other agents' foldseek jobs; timings are upper bounds)
Input: 78 canonical QfO 2022_02 (updZEBRAFISH) proteomes, N=988,778, 445 Maa.
- `diamond deepclust -d qfo.faa --threads 24 -M 100G` (v2.2.8 defaults): wall 21:02, RSS 9.2 GB, 328,252 clusters (237,024 singletons); max 2,871; sum g^2 = 5.0e7.
  `dcstats.py`: FastOMA root HOGs fully inside one cluster 36.8% (20.8% for HOGs>=10 genes); gene-weighted largest share 0.636.
- `famsizes.py` FastOMA root HOGs: 62,551 groups, max 3,527, sum g^2 5.1e7.
- `ccsizes.py` connected components of OrthoFinder QfO pairs: max 114,902 genes; sum g^2 1.35e10 (98% from the top component).
- Backend test on 7,951 FastOMA HOGs (85,032 proteins, 5.04e6 in-family ordered pairs; `mkpairs.py`):
  mmseqs tsv2db (custom prefilter) + `mmseqs align -e 1e-3 --alignment-mode 3`, 24 thr: 457 s, 3.92e6 hits (reference).
  DIAMOND, all 85k packed in one block, `-k 0 -e 1e-3 -b 4 -c 1 -p 24`: mid-sensitive 36 s (86.6% of mmseqs in-family pairs), sensitive 126 s (92.9%).
