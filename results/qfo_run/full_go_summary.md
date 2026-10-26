# Run summary: `work/qfo2022/full_go`

| stage | wall s | CPU s (self+children) | peak RSS MB (own / max child) | load start→end |
|---|---|---|---|---|
| merge | 347.9 | 345.9 | 8202.8 / 6900.3 | 28.5→31.1 |
| output | 35.53 | 35.3 | 5342.0 / 6900.3 | 31.1→32.4 |

- **families_S0**: 163092 families (30267 multi-gene, 132825 singletons = 13.4% of genes); largest [7464, 7410, 7389, 7364, 7344, 7306, 7288, 7258, 7251, 7203]; Σg²/N² = 0.002299
- **families_after_rescue**: 162721 families (29896 multi-gene, 132825 singletons = 13.4% of genes); largest [7464, 7410, 7389, 7364, 7344, 7317, 7306, 7298, 7288, 7258]; Σg²/N² = 0.002902
- **root HOGs**: 36479 top-level groups (740903 genes; 247875 singleton genes); largest 6528 genes/15 sp @Bilateria, 5549 genes/57 sp @cellular organisms, 4642 genes/53 sp @cellular organisms, 3909 genes/51 sp @cellular organisms, 3831 genes/49 sp @cellular organisms

| level | depth | #species | #HOGs (≥2 genes) | genes in HOGs | max HOG | orthologGroups with this TaxRange |
|---|---|---|---|---|---|---|
| cellular organisms | 0 | 78 | 36479 | 740903 | 6528 | 4687 |
| Archaea | 1 | 7 | 2259 | 11229 | 89 | 1230 |
| Bacteria | 1 | 23 | 9394 | 56322 | 102 | 4708 |
| Eukaryota | 1 | 48 | 32741 | 668085 | 6528 | 8933 |
| Bacillati | 2 | 7 | 4695 | 18198 | 54 | 2879 |
| Metamonada | 2 | 2 | 2755 | 13923 | 1058 | 1019 |
| Methanobacteriota | 2 | 4 | 1739 | 6203 | 61 | 1346 |
| Opisthokonta | 2 | 35 | 22811 | 483565 | 6528 | 5901 |
| Pseudomonadati | 2 | 13 | 5996 | 29099 | 88 | 4202 |
| Sar | 2 | 4 | 6963 | 38061 | 229 | 2728 |
| Thermoproteati | 2 | 3 | 1105 | 3559 | 41 | 661 |
| Thermotogati | 2 | 2 | 989 | 2482 | 18 | 684 |
| Viridiplantae | 2 | 5 | 13033 | 113074 | 1711 | 4997 |
| Actinomycetes | 3 | 2 | 2482 | 6555 | 36 | 1844 |
| Alveolata | 3 | 2 | 4545 | 17845 | 229 | 1593 |
| Cyanophyceae | 3 | 2 | 1951 | 4484 | 21 | 1687 |
| Embryophyta | 3 | 4 | 15667 | 101516 | 296 | 7674 |
| Eumetazoa | 3 | 22 | 21301 | 398042 | 6528 | 8598 |
| Fungi | 3 | 12 | 8005 | 75183 | 382 | 3426 |
| PVC group | 3 | 2 | 831 | 2424 | 31 | 420 |
| Pseudomonadota | 3 | 4 | 3851 | 12729 | 42 | 2682 |
| Stenosarchaea group | 3 | 2 | 1113 | 3006 | 32 | 825 |
| Stramenopiles | 3 | 2 | 3847 | 11191 | 65 | 2986 |
| Thermoproteota | 3 | 2 | 800 | 1875 | 29 | 684 |
| Bilateria | 4 | 21 | 21150 | 380631 | 6528 | 8075 |
| Dikarya | 4 | 11 | 8565 | 69588 | 245 | 4345 |
| Gammaproteobacteria | 4 | 2 | 2421 | 5687 | 25 | 1997 |
| Mesangiospermae | 4 | 3 | 18157 | 81002 | 170 | 11425 |
| Ascomycota | 5 | 8 | 7613 | 52118 | 243 | 3552 |
| Basidiomycota | 5 | 3 | 5142 | 15917 | 90 | 4660 |
| Chordata | 5 | 15 | 20007 | 304958 | 6019 | 10617 |
| Poaceae | 5 | 2 | 18312 | 44420 | 63 | 18312 |
| Protostomia | 5 | 6 | 10303 | 68158 | 710 | 5897 |
| Ecdysozoa | 6 | 5 | 9460 | 54310 | 707 | 4899 |
| Euteleostomi | 6 | 13 | 24128 | 274946 | 2690 | 15403 |
| saccharomyceta | 6 | 7 | 8524 | 47225 | 132 | 4195 |
| Arthropoda | 7 | 4 | 9682 | 43793 | 182 | 6115 |
| Neopterygii | 7 | 3 | 18605 | 57452 | 275 | 15418 |
| Saccharomycotina | 7 | 3 | 4482 | 13704 | 47 | 4295 |
| Tetrapoda | 7 | 10 | 20577 | 208759 | 1551 | 14252 |
| leotiomyceta | 7 | 4 | 8819 | 31301 | 130 | 7965 |
| Amniota | 8 | 9 | 21802 | 182082 | 840 | 13943 |
| Aparaglossata | 8 | 3 | 9298 | 31794 | 67 | 7713 |
| Clupeocephala | 8 | 2 | 16672 | 37353 | 77 | 16672 |
| sordariomyceta | 8 | 2 | 6428 | 13897 | 50 | 6063 |
| Diptera | 9 | 2 | 8407 | 20154 | 58 | 7787 |
| Theria | 9 | 8 | 22303 | 165042 | 229 | 16186 |
| Boreoeutheria | 10 | 7 | 22525 | 145412 | 120 | 18784 |
| Euarchontoglires | 11 | 5 | 22037 | 103710 | 109 | 17843 |
| Laurasiatheria | 11 | 2 | 17136 | 37819 | 27 | 17136 |
| Homininae | 12 | 3 | 20849 | 61108 | 23 | 20849 |
| Murinae | 12 | 2 | 18697 | 39902 | 92 | 18697 |
