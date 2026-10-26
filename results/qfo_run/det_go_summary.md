# Run summary: `work/qfo2022_det/det_go`

| stage | wall s | CPU s (self+children) | peak RSS MB (own / max child) | load start→end |
|---|---|---|---|---|
| merge | 388.75 | 369.8 | 9440.4 / 7753.7 | 23.3→39.6 |
| output | 36.49 | 34.0 | 5622.8 / 7753.7 | 39.6→36.9 |

- **families_S0**: 163187 families (30300 multi-gene, 132887 singletons = 13.4% of genes); largest [7454, 7433, 7382, 7339, 7337, 7321, 7232, 7196, 7112, 7087]; Σg²/N² = 0.002282
- **families_after_rescue**: 162829 families (29942 multi-gene, 132887 singletons = 13.4% of genes); largest [7454, 7442, 7433, 7382, 7381, 7339, 7337, 7321, 7302, 7232]; Σg²/N² = 0.002828
- **root HOGs**: 36628 top-level groups (740445 genes; 248333 singleton genes); largest 6154 genes/57 sp @cellular organisms, 6026 genes/14 sp @Chordata, 4668 genes/49 sp @cellular organisms, 3824 genes/52 sp @cellular organisms, 3148 genes/35 sp @cellular organisms

| level | depth | #species | #HOGs (≥2 genes) | genes in HOGs | max HOG | orthologGroups with this TaxRange |
|---|---|---|---|---|---|---|
| cellular organisms | 0 | 78 | 36628 | 740445 | 6154 | 4621 |
| Archaea | 1 | 7 | 2244 | 11229 | 89 | 1233 |
| Bacteria | 1 | 23 | 9436 | 56376 | 105 | 4710 |
| Eukaryota | 1 | 48 | 32761 | 667772 | 6026 | 9000 |
| Bacillati | 2 | 7 | 4724 | 18227 | 54 | 2877 |
| Metamonada | 2 | 2 | 2724 | 13786 | 1058 | 1021 |
| Methanobacteriota | 2 | 4 | 1745 | 6205 | 68 | 1352 |
| Opisthokonta | 2 | 35 | 22596 | 483267 | 6026 | 5931 |
| Pseudomonadati | 2 | 13 | 6001 | 29153 | 94 | 4194 |
| Sar | 2 | 4 | 7048 | 38021 | 233 | 2716 |
| Thermoproteati | 2 | 3 | 1116 | 3568 | 41 | 660 |
| Thermotogati | 2 | 2 | 990 | 2476 | 18 | 683 |
| Viridiplantae | 2 | 5 | 13066 | 113272 | 1695 | 5006 |
| Actinomycetes | 3 | 2 | 2484 | 6592 | 36 | 1844 |
| Alveolata | 3 | 2 | 4566 | 17671 | 111 | 1582 |
| Cyanophyceae | 3 | 2 | 1955 | 4483 | 20 | 1695 |
| Embryophyta | 3 | 4 | 15683 | 101620 | 298 | 7750 |
| Eumetazoa | 3 | 22 | 21328 | 397522 | 6026 | 8639 |
| Fungi | 3 | 12 | 7965 | 75134 | 389 | 3427 |
| PVC group | 3 | 2 | 842 | 2451 | 31 | 424 |
| Pseudomonadota | 3 | 4 | 3837 | 12729 | 48 | 2676 |
| Stenosarchaea group | 3 | 2 | 1111 | 3018 | 32 | 824 |
| Stramenopiles | 3 | 2 | 3913 | 11202 | 65 | 3018 |
| Thermoproteota | 3 | 2 | 816 | 1912 | 29 | 696 |
| Bilateria | 4 | 21 | 21067 | 380181 | 6026 | 8119 |
| Dikarya | 4 | 11 | 8524 | 69593 | 234 | 4359 |
| Gammaproteobacteria | 4 | 2 | 2426 | 5719 | 25 | 1999 |
| Mesangiospermae | 4 | 3 | 18191 | 81046 | 170 | 11432 |
| Ascomycota | 5 | 8 | 7553 | 52174 | 229 | 3560 |
| Basidiomycota | 5 | 3 | 5132 | 15915 | 90 | 4662 |
| Chordata | 5 | 15 | 20163 | 304355 | 6026 | 10595 |
| Poaceae | 5 | 2 | 18323 | 44472 | 34 | 18323 |
| Protostomia | 5 | 6 | 10239 | 67929 | 727 | 5923 |
| Ecdysozoa | 6 | 5 | 9443 | 54398 | 724 | 4921 |
| Euteleostomi | 6 | 13 | 24165 | 274639 | 2230 | 15417 |
| saccharomyceta | 6 | 7 | 8462 | 47260 | 138 | 4195 |
| Arthropoda | 7 | 4 | 9717 | 43811 | 141 | 6144 |
| Neopterygii | 7 | 3 | 18604 | 57470 | 275 | 15427 |
| Saccharomycotina | 7 | 3 | 4487 | 13717 | 47 | 4302 |
| Tetrapoda | 7 | 10 | 20573 | 208626 | 1277 | 14258 |
| leotiomyceta | 7 | 4 | 8823 | 31321 | 138 | 7964 |
| Amniota | 8 | 9 | 21752 | 182104 | 1147 | 13925 |
| Aparaglossata | 8 | 3 | 9277 | 31791 | 60 | 7691 |
| Clupeocephala | 8 | 2 | 16666 | 37387 | 113 | 16666 |
| sordariomyceta | 8 | 2 | 6435 | 13873 | 51 | 6064 |
| Diptera | 9 | 2 | 8411 | 20246 | 56 | 7784 |
| Theria | 9 | 8 | 22236 | 165147 | 446 | 16195 |
| Boreoeutheria | 10 | 7 | 22454 | 145493 | 445 | 18762 |
| Euarchontoglires | 11 | 5 | 21993 | 103723 | 140 | 17842 |
| Laurasiatheria | 11 | 2 | 17160 | 37868 | 31 | 17160 |
| Homininae | 12 | 3 | 20845 | 61147 | 37 | 20845 |
| Murinae | 12 | 2 | 18679 | 39882 | 82 | 18679 |
