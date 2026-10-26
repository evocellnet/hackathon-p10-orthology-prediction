# Run summary: `work/qfo2022_det/det_nogo`

| stage | wall s | CPU s (self+children) | peak RSS MB (own / max child) | load start→end |
|---|---|---|---|---|
| input | 0.58 | 0.6 | 607.7 / 639.7 | 30.9→30.9 |
| families | 1533.68 | 10786.4 | 6411.2 / 6781.9 | 30.9→21.7 |
| scores | 3886.42 | 31978.8 | 29515.2 / 29515.2 | 21.7→41.2 |
| merge | 144.17 | 142.5 | 8417.1 / 7215.7 | 37.6→40.8 |
| output | 36.42 | 35.7 | 5466.4 / 7215.7 | 40.8→39.8 |

- **families_S0**: 163187 families (30300 multi-gene, 132887 singletons = 13.4% of genes); largest [7454, 7433, 7382, 7339, 7337, 7321, 7232, 7196, 7112, 7087]; Σg²/N² = 0.002282
- **families_after_rescue**: 162829 families (29942 multi-gene, 132887 singletons = 13.4% of genes); largest [7454, 7442, 7433, 7382, 7381, 7339, 7337, 7321, 7302, 7232]; Σg²/N² = 0.002828
- **root HOGs**: 36399 top-level groups (740516 genes; 248262 singleton genes); largest 6183 genes/57 sp @cellular organisms, 5889 genes/14 sp @Chordata, 4399 genes/49 sp @cellular organisms, 3832 genes/52 sp @cellular organisms, 3154 genes/35 sp @cellular organisms

| level | depth | #species | #HOGs (≥2 genes) | genes in HOGs | max HOG | orthologGroups with this TaxRange |
|---|---|---|---|---|---|---|
| cellular organisms | 0 | 78 | 36399 | 740516 | 6183 | 4602 |
| Archaea | 1 | 7 | 2243 | 11230 | 89 | 1231 |
| Bacteria | 1 | 23 | 9432 | 56393 | 105 | 4710 |
| Eukaryota | 1 | 48 | 32587 | 668027 | 5889 | 8998 |
| Bacillati | 2 | 7 | 4718 | 18230 | 54 | 2871 |
| Metamonada | 2 | 2 | 2724 | 13786 | 1058 | 1021 |
| Methanobacteriota | 2 | 4 | 1746 | 6206 | 68 | 1353 |
| Opisthokonta | 2 | 35 | 22494 | 483390 | 5889 | 5939 |
| Pseudomonadati | 2 | 13 | 5992 | 29144 | 94 | 4191 |
| Sar | 2 | 4 | 7052 | 37998 | 233 | 2714 |
| Thermoproteati | 2 | 3 | 1112 | 3568 | 41 | 660 |
| Thermotogati | 2 | 2 | 990 | 2477 | 18 | 683 |
| Viridiplantae | 2 | 5 | 13019 | 113361 | 1675 | 5003 |
| Actinomycetes | 3 | 2 | 2483 | 6594 | 36 | 1844 |
| Alveolata | 3 | 2 | 4566 | 17671 | 111 | 1582 |
| Cyanophyceae | 3 | 2 | 1955 | 4476 | 20 | 1695 |
| Embryophyta | 3 | 4 | 15530 | 101735 | 298 | 7713 |
| Eumetazoa | 3 | 22 | 21094 | 397788 | 5889 | 8640 |
| Fungi | 3 | 12 | 7960 | 75130 | 390 | 3427 |
| PVC group | 3 | 2 | 842 | 2451 | 31 | 424 |
| Pseudomonadota | 3 | 4 | 3832 | 12720 | 48 | 2675 |
| Stenosarchaea group | 3 | 2 | 1112 | 3022 | 32 | 825 |
| Stramenopiles | 3 | 2 | 3913 | 11202 | 65 | 3018 |
| Thermoproteota | 3 | 2 | 816 | 1912 | 29 | 696 |
| Bilateria | 4 | 21 | 20886 | 380354 | 5889 | 8095 |
| Dikarya | 4 | 11 | 8527 | 69591 | 234 | 4357 |
| Gammaproteobacteria | 4 | 2 | 2426 | 5716 | 25 | 1998 |
| Mesangiospermae | 4 | 3 | 18101 | 81080 | 170 | 11353 |
| Ascomycota | 5 | 8 | 7554 | 52175 | 229 | 3559 |
| Basidiomycota | 5 | 3 | 5132 | 15915 | 90 | 4662 |
| Chordata | 5 | 15 | 20067 | 304557 | 5889 | 10584 |
| Poaceae | 5 | 2 | 17901 | 43995 | 34 | 17901 |
| Protostomia | 5 | 6 | 10243 | 67937 | 724 | 5921 |
| Ecdysozoa | 6 | 5 | 9450 | 54406 | 723 | 4920 |
| Euteleostomi | 6 | 13 | 24030 | 274738 | 2282 | 15425 |
| saccharomyceta | 6 | 7 | 8462 | 47262 | 138 | 4195 |
| Arthropoda | 7 | 4 | 9707 | 43852 | 144 | 6144 |
| Neopterygii | 7 | 3 | 18317 | 56984 | 278 | 15278 |
| Saccharomycotina | 7 | 3 | 4487 | 13718 | 47 | 4302 |
| Tetrapoda | 7 | 10 | 20589 | 208751 | 1629 | 14232 |
| leotiomyceta | 7 | 4 | 8820 | 31321 | 138 | 7963 |
| Amniota | 8 | 9 | 21741 | 181850 | 1157 | 13692 |
| Aparaglossata | 8 | 3 | 9277 | 31778 | 60 | 7689 |
| Clupeocephala | 8 | 2 | 16292 | 36997 | 113 | 16292 |
| sordariomyceta | 8 | 2 | 6433 | 13874 | 51 | 6063 |
| Diptera | 9 | 2 | 8407 | 20236 | 56 | 7781 |
| Theria | 9 | 8 | 22162 | 165021 | 509 | 15998 |
| Boreoeutheria | 10 | 7 | 22387 | 145285 | 468 | 18682 |
| Euarchontoglires | 11 | 5 | 21907 | 103730 | 146 | 17763 |
| Laurasiatheria | 11 | 2 | 16908 | 37635 | 43 | 16908 |
| Homininae | 12 | 3 | 20769 | 61099 | 37 | 20769 |
| Murinae | 12 | 2 | 18568 | 39939 | 103 | 18568 |
