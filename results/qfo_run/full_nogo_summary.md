# Run summary: `work/qfo2022/full_nogo`

| stage | wall s | CPU s (self+children) | peak RSS MB (own / max child) | load start→end |
|---|---|---|---|---|
| input | 0.62 | 0.6 | 607.5 / 639.4 | 43.7→43.0 |
| families | 2136.7 | 11565.3 | 6201.9 / 6853.6 | 43.0→42.7 |
| scores | 5404.08 | 31804.9 | 29858.9 / 29858.9 | 42.7→39.4 |
| merge | 131.16 | 130.3 | 7051.6 / 6393.3 | 37.1→30.7 |
| output | 34.35 | 34.1 | 5218.5 / 6393.3 | 30.7→27.7 |

- **families_S0**: 163092 families (30267 multi-gene, 132825 singletons = 13.4% of genes); largest [7464, 7410, 7389, 7364, 7344, 7306, 7288, 7258, 7251, 7203]; Σg²/N² = 0.002299
- **families_after_rescue**: 162721 families (29896 multi-gene, 132825 singletons = 13.4% of genes); largest [7464, 7410, 7389, 7364, 7344, 7317, 7306, 7298, 7288, 7258]; Σg²/N² = 0.002902
- **root HOGs**: 36315 top-level groups (741173 genes; 247605 singleton genes); largest 7029 genes/15 sp @Bilateria, 5551 genes/57 sp @cellular organisms, 4799 genes/53 sp @cellular organisms, 3924 genes/49 sp @cellular organisms, 3914 genes/51 sp @cellular organisms

| level | depth | #species | #HOGs (≥2 genes) | genes in HOGs | max HOG | orthologGroups with this TaxRange |
|---|---|---|---|---|---|---|
| cellular organisms | 0 | 78 | 36315 | 741173 | 7029 | 4647 |
| Archaea | 1 | 7 | 2257 | 11229 | 89 | 1229 |
| Bacteria | 1 | 23 | 9391 | 56335 | 102 | 4703 |
| Eukaryota | 1 | 48 | 32685 | 668294 | 7029 | 8944 |
| Bacillati | 2 | 7 | 4691 | 18200 | 54 | 2876 |
| Metamonada | 2 | 2 | 2755 | 13923 | 1058 | 1019 |
| Methanobacteriota | 2 | 4 | 1738 | 6203 | 61 | 1345 |
| Opisthokonta | 2 | 35 | 22663 | 483710 | 7029 | 5901 |
| Pseudomonadati | 2 | 13 | 5989 | 29091 | 88 | 4201 |
| Sar | 2 | 4 | 6965 | 38036 | 229 | 2727 |
| Thermoproteati | 2 | 3 | 1101 | 3559 | 41 | 661 |
| Thermotogati | 2 | 2 | 989 | 2483 | 18 | 684 |
| Viridiplantae | 2 | 5 | 13029 | 113161 | 1696 | 4995 |
| Actinomycetes | 3 | 2 | 2481 | 6557 | 36 | 1844 |
| Alveolata | 3 | 2 | 4545 | 17845 | 229 | 1593 |
| Cyanophyceae | 3 | 2 | 1951 | 4477 | 21 | 1687 |
| Embryophyta | 3 | 4 | 15560 | 101576 | 296 | 7654 |
| Eumetazoa | 3 | 22 | 21136 | 398196 | 7029 | 8569 |
| Fungi | 3 | 12 | 8000 | 75180 | 382 | 3428 |
| PVC group | 3 | 2 | 831 | 2424 | 31 | 420 |
| Pseudomonadota | 3 | 4 | 3846 | 12721 | 42 | 2680 |
| Stenosarchaea group | 3 | 2 | 1113 | 3010 | 32 | 826 |
| Stramenopiles | 3 | 2 | 3848 | 11191 | 65 | 2987 |
| Thermoproteota | 3 | 2 | 800 | 1875 | 29 | 684 |
| Bilateria | 4 | 21 | 21035 | 380783 | 7029 | 8086 |
| Dikarya | 4 | 11 | 8562 | 69589 | 245 | 4345 |
| Gammaproteobacteria | 4 | 2 | 2420 | 5685 | 25 | 1994 |
| Mesangiospermae | 4 | 3 | 18052 | 81032 | 170 | 11350 |
| Ascomycota | 5 | 8 | 7618 | 52120 | 243 | 3551 |
| Basidiomycota | 5 | 3 | 5142 | 15917 | 90 | 4660 |
| Chordata | 5 | 15 | 19932 | 305039 | 5890 | 10577 |
| Poaceae | 5 | 2 | 17882 | 43916 | 63 | 17882 |
| Protostomia | 5 | 6 | 10306 | 68175 | 708 | 5896 |
| Ecdysozoa | 6 | 5 | 9467 | 54316 | 707 | 4900 |
| Euteleostomi | 6 | 13 | 23961 | 275013 | 2282 | 15429 |
| saccharomyceta | 6 | 7 | 8527 | 47227 | 132 | 4195 |
| Arthropoda | 7 | 4 | 9671 | 43832 | 182 | 6116 |
| Neopterygii | 7 | 3 | 18325 | 56971 | 278 | 15273 |
| Saccharomycotina | 7 | 3 | 4482 | 13705 | 47 | 4295 |
| Tetrapoda | 7 | 10 | 20589 | 208837 | 1629 | 14224 |
| leotiomyceta | 7 | 4 | 8816 | 31302 | 130 | 7964 |
| Amniota | 8 | 9 | 21730 | 181767 | 849 | 13715 |
| Aparaglossata | 8 | 3 | 9298 | 31773 | 67 | 7710 |
| Clupeocephala | 8 | 2 | 16301 | 36960 | 113 | 16301 |
| sordariomyceta | 8 | 2 | 6426 | 13898 | 50 | 6062 |
| Diptera | 9 | 2 | 8403 | 20139 | 58 | 7785 |
| Theria | 9 | 8 | 22159 | 164908 | 509 | 15969 |
| Boreoeutheria | 10 | 7 | 22402 | 145132 | 185 | 18680 |
| Euarchontoglires | 11 | 5 | 21892 | 103696 | 185 | 17759 |
| Laurasiatheria | 11 | 2 | 16881 | 37558 | 36 | 16881 |
| Homininae | 12 | 3 | 20770 | 61054 | 33 | 20770 |
| Murinae | 12 | 2 | 18567 | 39890 | 92 | 18567 |
