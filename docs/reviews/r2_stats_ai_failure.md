# Review R2: statistics and AI-agent failure modes (DESIGN.md v1)

**Clock.** `scripts/clock.sh` gives 15:17 (T+0h37), but DESIGN v1 is stamped 15:20 (T+0h40). That is a future stamp again. Use the script's output.

**Round-1 items addressed.**
- VGNC removed from dev.
- Peek budget.
- Leak test.
- S3 truth taken from RefOG trees.
- Primary export pre-registered.
- Shuffled controls.
- Runtime protocol.
- Gates.
- Manifests and CLAIMS.md.

The items below are still open or new.

## Must
1. **The 35/35 split leaks through homology.** I ran DIAMOND sensitive (e < 1e-5) on `refog_sequences/`. Several RefOGs are paralogous:
   - {005, 006}
   - {016, 026, 027, 032, 037, 039, 058}
   - {038, 052}
   - {015, 023}
   - {021, 068}

   If the merge is tuned to keep a dev RefOG apart from its test paralog, that tuning sets the test errors. **Split by homology component (60 units), stratified by size (3–125 genes).** Commit the components and the split hash. Pre-register which tree version S3 reads when several exist (v1/v2, `trees_tight`). Only dev trees may be read.

2. **35 test RefOGs is underpowered, and the success rule rewards that.** I simulated with the real RefOG sizes and benchmark.py's even weighting, with methods that differ only in split rate.
   - The paired-bootstrap SE of ΔF is about 0.04 at n = 35 and about 0.03 at n = 70.
   - A true ΔF of 0.03 is detected in fewer than 15% of runs.
   - Real methods have correlated errors, which lowers the SE, but probably only to about 0.02–0.03.

   The rule "no comparator dominates beyond the CI" is almost automatically met when CIs are wide. Replace it with **non-inferiority**: the lower 95% bound of Δ(ours − comparator) must be ≥ −δ. Pre-register δ = 0.03 for OrthoBench F and for SwissTree/TreeFam-A TPR, with PPV no worse than −0.02. Point-estimate tests, and tests with CI width > 2δ, are reported as "inconclusive", not "comparable". REPORT.md states the minimum detectable effect.

3. **Name the primary endpoints.** "Every QfO test" × 5 comparators means more than 40 checks. Re-adopt the R1 list and treat everything else as exploratory:
   - OrthoBench-test F;
   - SwissTree;
   - TreeFam-A;
   - G-STD2-Euk RF together with n;
   - EC with n.

4. **The T+20h gate is decided on test data and is undefined.** SwissTree and TreeFam-A are both test sets, and the dev data is Bilateria-only with no HGT.
   - Adopt R1's 30% TreeFam-A dev split by family-ID hash, plus a dev QfO clade that includes bacteria. Run the gate there.
   - Define "dominated" as follows, using a paired bootstrap over families: lower TPR at ≤ PPV, by more than 2 SE of Δ.
   - The T+12h gate needs OrthoFinder on OrthoBench-dev before T+12h.
   - State units: F as a fraction. `orthobench_eval.py` prints percentages.

5. **Per-RefOG counts.** `orthobench_eval.py` outputs only pooled values. Make it output per-RefOG TP/FP/FN, with a test that their sum reproduces benchmark.py's F exactly. Report macro-F as well.

## Should
6. **The layer keep rule is lenient.** "Killed if Δ ≤ 0" keeps noise. A layer enters the frozen config only if dev Δ(real − shuffled) > 0 **and** Δ(real − off) exceeds the seed spread. Before the first synteny run, fix the S5 ambiguous stratum (top-2 score ratio) and the WGD target set (teleost, yeast).

7. **The leak test is a string denylist, which agents route around.**
   - Add a runtime allowlist: `sys.addaudithook` in the pipeline entry point, failing on any `open()` outside `data/inputs/` and the run directory.
   - Extend the denylist: `controls`, `oracle`, `vgnc`, `ReconciledTrees`, `Data_for_RefOGs`, `participants_2022`, `reference_results`.
   - Competitor predictions never feed the method, including as "consensus" pseudo-labels.

8. **Simulator circularity.** A different agent from the method author owns the toy simulator and its truth. Zombi is mandatory for S3b, not "if installable".

9. **No fallback for the 0.005 reproduction check.** Pre-state what happens if it fails: either score all methods with the reimplemented tests, or report no comparison.

## Could
10. The final evaluation script refuses any config whose hash is not in `results/freeze.txt`. Peek outputs are logged with their config hash.

No further changes are needed elsewhere in v1.
