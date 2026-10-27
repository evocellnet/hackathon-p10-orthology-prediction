# Protocol deviations and disclosures (must appear in REPORT.md)
Collected at M3 (T+14h15) by an independent retrospective agent; maintained by the PI.

1. **The TreeFam-A dev split was created after QfO peek 1** (peek 1 cafe3f5 → split f65d5b7). Peek 1 had already shown aggregate TreeFam-A and SwissTree recall, and that motivated the recall work. So dev/test is not blind in the DESIGN §4 sense. Test-only TreeFam-A families were never scored separately before the freeze (verify at freeze: results/test_log.tsv).
2. **Heavy dev reuse:** OrthoBench dev has 243 runs on 34 RefOGs, and synteny parameters, local trees, FastTree, the presets and deep absorption were all selected on it. Dev numbers (79.7 vs OrthoFinder 70.5) are optimistic; the comparator was not tuned. Claims rely on OrthoBench **test** only.
3. **Deep-node threshold** (median seed −ln S ≥ 1.5) was defined and tuned on QfO TreeFam-A dev. No OrthoBench node qualifies, so it is QfO-specific and untested on OrthoBench.
4. **The S0 completeness target (≥95%) was missed:** dev 0.794, test 0.722, against a sequence-similarity ceiling ≈0.82. Recorded as a deviation; test completeness was logged.
5. **Non-inferiority on recall failed at peek 1** (TreeFam-A TPR .55 vs OrthoFinder .74); the final result must be reported against the pre-registered rule as it stands.
6. **Comparator rows** are published OpenEBench values, produced under other conditions. FastOMA benefits from the OMA database. Our local FastOMA OrthoBench run (dev 29.6) is likely misconfigured and is not used.
7. **Design review round 3 was not run** (deferred, then superseded by data-driven reviews). The freeze moved from T+36h (DESIGN §7) to T+20h, to leave time for evaluation and review.
8. **Post-hoc redesigns.** Constrained synteny replaced the additive bonus after it failed. The micro → macro pair-F1 switch was made on dev. Shuffled-control spans are ±1–3 OG-F, so OG-level synteny gains of +4–7 are moderate evidence. VGNC PPV .982 suggests vertebrate over-merging.
9. **No claim was independently verified before M3.** Independent verification of CLAIMS.md is scheduled after the freeze. G-STD2/STD tests are run once, on the frozen config.
10. **Adoption rules were first checked only on the fast preset.** so_struct passed on fast and failed on accurate; the presets now differ in options (each adopted on its own base).
