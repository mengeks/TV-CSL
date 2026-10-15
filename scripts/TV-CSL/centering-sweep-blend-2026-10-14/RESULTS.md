# Centering sweep at eta_0 = blend-0.25 — results (2026-10-14)
eta_0 = 0.25 * hard + 0.75 * smooth (complex-lasso eta_0 error ~28% of SD), complex spec, n = 4000, ~49 reps per m
(missing: m0_i37, m1_i36, m1.5_i36, m2_i36). Reproduce: `Rscript scripts/TV-CSL/centering-sweep-blend-2026-10-14/summarise_sweep.R`.

| m | slope bias marg / tv / oracle | coverage (naive) marg / tv / oracle | MSE marg / tv / oracle |
|---|---|---|---|
| 0   | −0.004 / −0.006 / −0.004 | 0.88 / 0.88 / 0.88 | 0.013 / 0.013 / 0.013 |
| 0.5 | −0.007 / −0.009 / −0.003 | 0.93 / 0.91 / 0.93 | 0.013 / 0.014 / 0.013 |
| 1   | −0.018 / −0.012 / +0.001 | 0.91 / 0.91 / 0.93 | 0.015 / 0.016 / 0.015 |
| 1.5 | −0.033 / −0.022 / −0.000 | 0.89 / 0.91 / 0.93 | 0.018 / 0.019 / 0.018 |
| 2   | −0.048 / −0.031 / −0.002 | 0.83 / 0.90 / 0.91 | 0.025 / 0.024 / 0.022 |

Paired slope-bias reduction (tv − marginal, both biases negative): m = 1: +0.006 (z 2.7), m = 1.5: +0.011 (z 4.2), m = 2: +0.017 (z 4.9).
Paired MSE difference: none significant at m ≥ 1 (m = 2: −0.0009, z −0.5).

Conclusions: time-varying centering cuts the Gao–Hastie (marginal) bias by about 35%, significantly, and the gap grows with m as the theory predicts.
It keeps coverage at about 0.90 at m = 2, where marginal drops to 0.83. There is no MSE gain, because its variance is larger
(emp SD 0.078 vs 0.071 at m = 2). The null point m = 0 ties. Coverage is still about 0.88–0.93, not 0.95, for every centering.
