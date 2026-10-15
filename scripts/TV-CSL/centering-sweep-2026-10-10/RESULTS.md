# Centering sweep — results (2026-10-10)
Non-linear eta_0, complex spec, n = 4000, 50 reps per m, beta0 = (0, m, m, m), common random numbers across m.
Fixed DGP (branch fix-dgp-exact-event-times @ 42c0fb9). Reproduce: `Rscript scripts/TV-CSL/centering-sweep-2026-10-10/summarise_sweep.R`.
Slopes = beta_2..beta_4 (true m); coverage uses the naive SE (the sandwich SE is within 0.01 of it throughout).

| m | MSE marg / tv / oracle | slope bias marg / tv / oracle | coverage marg / tv / oracle |
|---|---|---|---|
| 0   | 0.017 / 0.017 / 0.017 | −0.006 / −0.008 / −0.006 | 0.90 / 0.91 / 0.90 |
| 0.5 | 0.020 / 0.021 / 0.020 | −0.021 / −0.023 / −0.003 | 0.90 / 0.87 / 0.89 |
| 1   | 0.025 / 0.027 / 0.024 | −0.011 / −0.037 / +0.003 | 0.91 / 0.88 / 0.90 |
| 1.5 | 0.032 / 0.041 / 0.028 | −0.034 / −0.072 / +0.002 | 0.84 / 0.78 / 0.92 |
| 2   | 0.069 / 0.073 / 0.035 | −0.107 / −0.120 / −0.002 | 0.57 / 0.65 / 0.91 |

Paired MSE difference vs marginal (z):
- time-varying: −0.2, +0.8, +1.1, +2.9, +0.7. Never better; worse at m = 1.5.
- oracle: +0.3, −0.9, −0.4, −1.7, −6.6.

Conclusions:
1. The theory holds with the true time-varying propensity. Oracle centering is unbiased at every m, with coverage 0.89–0.92. The marginal (Gao–Hastie) centering bias grows with m (−0.006 → −0.107) and its coverage falls from 0.90 to 0.57. At m = 2 the oracle halves the MSE.
2. The estimated time-varying centering (`cox-time-varying-prop`: GAM `W ~ s(tstop) + X1 + X2 + X3` on event rows) does not realise this gain. It is as biased as marginal, or more. The true logit, log r + tau + r(x) t − d(x) t² + log I(t, x), has strong time-by-covariate interactions that an additive model cannot fit.
3. The null point m = 0 ties exactly, as predicted.
4. Coverage is about 0.90, not 0.95, even for the oracle and at m = 0: the empirical SD is about 11% above the mean SE (0.0675 vs 0.0607 at m = 0). This is a separate SE under-estimate that affects all centerings.
