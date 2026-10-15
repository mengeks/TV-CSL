# True-eta check (m = 2, n = 4000, non-linear eta_0, complex spec, 20 reps; same datasets as the centering sweep)
Estimated-eta rows come from the centering sweep; true-eta rows use options(tvcsl.true_eta = TRUE).

| centering | eta | MSE | slope bias (SE) | coverage |
|---|---|---|---|---|
| marginal | estimated | 0.063 | −0.102 (0.012) | 0.57 |
| marginal | true | 0.017 | +0.014 (0.009) | 0.93 |
| time-varying (GAM) | estimated | 0.076 | −0.129 (0.011) | 0.57 |
| time-varying (GAM) | true | 0.021 | +0.004 (0.010) | 0.95 |
| oracle | estimated | 0.031 | −0.007 (0.012) | 0.90 |
| oracle | true | 0.023 | −0.001 (0.011) | 0.93 |

Conclusion: the bias comes from eta_hat (S-lasso, complex spec; eta_0 error about 65% of its SD).
With the true eta, every centering is unbiased with coverage 0.93–0.95. With estimated eta, only the exact a0 protects;
the GAM, though 3x closer to a0 than marginal, does not.
