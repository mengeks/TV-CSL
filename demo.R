# Run from the repository root: Rscript demo.R
suppressPackageStartupMessages({
  library(mgcv)
  source("R/datagen-helper.R")
  source("scripts/TV-CSL/time-varying-estimate.R")
})

# 1. Simulate survival data with tau(x) = x1 + x2 + x3.
dat <- generate_simulated_data(
  n = 500, eta_type = "linear", HTE_type = "linear", seed_value = 123
)
tv <- create_pseudo_dataset(dat)

# 2. Fit TV-CSL with five-fold cross-fitting.
invisible(capture.output(fit <- TV_CSL(
  train_data = tv, test_data = dat, train_data_original = dat,
  HTE_type = "linear", eta_type = "linear", K = 5,
  prop_score_spec = "cox-time-varying-prop",
  lasso_type = "S-lasso", regressor_spec = "linear",
  final_model_method = "lasso_coxph", HTE_spec = "linear", verbose = 0
)))

# 3. Compare estimated log-hazard-ratio coefficients with the truth.
results <- data.frame(
  term = c("Intercept", "X.1", "X.2", "X.3"),
  truth = c(0, 1, 1, 1), estimate = round(unname(fit$beta_HTE), 3)
)
print(results, row.names = FALSE)
