# Does TV-CSL's bias come from the first-stage outcome model eta_hat?
# Same datasets as the centering sweep at m = 2, n = 4000, but with the TRUE eta (tvcsl.true_eta).
# Compare with the sweep's estimated-eta results for the same reps.
# Usage (from project root): Rscript eta_rep.R <i> <prop_score_spec>
suppressMessages({ source("R/datagen-helper.R"); source("scripts/TV-CSL/time-varying-estimate.R") })
options(tvcsl.true_eta = TRUE)
D <- "scripts/TV-CSL/eta-check-2026-10-10"
a <- commandArgs(TRUE); i <- as.integer(a[1]); ps <- a[2]
m <- as.numeric(Sys.getenv("M", "2")); n <- as.integer(Sys.getenv("N", "4000"))
RESULTS_DIR <- paste0(tempfile("sidefiles"), "/")
gen <- function(seed) generate_simulated_data(n = n, eta_type = "non-linear", HTE_type = "linear",
                                              linear_HTE_multiplier = m, seed_value = seed)
train <- gen(1e6 + i); test <- gen(2e6 + i)          # same seeds as the centering sweep
set.seed(123 + 11 * i)
cfg <- list(prop_score_specs = ps, regressor_specs = "complex", lasso_types = "S-lasso",
            HTE_specs = "linear", final_model_methods = "lasso_coxph")
s0 <- Sys.time()
r <- run_TV_CSL_estimation(train, test, cfg, i = i, K = 5, HTE_type = "linear", eta_type = "non-linear",
                           temp_result_csv_file = tempfile())[[1]]
res <- data.frame(i = i, m = m, n = n, prop = ps, eta = "true", MSE = r$MSE, coef_idx = 1:4,
                  beta = r$beta_HTE, true = c(0, m, m, m), se_naive = r$se_naive,
                  se_sandwich = r$se_sandwich, fit_secs = as.numeric(difftime(Sys.time(), s0, units = "secs")))
write.csv(res, sprintf("%s/res/rep_n%d_m%s_%s_i%d.csv", D, n, format(m), ps, i), row.names = FALSE)
cat("DONE", i, ps, "\n")
