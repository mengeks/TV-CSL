# One rep of the centering sweep: TV-CSL marginal vs time-varying centering (+ oracle),
# non-linear eta_0, complex eta_0 spec, beta0 = (0, m, m, m).
# Usage (from project root): Rscript sweep_rep.R <i> <m> <n> [specs, comma-separated]
suppressMessages({ source("R/datagen-helper.R"); source("scripts/TV-CSL/time-varying-estimate.R") })
D <- "scripts/TV-CSL/centering-sweep-2026-10-10"
a <- commandArgs(TRUE); i <- as.integer(a[1]); m <- as.numeric(a[2]); n <- as.integer(a[3])
specs <- if (length(a) >= 4) strsplit(a[4], ",")[[1]] else
  c("cox-linear-censored-only", "cox-time-varying-prop", "cox-time-varying-oracle")
RESULTS_DIR <- paste0(tempfile("sidefiles"), "/")          # per-process; side CSVs are not used
beta0 <- c(0, m, m, m)
# Common random numbers: rep i uses the same seeds for every m.
gen <- function(seed) generate_simulated_data(n = n, eta_type = "non-linear", HTE_type = "linear",
                                              linear_HTE_multiplier = m, seed_value = seed)
t0 <- Sys.time(); train <- gen(1e6 + i); test <- gen(2e6 + i)
gen_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
out <- list()
for (ps in specs) {
  set.seed(123 + 11 * i)
  cfg <- list(prop_score_specs = ps, regressor_specs = "complex", lasso_types = "S-lasso",
              HTE_specs = "linear", final_model_methods = "lasso_coxph")
  s0 <- Sys.time()
  r <- run_TV_CSL_estimation(train, test, cfg, i = i, K = 5, HTE_type = "linear", eta_type = "non-linear",
                             temp_result_csv_file = tempfile())[[1]]
  out[[ps]] <- data.frame(i = i, m = m, n = n, prop = ps, MSE = r$MSE, coef_idx = 1:4, beta = r$beta_HTE,
                          true = beta0, se_naive = r$se_naive, se_sandwich = r$se_sandwich,
                          fit_secs = as.numeric(difftime(Sys.time(), s0, units = "secs")),
                          gen_secs = gen_secs, n_events = sum(train$Delta))
}
f <- sprintf("%s/res/rep_n%d_m%s_i%d.csv", D, n, format(m), i)
write.csv(do.call(rbind, out), f, row.names = FALSE)
cat("DONE", f, "fit secs:", round(sapply(out, function(o) o$fit_secs[1])), "\n")
