#!/usr/bin/env Rscript
# Diagnostic: oracle vs estimated risk-set-adjusted propensity.
# Q: Does using the TRUE τ₀(x) in the propensity give correct coverage?
# Run with n=200, B=10 for a fast check; B=50 for reliable inference.

suppressPackageStartupMessages(library(tidyverse))
source("R/old/TV-CSL.R")
source("R/old/data-handler.R")
source("R/datagen-helper.R")
source("scripts/TV-CSL/time-varying-estimate.R")

eta_type  <- "non-linear"
HTE_type  <- "linear"
n_sim     <- 200        # n=200 is fast enough for a quick check
K         <- 5
true_beta <- c(0, 1, 1, 1)

datagen_params <- list(
  is_time_varying        = TRUE,
  light_censoring        = FALSE,
  lambda_C               = 0.1,
  p                      = 3,
  X_distribution         = "normal",
  X_cov_type             = "identity",
  tx_difficulty          = "simple",
  linear_HTE_multiplier  = 1,
  linear_intercept       = 0,
  linear_slope_multiplier = 2.5,
  eta_type               = eta_type,
  HTE_type               = HTE_type
)

data_dir <- here::here("data", paste0(eta_type, "_", HTE_type))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

methods_TV_CSL <- list(
  prop_score_specs    = "PLACEHOLDER",
  regressor_specs     = "complex",
  lasso_types         = "S-lasso",
  final_model_methods = "lasso_coxph",
  HTE_specs           = "linear"
)

run_one <- function(iter, prop_spec) {
  set.seed(iter)
  generate_and_save_data(i=iter, n=n_sim,
                         path_for_sim_data=data_dir,
                         params=datagen_params, verbose=0)
  raw <- read_single_simulation_data(n=n_sim, i=iter,
                                     eta_type=eta_type, HTE_type=HTE_type)$data
  methods_TV_CSL$prop_score_specs <- prop_spec
  results <- tryCatch(
    run_TV_CSL_estimation(
      train_data_original  = raw,
      test_data            = raw,
      methods_TV_CSL       = methods_TV_CSL,
      i                    = iter,
      K                    = K,
      HTE_type             = HTE_type,
      eta_type             = eta_type,
      temp_result_csv_file = tempfile(fileext=".csv"),
      verbose              = 0
    ),
    error = function(e) { message("Error: ", e$message); NULL }
  )
  if (is.null(results) || length(results) == 0) return(NULL)
  ret <- results[[1]]
  se  <- ret$se_sandwich
  covered <- is.finite(se) &
    (ret$beta_HTE - 1.96*se <= true_beta) &
    (true_beta <= ret$beta_HTE + 1.96*se)
  data.frame(iter=iter, spec=prop_spec, coef=1:4,
             beta=ret$beta_HTE, se=se, covered=covered,
             bias=ret$beta_HTE - true_beta)
}

B <- 10
cat("Running", B, "reps per spec, n =", n_sim, ", eta = non-linear\n")
cat("Specs: estimated τ̂ (S-lasso complex) vs oracle τ₀(x)\n\n")

rows <- list()
for (b in seq_len(B)) {
  t0 <- proc.time()["elapsed"]
  for (spec in c("cox-risk-set-adjusted", "cox-risk-set-adjusted-oracle")) {
    r <- run_one(iter = 700 + b, prop_spec = spec)
    if (!is.null(r)) rows <- c(rows, list(r))
  }
  cat("  rep", b, "/", B, "  (", round(proc.time()["elapsed"] - t0, 1), "s)\n")
}

if (length(rows) == 0) stop("No results collected")
mc <- do.call(rbind, rows)

cat("\n=== Coverage (nominal 95%) and mean bias ===\n")
cat("n =", n_sim, ", eta = non-linear, β_true =", true_beta, "\n")
cat("(", B, "reps; note: n=200 has high variance so intervals are wide)\n\n")
summary_tbl <- mc %>%
  group_by(spec, coef) %>%
  summarise(coverage  = mean(covered,   na.rm=TRUE),
            mean_bias = mean(bias,      na.rm=TRUE),
            n_runs    = n(),
            .groups   = "drop")
print(as.data.frame(summary_tbl))
