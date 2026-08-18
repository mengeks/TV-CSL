#!/usr/bin/env Rscript
# Local timing and memory test before cluster submission.
#
# Run from project root:
#   cd /homes2/xmeng/TV-CSL
#   Rscript tests/timing-baby-example.R
#
# Tests all four methods (S-Cox, S-lasso, TV-CSL correct/misspec prop score)
# with n=200, K=2, eta=non-linear, HTE=linear, iter=1.
# Generates the two required datasets (train + test) on the fly if missing.

library(here)
library(jsonlite)

# ---- 0. Generate data if not already on disk --------------------------------
source(here::here("R/datagen-helper.R"))

n        <- 200
eta_type <- "non-linear"
HTE_type <- "linear"

datagen_params <- list(
  is_time_varying         = TRUE,
  light_censoring         = FALSE,
  lambda_C                = 0.1,
  p                       = 3,
  eta_type                = eta_type,
  HTE_type                = HTE_type,
  X_distribution          = "normal",
  X_cov_type              = "identity",
  tx_difficulty           = "simple",
  linear_HTE_multiplier   = 1,
  linear_intercept        = 0,
  linear_slope_multiplier = 2.5
)

data_dir <- here::here("data", paste0(eta_type, "_", HTE_type))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

for (idx in c(1, 101)) {   # iter 1 = train; iter 101 = held-out test (i + 100)
  fpath <- file.path(data_dir, paste0("sim_data_", idx, "_n_", n, ".rds"))
  if (!file.exists(fpath)) {
    cat(sprintf("Generating data: iter=%d, n=%d ...\n", idx, n))
    generate_and_save_data(
      i                 = idx,
      n                 = n,
      path_for_sim_data = data_dir,
      params            = datagen_params,
      verbose           = 0
    )
  } else {
    cat(sprintf("Data already exists: %s\n", basename(fpath)))
  }
}

# ---- 1. Boot the runner (sources time-varying-estimate + TV-CSL) ------------
source(here::here("scripts/TV-CSL/TV-CSL-runner.R"))

# ---- 2. Write a temporary JSON with K=2 for faster local testing ------------
config   <- fromJSON(here::here("scripts/TV-CSL/params-cluster.json"))
config$K <- 2
tmp_json <- tempfile(fileext = ".json")
write_json(config, tmp_json, auto_unbox = TRUE, pretty = TRUE)

# ---- 3. Run experiment and time it ------------------------------------------
cat("\n=================================================================\n")
cat(" TV-CSL Baby Example\n")
cat(sprintf(" n=%d | eta=%s | HTE=%s | iter=1 | K=2\n", n, eta_type, HTE_type))
cat("=================================================================\n\n")

t <- system.time({
  run_experiment_iteration(
    i         = 1,
    json_file  = tmp_json,
    eta_type   = eta_type,
    HTE_type   = HTE_type,
    n          = n,
    verbose    = 1
  )
})

cat("\n=================================================================\n")
cat(" Timing summary\n")
cat(sprintf(" Wall clock : %.1f s\n", t["elapsed"]))
cat(sprintf(" User time  : %.1f s\n", t["user.self"]))
cat("=================================================================\n")
cat("Memory (after gc):\n")
print(gc())
