#!/usr/bin/env Rscript
# Baby example for local timing and memory profiling before cluster submission.
#
# Run from the project root:
#   cd /homes2/xmeng/TV-CSL
#   Rscript scripts/TV-CSL/test_baby_example.R
#
# What it tests:
#   - All four methods: S-Cox, S-lasso, TV-CSL (correct prop), TV-CSL (intercept-only prop)
#   - eta_type = "non-linear"  (harder DGP)
#   - HTE_type = "linear"
#   - n = 200, iteration = 1, K = 2 folds (to keep runtime short)

source("scripts/TV-CSL/TV-CSL-runner.R")

json_file <- "scripts/TV-CSL/params-cluster.json"
n         <- 200
eta_type  <- "non-linear"
HTE_type  <- "linear"
i         <- 1

# Override K to 2 just for this local test (keep params-cluster.json at K=5 for cluster)
config       <- jsonlite::fromJSON(json_file)
config$K     <- 2
tmp_json     <- tempfile(fileext = ".json")
jsonlite::write_json(config, tmp_json, auto_unbox = TRUE, pretty = TRUE)

cat("=================================================================\n")
cat(" TV-CSL Baby Example\n")
cat(sprintf(" n=%d | eta=%s | HTE=%s | iter=%d | K=2\n", n, eta_type, HTE_type, i))
cat("=================================================================\n\n")

t <- system.time({
  run_experiment_iteration(
    i         = i,
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
cat("Memory in use (peak approximate):\n")
print(gc())
