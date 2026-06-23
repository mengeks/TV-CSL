# =============================================================
# 04b_repair_oracle_lcsf.R
# Re-runs oracle_lcsf on completed MC iterations where it
# returned NA (due to a scoping bug now fixed).
# Uses the same seeds so data is identical.
# =============================================================

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(parallel)
})

.args <- commandArgs(trailingOnly = FALSE)
.fa   <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.fa) && nchar(.fa)) {
  normalizePath(dirname(.fa), mustWork = FALSE)
} else {
  tryCatch(normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
           error = function(e) getwd())
}
setwd(SCRIPT_DIR)

source("01_utils_dgp.R")
source("02_utils_oracle.R")
source("03_utils_methods.R")

SIM <- list(
  n_train   = 1000,
  n_test    = 500,
  h         = DGP_PARAMS$h,
  seed_base = 2024
)
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")

# Find result files where oracle_lcsf is NA
rds_files <- sort(list.files(RESULTS_DIR, "mc_iter_.*\\.rds$", full.names = TRUE))
needs_repair <- vapply(rds_files, function(f) {
  r <- readRDS(f)
  is.null(r$oracle_lcsf) || all(is.na(r$oracle_lcsf))
}, logical(1))

cat(sprintf("%d of %d files need oracle_lcsf repair\n", sum(needs_repair), length(rds_files)))
repair_files <- rds_files[needs_repair]

if (length(repair_files) == 0) {
  cat("Nothing to repair. Exiting.\n")
  quit(status = 0)
}

repair_one <- function(fpath) {
  r  <- readRDS(fpath)
  mc <- r$mc

  set.seed(SIM$seed_base + mc)
  d_tr <- simulate_data(SIM$n_train, DGP_PARAMS)
  d_te <- simulate_data(SIM$n_test,  DGP_PARAMS)

  X_tr <- d_tr$X; U_tr <- d_tr$U; D_tr <- d_tr$Delta; A_tr <- d_tr$A
  X_te <- d_te$X
  tau_true <- r$tau_true   # use saved ground truth

  t1 <- proc.time()["elapsed"]
  tryCatch({
    fit <- fit_linear_csf_oracle(X_tr, U_tr, D_tr, A_tr, X_te, SIM$h, DGP_PARAMS)
    r$oracle_lcsf_hat  <- fit$tau_hat
    r$oracle_lcsf      <- evaluate_metrics(fit$tau_hat, tau_true)
    r$oracle_lcsf_beta <- fit$beta
    cat(sprintf("[MC %03d] oracle_lcsf repaired: MSE=%.5f Cor=%.3f  beta=[%.4f %.4f]\n",
        mc, r$oracle_lcsf["MSE"], r$oracle_lcsf["Cor"],
        fit$beta[1], fit$beta[2]))
  }, error = function(e) {
    cat(sprintf("[MC %03d] oracle_lcsf ERROR: %s\n", mc, conditionMessage(e)))
  })
  r$timing$oracle_lcsf <- round(proc.time()["elapsed"] - t1, 1)

  saveRDS(r, fpath)
  invisible(r)
}

# Run in parallel
n_cores <- min(24L, max(1L, detectCores(logical = FALSE) - 2L), length(repair_files))
cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

clusterExport(cl, c("SCRIPT_DIR", "RESULTS_DIR", "SIM", "DGP_PARAMS"))
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(grf); library(randomForestSRC); library(survival)})
  source(file.path(SCRIPT_DIR, "01_utils_dgp.R"))
  source(file.path(SCRIPT_DIR, "02_utils_oracle.R"))
  source(file.path(SCRIPT_DIR, "03_utils_methods.R"))
  invisible(NULL)
}))

t_start <- proc.time()["elapsed"]
foreach(fpath = repair_files, .errorhandling = "pass") %dopar% {
  repair_one(fpath)
}
stopCluster(cl)

elapsed <- round(proc.time()["elapsed"] - t_start, 0)
cat(sprintf("Repair complete in %.0f seconds.\n", elapsed))

# Re-aggregate
all_res <- lapply(sort(list.files(RESULTS_DIR, "mc_iter_.*\\.rds$", full.names = TRUE)), readRDS)
saveRDS(all_res, file.path(RESULTS_DIR, "all_results.rds"))
cat("Updated all_results.rds\n")
