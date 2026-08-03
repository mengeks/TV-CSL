#!/usr/bin/env Rscript
# Monte Carlo simulation: Semenova-Gao-Hastie risk-set orthogonalized Cox PL
#
# Usage:
#   Rscript sim.R                            # defaults: B=200, n=1000, oracle
#   Rscript sim.R --n 2000 --B 100 --type crossfit --folds 2 --cores 4

suppressPackageStartupMessages(library(survival))

script_dir <- local({
  cli_args  <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", cli_args, value = TRUE)
  if (length(file_flag)) {
    dirname(normalizePath(sub("^--file=", "", file_flag[1])))
  } else {
    f <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
    if (!is.null(f)) dirname(f) else "."
  }
})
for (f in c("dgp.R", "nuisance.R", "estimator.R", "inference.R"))
  source(file.path(script_dir, f))

# ── CLI args ───────────────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) && length(args) >= i + 1L) return(args[i + 1L])
  default
}
n_obs         <- as.integer(get_arg("--n",     1000))
B             <- as.integer(get_arg("--B",      200))
nuisance_type <- get_arg("--type",  "oracle")
n_folds       <- as.integer(get_arg("--folds",    2))
n_cores       <- as.integer(get_arg("--cores",    1))

beta_true <- DGP_PARAMS$beta_true
d         <- length(beta_true)

cat(sprintf(
  "\n=== Semenova-Gao-Hastie Cox PL ===\n  n=%d  B=%d  nuisances=%s  folds=%d\n\n",
  n_obs, B, nuisance_type, n_folds
))

# ── one replication (importable by mclapply workers) ──────────────────────────
run_rep <- function(i) {
  for (f in c("dgp.R", "nuisance.R", "estimator.R", "inference.R"))
    source(file.path(script_dir, f), local = FALSE)

  set.seed(1000L + i)
  dat   <- simulate_data(n_obs)
  basis <- make_basis(dat$X)

  nuis  <- if (nuisance_type == "oracle") {
    oracle_nuisances(dat)
  } else {
    estimate_nuisances_crossfit(dat, n_folds = n_folds)
  }

  res  <- fit_rso_coxph(dat, basis, nuis)
  beta <- res$beta_hat
  svar <- compute_sandwich(dat, basis, nuis, beta, res$precomp)

  se      <- sqrt(diag(svar$Omega) / n_obs)
  covered <- abs(beta - beta_true) < 1.96 * se

  list(beta_hat = beta, se = se, covered = covered,
       event_rate = dat$event_rate, elapsed = proc.time()[["elapsed"]])
}


# ── run ────────────────────────────────────────────────────────────────────────
t_start <- proc.time()[["elapsed"]]

if (n_cores > 1L) {
  library(parallel)
  results <- mclapply(seq_len(B), run_rep, mc.cores = n_cores)
} else {
  results <- vector("list", B)
  for (i in seq_len(B)) {
    if (i == 1L || i %% 25L == 0L) cat(sprintf("  rep %3d / %d\n", i, B))
    results[[i]] <- run_rep(i)
  }
}

t_total <- proc.time()[["elapsed"]] - t_start

# ── summarise ─────────────────────────────────────────────────────────────────
beta_mat  <- do.call(rbind, lapply(results, `[[`, "beta_hat"))
se_mat    <- do.call(rbind, lapply(results, `[[`, "se"))
cover_mat <- do.call(rbind, lapply(results, `[[`, "covered"))
ev_rate   <- mean(sapply(results, `[[`, "event_rate"))

coef_names <- sprintf("beta_%d (true=%.1f)", seq_len(d) - 1L, beta_true)
cat(sprintf("\n--- Results [event rate: %.3f] ---\n", ev_rate))
cat(sprintf("%-18s  %8s  %8s  %8s  %8s\n",
            "Coef", "Bias", "EmpSD", "MeanSE", "Cover"))
cat(strrep("-", 58), "\n")
for (j in seq_len(d)) {
  cat(sprintf("%-18s  %+8.4f  %8.4f  %8.4f  %8.3f\n",
              coef_names[j],
              mean(beta_mat[, j]) - beta_true[j],
              sd(beta_mat[, j]),
              mean(se_mat[, j]),
              mean(cover_mat[, j])))
}
elapsed_v <- sapply(results, `[[`, "elapsed")
cat(sprintf(
  "\n--- Timing [%d reps, n=%d] ---\n  total %.1fs  |  mean %.2fs/rep  |  median %.2fs/rep\n",
  B, n_obs, t_total, mean(diff(c(t_start, elapsed_v))),
  median(diff(c(t_start, elapsed_v)))
))
cat("\nDone.\n")
