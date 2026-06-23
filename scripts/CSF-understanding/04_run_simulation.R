# =============================================================
# 04_run_simulation.R
# Main simulation: CSF vs baselines (100 MC iterations)
#
# Parallelises over MC iterations via foreach + doParallel.
# Each iteration is saved individually as results/mc_iter_NNN.rds
# so progress is preserved even if the job is interrupted.
# =============================================================

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(parallel)
  library(grf)
  library(randomForestSRC)
  library(survival)
})

# ---- Working directory and source utilities ----
# Works both with Rscript and source()
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.file_arg) && nchar(.file_arg)) {
  normalizePath(dirname(.file_arg), mustWork = FALSE)
} else {
  tryCatch(normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
           error = function(e) getwd())
}
setwd(SCRIPT_DIR)  # ensure relative paths work

source(file.path(SCRIPT_DIR, "01_utils_dgp.R"))
source(file.path(SCRIPT_DIR, "02_utils_oracle.R"))
source(file.path(SCRIPT_DIR, "03_utils_methods.R"))

# ---- DGP choice: "dgp1" (Weibull log-HR) or "dgp2" (linear RMST-CATE) ----
DGP_CHOICE <- "dgp2"   # <-- change this line to switch DGPs

# ---- Simulation parameters ----
SIM <- list(
  n_train        = 1000,
  n_test         = 500,
  n_mc           = 100,
  h              = DGP_PARAMS$h,
  num_trees_csf  = 1000,
  num_trees_rsf  = 500,
  seed_base      = 2024
)

RESULTS_DIR <- file.path(SCRIPT_DIR, paste0("results_", DGP_CHOICE))
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Single MC iteration ----
run_one_mc <- function(mc, sim = SIM, dgp = DGP_PARAMS, dgp_choice = DGP_CHOICE,
                       results_dir = RESULTS_DIR) {
  out_file <- file.path(results_dir, sprintf("mc_iter_%03d.rds", mc))

  # Skip if already done (allows restarts)
  if (file.exists(out_file)) {
    message(sprintf("[MC %03d] already done, skipping.", mc))
    return(readRDS(out_file))
  }

  set.seed(sim$seed_base + mc)
  t0 <- proc.time()["elapsed"]
  msg <- function(...) message(sprintf("[MC %03d] ", mc), ..., appendLF = TRUE)

  msg("Generating data (", dgp_choice, ")")
  if (dgp_choice == "dgp2") {
    d_tr <- simulate_data_dgp2(sim$n_train, dgp)
    d_te <- simulate_data_dgp2(sim$n_test,  dgp)
  } else {
    d_tr <- simulate_data(sim$n_train, dgp)
    d_te <- simulate_data(sim$n_test,  dgp)
  }

  X_tr <- d_tr$X;  U_tr <- d_tr$U
  A_tr <- d_tr$A;  D_tr <- d_tr$Delta
  X_te <- d_te$X

  msg("Computing true tau on test set")
  tau_true <- if (dgp_choice == "dgp2") {
    compute_true_tau_dgp2(X_te, dgp)
  } else {
    compute_true_tau(X_te, sim$h, dgp)
  }

  res <- list(mc = mc, tau_true = tau_true, timing = list())

  run_method <- function(name, expr) {
    t1 <- proc.time()["elapsed"]
    msg("Fitting ", name)
    tryCatch({
      fit <- eval(expr)
      res[[paste0(name, "_hat")]]     <<- fit$tau_hat
      res[[name]]                     <<- evaluate_metrics(fit$tau_hat, tau_true)
    }, error = function(e) {
      message(sprintf("[MC %03d] %s FAILED: %s", mc, name, conditionMessage(e)))
      res[[paste0(name, "_hat")]] <<- rep(NA_real_, sim$n_test)
      res[[name]]                 <<- setNames(rep(NA_real_, 5),
                                               c("MSE","RMSE","Bias","Cor","SignError"))
    })
    res$timing[[name]] <<- round(proc.time()["elapsed"] - t1, 1)
    invisible(NULL)
  }

  run_method("csf",
    quote(fit_csf(X_tr, U_tr, D_tr, A_tr, X_te, sim$h, sim$num_trees_csf, num.threads = 1L)))

  run_method("tlearner",
    quote(fit_tlearner(X_tr, U_tr, D_tr, A_tr, X_te, sim$h, sim$num_trees_rsf, ncore = 1L)))

  run_method("src1",
    quote(fit_slearner(X_tr, U_tr, D_tr, A_tr, X_te, sim$h, sim$num_trees_rsf, ncore = 1L)))

  run_method("src2",
    quote(fit_src2(X_tr, U_tr, D_tr, A_tr, X_te, sim$h, sim$num_trees_rsf, ncore = 1L)))

  run_method("cox",
    quote(fit_cox_interaction(X_tr, U_tr, D_tr, A_tr, X_te, sim$h)))

  run_method("ipcw_cf",
    quote(fit_ipcw_cf(X_tr, U_tr, D_tr, A_tr, X_te, sim$h, sim$num_trees_csf, num.threads = 1L)))

  # Oracle linear CSF (DGP-specific)
  {
    t1_olcsf <- proc.time()["elapsed"]
    msg("Fitting oracle_lcsf")
    tryCatch({
      fit_olcsf <- if (dgp_choice == "dgp2") {
        fit_linear_csf_oracle_dgp2(X_tr, U_tr, D_tr, A_tr, X_te,
                                   d_tr$r0, d_tr$r1, sim$h, dgp)
      } else {
        fit_linear_csf_oracle(X_tr, U_tr, D_tr, A_tr, X_te, sim$h, dgp)
      }
      res$oracle_lcsf_hat  <- fit_olcsf$tau_hat
      res$oracle_lcsf      <- evaluate_metrics(fit_olcsf$tau_hat, tau_true)
      res$oracle_lcsf_beta <- fit_olcsf$beta
    }, error = function(e) {
      message(sprintf("[MC %03d] oracle_lcsf FAILED: %s", mc, conditionMessage(e)))
      res$oracle_lcsf_hat <<- rep(NA_real_, sim$n_test)
      res$oracle_lcsf     <<- setNames(rep(NA_real_, 5),
                                       c("MSE","RMSE","Bias","Cor","SignError"))
    })
    res$timing$oracle_lcsf <- round(proc.time()["elapsed"] - t1_olcsf, 1)
  }

  res$total_time <- round(proc.time()["elapsed"] - t0, 1)
  msg(sprintf("Done (%.0fs)", res$total_time))

  saveRDS(res, out_file)
  res
}

# ---- Parallel setup ----
# Cap at 24 workers (1 thread per worker = no oversubscription on 64-core machine)
n_cores <- min(24L, max(1L, detectCores(logical = FALSE) - 2L), SIM$n_mc)
message(sprintf("Launching %d parallel workers for %d MC iterations", n_cores, SIM$n_mc))

cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

# Source utilities on each worker; disable OpenMP sub-threads
clusterExport(cl, c("SCRIPT_DIR", "RESULTS_DIR", "SIM", "DGP_PARAMS", "DGP2_PARAMS", "DGP_CHOICE"))
invisible(clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  suppressPackageStartupMessages({
    library(grf)
    library(randomForestSRC)
    library(survival)
  })
  source(file.path(SCRIPT_DIR, "01_utils_dgp.R"))
  source(file.path(SCRIPT_DIR, "02_utils_oracle.R"))
  source(file.path(SCRIPT_DIR, "03_utils_methods.R"))
  invisible(NULL)
}))

# ---- Run simulation ----
t_start <- proc.time()["elapsed"]

all_results <- foreach(
  mc      = seq_len(SIM$n_mc),
  .errorhandling = "pass",
  .verbose = FALSE
) %dopar% {
  dgp <- if (DGP_CHOICE == "dgp2") DGP2_PARAMS else DGP_PARAMS
  run_one_mc(mc, SIM, dgp, DGP_CHOICE, RESULTS_DIR)
}

stopCluster(cl)

elapsed <- round(proc.time()["elapsed"] - t_start, 0)
message(sprintf("All %d MC iterations completed in %.0f minutes.",
                SIM$n_mc, elapsed / 60))

# Save combined results
saveRDS(all_results, file.path(RESULTS_DIR, "all_results.rds"))
message("Saved: ", file.path(RESULTS_DIR, "all_results.rds"))
