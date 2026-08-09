#!/usr/bin/env Rscript
# One Monte Carlo replication

`%||%` <- function(a, b) if (!is.null(a)) a else b

one_rep <- function(n = 1000, nuisance_type = c("oracle", "crossfit"),
                    n_folds = 2, method = c("nr", "coxph"), seed = NULL) {
  nuisance_type <- match.arg(nuisance_type)
  method        <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  dat   <- simulate_data(n)
  basis <- make_basis(dat$X)

  nuis  <- if (nuisance_type == "oracle") {
    oracle_nuisances(dat)
  } else {
    estimate_nuisances_crossfit(dat, n_folds = n_folds)
  }

  res   <- fit_rso_coxph(dat, basis, nuis, method = method)
  beta  <- res$beta_hat

  # Always use compute_sandwich for a consistent SE across methods
  svar  <- compute_sandwich(dat, basis, nuis, beta, res$precomp)
  se    <- sqrt(diag(svar$Omega) / n)

  list(
    beta_hat   = beta,
    se         = se,
    covered    = abs(beta - DGP_PARAMS$beta_true) < 1.96 * se,
    event_rate = dat$event_rate
  )
}
