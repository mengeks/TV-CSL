#!/usr/bin/env Rscript
# One Monte Carlo replication

`%||%` <- function(a, b) if (!is.null(a)) a else b

one_rep <- function(n = 1000, nuisance_type = c("oracle", "crossfit"),
                    n_folds = 2, seed = NULL) {
  nuisance_type <- match.arg(nuisance_type)
  if (!is.null(seed)) set.seed(seed)

  dat   <- simulate_data(n)
  basis <- make_basis(dat$X)

  nuis  <- if (nuisance_type == "oracle") {
    oracle_nuisances(dat)
  } else {
    estimate_nuisances_crossfit(dat, n_folds = n_folds)
  }

  res   <- fit_rso_coxph(dat, basis, nuis)
  beta  <- res$beta_hat

  svar  <- compute_sandwich(dat, basis, nuis, beta, res$precomp)

  se       <- sqrt(diag(svar$Omega) / n)
  covered  <- abs(beta - DGP_PARAMS$beta_true) < 1.96 * se

  list(
    beta_hat   = beta,
    se         = se,
    covered    = covered,
    event_rate = dat$event_rate
  )
}
