#!/usr/bin/env Rscript
# Risk-set orthogonalized Cox partial-likelihood estimator
#
# Solves the orthogonal score equation via Newton-Raphson directly,
# without expanding to counting-process format (avoids coxph interval issues).
#
# Score (Jacobian-based stabilized):
#   G(beta) = (1/n) sum_{k: events at t_k} [ q_i(t_k) - q_bar(t_k; beta) ]
#
# where:
#   q_j(t)    = (W_j - a_t(X_j)) p(X_j)
#   theta_j(t) = nu_t(X_j) + q_j(t)' beta
#   q_bar(t)  = sum_j Y_j(t) exp(theta_j) q_j / sum_j Y_j(t) exp(theta_j)

# Pre-compute time-varying nuisances at all event times.
# Calls R_fn_vec once per arm (2 calls total) instead of once per event time.
precompute_nuis_at <- function(dat, nuis) {
  event_times <- sort(unique(dat$Y[dat$Del == 1]))
  K           <- length(event_times)

  R0_mat <- nuis$R_fn_vec(event_times, 0)   # n x K
  R1_mat <- nuis$R_fn_vec(event_times, 1)   # n x K

  e        <- nuis$e
  eta0     <- nuis$eta0
  eta1     <- nuis$eta1
  exp_eta0 <- exp(eta0)
  exp_eta1 <- exp(eta1)

  nuis_at <- vector("list", K)
  for (k in seq_len(K)) {
    num  <- e * R1_mat[, k] * exp_eta1
    den  <- num + (1 - e) * R0_mat[, k] * exp_eta0
    a_t  <- pmin(pmax(num / den, 1e-6), 1 - 1e-6)
    nuis_at[[k]] <- list(
      a_t  = a_t,
      nu_t = a_t * eta1 + (1 - a_t) * eta0
    )
  }

  list(event_times = event_times, nuis_at = nuis_at)
}

# Score G(beta) — d-vector
compute_score <- function(beta, dat, basis, precomp) {
  Y   <- dat$Y;  Del <- dat$Del;  W <- dat$W
  d   <- ncol(basis)
  G   <- numeric(d)

  for (k in seq_along(precomp$event_times)) {
    t       <- precomp$event_times[k]
    rs      <- precomp$nuis_at[[k]]
    in_risk <- Y >= t
    evt_now <- Del == 1 & Y == t
    if (!any(evt_now)) next

    q_mat <- (W[in_risk] - rs$a_t[in_risk]) * basis[in_risk, , drop = FALSE]
    theta  <- rs$nu_t[in_risk] + as.numeric(q_mat %*% beta)
    ew     <- exp(theta - max(theta))    # log-sum-exp stabilization
    S0     <- sum(ew)
    q_bar  <- colSums(ew * q_mat) / S0

    n_evt <- sum(evt_now[in_risk])
    G     <- G + colSums(q_mat[evt_now[in_risk], , drop = FALSE]) - n_evt * q_bar
  }
  G / nrow(basis)
}

# Jacobian dG/dbeta — d x d matrix (negative information)
compute_jacobian <- function(beta, dat, basis, precomp) {
  Y   <- dat$Y;  Del <- dat$Del;  W <- dat$W
  d   <- ncol(basis)
  J   <- matrix(0, d, d)

  for (k in seq_along(precomp$event_times)) {
    t       <- precomp$event_times[k]
    rs      <- precomp$nuis_at[[k]]
    in_risk <- Y >= t
    evt_now <- Del == 1 & Y == t
    if (!any(evt_now)) next

    q_mat <- (W[in_risk] - rs$a_t[in_risk]) * basis[in_risk, , drop = FALSE]
    theta  <- rs$nu_t[in_risk] + as.numeric(q_mat %*% beta)
    ew     <- exp(theta - max(theta))
    S0     <- sum(ew)
    q_bar  <- colSums(ew * q_mat) / S0

    q_c <- sweep(q_mat, 2, q_bar, "-")
    V_q <- crossprod(sqrt(ew) * q_c) / S0

    n_evt <- sum(evt_now[in_risk])
    J     <- J - n_evt * V_q
  }
  J / nrow(basis)
}

# Newton-Raphson solver
fit_rso_coxph <- function(dat, basis, nuis, max_iter = 100, tol = 1e-10) {
  precomp <- precompute_nuis_at(dat, nuis)
  beta    <- rep(0, ncol(basis))

  for (iter in seq_len(max_iter)) {
    g    <- compute_score(beta, dat, basis, precomp)
    if (max(abs(g)) < tol) break
    J    <- compute_jacobian(beta, dat, basis, precomp)
    step <- tryCatch(solve(J, g), error = function(e) rep(NA_real_, ncol(basis)))
    if (anyNA(step)) break
    if (max(abs(step)) > 2) step <- step * 2 / max(abs(step))
    beta <- beta - step
  }

  list(beta_hat = beta, precomp = precomp)
}
