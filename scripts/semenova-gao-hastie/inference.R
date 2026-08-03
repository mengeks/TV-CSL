#!/usr/bin/env Rscript
# Sandwich variance estimator for the risk-set orthogonalized Cox PL
#
# Uses pre-computed time-varying nuisances from fit_rso_coxph() output.
#
# At each event time t_k:
#   q_bar(t)    = sum_j Y_j(t) exp(theta_j) q_j / sum_j Y_j(t) exp(theta_j)
#   V_q(t)      = variance of q under the risk-set measure
#   dLambda0(t) = dN.(t) / sum_j Y_j(t) exp(theta_j)
#   dM_i(t)     = dN_i(t) - Y_i(t) exp(theta_i(t)) dLambda0(t)
#   psi_i       = sum_t (q_i(t) - q_bar(t)) dM_i(t)
#
#   Q_hat  = (1/n) sum_k V_q(t_k) * dN.(t_k)
#   Sigma  = (1/n) sum_i psi_i psi_i'
#   Omega  = Q_hat^{-1} Sigma Q_hat^{-1}
#   se(tau(x)) = sqrt( p(x)' Omega p(x) / n )

compute_sandwich <- function(dat, basis, nuis, beta_hat, precomp) {
  n   <- nrow(basis)
  d   <- ncol(basis)
  Y   <- dat$Y;  Del <- dat$Del;  W <- dat$W

  Q_sum   <- matrix(0, d, d)
  psi_mat <- matrix(0, n, d)

  for (k in seq_along(precomp$event_times)) {
    t       <- precomp$event_times[k]
    rs      <- precomp$nuis_at[[k]]
    in_risk <- Y >= t
    evt_now <- Del == 1 & Y == t
    idx     <- which(in_risk)

    q_mat   <- (W[in_risk] - rs$a_t[in_risk]) * basis[in_risk, , drop = FALSE]
    theta   <- rs$nu_t[in_risk] + as.numeric(q_mat %*% beta_hat)
    ew      <- exp(theta - max(theta))   # stabilized weights
    S0      <- sum(ew)
    if (S0 == 0) next

    q_bar   <- colSums(ew * q_mat) / S0
    q_c     <- sweep(q_mat, 2, q_bar, "-")
    V_q     <- crossprod(sqrt(ew) * q_c) / S0

    dN_t    <- sum(evt_now)
    Q_sum   <- Q_sum + V_q * dN_t

    # Breslow martingale increment: exp(-max) cancels in ew/S0, so use stabilized weights
    dM      <- as.integer(evt_now[in_risk]) - ew * dN_t / S0
    psi_mat[idx, ] <- psi_mat[idx, ] + sweep(q_c, 1, dM, "*")
  }

  Q_hat  <- Q_sum / n
  Sigma  <- crossprod(psi_mat) / n
  Q_inv  <- solve(Q_hat)
  Omega  <- Q_inv %*% Sigma %*% Q_inv

  list(Q_hat = Q_hat, Sigma = Sigma, Omega = Omega, psi_mat = psi_mat)
}

# Pointwise SE and 95% CI for tau(x) = p(x)' beta
pointwise_inference <- function(px, beta_hat, svar, n) {
  tau_hat <- as.numeric(px %*% beta_hat)
  se      <- sqrt(as.numeric(t(px) %*% svar$Omega %*% px) / n)
  list(
    tau_hat = tau_hat,
    se      = se,
    ci_lo   = tau_hat - 1.96 * se,
    ci_hi   = tau_hat + 1.96 * se
  )
}
