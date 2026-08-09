#!/usr/bin/env Rscript
# Risk-set orthogonalized Cox partial-likelihood estimator
#
# Two equivalent implementations selectable via method=:
#
#   "nr"    — direct Newton-Raphson on the orthogonal PL score (default)
#   "coxph" — expand to counting-process (start,stop,event) format and call
#              survival::coxph(..., ties="breslow", robust=TRUE)
#
# Both return list(beta_hat, precomp, coxph_fit) where coxph_fit is NULL
# for method="nr".  SE/inference always use compute_sandwich() so that
# results are directly comparable across methods.

suppressPackageStartupMessages(library(survival))

# ── Pre-compute time-varying nuisances at all event times ─────────────────────
# Calls R_fn_vec once per arm (2 total) instead of once per event time.
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

# ── Score G(beta) — d-vector ──────────────────────────────────────────────────
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
    ew     <- exp(theta - max(theta))
    S0     <- sum(ew)
    q_bar  <- colSums(ew * q_mat) / S0

    n_evt <- sum(evt_now[in_risk])
    G     <- G + colSums(q_mat[evt_now[in_risk], , drop = FALSE]) - n_evt * q_bar
  }
  G / nrow(basis)
}

# ── Jacobian dG/dbeta — d x d matrix ─────────────────────────────────────────
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

# ── Expand to counting-process format (for method="coxph") ───────────────────
# Returns one row per (subject, event-time interval) the subject is at risk.
# Intervals are (t_{k-1}, t_k] with consecutive event times, so coxph sees
# non-overlapping intervals per subject and handles ties="breslow" correctly.
make_expanded <- function(dat, basis, precomp) {
  d    <- ncol(basis)
  K    <- length(precomp$event_times)
  rows <- vector("list", K)

  for (k in seq_len(K)) {
    t   <- precomp$event_times[k]
    t0  <- if (k == 1L) 0 else precomp$event_times[k - 1L]
    rs  <- precomp$nuis_at[[k]]
    idx <- which(dat$Y >= t)

    q     <- (dat$W[idx] - rs$a_t[idx]) * basis[idx, , drop = FALSE]
    df_k  <- as.data.frame(q)
    names(df_k) <- paste0("q", seq_len(d))
    df_k$id    <- idx
    df_k$start <- t0
    df_k$stop  <- t
    df_k$event <- as.integer(dat$Del[idx] == 1L & dat$Y[idx] == t)
    df_k$nut   <- rs$nu_t[idx]
    rows[[k]]  <- df_k
  }
  do.call(rbind, rows)
}

# ── Main estimator ────────────────────────────────────────────────────────────
fit_rso_coxph <- function(dat, basis, nuis,
                          method = c("nr", "coxph"),
                          max_iter = 100, tol = 1e-10) {
  method  <- match.arg(method)
  precomp <- precompute_nuis_at(dat, nuis)
  d       <- ncol(basis)

  if (method == "nr") {
    beta <- rep(0, d)
    for (iter in seq_len(max_iter)) {
      g    <- compute_score(beta, dat, basis, precomp)
      if (max(abs(g)) < tol) break
      J    <- compute_jacobian(beta, dat, basis, precomp)
      step <- tryCatch(solve(J, g), error = function(e) rep(NA_real_, d))
      if (anyNA(step)) break
      if (max(abs(step)) > 2) step <- step * 2 / max(abs(step))
      beta <- beta - step
    }
    list(beta_hat = beta, precomp = precomp, coxph_fit = NULL)

  } else {
    expanded <- make_expanded(dat, basis, precomp)
    q_terms  <- paste(paste0("q", seq_len(d)), collapse = " + ")
    fml      <- as.formula(
      paste0("Surv(start, stop, event) ~ ", q_terms,
             " + offset(nut) + cluster(id)")
    )
    coxfit   <- coxph(fml, data = expanded, ties = "breslow")
    list(beta_hat = coef(coxfit), precomp = precomp, coxph_fit = coxfit)
  }
}
