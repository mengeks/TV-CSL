#!/usr/bin/env Rscript
# Nuisance estimation: propensity, outcome LP, censoring survival
#
# Exported:
#   estimate_nuisances_crossfit(dat, n_folds)
#     -> list(e, eta0, eta1, R_fn_vec)
#
# Nuisance models: pooled Cox with W as binary covariate (no interaction).
# This uses a SHARED baseline hazard for both arms, so the LP predictions
# at W=0 and W=1 are on the same scale — critical for a_t(X) calibration.
# It also doubles the effective training size vs arm-specific Cox models.
#
# R_fn_vec(t_vec, w): n x length(t_vec) matrix of R_w(t_k, X_i).
# Called ONCE per arm from precompute_nuis_at().

suppressPackageStartupMessages(library(survival))

estimate_nuisances_crossfit <- function(dat, n_folds = 2) {
  n   <- nrow(dat$X)
  X   <- dat$X
  W   <- dat$W
  Y   <- dat$Y
  Del <- dat$Del

  fold <- sample(rep(seq_len(n_folds), length.out = n))

  e    <- numeric(n)
  eta0 <- numeric(n)
  eta1 <- numeric(n)

  fold_cache <- vector("list", n_folds)

  for (k in seq_len(n_folds)) {
    tr <- fold != k
    ev <- fold == k

    df_tr   <- data.frame(W = W[tr], X1 = X[tr, 1], X2 = X[tr, 2],
                          Y = Y[tr],  Del = Del[tr])
    # Evaluation frames at W=0 and W=1 for the held-out fold
    df_ev0  <- data.frame(W = 0L, X1 = X[ev, 1], X2 = X[ev, 2])
    df_ev1  <- data.frame(W = 1L, X1 = X[ev, 1], X2 = X[ev, 2])

    # Propensity
    fit_e   <- glm(W ~ X1 + X2, family = binomial, data = df_tr)
    e[ev]   <- predict(fit_e, newdata = df_ev0, type = "response")

    # Pooled outcome Cox: shared baseline -> LP at W=0 and W=1 are comparable
    fit_cox  <- coxph(Surv(Y, Del) ~ X1 + X2 + W, data = df_tr)
    eta0[ev] <- predict(fit_cox, newdata = df_ev0, type = "lp")
    eta1[ev] <- predict(fit_cox, newdata = df_ev1, type = "lp")

    # Pooled censoring Cox (flip Delta)
    df_cens  <- df_tr; df_cens$Del <- 1L - df_cens$Del
    fit_cens <- coxph(Surv(Y, Del) ~ X1 + X2 + W, data = df_cens)

    fold_cache[[k]] <- list(
      ev    = ev,
      # LP at W=0 and W=1 from pooled models (computed once per fold)
      lp_s0 = predict(fit_cox,  newdata = df_ev0, type = "lp"),
      lp_s1 = predict(fit_cox,  newdata = df_ev1, type = "lp"),
      lp_g0 = predict(fit_cens, newdata = df_ev0, type = "lp"),
      lp_g1 = predict(fit_cens, newdata = df_ev1, type = "lp"),
      # Shared baseline cumulative hazard (centered=FALSE -> baseline at W=0)
      bh_s  = basehaz(fit_cox,  centered = FALSE),
      bh_g  = basehaz(fit_cens, centered = FALSE)
    )
  }

  e <- pmin(pmax(e, 0.02), 0.98)

  # R_fn_vec(t_vec, w): n x K matrix, no repeated model/basehaz calls
  R_fn_vec <- function(t_vec, w) {
    R_mat <- matrix(0, n, length(t_vec))
    for (k in seq_len(n_folds)) {
      fc   <- fold_cache[[k]]
      ev   <- fc$ev
      lp_s <- if (w == 0L) fc$lp_s0 else fc$lp_s1
      lp_g <- if (w == 0L) fc$lp_g0 else fc$lp_g1

      H0s <- approx_cumhaz_vec(fc$bh_s, t_vec)   # length-K, shared baseline
      H0g <- approx_cumhaz_vec(fc$bh_g, t_vec)

      S_w <- exp(-outer(exp(lp_s), H0s))           # n_ev x K
      G_w <- exp(-outer(exp(lp_g), H0g))
      R_mat[ev, ] <- S_w * G_w
    }
    R_mat
  }

  list(e = e, eta0 = eta0, eta1 = eta1, R_fn_vec = R_fn_vec)
}

# Interpolate cumulative hazard at a vector of times from basehaz() output.
# Prepend (0, 0) so times before the first event return 0.
approx_cumhaz_vec <- function(bh, t_vec) {
  approx(c(0, bh$time), c(0, bh$hazard), xout = t_vec, rule = 2)$y
}
