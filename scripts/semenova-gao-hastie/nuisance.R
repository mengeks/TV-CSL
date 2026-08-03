#!/usr/bin/env Rscript
# Nuisance estimation: propensity, outcome LP, censoring survival
#
# Exported:
#   estimate_nuisances_crossfit(dat, n_folds)
#     -> list(e, eta0, eta1, R_fn_vec)
#
#   R_fn_vec(t_vec, w) returns an n x length(t_vec) matrix of
#   R_w(t_k, X_i) = S_w(t_k|X_i) * G_w(t_k|X_i), cross-fit over folds.
#   Called ONCE per arm from precompute_nuis_at() rather than per event time.

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

  # Per-fold: store evaluation indices, linear predictors, basehaz tables
  fold_cache <- vector("list", n_folds)

  for (k in seq_len(n_folds)) {
    tr <- fold != k
    ev <- fold == k

    df_tr  <- data.frame(W = W[tr], X1 = X[tr, 1], X2 = X[tr, 2],
                         Y = Y[tr], Del = Del[tr])
    df_ev  <- data.frame(X1 = X[ev, 1], X2 = X[ev, 2])

    # Propensity
    fit_e   <- glm(W ~ X1 + X2, family = binomial, data = df_tr)
    e[ev]   <- predict(fit_e, newdata = df_ev, type = "response")

    # Outcome LP by arm
    fit_cox0 <- coxph(Surv(Y, Del) ~ X1 + X2, data = df_tr[df_tr$W == 0, ])
    fit_cox1 <- coxph(Surv(Y, Del) ~ X1 + X2, data = df_tr[df_tr$W == 1, ])
    eta0[ev] <- predict(fit_cox0, newdata = df_ev, type = "lp")
    eta1[ev] <- predict(fit_cox1, newdata = df_ev, type = "lp")

    # Censoring (flip Delta) by arm
    df_cens <- df_tr; df_cens$Del <- 1L - df_cens$Del
    fit_cen0 <- coxph(Surv(Y, Del) ~ X1 + X2, data = df_cens[df_cens$W == 0, ])
    fit_cen1 <- coxph(Surv(Y, Del) ~ X1 + X2, data = df_cens[df_cens$W == 1, ])

    fold_cache[[k]] <- list(
      ev    = ev,
      # linear predictors (computed once per fold, reused across all event times)
      lp_s0 = predict(fit_cox0, newdata = df_ev, type = "lp"),
      lp_s1 = predict(fit_cox1, newdata = df_ev, type = "lp"),
      lp_g0 = predict(fit_cen0, newdata = df_ev, type = "lp"),
      lp_g1 = predict(fit_cen1, newdata = df_ev, type = "lp"),
      # baseline cumulative hazard tables (computed once per fold)
      bh_s0 = basehaz(fit_cox0, centered = FALSE),
      bh_s1 = basehaz(fit_cox1, centered = FALSE),
      bh_g0 = basehaz(fit_cen0, centered = FALSE),
      bh_g1 = basehaz(fit_cen1, centered = FALSE)
    )
  }

  e <- pmin(pmax(e, 0.02), 0.98)

  # R_fn_vec(t_vec, w): n x K matrix, no repeated model calls
  R_fn_vec <- function(t_vec, w) {
    R_mat <- matrix(0, n, length(t_vec))
    for (k in seq_len(n_folds)) {
      fc  <- fold_cache[[k]]
      ev  <- fc$ev
      lp_s <- if (w == 0) fc$lp_s0 else fc$lp_s1
      lp_g <- if (w == 0) fc$lp_g0 else fc$lp_g1
      bh_s <- if (w == 0) fc$bh_s0 else fc$bh_s1
      bh_g <- if (w == 0) fc$bh_g0 else fc$bh_g1

      H0s <- approx_cumhaz_vec(bh_s, t_vec)   # length-K vector
      H0g <- approx_cumhaz_vec(bh_g, t_vec)

      # outer products: n_ev x K, then element-wise product
      S_w <- exp(-outer(exp(lp_s), H0s))
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
