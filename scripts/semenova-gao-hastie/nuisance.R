#!/usr/bin/env Rscript
# Nuisance estimation: propensity, outcome LP, censoring survival
#
# Exported functions:
#   estimate_nuisances_crossfit(dat, n_folds)  -> list(e, eta0, eta1, R_fn)
#   nuisances_to_riskset(t, nuis)              -> list(a_t, nu_t)  per-row at time t

suppressPackageStartupMessages(library(survival))

# Cross-fitted nuisance estimation using parametric models
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

  # Store per-fold survival model objects to compute R_fn later
  fold_models <- vector("list", n_folds)

  for (k in seq_len(n_folds)) {
    tr <- fold != k
    ev <- fold == k

    df_tr  <- data.frame(W = W[tr],  X1 = X[tr, 1], X2 = X[tr, 2],
                         Y = Y[tr],  Del = Del[tr])
    df_ev  <- data.frame(X1 = X[ev, 1], X2 = X[ev, 2], W = W[ev])

    # Propensity
    fit_e     <- glm(W ~ X1 + X2, family = binomial, data = df_tr)
    e[ev]     <- predict(fit_e, newdata = df_ev, type = "response")

    # Outcome LP: separate Cox models by arm to get eta0, eta1
    df_tr0    <- df_tr[df_tr$W == 0, ]
    df_tr1    <- df_tr[df_tr$W == 1, ]
    fit_cox0  <- coxph(Surv(Y, Del) ~ X1 + X2, data = df_tr0, x = TRUE)
    fit_cox1  <- coxph(Surv(Y, Del) ~ X1 + X2, data = df_tr1, x = TRUE)
    eta0[ev]  <- predict(fit_cox0, newdata = df_ev, type = "lp")
    eta1[ev]  <- predict(fit_cox1, newdata = df_ev, type = "lp")

    # Censoring survival: model 1-Del (censoring indicator) by arm
    df_tr0c   <- df_tr; df_tr0c$Del <- 1L - df_tr0c$Del
    fit_cens0 <- coxph(Surv(Y, Del) ~ X1 + X2,
                       data = df_tr0c[df_tr0c$W == 0, ], x = TRUE)
    fit_cens1 <- coxph(Surv(Y, Del) ~ X1 + X2,
                       data = df_tr0c[df_tr0c$W == 1, ], x = TRUE)

    fold_models[[k]] <- list(
      ev = ev,
      fit_cox0 = fit_cox0, fit_cox1 = fit_cox1,
      fit_cens0 = fit_cens0, fit_cens1 = fit_cens1
    )
  }

  e <- pmin(pmax(e, 0.02), 0.98)

  # Build R_fn: R_w(t, X) = S_w(t|X) * G_w(t|X), evaluated at all n rows
  # Uses fold-specific models to avoid data leakage
  R_fn <- function(t, w) {
    R_vals <- numeric(n)
    for (k in seq_len(n_folds)) {
      fm  <- fold_models[[k]]
      ev  <- fm$ev
      df_ev <- data.frame(X1 = X[ev, 1], X2 = X[ev, 2])

      fit_s <- if (w == 0) fm$fit_cox0 else fm$fit_cox1
      fit_g <- if (w == 0) fm$fit_cens0 else fm$fit_cens1

      lp_s <- predict(fit_s, newdata = df_ev, type = "lp")
      lp_g <- predict(fit_g, newdata = df_ev, type = "lp")

      # Breslow cumulative baseline hazard at time t for each model
      bh_s <- basehaz(fit_s, centered = FALSE)
      bh_g <- basehaz(fit_g, centered = FALSE)
      H0s_t <- approx_cumhaz(bh_s, t)
      H0g_t <- approx_cumhaz(bh_g, t)

      S_w <- exp(-H0s_t * exp(lp_s))
      G_w <- exp(-H0g_t * exp(lp_g))
      R_vals[ev] <- S_w * G_w
    }
    R_vals
  }

  list(e = e, eta0 = eta0, eta1 = eta1, R_fn = R_fn)
}

# Linearly interpolate cumulative hazard at time t from basehaz() output
approx_cumhaz <- function(bh, t) {
  if (t <= bh$time[1]) return(0)
  if (t >= bh$time[length(bh$time)]) return(bh$hazard[length(bh$hazard)])
  approx(bh$time, bh$hazard, xout = t, rule = 2)$y
}

# Compute risk-set modified propensity a_t(X) and offset nu_t(X) at time t
nuisances_to_riskset <- function(t, nuis) {
  e    <- nuis$e
  eta0 <- nuis$eta0
  eta1 <- nuis$eta1
  R0   <- nuis$R_fn(t, 0)
  R1   <- nuis$R_fn(t, 1)

  # Numerator / denominator of a_t
  num  <- e * R1 * exp(eta1)
  den  <- num + (1 - e) * R0 * exp(eta0)
  a_t  <- num / den

  # Clamp to avoid 0/1
  a_t  <- pmin(pmax(a_t, 1e-6), 1 - 1e-6)

  nu_t <- a_t * eta1 + (1 - a_t) * eta0

  list(a_t = a_t, nu_t = nu_t)
}
