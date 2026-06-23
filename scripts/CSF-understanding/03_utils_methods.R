# =============================================================
# 03_utils_methods.R
# All estimators for CSF simulation study.
# Assumes 01_utils_dgp.R and 02_utils_oracle.R are sourced.
# =============================================================

library(grf)
library(randomForestSRC)
library(survival)

# ---- Helper: RMST from step-function survival curves ----
# surv_mat : n_test x n_times matrix of survival probabilities
# times    : corresponding event times (from training set)
# h        : horizon
rmst_from_stepfun <- function(surv_mat, times, h) {
  idx_h <- times <= h
  t_use <- c(0, times[idx_h], h)

  apply(surv_mat, 1, function(s) {
    s_vals  <- c(1, s[idx_h])
    # extrapolate: S(h) = last observed value at or before h
    s_at_h  <- if (any(idx_h)) s[max(which(idx_h))] else 1
    s_vals  <- c(s_vals, s_at_h)
    # left Riemann sum (step function)
    sum(diff(t_use) * head(s_vals, -1))
  })
}

# ---- Method 1: CSF (grf::causal_survival_forest) ----

fit_csf <- function(X_train, U_train, Delta_train, A_train, X_test,
                    h, num.trees = 1000, num.threads = 1) {
  csf <- causal_survival_forest(
    X           = X_train,
    Y           = U_train,
    W           = A_train,
    D           = Delta_train,
    horizon     = h,
    target      = "RMST",
    num.trees   = num.trees,
    num.threads = num.threads
  )
  list(tau_hat = predict(csf, X_test, num.threads = num.threads)$predictions)
}

# ---- Method 2: T-learner / VT (two survival forests) ----

fit_tlearner <- function(X_train, U_train, Delta_train, A_train, X_test,
                         h, num.trees = 500, ncore = 1) {
  idx0 <- A_train == 0; idx1 <- A_train == 1

  dat0 <- data.frame(time = U_train[idx0], status = Delta_train[idx0],
                     X_train[idx0, , drop = FALSE])
  dat1 <- data.frame(time = U_train[idx1], status = Delta_train[idx1],
                     X_train[idx1, , drop = FALSE])

  rsf0 <- rfsrc(Surv(time, status) ~ ., data = dat0, ntree = num.trees, nodesize = 15, ncore = ncore)
  rsf1 <- rfsrc(Surv(time, status) ~ ., data = dat1, ntree = num.trees, nodesize = 15, ncore = ncore)

  X_test_df <- as.data.frame(X_test)
  pred0 <- predict(rsf0, X_test_df)
  pred1 <- predict(rsf1, X_test_df)

  rmst0 <- rmst_from_stepfun(pred0$survival, pred0$time.interest, h)
  rmst1 <- rmst_from_stepfun(pred1$survival, pred1$time.interest, h)

  list(tau_hat = rmst1 - rmst0)
}

# ---- Method 3: S-learner (SRC1: one RSF with A, X) ----

fit_slearner <- function(X_train, U_train, Delta_train, A_train, X_test,
                         h, num.trees = 500, ncore = 1) {
  datS  <- data.frame(time = U_train, status = Delta_train,
                      A = A_train, X_train)
  rsfS  <- rfsrc(Surv(time, status) ~ ., data = datS, ntree = num.trees, nodesize = 15, ncore = ncore)

  test0 <- data.frame(A = 0, X_test)
  test1 <- data.frame(A = 1, X_test)

  pred0 <- predict(rsfS, test0)
  pred1 <- predict(rsfS, test1)

  rmst0 <- rmst_from_stepfun(pred0$survival, pred0$time.interest, h)
  rmst1 <- rmst_from_stepfun(pred1$survival, pred1$time.interest, h)

  list(tau_hat = rmst1 - rmst0)
}

# ---- Method 4: SRC2 (S-learner with A*X interaction features) ----

fit_src2 <- function(X_train, U_train, Delta_train, A_train, X_test,
                     h, num.trees = 500, ncore = 1) {
  AX_tr <- A_train * X_train
  colnames(AX_tr) <- paste0("AX", seq_len(ncol(X_train)))

  datS2 <- data.frame(time = U_train, status = Delta_train,
                       A = A_train, X_train, AX_tr)
  rsfS2 <- rfsrc(Surv(time, status) ~ ., data = datS2, ntree = num.trees, nodesize = 15, ncore = ncore)

  AX_te0 <- matrix(0, nrow(X_test), ncol(X_test), dimnames = list(NULL, paste0("AX", seq_len(ncol(X_test)))))
  AX_te1 <- X_test; colnames(AX_te1) <- paste0("AX", seq_len(ncol(X_test)))

  pred0 <- predict(rsfS2, data.frame(A = 0, X_test, AX_te0))
  pred1 <- predict(rsfS2, data.frame(A = 1, X_test, AX_te1))

  rmst0 <- rmst_from_stepfun(pred0$survival, pred0$time.interest, h)
  rmst1 <- rmst_from_stepfun(pred1$survival, pred1$time.interest, h)

  list(tau_hat = rmst1 - rmst0)
}

# ---- Method 5: Cox interaction model ----
# Fits Cox with A*X1 interaction; converts to RMST-CATE via numerical integral

fit_cox_interaction <- function(X_train, U_train, Delta_train, A_train, X_test, h) {
  p        <- ncol(X_train)
  Xdf_tr   <- as.data.frame(X_train)
  colnames(Xdf_tr) <- paste0("X", seq_len(p))

  dat_cox  <- data.frame(time = U_train, status = Delta_train, A = A_train, Xdf_tr)
  fmla     <- as.formula(
    paste("Surv(time, status) ~",
          paste(paste0("X", seq_len(p)), collapse = " + "),
          "+ A + A:X1")
  )
  cox_fit  <- coxph(fmla, data = dat_cox)

  Xdf_te   <- as.data.frame(X_test)
  colnames(Xdf_te) <- paste0("X", seq_len(p))
  newdat0  <- data.frame(A = 0, Xdf_te)
  newdat1  <- data.frame(A = 1, Xdf_te)

  lp0 <- predict(cox_fit, newdata = newdat0, type = "lp")
  lp1 <- predict(cox_fit, newdata = newdat1, type = "lp")

  # centered=TRUE: H0 is at column means, consistent with predict(type="lp")
  bh      <- basehaz(cox_fit, centered = TRUE)
  t_b     <- bh$time
  H_b     <- bh$hazard
  idx_h   <- t_b <= h
  t_grid  <- c(0, t_b[idx_h], h)
  H_grid  <- c(0, H_b[idx_h], if (any(idx_h)) tail(H_b[idx_h], 1) else 0)

  n_test  <- nrow(X_test)
  rmst0   <- rmst1 <- numeric(n_test)
  for (i in seq_len(n_test)) {
    s0 <- exp(-H_grid * exp(lp0[i]))
    s1 <- exp(-H_grid * exp(lp1[i]))
    rmst0[i] <- sum(diff(t_grid) * head(s0, -1))
    rmst1[i] <- sum(diff(t_grid) * head(s1, -1))
  }

  list(tau_hat = rmst1 - rmst0)
}

# ---- Method 6: IPCW Causal Forest ----
# Estimates censoring survival via Cox, weights complete cases,
# then runs grf::causal_forest on truncated outcomes

fit_ipcw_cf <- function(X_train, U_train, Delta_train, A_train, X_test,
                        h, num.trees = 1000, num.threads = 1) {
  n        <- length(U_train)
  Delta_h  <- as.integer(Delta_train == 1 | U_train >= h)

  # Fit censoring Cox model (event = 1 - Delta, i.e. censored)
  p        <- ncol(X_train)
  Xdf_tr   <- as.data.frame(X_train)
  colnames(Xdf_tr) <- paste0("X", seq_len(p))
  dat_cens <- data.frame(time = U_train, status = 1L - Delta_train,
                         A = A_train, Xdf_tr)
  fmla_c   <- as.formula(
    paste("Surv(time, status) ~",
          paste(paste0("X", seq_len(p)), collapse = " + "),
          "+ A")
  )
  cox_cens <- coxph(fmla_c, data = dat_cens)

  newdat_tr <- data.frame(A = A_train, Xdf_tr)
  lp_c      <- predict(cox_cens, newdata = newdat_tr, type = "lp")
  bh_c      <- basehaz(cox_cens, centered = FALSE)
  t_c       <- c(0, bh_c$time)
  H_c       <- c(0, bh_c$hazard)

  # S_C(min(U_i, h) | X_i, A_i)
  u_h_vec <- pmin(U_train, h)
  sc_vec  <- sapply(seq_len(n), function(i) {
    H_ui <- approx(t_c, H_c, xout = u_h_vec[i], method = "linear",
                   rule = 2, ties = "ordered")$y
    exp(-H_ui * exp(lp_c[i]))
  })
  ipcw_w  <- 1 / pmax(sc_vec, 0.02)  # truncate to cap extreme weights

  idx_cc  <- Delta_h == 1
  if (sum(idx_cc) < 30) {
    warning("Too few complete cases for IPCW-CF")
    return(list(tau_hat = rep(NA_real_, nrow(X_test))))
  }

  cf <- causal_forest(
    X              = X_train[idx_cc, ],
    Y              = pmin(U_train[idx_cc], h),
    W              = A_train[idx_cc],
    sample.weights = ipcw_w[idx_cc],
    num.trees      = num.trees,
    num.threads    = num.threads
  )
  list(tau_hat = predict(cf, X_test, num.threads = num.threads)$predictions)
}

# ---- Method 7: Oracle Linear CSF ----
# Solves the linear CSF estimating equation with oracle nuisances.

fit_linear_csf_oracle <- function(X_train, U_train, Delta_train, A_train, X_test,
                                  h = DGP_PARAMS$h, params = DGP_PARAMS,
                                  verbose = FALSE) {
  n    <- nrow(X_train)
  BH   <- compute_BH_oracle(X_train, U_train, Delta_train, A_train, h, params, verbose)
  B    <- BH[, "B"]
  H    <- BH[, "H"]

  # Use only rows with valid BH
  ok   <- is.finite(B) & is.finite(H)
  if (sum(ok) < 10) {
    warning("Too few valid BH rows in oracle linear CSF")
    return(list(tau_hat = rep(NA_real_, nrow(X_test)), beta = c(NA, NA)))
  }

  e_hat <- efun(X_train[ok, , drop = FALSE])
  r     <- A_train[ok] - e_hat
  P     <- cbind(1, X_train[ok, 1])

  M <- t(P) %*% (P * (r^2 * H[ok])) / sum(ok)
  v <- colMeans(P * (r * B[ok]))

  beta_hat <- tryCatch(solve(M, v), error = function(e) c(NA_real_, NA_real_))

  P_test   <- cbind(1, X_test[, 1])
  list(tau_hat = as.numeric(P_test %*% beta_hat), beta = beta_hat)
}

# ---- Method 7b: Oracle Linear CSF for DGP2 ----
# Uses compute_BH_oracle_dgp2 (analytical Q_w via exponential model).

fit_linear_csf_oracle_dgp2 <- function(X_train, U_train, Delta_train, A_train, X_test,
                                        r0_vec, r1_vec,
                                        h, params, verbose = FALSE) {
  BH <- compute_BH_oracle_dgp2(X_train, U_train, Delta_train, A_train,
                                r0_vec, r1_vec, h, params, verbose)
  B  <- BH[, "B"]
  H  <- BH[, "H"]

  ok <- is.finite(B) & is.finite(H)
  if (sum(ok) < 10) {
    warning("Too few valid BH rows in oracle linear CSF (DGP2)")
    return(list(tau_hat = rep(NA_real_, nrow(X_test)), beta = c(NA, NA)))
  }

  e_hat <- efun(X_train[ok, , drop = FALSE])
  r     <- A_train[ok] - e_hat
  P     <- cbind(1, X_train[ok, 1])

  M <- t(P) %*% (P * (r^2 * H[ok])) / sum(ok)
  v <- colMeans(P * (r * B[ok]))

  beta_hat <- tryCatch(solve(M, v), error = function(e) c(NA_real_, NA_real_))
  P_test   <- cbind(1, X_test[, 1])
  list(tau_hat = as.numeric(P_test %*% beta_hat), beta = beta_hat)
}

# ---- Evaluation metrics ----

evaluate_metrics <- function(tau_hat, tau_true) {
  if (any(is.na(tau_hat))) {
    nm <- c("MSE", "RMSE", "Bias", "Cor", "SignError")
    return(setNames(rep(NA_real_, 5), nm))
  }
  c(
    MSE       = mean((tau_hat - tau_true)^2),
    RMSE      = sqrt(mean((tau_hat - tau_true)^2)),
    Bias      = mean(tau_hat - tau_true),
    Cor       = cor(tau_hat, tau_true),
    SignError = mean(sign(tau_hat) != sign(tau_true))
  )
}
