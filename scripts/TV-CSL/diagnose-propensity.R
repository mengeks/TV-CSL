# diagnose-propensity.R
#
# Compares the GAM-estimated event propensity a_hat_t(x) against the true
# analytical a_t(x) derived from the simulation DGP.  Answers three questions:
#
#   1. How bad is the GAM estimate?  (MSE, MAE, correlation on the prob scale)
#   2. Where does it fail?  (residual patterns vs. t and covariates)
#   3. How complex is the truth?  (R^2 of additive vs. interaction models of true logit)
#
# DGP: A|X ~ Exp(exp(x2+x3)),  h0(t)=t,  tau(x)=x1+x2+x3,
#      logit a_t(x) = x1+2x2+2x3 + r*t - d*t^2 + log I(t,x)
#      where r=exp(x2+x3), d=0.5*exp(eta0)*(exp(tau)-1), I = integral_0^t exp(d*s^2-r*s) ds.

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
})

source("R/datagen-helper.R")
source("R/old/data-handler.R")

set.seed(2024)
N_DIAG <- 2000  # subjects for the diagnosis dataset

# ---- True a_t(x) -------------------------------------------------------
compute_true_a_t <- function(t, x1, x2, x3, eta0) {
  if (is.na(t) || t <= 0) return(NA_real_)
  r   <- exp(x2 + x3)
  tau <- x1 + x2 + x3
  d   <- 0.5 * exp(eta0) * (exp(tau) - 1)
  M   <- max(0, d * t^2 - r * t)
  val <- tryCatch(
    integrate(function(s) exp(d * s^2 - r * s - M), 0, t, rel.tol = 1e-6)$value,
    error = function(e) NA_real_
  )
  if (is.na(val) || val <= 0) return(NA_real_)
  plogis(x1 + 2*x2 + 2*x3 + r*t - d*t^2 + M + log(val))
}

# ---- Run diagnosis for each DGP ----------------------------------------
for (eta_type in c("linear", "non-linear")) {

  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("  eta_type = %s\n", eta_type))
  cat(sprintf("%s\n", strrep("=", 60)))

  dat <- generate_simulated_data(
    n                     = N_DIAG,
    is_time_varying       = TRUE,
    light_censoring       = FALSE,
    lambda_C              = 0.1,
    p                     = 3,
    baseline_type         = "linear",
    eta_type              = eta_type,
    X_distribution        = "normal",
    X_cov_type            = "identity",
    tx_difficulty         = "simple",
    HTE_type              = "linear",
    linear_intercept      = 0,
    linear_slope_multiplier = 2.5,
    linear_HTE_multiplier = 1,
    seed_value            = 2024,
    verbose               = 0
  )

  # Treatment status at observed time: W(U_i) = 1{A_i < U_i}
  dat$W_at_event <- as.integer(dat$A < dat$U)

  n_events <- sum(dat$Delta)
  cat(sprintf("n=%d  events=%d (%.0f%%)  treated-at-event=%d (%.0f%% of events)\n",
              N_DIAG, n_events, 100*n_events/N_DIAG,
              sum(dat$Delta & dat$W_at_event),
              100 * mean(dat$W_at_event[dat$Delta == 1])))

  # ---- 1. Fit GAM exactly as in cox-time-varying-prop --------------------
  events <- dat %>% filter(Delta == 1) %>%
    rename(tstop = U, W = W_at_event)

  gam_mod <- tryCatch(
    gam(W ~ s(tstop, k = 5) + X.1 + X.2 + X.3,
        data = events, family = binomial, method = "REML"),
    error = function(e) {
      message("GAM failed, using GLM fallback"); NULL
    }
  )
  if (is.null(gam_mod))
    gam_mod <- glm(W ~ tstop + X.1 + X.2 + X.3, data = events, family = binomial)

  cat("\nGAM coefficients (parametric part):\n")
  ptab <- summary(gam_mod)$p.table
  print(round(ptab, 4))
  if (inherits(gam_mod, "gam"))
    cat(sprintf("Deviance explained by GAM: %.1f%%\n", 100*summary(gam_mod)$dev.expl))

  # ---- 2. Compute true a_t(x) for event subjects ------------------------
  # Events are the natural evaluation domain: the propensity is used only when dN(t)=1.
  eval <- events  # already filtered to events (Delta==1)

  cat(sprintf("\nComputing true a_t(x) for %d event subjects...\n", nrow(eval)))
  a_true <- mapply(compute_true_a_t,
                   t    = eval$tstop,
                   x1   = eval$X.1,
                   x2   = eval$X.2,
                   x3   = eval$X.3,
                   eta0 = eval$eta_0)

  a_hat <- as.vector(plogis(predict(gam_mod, newdata = eval, type = "link")))

  ok <- is.finite(a_hat) & is.finite(a_true) &
        a_true > 1e-6 & a_true < 1 - 1e-6
  cat(sprintf("Valid pairs: %d / %d\n", sum(ok), length(ok)))

  # ---- 3. Quantify error -------------------------------------------------
  err_p <- a_hat[ok] - a_true[ok]           # probability scale
  err_l <- qlogis(a_hat[ok]) - qlogis(a_true[ok])  # logit scale

  cat(sprintf("\n--- Propensity estimation error (a_hat vs a_true) ---\n"))
  cat(sprintf("  MSE  (prob  scale): %.5f\n",  mean(err_p^2)))
  cat(sprintf("  MAE  (prob  scale): %.5f\n",  mean(abs(err_p)))  )
  cat(sprintf("  Bias (prob  scale): %+.5f\n", mean(err_p)))
  cat(sprintf("  MSE  (logit scale): %.5f\n",  mean(err_l^2)))
  cat(sprintf("  Bias (logit scale): %+.5f\n", mean(err_l)))
  cat(sprintf("  Pearson corr(a_hat, a_true): %.4f\n", cor(a_hat[ok], a_true[ok])))
  cat(sprintf("  Var(true logit): %.4f  |  Var(fitted logit): %.4f\n",
              var(qlogis(a_true[ok])), var(qlogis(a_hat[ok]))))

  # ---- 4. Residual patterns -----------------------------------------------
  df <- data.frame(
    t    = eval$tstop[ok],
    x1   = eval$X.1[ok],
    x2   = eval$X.2[ok],
    x3   = eval$X.3[ok],
    err  = err_l,
    a_hat  = a_hat[ok],
    a_true = a_true[ok]
  )

  tq <- quantile(df$t, c(0, .25, .5, .75, 1))
  df$t_quartile <- cut(df$t, tq, labels = paste0("Q", 1:4), include.lowest = TRUE)

  cat("\nMean logit-error (GAM minus truth) by event-time quartile:\n")
  cat("  (Q1=early events, Q4=late events; systematic pattern → time trend misfit)\n")
  err_by_t <- tapply(df$err, df$t_quartile, function(x)
    sprintf("%+.3f (n=%d)", mean(x), length(x)))
  print(err_by_t)

  cat("\nMean logit-error by X.1 quantile:\n")
  x1q <- quantile(df$x1, c(0, .33, .67, 1))
  df$x1_tile <- cut(df$x1, x1q, labels = c("lo","mid","hi"), include.lowest = TRUE)
  print(tapply(df$err, df$x1_tile, function(x) sprintf("%+.3f", mean(x))))

  cat("\nMean logit-error by (X.2+X.3) quantile (drives r(x)=exp(x2+x3)):\n")
  r_sum <- df$x2 + df$x3
  rq <- quantile(r_sum, c(0, .33, .67, 1))
  df$r_tile <- cut(r_sum, rq, labels = c("lo","mid","hi"), include.lowest = TRUE)
  print(tapply(df$err, df$r_tile, function(x) sprintf("%+.3f", mean(x))))

  # ---- 5. Complexity of the truth: how much do interactions matter? -------
  cat("\n--- R^2 of OLS models for logit(a_true): measures interaction complexity ---\n")

  lm_add <- lm(qlogis(a_true[ok]) ~ t + x1 + x2 + x3, data = df)
  cat(sprintf("  Additive  t + x1 + x2 + x3                     R^2 = %.4f\n",
              summary(lm_add)$r.squared))

  lm_t2  <- lm(qlogis(a_true[ok]) ~ poly(t,2) + x1 + x2 + x3, data = df)
  cat(sprintf("  + t^2 term                                       R^2 = %.4f\n",
              summary(lm_t2)$r.squared))

  lm_int <- lm(qlogis(a_true[ok]) ~ poly(t,2) + x1 + x2 + x3 +
                 t:x1 + t:x2 + t:x3, data = df)
  cat(sprintf("  + t:x1 + t:x2 + t:x3 interactions               R^2 = %.4f\n",
              summary(lm_int)$r.squared))

  lm_full <- lm(qlogis(a_true[ok]) ~ poly(t,2) + x1 + x2 + x3 +
                  t:x1 + t:x2 + t:x3 + I(x2+x3) + I((x2+x3)^2) +
                  t:I(x2+x3) + I(t^2):I(x2+x3), data = df)
  cat(sprintf("  + exp(x2+x3) terms (key driver of r*t - d*t^2)  R^2 = %.4f\n",
              summary(lm_full)$r.squared))

  cat(sprintf("\nGap: additive → full interaction model: delta-R^2 = %.4f\n",
              summary(lm_full)$r.squared - summary(lm_add)$r.squared))
  cat("  → This gap is the information the additive GAM cannot capture.\n")

  # ---- 6. GAM coefficient on X.1 vs. true marginal effect ----------------
  # In the truth: logit a_t(x) has coefficient 1 on x1 (from x1+2x2+2x3 term)
  # plus an indirect contribution via d(x)*t^2 which depends on eta_0(x).
  # The GAM fits a single linear coefficient for X.1.
  cat("\nGAM coefficient on X.1:", round(coef(gam_mod)["X.1"], 4),
      "  (truth: +1 from logit formula, modified by d(x)*t^2 interaction)\n")

  # ---- 7. Where does a_t(x) vary most? -----------------------------------
  # Range of true a_t(x): shows how much propensity varies across evaluation points.
  cat(sprintf("\nTrue a_t(x) range: [%.3f, %.3f]  mean=%.3f  sd=%.3f\n",
              min(a_true[ok]), max(a_true[ok]),
              mean(a_true[ok]), sd(a_true[ok])))
  cat(sprintf("Fitted a_hat range: [%.3f, %.3f]  mean=%.3f  sd=%.3f\n",
              min(a_hat[ok]), max(a_hat[ok]),
              mean(a_hat[ok]), sd(a_hat[ok])))
}

cat(sprintf("\n%s\n", strrep("=", 60)))
cat("Proposal for better a_t(x) estimation:\n")
cat(sprintf("%s\n", strrep("-", 60)))
cat("The additive GAM fails for two structural reasons:
  1. The truth contains d(x)*t^2 where d(x)=0.5*exp(eta_0(x))*(exp(tau(x))-1),
     which is a product of a nonlinear covariate function and t^2.
  2. The integral I(t,x) = int_0^t exp(d*s^2 - r*s)ds introduces further
     t-covariate interactions through r(x)=exp(x2+x3).

Better estimation classes (in order of difficulty):
  a. Tensor-product GAM:  W ~ te(tstop, X.2, X.3, k=c(5,3,3)) + X.1
     Captures the r*t and d*t^2 structure if the smooth is flexible enough.
     Still parametric in X.1 (coefficient +1 is exact in the logit formula).
  b. Gradient-boosted trees (XGBoost/LightGBM):
     Naturally capture multiplicative interactions through tree splits.
     No assumption about the functional form of interactions.
  c. Neural network / deep hazard model:
     Can approximate arbitrary continuous functions of (t,x).
  d. Model-assisted: pre-specify the structural form
       logit a_t(x) = alpha0 + alpha1*(x1+x2+x3) + alpha2*(x2+x3)
                    + alpha3*exp(x2+x3)*t + alpha4*d_hat(x)*t^2
     and estimate (alpha0,...,alpha4) by logistic regression,
     where d_hat(x) uses a first-stage estimate of eta_0(x).
     This is semi-parametric and should be near-oracle if eta_0 is well estimated.
")
