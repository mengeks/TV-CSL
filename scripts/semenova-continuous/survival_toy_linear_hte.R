#!/usr/bin/env Rscript
# Survival analog of original_semenova_2020_toy_linear_hte.R
# (Section 5 of TV-CSL-semenova-inference-MRE.tex)
#
# DGP: Cox PH with tau(x) = 1 + x1, eta0(x) = 0.4*x2, lam0 = 0.05
#      Censoring: Exp(0.03)
#      Gamma_surv = (W - ehat) * dM / kappa  (one-step martingale pseudo-outcome)
#
# beta_true = plim(bhat) is determined empirically via a large-n pilot
# (analytically it equals c*(1,1,0) for some c; the pilot measures c).
#
# Outputs same columns as the continuous script:
#   bias, true SE, mean estimated SE, Wald coverage, bootstrap coverage
#   for all three basis coefficients + timing.

suppressPackageStartupMessages(library(survival))

# ── one MC replication ─────────────────────────────────────────────────────────
one_rep_surv <- function(n = 4000, B_boot = 500, seed = NULL, run_boot = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  ## DGP (matches LaTeX Section 5)
  X    <- matrix(rnorm(n * 2), n, 2)
  e    <- plogis(0.7 * X[, 1])
  W    <- rbinom(n, 1, e)
  tau  <- 1 + X[, 1]
  eta0 <- 0.4 * X[, 2]
  lam0 <- 0.05
  Tevt <- rexp(n, rate = lam0 * exp(eta0 + W * tau))
  Cens <- rexp(n, rate = 0.03)
  U    <- pmin(Tevt, Cens)
  Del  <- as.integer(Tevt <= Cens)

  ## Step A: 2-fold cross-fit propensity + Cox (no W) -> etahat
  fold   <- sample(rep(1:2, length.out = n))
  ehat   <- numeric(n)
  etahat <- numeric(n)
  for (k in 1:2) {
    tr   <- fold != k
    ev   <- fold == k
    dat_tr <- data.frame(W  = W[tr],  X1 = X[tr, 1], X2 = X[tr, 2],
                         U  = U[tr],  Del = Del[tr])
    dat_ev <- data.frame(X1 = X[ev, 1], X2 = X[ev, 2])

    fit_e       <- glm(W ~ X1 + X2, family = binomial, data = dat_tr)
    ehat[ev]    <- predict(fit_e, newdata = dat_ev, type = "response")

    fit_cox     <- coxph(Surv(U, Del) ~ X1 + X2, data = dat_tr)
    etahat[ev]  <- predict(fit_cox, newdata = dat_ev, type = "lp")
  }
  ehat <- pmin(pmax(ehat, 0.02), 0.98)

  ## Vectorised Breslow baseline hazard (model without W, full cross-fitted etahat)
  risk_score0   <- exp(etahat)
  ord           <- order(U)
  rs_ord        <- risk_score0[ord]
  Del_ord       <- Del[ord]
  risk_set_sum  <- rev(cumsum(rev(rs_ord)))          # sum from position j to n
  dLambda0_ord  <- ifelse(Del_ord == 1, 1 / risk_set_sum, 0)
  dLambda0      <- numeric(n)
  dLambda0[ord] <- dLambda0_ord

  ## Survival pseudo-outcome
  dM         <- Del - risk_score0 * dLambda0
  Ares       <- W - ehat
  kappa      <- mean(Ares^2 * Del)
  Gamma_surv <- Ares * dM / kappa

  ## Step B: OLS of Gamma_surv on p(X) = [1, x1, x2]
  P    <- cbind(1, X)
  bhat <- solve(crossprod(P), crossprod(P, Gamma_surv))
  res  <- as.numeric(Gamma_surv - P %*% bhat)
  Q    <- crossprod(P) / n
  Qinv <- solve(Q)

  ## Step C-i: Wald sandwich SE (timed)
  t0_wald <- proc.time()["elapsed"]
  Sig     <- crossprod(P * res) / n
  Omega   <- Qinv %*% Sig %*% Qinv / n
  se_wald <- sqrt(diag(Omega))
  t_wald  <- proc.time()["elapsed"] - t0_wald

  ## Step C-ii: multiplier bootstrap (vectorised, timed)
  if (run_boot) {
    t0_boot    <- proc.time()["elapsed"]
    XI         <- matrix(rnorm(n * B_boot), n, B_boot)
    boot_draws <- t(Qinv %*% (t(P) %*% (XI * res) / n))   # B_boot x 3
    ci_boot    <- apply(abs(boot_draws), 2, quantile, 0.95)
    t_boot     <- proc.time()["elapsed"] - t0_boot
  } else {
    ci_boot <- rep(NA_real_, 3)
    t_boot  <- 0
  }

  list(
    bhat       = as.numeric(bhat),
    se_wald    = se_wald,
    ci_boot    = ci_boot,
    t_wald     = t_wald,
    t_boot     = t_boot,
    event_rate = mean(Del)
  )
}

# ── Pilot: large n to pin down beta_true = c*(1,1,0) ──────────────────────────
cat("\n========== PILOT: n=100 000, seed=99 (no bootstrap) ==========\n")
cat("(Determines beta_true = plim(bhat) for the simplified martingale pseudo-outcome)\n\n")
pilot    <- one_rep_surv(n = 100000, B_boot = 0, seed = 99, run_boot = FALSE)
beta_true_surv <- pilot$bhat
cat(sprintf("Pilot bhat:    (%.5f,  %.5f,  %.5f)\n",
            beta_true_surv[1], beta_true_surv[2], beta_true_surv[3]))
cat(sprintf("Event rate:    %.3f\n", pilot$event_rate))
cat(sprintf("bhat[1]/bhat[0] = %.4f  (expect ~1.0 if direction is (1,1,0))\n",
            beta_true_surv[2] / beta_true_surv[1]))
cat(sprintf("bhat[2]/bhat[0] = %.4f  (expect ~0.0)\n",
            beta_true_surv[3] / beta_true_surv[1]))

# ── Smoke test (10 runs) ───────────────────────────────────────────────────────
cat("\n========== SMOKE TEST (10 runs, n=4000, B=500) ==========\n")
cat(sprintf("beta_true (pilot) = (%.4f, %.4f, %.4f)\n\n",
            beta_true_surv[1], beta_true_surv[2], beta_true_surv[3]))
cat(sprintf("%-5s  %8s  %8s  %8s  %12s  %12s\n",
            "run", "bhat[0]", "bhat[1]", "bhat[2]", "t_wald(s)", "t_boot(s)"))
cat(strrep("-", 62), "\n")
for (i in seq_len(10)) {
  r <- one_rep_surv(n = 4000, B_boot = 500, seed = i, run_boot = TRUE)
  cat(sprintf("%-5d  %8.4f  %8.4f  %8.4f  %12.5f  %12.5f\n",
              i, r$bhat[1], r$bhat[2], r$bhat[3], r$t_wald, r$t_boot))
}

# ── MC simulation (300 reps) ───────────────────────────────────────────────────
nmc <- 300
cat(sprintf("\n========== MC SIMULATION (%d reps, n=4000, B=500) ==========\n", nmc))
mc_results <- vector("list", nmc)
for (i in seq_len(nmc)) {
  if (i == 1 || i %% 50 == 0) cat(sprintf("  rep %3d / %d\n", i, nmc))
  mc_results[[i]] <- one_rep_surv(n = 4000, B_boot = 500,
                                  seed = 2000 + i, run_boot = TRUE)
}

## Collect
bhat_mat   <- do.call(rbind, lapply(mc_results, `[[`, "bhat"))
se_wald_m  <- do.call(rbind, lapply(mc_results, `[[`, "se_wald"))
ci_boot_m  <- do.call(rbind, lapply(mc_results, `[[`, "ci_boot"))
t_wald_v   <- sapply(mc_results, `[[`, "t_wald")
t_boot_v   <- sapply(mc_results, `[[`, "t_boot")
ev_rate_v  <- sapply(mc_results, `[[`, "event_rate")

## Coverage relative to pilot beta_true
cover_wald <- abs(bhat_mat - matrix(beta_true_surv, nmc, 3, byrow = TRUE)) <
              1.96 * se_wald_m
cover_boot <- abs(bhat_mat - matrix(beta_true_surv, nmc, 3, byrow = TRUE)) <
              ci_boot_m

## Print
coef_names <- c(
  sprintf("beta_0 (true=%.4f)", beta_true_surv[1]),
  sprintf("beta_1 (true=%.4f)", beta_true_surv[2]),
  sprintf("beta_2 (true=%.4f)", beta_true_surv[3])
)
cat(sprintf("\n--- Per-coefficient summary  [event rate: mean=%.3f] ---\n",
            mean(ev_rate_v)))
cat(sprintf("%-32s  %8s  %8s  %10s  %10s  %10s  %10s\n",
            "Coefficient", "Bias", "True SE", "MeanSE_wald",
            "Cov_wald", "CIhalf_boot", "Cov_boot"))
cat(strrep("-", 96), "\n")
for (j in 1:3) {
  bias      <- mean(bhat_mat[, j]) - beta_true_surv[j]
  true_se   <- sd(bhat_mat[, j])
  mean_se_w <- mean(se_wald_m[, j])
  mean_ci_b <- mean(ci_boot_m[, j])
  cov_w     <- mean(cover_wald[, j])
  cov_b     <- mean(cover_boot[, j])
  cat(sprintf("%-32s  %+8.4f  %8.4f  %10.4f  %10.3f  %10.4f  %10.3f\n",
              coef_names[j], bias, true_se, mean_se_w, cov_w, mean_ci_b, cov_b))
}

cat(sprintf("\n--- Timing over %d reps ---\n", nmc))
cat(sprintf("  Wald sandwich:  mean %8.5f s  median %8.5f s  total %7.2f s\n",
            mean(t_wald_v), median(t_wald_v), sum(t_wald_v)))
cat(sprintf("  Mplier boot:    mean %8.5f s  median %8.5f s  total %7.2f s\n",
            mean(t_boot_v), median(t_boot_v), sum(t_boot_v)))

cat("\nDone.\n")
