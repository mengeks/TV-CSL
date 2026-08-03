#!/usr/bin/env Rscript
# MC coverage study for the continuous-case Semenova example
# (Section 4 of TV-CSL-semenova-inference-MRE.tex)
#
# DGP : tau(x) = 1 + x1,  basis p(x) = (1, x1, x2)
#        => true beta = (1, 1, 0)          <-- this answers the "0? or 1?" question
#
# Outputs:
#   - Smoke test (5-10 runs): bhat and per-method timing
#   - MC simulation (300 reps, n=4000, B=500 bootstrap):
#       bias, true SE, mean estimated SE, Wald coverage, bootstrap coverage
#       for intercept (beta_0=1) and slope (beta_1=1)

# ── one MC replication ─────────────────────────────────────────────────────────
one_rep <- function(n = 4000, B_boot = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  ## DGP
  X   <- matrix(rnorm(n * 2), n, 2)
  e   <- plogis(0.7 * X[, 1])
  W   <- rbinom(n, 1, e)
  tau <- 1 + X[, 1]             # true CATE  => true beta = (1, 1, 0)
  mu0 <- X[, 2]^2
  Y   <- mu0 + W * tau + rnorm(n)

  ## Step A: 2-fold cross-fit nuisances -> AIPW signal Gamma
  fold  <- sample(rep(1:2, length.out = n))
  ehat  <- mhat1 <- mhat0 <- numeric(n)
  for (k in 1:2) {
    tr <- fold != k
    ev <- fold == k
    Xev <- X[ev, ]

    fit_e <- glm(W[tr] ~ X[tr, ], family = binomial)
    ehat[ev] <- as.numeric(cbind(1, Xev) %*% coef(fit_e))
    ehat[ev] <- plogis(ehat[ev])   # apply logistic link

    idx1 <- which(tr & W == 1)
    fit1 <- lm(Y[idx1] ~ X[idx1, ])
    mhat1[ev] <- as.numeric(cbind(1, Xev) %*% coef(fit1))

    idx0 <- which(tr & W == 0)
    fit0 <- lm(Y[idx0] ~ X[idx0, ])
    mhat0[ev] <- as.numeric(cbind(1, Xev) %*% coef(fit0))
  }
  ehat  <- pmin(pmax(ehat, 0.01), 0.99)
  Gamma <- mhat1 - mhat0 +
           W * (Y - mhat1) / ehat -
           (1 - W) * (Y - mhat0) / (1 - ehat)

  ## Step B: OLS of Gamma on p(X) = [1, x1, x2]
  P    <- cbind(1, X)
  bhat <- solve(crossprod(P), crossprod(P, Gamma))
  res  <- as.numeric(Gamma - P %*% bhat)
  Q    <- crossprod(P) / n
  Qinv <- solve(Q)

  ## Step C-i: Wald sandwich SE (timed)
  t0_wald <- proc.time()["elapsed"]
  Sig     <- crossprod(P * res) / n
  Omega   <- Qinv %*% Sig %*% Qinv / n
  se_wald <- sqrt(diag(Omega))
  t_wald  <- proc.time()["elapsed"] - t0_wald

  ## Step C-ii: multiplier bootstrap (vectorized, timed)
  ## boot_draws[b, ] ~ bhat - beta asymptotically
  t0_boot    <- proc.time()["elapsed"]
  XI         <- matrix(rnorm(n * B_boot), n, B_boot)   # n x B_boot
  boot_draws <- t(Qinv %*% (t(P) %*% (XI * res) / n))  # B_boot x 3
  ci_boot    <- apply(abs(boot_draws), 2, quantile, 0.95)  # 95th pct margin
  t_boot     <- proc.time()["elapsed"] - t0_boot

  list(
    bhat     = as.numeric(bhat),
    se_wald  = se_wald,
    ci_boot  = ci_boot,    # half-width of 95% bootstrap CI for each coef
    t_wald   = t_wald,
    t_boot   = t_boot
  )
}

# ── Smoke test (10 runs) ───────────────────────────────────────────────────────
cat("\n========== SMOKE TEST (10 runs, n=4000, B=500) ==========\n")
cat(sprintf("True beta = (1, 1, 0) — recovering tau(x) = 1 + x1\n\n"))
cat(sprintf("%-5s  %8s  %8s  %8s  %12s  %12s\n",
            "run", "bhat[0]", "bhat[1]", "bhat[2]", "t_wald(s)", "t_boot(s)"))
cat(strrep("-", 62), "\n")
for (i in seq_len(10)) {
  r <- one_rep(n = 4000, B_boot = 500, seed = i)
  cat(sprintf("%-5d  %8.4f  %8.4f  %8.4f  %12.5f  %12.5f\n",
              i, r$bhat[1], r$bhat[2], r$bhat[3], r$t_wald, r$t_boot))
}

# ── MC simulation (300 reps) ───────────────────────────────────────────────────
nmc       <- 300
beta_true <- c(1, 1, 0)

cat(sprintf("\n========== MC SIMULATION (%d reps, n=4000, B=500) ==========\n", nmc))
mc_results <- vector("list", nmc)
for (i in seq_len(nmc)) {
  if (i == 1 || i %% 50 == 0) cat(sprintf("  rep %3d / %d\n", i, nmc))
  mc_results[[i]] <- one_rep(n = 4000, B_boot = 500, seed = 1000 + i)
}

## Collect
bhat_mat   <- do.call(rbind, lapply(mc_results, `[[`, "bhat"))
se_wald_m  <- do.call(rbind, lapply(mc_results, `[[`, "se_wald"))
ci_boot_m  <- do.call(rbind, lapply(mc_results, `[[`, "ci_boot"))
t_wald_v   <- sapply(mc_results, `[[`, "t_wald")
t_boot_v   <- sapply(mc_results, `[[`, "t_boot")

## Coverage
#   Wald:       CI = bhat_j ± 1.96 * se_wald_j
#   Bootstrap:  CI = bhat_j ± ci_boot_j  (95th pct of |boot draw_j|)
cover_wald <- abs(bhat_mat - matrix(beta_true, nmc, 3, byrow = TRUE)) <
              1.96 * se_wald_m
cover_boot <- abs(bhat_mat - matrix(beta_true, nmc, 3, byrow = TRUE)) <
              ci_boot_m

## Print results
coef_names <- c("Intercept (beta_0 = 1)", "x1 slope (beta_1 = 1)", "x2 slope (beta_2 = 0)")
cat("\n--- Per-coefficient summary ---\n")
cat(sprintf("%-26s  %8s  %8s  %10s  %10s  %10s  %10s\n",
            "Coefficient", "Bias", "True SE", "MeanSE_wald",
            "Cov_wald", "CIhalf_boot", "Cov_boot"))
cat(strrep("-", 90), "\n")
for (j in 1:3) {
  bias        <- mean(bhat_mat[, j]) - beta_true[j]
  true_se     <- sd(bhat_mat[, j])
  mean_se_w   <- mean(se_wald_m[, j])
  mean_ci_b   <- mean(ci_boot_m[, j])
  cov_w       <- mean(cover_wald[, j])
  cov_b       <- mean(cover_boot[, j])
  cat(sprintf("%-26s  %+8.4f  %8.4f  %10.4f  %10.3f  %10.4f  %10.3f\n",
              coef_names[j], bias, true_se, mean_se_w, cov_w, mean_ci_b, cov_b))
}

## Timing
cat(sprintf("\n--- Timing over %d reps ---\n", nmc))
cat(sprintf("  Wald sandwich:  mean %8.5f s  median %8.5f s  total %7.2f s\n",
            mean(t_wald_v), median(t_wald_v), sum(t_wald_v)))
cat(sprintf("  Mplier boot:    mean %8.5f s  median %8.5f s  total %7.2f s\n",
            mean(t_boot_v), median(t_boot_v), sum(t_boot_v)))

cat("\nDone.\n")
