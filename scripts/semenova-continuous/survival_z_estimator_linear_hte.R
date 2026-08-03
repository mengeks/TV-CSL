#!/usr/bin/env Rscript
# Semenova-style orthogonal Z-estimation for heterogeneous proportional-hazard effects
#
# DGP: Cox PH with baseline (not time-varying) treatment
#   lambda_i(t) = lambda0(t) * exp(eta0(X_i) + W_i * tau(X_i))
#   tau(X) = p(X)' beta_0,   p(X) = (1, X1, X2)
#   true beta_0 = (1, 1, 0)  =>  tau(X) = 1 + X1
#   eta0(X) = 0.4 * X2,  lambda0 = 0.05 (constant baseline)
#   W | X ~ Bernoulli(logistic(0.7 * X1))
#   C ~ Exp(0.03)
#
# Method: orthogonal series Z-estimator (Semenova PH analog)
#   Solve: G(beta) = (1/n) sum_i (W_i - ehat_i) p(X_i) dMhat_i(beta) = 0
#   dMhat_i(beta) = Delta_i - exp(etahat_i + W_i p(X_i)'beta) * Lambda0hat(U_i; beta)
#   Lambda0hat = Breslow baseline cumulative hazard given beta.
#
# Nuisance estimation (2-fold cross-fit):
#   ehat(X)   <- logistic regression of W on X   [uses all training data]
#   etahat(X) <- Cox LP from W=0 arm only        [identifies eta0 without W confounding]
#
# Semenova analogy:
#   continuous: regress orthogonal signal Gamma on p(X)
#   PH survival: solve orthogonal martingale score projected on p(X)
#
# Inference:
#   Wald:      Var(betahat) ~ J^{-1} Sigma J^{-T} / n,  SE = sqrt(diag(Var))
#   Uniform:   multiplier bootstrap on phi_i = -J^{-1} psi_i

suppressPackageStartupMessages(library(survival))

# ── Breslow cumulative baseline hazard ─────────────────────────────────────────
# Returns Lambda0(U_i; beta) for each i.
# Risk: r_i(beta) = exp(etahat_i + W_i * P_i' beta).
# Breslow: dLambda0(tj) = d_j / sum_{m: U_m >= tj} r_m(beta).
compute_Lambda0 <- function(beta, U, Del, etahat, W, P) {
  risk     <- exp(etahat + W * as.numeric(P %*% beta))
  ord      <- order(U)
  risk_set <- rev(cumsum(rev(risk[ord])))
  dL_ord   <- Del[ord] / risk_set
  L0_ord   <- cumsum(dL_ord)
  L0_ord[order(ord)]
}

# ── Estimating equation G(beta) ────────────────────────────────────────────────
# G(beta) = (1/n) sum_i (W_i - ehat_i) p(X_i) [Delta_i - r_i * Lambda0_i]
# At beta_0 with consistent nuisances: E[G(beta_0)] = 0 (orthogonal martingale score).
compute_G <- function(beta, U, Del, W, P, ehat, etahat) {
  risk    <- exp(etahat + W * as.numeric(P %*% beta))
  Lambda0 <- compute_Lambda0(beta, U, Del, etahat, W, P)
  M_res   <- Del - risk * Lambda0
  colMeans((W - ehat) * P * M_res)
}

# ── Numerical Jacobian (central differences) ───────────────────────────────────
compute_Jacobian <- function(beta, G_fn, eps = 1e-5) {
  q <- length(beta)
  J <- matrix(0, q, q)
  for (j in seq_len(q)) {
    bp <- bm <- beta
    bp[j] <- bp[j] + eps
    bm[j] <- bm[j] - eps
    J[, j] <- (G_fn(bp) - G_fn(bm)) / (2 * eps)
  }
  J
}

# ── One MC replication ─────────────────────────────────────────────────────────
one_rep_surv <- function(n = 4000, B_boot = 500, seed = NULL, run_boot = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  ## ── DGP ────────────────────────────────────────────────────────────────────
  X    <- matrix(rnorm(n * 2), n, 2)
  e    <- plogis(0.7 * X[, 1])
  W    <- rbinom(n, 1, e)
  tau  <- 1 + X[, 1]         # true CATE (log-HR scale): beta_0 = (1, 1, 0)
  eta0 <- 0.4 * X[, 2]
  lam0 <- 0.05
  Tevt <- rexp(n, rate = lam0 * exp(eta0 + W * tau))
  Cens <- rexp(n, rate = 0.03)
  U    <- pmin(Tevt, Cens)
  Del  <- as.integer(Tevt <= Cens)

  P <- cbind(1, X)   # basis p(X) = (1, X1, X2)
  q <- ncol(P)

  ## ── Step A: 2-fold cross-fit nuisances ────────────────────────────────────
  # ehat  <- logistic regression on all training data
  # etahat <- Cox LP fitted on W=0 arm only (identifies eta0(X) without treatment confounding)
  fold   <- sample(rep(1:2, length.out = n))
  ehat   <- numeric(n)
  etahat <- numeric(n)
  for (k in 1:2) {
    tr  <- fold != k
    ev  <- fold == k
    dat_tr_all <- data.frame(W   = W[tr],  X1 = X[tr, 1], X2 = X[tr, 2],
                             U   = U[tr],  Del = Del[tr])
    dat_ev     <- data.frame(X1  = X[ev, 1], X2 = X[ev, 2])

    # Propensity (all training data)
    fit_e    <- glm(W ~ X1 + X2, family = binomial, data = dat_tr_all)
    ehat[ev] <- predict(fit_e, newdata = dat_ev, type = "response")

    # eta0: Cox on W=0 arm only -> LP identifies eta0(X) = 0.4*X2
    ctrl_tr    <- tr & W == 0
    dat_ctrl   <- data.frame(X1  = X[ctrl_tr, 1], X2 = X[ctrl_tr, 2],
                             U   = U[ctrl_tr],     Del = Del[ctrl_tr])
    fit_cox    <- coxph(Surv(U, Del) ~ X1 + X2, data = dat_ctrl)
    etahat[ev] <- predict(fit_cox, newdata = dat_ev, type = "lp")
  }
  ehat <- pmin(pmax(ehat, 0.02), 0.98)

  G_fn <- function(beta) compute_G(beta, U, Del, W, P, ehat, etahat)

  ## ── Step B: damped Newton-Raphson to solve G(beta) = 0 ────────────────────
  beta   <- rep(0, q)
  n_iter <- 0L
  for (iter in seq_len(100)) {
    g <- G_fn(beta)
    if (max(abs(g)) < 1e-10) break
    J    <- compute_Jacobian(beta, G_fn)
    step <- tryCatch(solve(J, g), error = function(e) rep(NA_real_, q))
    if (anyNA(step)) break
    step_max <- max(abs(step))
    if (step_max > 1) step <- step / step_max
    beta   <- beta - step
    n_iter <- iter
  }
  beta_hat <- beta

  ## ── Step C: score contributions and Jacobian at beta_hat ──────────────────
  risk_hat    <- exp(etahat + W * as.numeric(P %*% beta_hat))
  Lambda0_hat <- compute_Lambda0(beta_hat, U, Del, etahat, W, P)
  M_res_hat   <- Del - risk_hat * Lambda0_hat
  psi_mat     <- (W - ehat) * P * M_res_hat    # n x q: psi_i in rows

  J_hat <- compute_Jacobian(beta_hat, G_fn)
  Jinv  <- solve(J_hat)

  ## ── Step C-i: Wald sandwich SE ─────────────────────────────────────────────
  t0_wald <- proc.time()["elapsed"]
  Sigma   <- crossprod(psi_mat) / n
  Omega   <- Jinv %*% Sigma %*% t(Jinv) / n    # Var(betahat)
  se_wald <- sqrt(diag(Omega))
  t_wald  <- proc.time()["elapsed"] - t0_wald

  ## ── Step C-ii: multiplier bootstrap ────────────────────────────────────────
  # Influence function: phi_i = -J^{-1} psi_i
  # Bootstrap draw:    (1/n) sum_i xi_i phi_i ~ beta_hat - beta_0
  if (run_boot) {
    t0_boot     <- proc.time()["elapsed"]
    XI          <- matrix(rnorm(n * B_boot), n, B_boot)
    boot_scores <- t(psi_mat) %*% XI / n       # q x B_boot
    boot_draws  <- t(-Jinv %*% boot_scores)    # B_boot x q
    ci_boot     <- apply(abs(boot_draws), 2, quantile, 0.95)
    t_boot      <- proc.time()["elapsed"] - t0_boot
  } else {
    ci_boot <- rep(NA_real_, q)
    t_boot  <- 0
  }

  list(
    beta_hat   = beta_hat,
    se_wald    = se_wald,
    ci_boot    = ci_boot,
    t_wald     = t_wald,
    t_boot     = t_boot,
    event_rate = mean(Del),
    n_iter     = n_iter
  )
}

# ── Pilot: n = 100,000 — verify Z-estimator recovers beta_0 = (1, 1, 0) ───────
cat("\n========== PILOT: n=100 000, seed=99 (no bootstrap) ==========\n")
cat("True beta_0 = (1, 1, 0).  Z-estimator should recover this directly.\n\n")
pilot <- one_rep_surv(n = 100000, run_boot = FALSE, seed = 99)
cat(sprintf("beta_hat:        (%.5f,  %.5f,  %.5f)\n",
            pilot$beta_hat[1], pilot$beta_hat[2], pilot$beta_hat[3]))
cat(sprintf("Event rate:       %.3f\n", pilot$event_rate))
cat(sprintf("Newton iters:     %d\n",   pilot$n_iter))
cat(sprintf("beta[1]/beta[0] = %.4f  (expect 1.000)\n",
            pilot$beta_hat[2] / pilot$beta_hat[1]))
cat(sprintf("beta[2]/beta[0] = %.4f  (expect 0.000)\n",
            pilot$beta_hat[3] / pilot$beta_hat[1]))

# ── Smoke test: 10 runs, n = 4000 ─────────────────────────────────────────────
cat("\n========== SMOKE TEST (10 runs, n=4000, B=500) ==========\n")
cat("True beta_0 = (1, 1, 0)\n\n")
cat(sprintf("%-5s  %8s  %8s  %8s  %5s  %10s  %10s\n",
            "run", "bhat[0]", "bhat[1]", "bhat[2]", "iter",
            "t_wald(s)", "t_boot(s)"))
cat(strrep("-", 68), "\n")
for (i in seq_len(10)) {
  r <- one_rep_surv(n = 4000, B_boot = 500, seed = i, run_boot = TRUE)
  cat(sprintf("%-5d  %8.4f  %8.4f  %8.4f  %5d  %10.5f  %10.5f\n",
              i, r$beta_hat[1], r$beta_hat[2], r$beta_hat[3],
              r$n_iter, r$t_wald, r$t_boot))
}

# ── MC simulation: 300 reps, n = 4000 ─────────────────────────────────────────
nmc       <- 300
beta_true <- c(1, 1, 0)

cat(sprintf("\n========== MC SIMULATION (%d reps, n=4000, B=500) ==========\n", nmc))
mc_results <- vector("list", nmc)
for (i in seq_len(nmc)) {
  if (i == 1 || i %% 50 == 0) cat(sprintf("  rep %3d / %d\n", i, nmc))
  mc_results[[i]] <- one_rep_surv(n = 4000, B_boot = 500,
                                  seed = 2000 + i, run_boot = TRUE)
}

## Collect
beta_mat  <- do.call(rbind, lapply(mc_results, `[[`, "beta_hat"))
se_wald_m <- do.call(rbind, lapply(mc_results, `[[`, "se_wald"))
ci_boot_m <- do.call(rbind, lapply(mc_results, `[[`, "ci_boot"))
t_wald_v  <- sapply(mc_results, `[[`, "t_wald")
t_boot_v  <- sapply(mc_results, `[[`, "t_boot")
ev_rate_v <- sapply(mc_results, `[[`, "event_rate")
n_iter_v  <- sapply(mc_results, `[[`, "n_iter")

## Coverage against beta_0 = (1, 1, 0)
cover_wald <- abs(beta_mat - matrix(beta_true, nmc, 3, byrow = TRUE)) <
              1.96 * se_wald_m
cover_boot <- abs(beta_mat - matrix(beta_true, nmc, 3, byrow = TRUE)) <
              ci_boot_m

## Summary table
coef_names <- c("beta_0  (true = 1)", "beta_1  (true = 1)", "beta_2  (true = 0)")
cat(sprintf(
  "\n--- Per-coefficient summary  [event rate: %.3f | mean Newton iters: %.1f] ---\n",
  mean(ev_rate_v), mean(n_iter_v)))
cat(sprintf("%-22s  %8s  %8s  %10s  %10s  %10s  %10s\n",
            "Coefficient", "Bias", "True SE", "MeanSE_wald",
            "Cov_wald", "CIhalf_boot", "Cov_boot"))
cat(strrep("-", 92), "\n")
for (j in 1:3) {
  bias      <- mean(beta_mat[, j]) - beta_true[j]
  true_se   <- sd(beta_mat[, j])
  mean_se_w <- mean(se_wald_m[, j])
  mean_ci_b <- mean(ci_boot_m[, j])
  cov_w     <- mean(cover_wald[, j])
  cov_b     <- mean(cover_boot[, j])
  cat(sprintf("%-22s  %+8.4f  %8.4f  %10.4f  %10.3f  %10.4f  %10.3f\n",
              coef_names[j], bias, true_se, mean_se_w, cov_w, mean_ci_b, cov_b))
}

cat(sprintf("\n--- Timing over %d reps ---\n", nmc))
cat(sprintf("  Wald sandwich:  mean %8.5f s  median %8.5f s  total %7.2f s\n",
            mean(t_wald_v), median(t_wald_v), sum(t_wald_v)))
cat(sprintf("  Mplier boot:    mean %8.5f s  median %8.5f s  total %7.2f s\n",
            mean(t_boot_v), median(t_boot_v), sum(t_boot_v)))

cat("\nDone.\n")
