# Distance between estimated and oracle propensity a0(t, x) = P{W(t)=1 | event at t, X=x}.
# No TV-CSL fit: each model is fit on a training set and evaluated on an independent test set, at
#   "events":  (U_i, X_i) for test events, and
#   "riskset": (t_k, X_j) for every test subject j at risk at t_k, for 50 sampled event times t_k
#              (TV-CSL uses a at every risk-set row); averaged within t_k, then over t_k. Same DGP and seeds as the centering sweep (non-linear eta_0, beta0 = (0,m,m,m)).
# Usage (from project root): Rscript dist_task.R <m> [reps = 20] [n = 4000]
suppressMessages({ source("R/datagen-helper.R"); library(survival); library(mgcv) })
D <- "scripts/TV-CSL/propensity-distance-2026-10-10"
a <- commandArgs(TRUE); m <- as.numeric(a[1])
reps <- if (length(a) >= 2) as.integer(a[2]) else 20L
n    <- if (length(a) >= 3) as.integer(a[3]) else 4000L
gen <- function(seed) generate_simulated_data(n = n, eta_type = "non-linear", HTE_type = "linear",
                                              linear_HTE_multiplier = m, seed_value = seed)

# Oracle, as in R/TV-CSL.R (cox-time-varying-oracle):
# logit a0 = log r + tau + r t - d t^2 + log int_0^t exp(d s^2 - r s) ds, r = exp(x2 + x3),
# d = exp(eta_0) (exp(tau) - 1) / 2.
a0_fun <- function(t, x2, x3, eta0, tau) {
  r <- exp(x2 + x3); d <- 0.5 * exp(eta0) * (exp(tau) - 1); M <- max(0, d * t^2 - r * t)
  I <- integrate(function(s) exp(d * s^2 - r * s - M), 0, t, rel.tol = 1e-8)$value
  plogis(log(r) + tau + r * t - d * t^2 + M + log(I))
}

fit_models <- function(tr) {
  # Event rows of the (tstart, tstop] data: tstop = U, W = 1(A < U).
  ev <- data.frame(W = as.numeric(tr$A < tr$U), tstop = tr$U, tr[, c("X.1", "X.2", "X.3")])[tr$Delta == 1, ]
  # 1. marginal (cox-linear-censored-only): Cox for A on event subjects, a = 1 - exp(-t exp(X alpha))
  dA <- data.frame(U_A = pmin(tr$A, tr$U), Delta_A = tr$A <= tr$U, tr[, c("X.1", "X.2", "X.3")])[tr$Delta == 1, ]
  alpha <- coef(coxph(Surv(U_A, Delta_A) ~ X.1 + X.2 + X.3, data = dA, ties = "breslow"))
  # 2. original (cox-time-varying-prop): additive GAM on event rows
  g0 <- gam(W ~ s(tstop, k = 5) + X.1 + X.2 + X.3, family = binomial, data = ev)
  # 3. new: GAM with time x covariate interactions
  g1 <- gam(W ~ s(tstop, k = 5) + X.1 + X.2 + X.3 + ti(tstop, X.1) + ti(tstop, X.2) + ti(tstop, X.3),
            family = binomial, data = ev, method = "REML")
  list(alpha = alpha, g0 = g0, g1 = g1)
}

out <- list()
for (i in seq_len(reps)) {
  t0 <- Sys.time()
  tr <- gen(1e6 + i); te <- gen(2e6 + i)           # same seeds as the centering sweep
  f <- fit_models(tr)
  predict_all <- function(t, X) list(
    marginal          = 1 - exp(-t * exp(as.vector(as.matrix(X) %*% f$alpha))),
    "gam-additive"    = as.vector(predict(f$g0, data.frame(tstop = t, X), type = "response")),
    "gam-interaction" = as.vector(predict(f$g1, data.frame(tstop = t, X), type = "response")))
  Xc <- c("X.1", "X.2", "X.3")
  # events
  ev <- te[te$Delta == 1, ]
  a0_ev  <- with(ev, mapply(a0_fun, U, X.2, X.3, eta_0, HTE))
  est_ev <- predict_all(ev$U, ev[, Xc])
  # risk sets at 50 sampled event times
  set.seed(3e6 + i)
  tk <- sort(sample(ev$U, 50))
  rs <- do.call(rbind, lapply(seq_along(tk), function(k) { j <- te$U >= tk[k]
    data.frame(k = k, t = tk[k], te[j, c(Xc, "eta_0", "HTE")]) }))
  a0_rs  <- with(rs, mapply(a0_fun, t, X.2, X.3, eta_0, HTE))
  est_rs <- predict_all(rs$t, rs[, Xc])
  per_t  <- function(x) mean(tapply(x, rs$k, mean))          # equal weight per event time
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out[[i]] <- do.call(rbind, lapply(names(est_ev), function(k) rbind(
    data.frame(m = m, i = i, n = n, model = k, where = "events", rows = nrow(ev),
               mae = mean(abs(est_ev[[k]] - a0_ev)), rmse = sqrt(mean((est_ev[[k]] - a0_ev)^2)),
               bias = mean(est_ev[[k]] - a0_ev), secs = secs),
    data.frame(m = m, i = i, n = n, model = k, where = "riskset", rows = nrow(rs),
               mae = per_t(abs(est_rs[[k]] - a0_rs)), rmse = sqrt(per_t((est_rs[[k]] - a0_rs)^2)),
               bias = per_t(est_rs[[k]] - a0_rs), secs = secs))))
  cat(sprintf("m=%s rep %d: %.1f s\n", format(m), i, secs))
}
write.csv(do.call(rbind, out), sprintf("%s/res/dist_n%d_m%s.csv", D, n, format(m)), row.names = FALSE)
cat("DONE m =", m, "\n")
