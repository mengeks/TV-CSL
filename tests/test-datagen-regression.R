# tests/test-datagen-regression.R
# Regression tests for data-generation defects found in Oct 2026; each check failed on the old code.
# Run from project root: Rscript tests/test-datagen-regression.R

suppressPackageStartupMessages({
  library(survival)
  library(tidyverse)
})
source("R/datagen-helper.R")
source("scripts/TV-CSL/time-varying-estimate.R")

pass <- function(msg) cat(sprintf("PASS: %s\n", msg))
fail <- function(msg) stop(sprintf("FAIL: %s", msg))

# Exact cumulative hazard of calculate_hazard() at t (linear baseline h0(t) = t).
H_exact <- function(t, d) {
  exp(d$eta_0) * (pmin(t, d$A)^2 / 2 + (t > d$A) * exp(d$HTE) * (t^2 - pmin(t, d$A)^2) / 2)
}

# ---- 1. Event times solve H(T) = E exactly ------------------------------------
# Old: simsurv's 15-node quadrature across the hazard jump at A gave errors in H(T)
# of up to ~70% for subjects treated before T.
{
  covs <- generate_covariates(n = 2000, p = 3, X_distribution = "normal", eta_type = "linear",
                              HTE_type = "linear", linear_HTE_multiplier = 2)
  set.seed(1); E <- rexp(nrow(covs))
  set.seed(1); sim <- simulate_event_times(covs, is_time_varying = TRUE, baseline_type = "linear", maxt = 20)
  ev <- sim$status == 1
  treated <- ev & covs$A < sim$eventtime
  if (sum(treated) < 100) fail("too few subjects treated before T to test the jump at A")
  rel_err <- abs(H_exact(sim$eventtime, covs)[ev] / E[ev] - 1)
  if (max(rel_err) > 1e-8) fail(sprintf("H(T) != E (max relative error %.2g)", max(rel_err)))
  pass("event times invert the exact cumulative hazard, including across the jump at A")
}

# ---- 2. Survivors at maxt are censored, not events ------------------------------
# Old: Delta = (T <= C) with T truncated to maxt = 20 and C = min(20, C_gen) turned every
# survivor with C_gen >= 20 into a tied "event" at t = 20.
{
  dat <- generate_simulated_data(n = 2000, light_censoring = TRUE, eta_type = "linear",
                                 HTE_type = "linear", seed_value = 7)   # C = 20 for everyone
  survivors <- dat$status == 0
  if (sum(survivors) == 0) fail("no subjects survived to maxt; test is vacuous")
  if (any(dat$Delta[survivors] == 1)) fail("subjects alive at maxt were coded as events")
  if (sum(duplicated(dat$U[dat$Delta == 1])) > 0) fail("tied event times found")
  pass(sprintf("%d survivors at maxt censored; no tied event times", sum(survivors)))
}

# ---- 3. The (tstart, tstop] split keeps every event -----------------------------
# Old: create_pseudo_dataset() dropped intervals shorter than 0.001, deleting events that
# occur right after adoption (~7% of treated events when HTE is doubled).
{
  toy <- data.frame(id = 1:4, U = c(1, 1, 1, 2), Delta = c(1, 1, 1, 0),
                    A = c(1 - 1e-4, 1 - 1e-12, 0.5, Inf), X.1 = c(0.1, -0.2, 0.3, 0.4))
  p <- create_pseudo_dataset(toy)
  if (sum(p$Delta) != sum(toy$Delta)) fail("events lost in the split")
  ev_rows <- p[p$Delta == 1, ]
  if (!all(ev_rows$W[ev_rows$id %in% 1:3] == 1)) fail("post-adoption events not on a treated row")
  if (any(p$tstop - p$tstart <= 0)) fail("non-positive interval length")
  fit <- tryCatch(coxph(Surv(tstart, tstop, Delta) ~ W, data = p), error = function(e) e)
  if (inherits(fit, "error")) fail(paste("coxph rejects the split:", conditionMessage(fit)))

  dat <- generate_simulated_data(n = 500, eta_type = "linear", HTE_type = "linear",
                                 linear_HTE_multiplier = 2, seed_value = 3)
  if (sum(create_pseudo_dataset(dat)$Delta) != sum(dat$Delta)) fail("events lost in the split (simulated data)")
  pass("split keeps every event and coxph accepts it")
}

# ---- 4. generate_and_save_data() is reproducible from its seed -------------------
# Old: it passed params$seed_value (NULL) instead of params$seed, so set.seed(NULL)
# re-randomised every dataset.
params <- list(is_time_varying = TRUE, light_censoring = FALSE, lambda_C = 0.1, p = 3,
               X_distribution = "normal", X_cov_type = "identity", tx_difficulty = "simple",
               eta_type = "linear", HTE_type = "linear")
gen_dir <- tempfile(); dir.create(gen_dir)
gen <- function(i, prm = params) {
  invisible(capture.output(generate_and_save_data(i, 200, gen_dir, prm)))
  readRDS(file.path(gen_dir, paste0("sim_data_", i, "_n_200.rds")))$data
}
{
  if (!identical(gen(1), gen(1))) fail("same i gave different datasets")
  if (identical(gen(1), gen(2))) fail("different i gave identical datasets")
  pass("generate_and_save_data reproducible from i")
}

# ---- 5. JSON datagen params reach generate_simulated_data() ----------------------
# Old: linear_HTE_multiplier (and the eta_0 linear parameters) were silently ignored.
{
  d <- gen(3, modifyList(params, list(linear_HTE_multiplier = 2)))
  if (!isTRUE(all.equal(d$HTE, 2 * (d$X.1 + d$X.2 + d$X.3)))) fail("linear_HTE_multiplier ignored")
  pass("linear_HTE_multiplier passed through from params")
}

cat("All datagen regression tests passed.\n")
