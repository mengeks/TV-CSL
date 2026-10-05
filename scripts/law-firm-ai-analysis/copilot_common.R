## copilot_common.R
## Cohort construction, data preparation and estimator wrappers shared by
## copilot_analysis_tvcsl_lib.R and copilot_data_quality.R.
##
## Cohorts
##   "new_user": no recorded CC Assist use before time zero (cca_pre_washout == 0)
##   "all"     : every lawyer not already a regular user at entry
## Both drop lawyers whose follow-up has zero length (time_to_regular ==
## entry_week, i.e. regular user in their entry week).
##
## Tie conventions for adoption in the event week (A == U)
##   "lt": W(t) = 1(A < t)  -> adoption in the event week counts as after it
##   "le": W(t) = 1(A <= t) -> adoption in the event week counts as before it
## "le" is implemented by moving adoption to A - 0.5 (weeks are integers).

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

COHORT_LABELS <- c(new_user = "New users", all = "All lawyers")
TIES_LABELS   <- c(lt = "$W(t) = 1(A < t)$", le = "$W(t) = 1(A \\le t)$")

BASE_COVARIATES <- c("rank", "office", "practice_group",
                     "baseline_weekly_hours", "cca_pre_washout")

cohort_covariates <- function(cohort) {
  # cca_pre_washout is identically 0 in the new-user cohort
  if (cohort == "new_user") setdiff(BASE_COVARIATES, "cca_pre_washout")
  else BASE_COVARIATES
}

# Rows of df_raw in the cohort, before the zero-follow-up exclusion
cohort_rows <- function(df_raw, cohort) {
  keep <- df_raw$already_regular_at_entry == 0
  if (cohort == "new_user") keep <- keep & df_raw$cca_pre_washout == 0
  df_raw[keep, , drop = FALSE]
}

# Person-level analysis dataset for one cohort
prepare_cohort <- function(df_raw, cohort) {
  cohort_rows(df_raw, cohort) %>%
    filter(time_to_regular > entry_week) %>%   # exclude zero-length follow-up
    mutate(
      id        = row_number(),
      is_junior = as.integer(seniority == "Junior"),
      baseline_hours_missing = as.integer(is.na(baseline_weekly_hours)),
      baseline_weekly_hours = ifelse(
        is.na(baseline_weekly_hours),
        mean(baseline_weekly_hours, na.rm = TRUE),
        baseline_weekly_hours
      ),
      rank           = relevel(factor(rank), ref = "Trainees"),
      office         = droplevels(factor(office)),
      practice_group = droplevels(factor(practice_group)),
      U     = time_to_regular,
      A     = ifelse(is.na(first_cca_use_week), Inf, first_cca_use_week),
      Delta = event_regular
    )
}

# Sample-flow counts for both cohorts, plus descriptors of the analysed sets
cohort_counts <- function(df_raw) {
  n_all     <- nrow(df_raw)
  n_reg     <- sum(df_raw$already_regular_at_entry == 1)
  flow <- lapply(names(COHORT_LABELS), function(ch) {
    pre   <- cohort_rows(df_raw, ch)
    zero  <- sum(pre$time_to_regular <= pre$entry_week)
    d     <- prepare_cohort(df_raw, ch)
    tibble::tibble(
      cohort = ch,
      n_in_file                    = n_all,
      excl_already_regular         = n_reg,
      excl_prior_use               = if (ch == "new_user")
                                       sum(df_raw$already_regular_at_entry == 0 &
                                           df_raw$cca_pre_washout == 1) else 0L,
      n_cohort_before_zero_fu      = nrow(pre),
      excl_first_week_regulars     = zero,
      n_analysed                   = nrow(d),
      events                       = sum(d$Delta),
      adopted_by_entry             = sum(d$A <= d$entry_week),
      adopted_during_followup      = sum(d$A > d$entry_week & d$A < d$U),
      adopted_in_event_week        = sum(d$A == d$U & d$Delta == 1),
      adopted_after_followup_or_never = sum(d$A > d$U | !is.finite(d$A)),
      baseline_hours_missing       = sum(d$baseline_hours_missing)
    )
  })
  bind_rows(flow)
}

# Dummy-code covariates; columns renamed X.1, ..., X.p (TV-CSL.R convention).
# Reference levels: rank = Trainees, office/practice_group = first level.
add_design <- function(d, covariates) {
  mm <- model.matrix(reformulate(covariates), data = d)[, -1, drop = FALSE]
  xvars <- paste0("X.", seq_len(ncol(mm)))
  labels <- setNames(colnames(mm), xvars)
  colnames(mm) <- xvars
  list(data = cbind(d, mm), xvars = xvars, labels = labels)
}

# Adds U_A / Delta_A for the treatment-time model and builds the
# counting-process data with delayed entry, under tie convention `ties`.
make_tv <- function(d, ties = c("lt", "le")) {
  ties <- match.arg(ties)
  d <- d %>%
    mutate(
      A_eff   = if (ties == "lt") A else A - 0.5,
      U_A     = pmin(A_eff, U),
      Delta_A = as.integer(A_eff < U)
    )
  tv <- tmerge(d %>% select(-Delta), d, id = id,
               Delta = event(U, Delta), tstart = entry_week, tstop = U)
  tv <- tmerge(tv, d %>% filter(is.finite(A_eff)), id = id, W = tdc(A_eff))
  tv <- as.data.frame(tv) %>%
    mutate(W = as.integer(coalesce(W, 0L)), Delta = as.integer(Delta)) %>%
    filter(tstop > tstart) %>%
    arrange(id, tstart)
  attr(tv, "tcount") <- NULL
  attr(tv, "tm.retain") <- NULL
  list(orig = d, tv = tv)
}

# ── Estimators (all from R/TV-CSL.R; caller must source it) ──────────────────
# Each fit returns beta = (senior effect, junior - senior) and its 2x2 covariance.

METHODS <- c(
  fixed_cox  = "Cox, ignoring treatment time",
  tv_cox     = "Time-varying Cox",
  tvcsl_marg = "TV-CSL, marginal propensity",
  tvcsl_tvp  = "TV-CSL, time-varying propensity"
)
EFFECTS <- c(senior = "Used CC Assist (senior)",
             junior = "Used CC Assist (junior)",
             diff   = "Junior $-$ senior difference")

# TV-CSL.R functions print intermediate fits; keep the console readable.
quietly <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- expr))
  out
}

# S_lasso builds
#   Surv(tstart, tstop, Delta) ~ X.1 + ... + X.p + is_junior + W * (1 + is_junior)
# is_junior is collinear with the rank dummies (NA coefficient); listing it in
# outcome_vars keeps beta_HTE = (W, W:is_junior).
fit_cox_lib <- function(train_data, xvars) {
  fit <- quietly(S_lasso(
    train_data       = train_data,
    test_data        = train_data,
    regressor_spec   = "linear",
    HTE_spec         = "linear",
    outcome_vars     = c(xvars, "is_junior"),
    effect_modifiers = "is_junior"
  ))
  list(beta = fit$beta_HTE, V = fit$vcov_HTE)
}

# Nuisance nu: m-regression = Cox on X ignoring treatment (eta_0 = eta_1).
# Final model: coxph on (W - a(t,X)) * (1, is_junior) with offset(nu),
# cross-fitted over K folds; sandwich SE from stacked subject-level scores.
fit_tvcsl_lib <- function(dat, xvars, prop_score_spec, K) {
  fit <- quietly(TV_CSL(
    train_data          = dat$tv,
    test_data           = dat$orig,
    train_data_original = dat$orig,
    HTE_type            = NA,
    eta_type            = NA,
    K                   = K,
    prop_score_spec     = prop_score_spec,
    lasso_type          = "m-regression",
    regressor_spec      = "linear",
    final_model_method  = "lasso_coxph",
    HTE_spec            = "linear",
    lasso_warmstart     = 0,
    outcome_vars        = xvars,
    treatment_vars      = xvars,
    effect_modifiers    = "is_junior",
    verbose             = 0
  ))
  # Omega_hat is the asymptotic covariance of sqrt(n) * beta; n = subjects.
  list(beta = fit$beta_HTE, V = fit$Omega_hat / nrow(dat$orig))
}

# Long table of senior / junior / difference effects for one fit
effects_long <- function(fit) {
  contr <- list(senior = c(1, 0), junior = c(1, 1), diff = c(0, 1))
  bind_rows(lapply(names(contr), function(e) {
    cc  <- contr[[e]]
    est <- sum(cc * fit$beta)
    se  <- sqrt(drop(t(cc) %*% fit$V %*% cc))
    tibble::tibble(effect = e, est = est, se = se,
                   lo = est - 1.96 * se, hi = est + 1.96 * se,
                   p = 2 * pnorm(abs(est / se), lower.tail = FALSE))
  }))
}

# All four methods on a prepared person-level dataset `d` (from prepare_cohort)
run_methods <- function(d, covariates, ties, K = 5L) {
  des   <- add_design(d, covariates)
  dat   <- make_tv(des$data, ties)
  xvars <- des$xvars

  fits <- list(
    fixed_cox  = fit_cox_lib(
      dat$orig %>% mutate(tstart = entry_week, tstop = U, W = used_cca), xvars),
    tv_cox     = fit_cox_lib(dat$tv, xvars),
    tvcsl_marg = fit_tvcsl_lib(dat, xvars, "cox-linear-censored-only-breslow", K),
    tvcsl_tvp  = fit_tvcsl_lib(dat, xvars, "cox-time-varying-prop", K)
  )
  bind_rows(lapply(names(fits), function(m) {
    effects_long(fits[[m]]) %>% mutate(method = m, .before = 1)
  })) %>%
    mutate(ties = ties, n = nrow(dat$orig), .before = 1)
}
