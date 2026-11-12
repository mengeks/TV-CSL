## stanford_common.R
## Data preparation, fold checks and estimator wrappers shared by
## stanford_analysis.R (real data) and stanford_plasmode.R (semi-synthetic).
## Mirrors scripts/law-firm-ai-analysis/copilot_common.R; every estimator
## comes from R/TV-CSL.R, which the caller must source.
##
## Covariates (columns renamed X.1, ..., X.p, the TV-CSL.R convention):
##   X.1 age at acceptance (standardized)   adjustment + effect modifier
##   X.2 year of acceptance (standardized)  adjustment + effect modifier
##   X.3 prior surgery (0/1)                adjustment only
## Treatment W(t) = 1(A < t), A = days from acceptance to transplant.

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

STUDY_START <- as.Date("1967-10-01")
STUDY_CLOSE <- as.Date("1974-04-01")   # administrative censoring date

XVARS     <- c("X.1", "X.2", "X.3")
MODIFIERS <- c("X.1", "X.2")

METHODS <- c(
  fixed_cox  = "Cox, ignoring transplant time",
  tv_cox     = "Time-varying Cox",
  tvcsl_marg = "TV-CSL, marginal propensity",
  tvcsl_tvp  = "TV-CSL, time-varying propensity"
)
EFFECTS <- c(trt  = "Transplant",
             age  = "Age $\\times$ transplant",
             year = "Year $\\times$ transplant")

# Person-level data from survival::jasa. Time conventions follow
# process-stanford-RHC-data.R: follow-up of at least 0.5 days, and a
# transplant on the day of death or censoring is moved 0.5 days earlier.
load_stanford <- function() {
  survival::jasa %>%
    mutate(year = as.numeric(accept.dt - STUDY_START) / 365.25,
           fu   = as.numeric(fu.date - accept.dt),
           tx   = as.numeric(tx.date - accept.dt)) %>%
    transmute(
      id      = row_number(),
      age, year, surgery,
      U       = pmax(0.5, fu),
      A       = case_when(is.na(tx)            ~ Inf,
                          tx.date == fu.date   ~ tx - 0.5,
                          TRUE                 ~ tx),
      Delta   = as.integer(fustat),
      C_admin = as.numeric(STUDY_CLOSE - accept.dt),
      X.1     = as.numeric(scale(age)),
      X.2     = as.numeric(scale(year)),
      X.3     = as.integer(surgery)
    )
}

# Person-level data with U_A / Delta_A plus the counting-process data
make_tv <- function(d) {
  d <- d %>% mutate(U_A     = pmin(A, U),
                    Delta_A = as.integer(A < U))
  tv <- tmerge(d %>% select(-Delta), d, id = id, Delta = event(U, Delta))
  tv <- tmerge(tv, d %>% filter(is.finite(A)) %>% mutate(A0 = pmax(A, 0)),
               id = id, W = tdc(A0))
  tv <- as.data.frame(tv) %>%
    mutate(W = as.integer(coalesce(W, 0L)), Delta = as.integer(Delta)) %>%
    filter(tstop > tstart) %>%
    arrange(id, tstart)
  attr(tv, "tcount") <- NULL
  attr(tv, "tm.retain") <- NULL
  list(orig = d, tv = tv)
}

# TV_CSL() assigns folds as contiguous row blocks and jasa is sorted by
# acceptance date, so rows are shuffled first to avoid calendar-time folds.
shuffle_rows <- function(d) d[sample(nrow(d)), , drop = FALSE]

# Prior-surgery counts in each TV_CSL() training fold (rows in current order).
# The propensity models are fitted on training-fold events, so each needs
# surgery patients with events, both treated and untreated at the event,
# and (marginal propensity) with an observed transplant among those events.
fold_check <- function(d, K) {
  folds <- cut(seq_len(nrow(d)), breaks = K, labels = FALSE)
  treated_event <- d$Delta == 1 & d$A < d$U
  bind_rows(lapply(seq_len(K), function(k) {
    tr <- folds != k
    s  <- tr & d$X.3 == 1
    tibble::tibble(
      fold                      = k,
      n_train                   = sum(tr),
      surgery_train             = sum(s),
      surgery_heldout           = sum(!tr & d$X.3 == 1),
      events_train              = sum(tr & d$Delta == 1),
      surgery_events            = sum(s & d$Delta == 1),
      surgery_treated_events    = sum(s & treated_event),
      surgery_untreated_events  = sum(s & d$Delta == 1 & !treated_event)
    )
  })) %>%
    mutate(ok = surgery_treated_events > 0 & surgery_untreated_events > 0)
}

# ── Estimators ───────────────────────────────────────────────────────────────
# Each fit returns beta = (transplant, age x transplant, year x transplant)
# and covariance matrices.

quietly <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- expr))
  out
}

# S_lasso (linear/linear) = unpenalized coxph:
#   Surv(tstart, tstop, Delta) ~ X.1 + X.2 + X.3 + W * (1 + X.1 + X.2)
fit_cox_lib <- function(train_data) {
  fit <- quietly(S_lasso(
    train_data       = train_data,
    test_data        = train_data,
    regressor_spec   = "linear",
    HTE_spec         = "linear",
    outcome_vars     = XVARS,
    effect_modifiers = MODIFIERS
  ))
  list(beta = fit$beta_HTE, V = fit$vcov_HTE, V_info = fit$vcov_HTE)
}

# Nuisance nu: m-regression (Cox on X ignoring treatment); final model
# coxph on (W - a(t,X)) * (1, X.1, X.2) with offset(nu), cross-fitted.
fit_tvcsl_lib <- function(dat, prop_score_spec, K) {
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
    outcome_vars        = XVARS,
    treatment_vars      = XVARS,
    effect_modifiers    = MODIFIERS,
    verbose             = 0
  ))
  # Omega_hat is the asymptotic covariance of sqrt(n) * beta; n = subjects.
  list(beta = fit$beta_HTE, V = fit$Omega_hat / nrow(dat$orig),
       V_info = fit$vcov_naive)
}

effects_long <- function(fit) {
  tibble::tibble(
    effect  = names(EFFECTS),
    est     = unname(fit$beta),
    se      = sqrt(diag(fit$V)),
    se_info = sqrt(diag(fit$V_info))
  ) %>%
    mutate(lo = est - 1.96 * se, hi = est + 1.96 * se,
           p  = 2 * pnorm(abs(est / se), lower.tail = FALSE))
}

# All four methods on person-level data `d` (rows already in fold order).
# A method that fails returns NA rows rather than stopping the run.
run_methods <- function(d, K = 5L) {
  dat <- make_tv(d)
  fitters <- list(
    fixed_cox  = function() fit_cox_lib(
      dat$orig %>% mutate(tstart = 0, tstop = U, W = Delta_A)),
    tv_cox     = function() fit_cox_lib(dat$tv),
    tvcsl_marg = function() fit_tvcsl_lib(dat, "cox-linear-censored-only-breslow", K),
    tvcsl_tvp  = function() fit_tvcsl_lib(dat, "cox-time-varying-prop", K)
  )
  bind_rows(lapply(names(fitters), function(m) {
    res <- tryCatch(effects_long(fitters[[m]]()), error = function(e) {
      tibble::tibble(effect = names(EFFECTS), est = NA_real_, se = NA_real_,
                     se_info = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
    })
    res %>% mutate(method = m, .before = 1)
  }))
}

fmt_ci <- function(est, lo, hi) sprintf("%.3f (%.3f, %.3f)", est, lo, hi)
fmt_p  <- function(p) ifelse(p < 0.001, "$<$0.001", sprintf("%.3f", p))
