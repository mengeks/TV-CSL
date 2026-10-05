## copilot_analysis_tvcsl_lib.R
## Same analysis, input and output tables as copilot_analysis.R, but every
## estimator comes from R/TV-CSL.R instead of being written by hand:
##   - Cox fits (fixed and time-varying)  -> S_lasso(regressor_spec = "linear",
##                                           HTE_spec = "linear"), the unpenalized
##                                           coxph path with W * (1 + is_junior)
##   - TV-CSL (both propensity choices)   -> TV_CSL(lasso_type = "m-regression")
##       marginal propensity      = prop_score_spec "cox-linear-censored-only"
##       time-varying propensity  = prop_score_spec "cox-time-varying-prop"
##
## Produces:
##   <OUTPUT_DIR>/tab_copilot_ignore_time.tex
##   <OUTPUT_DIR>/tab_copilot_tvcsl.tex

# ── User configuration ────────────────────────────────────────────────────────

INPUT_FILE  <- here::here("scripts/law-firm-ai-analysis/survival_ccassist_SYNTHETIC.csv")
OUTPUT_DIR  <- here::here("scripts/law-firm-ai-analysis/tables-tvcsl-lib")
K_FOLDS     <- 5L

# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(survival)
  library(mgcv)
  library(here)
  source(here::here("R/TV-CSL.R"))
  source(here::here("R/data-handler.R"))
})

# TV-CSL.R functions print intermediate fits; keep the console readable.
quietly <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- expr))
  out
}

# ── 1. Load data ──────────────────────────────────────────────────────────────

df_raw <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

# ── 2. Prepare person-level dataset ──────────────────────────────────────────

df_base <- df_raw %>%
  filter(already_regular_at_entry == 0,
         time_to_regular > entry_week) %>%  # exclude zero-length follow-up
  mutate(
    id        = row_number(),
    is_junior = as.integer(seniority == "Junior"),
    baseline_weekly_hours = ifelse(
      is.na(baseline_weekly_hours),
      mean(baseline_weekly_hours, na.rm = TRUE),
      baseline_weekly_hours
    ),
    U     = time_to_regular,
    A     = ifelse(is.na(first_cca_use_week), Inf, first_cca_use_week),
    Delta = event_regular,
    U_A   = pmin(A, U),
    Delta_A = as.integer(A <= U & is.finite(A))
  )

# Adjustment covariates. A reference rank level is dropped so the design has
# no column aliased with the Cox baseline hazard: TV-CSL.R multiplies X by the
# coefficient vector directly, so an NA coefficient would turn nu/a(t,X) into
# NA. Trainees is the reference because that is the dummy coxph aliases (NA)
# in the `- 1` design of copilot_analysis.R. The choice matters for the
# marginal propensity a = 1 - exp(-t * exp(X alpha)), which has no baseline
# rate and so is not invariant to the reference level.
# Columns are renamed X.1, ..., X.p, the TV-CSL.R naming convention
# (m_regression() predicts on the test fold using the "X." columns).
mm <- model.matrix(
  ~ relevel(factor(rank), ref = "Trainees") + office + practice_group +
    baseline_weekly_hours + cca_pre_washout,
  data = df_base
)[, -1]
xvars <- paste0("X.", seq_len(ncol(mm)))
xlabels <- setNames(colnames(mm), xvars)
colnames(mm) <- xvars

df_orig <- cbind(df_base, mm)

# ── 3. Counting-process dataset (delayed entry handled) ───────────────────────
# W_i(t) = 1(A_i < t); convention: same-week event and treatment → untreated.
# create_pseudo_dataset() in TV-CSL.R starts every subject at t = 0, so the
# delayed-entry split uses survival::tmerge; TV_CSL() accepts any data that
# already has tstart/tstop.

df_tv <- tmerge(df_orig %>% select(-Delta), df_orig, id = id,
                Delta = event(U, Delta), tstart = entry_week, tstop = U)
df_tv <- tmerge(df_tv, df_orig %>% filter(is.finite(A)), id = id, W = tdc(A))
df_tv <- as.data.frame(df_tv) %>%
  mutate(W = as.integer(coalesce(W, 0L)), Delta = as.integer(Delta)) %>%
  filter(tstop > tstart) %>%
  arrange(id, tstart)
attr(df_tv, "tcount") <- NULL
attr(df_tv, "tm.retain") <- NULL

# ── 4–5. Cox fits via S_lasso (linear/linear → unpenalized coxph) ─────────────
# Formula built by S_lasso:
#   Surv(tstart, tstop, Delta) ~ X.1 + ... + X.p + is_junior + W * (1 + is_junior)
# is_junior is collinear with the rank dummies (NA coefficient); listing it in
# outcome_vars keeps beta_HTE = (W, W:is_junior).

fit_cox_lib <- function(train_data) {
  quietly(S_lasso(
    train_data       = train_data,
    test_data        = train_data,
    regressor_spec   = "linear",
    HTE_spec         = "linear",
    outcome_vars     = c(xvars, "is_junior"),
    effect_modifiers = "is_junior"
  ))
}

# Fixed treatment: used_cca is the ever-treated indicator over (entry_week, U]
fit_fixed <- fit_cox_lib(
  df_orig %>% mutate(tstart = entry_week, tstop = U, W = used_cca)
)
fit_tvcox <- fit_cox_lib(df_tv)

# ── 6. Format estimates ───────────────────────────────────────────────────────
# Each takes beta = (senior, junior - senior) and its 2x2 covariance.

fmt_est <- function(est, se_est, with_p = TRUE) {
  ci <- sprintf("%.3f (%.3f, %.3f)",
                est, est - 1.96 * se_est, est + 1.96 * se_est)
  if (!with_p) return(ci)
  pv <- 2 * pnorm(abs(est / se_est), lower.tail = FALSE)
  list(ci = ci, pval = if (pv < 0.001) "$<$0.001" else sprintf("%.3f", pv))
}

group_effects <- function(beta, V, with_p = TRUE) {
  c_jun <- c(1, 1)
  list(
    senior = fmt_est(beta[1], sqrt(V[1, 1]), with_p),
    inter  = fmt_est(beta[2], sqrt(V[2, 2]), with_p),
    junior = fmt_est(sum(beta), sqrt(drop(t(c_jun) %*% V %*% c_jun)), with_p)
  )
}

res_fixed <- group_effects(fit_fixed$beta_HTE, fit_fixed$vcov_HTE)
res_tvcox <- group_effects(fit_tvcox$beta_HTE, fit_tvcox$vcov_HTE)

# ── 7. TV-CSL via TV_CSL() ────────────────────────────────────────────────────
# Nuisance nu: m-regression = Cox on X ignoring treatment (eta_0 = eta_1).
# Final model: coxph on (W - a(t,X)) * (1, is_junior) with offset(nu),
# cross-fitted over K folds; sandwich SE from stacked subject-level scores.

run_tvcsl_lib <- function(prop_score_spec) {
  quietly(TV_CSL(
    train_data          = df_tv,
    test_data           = df_orig,
    train_data_original = df_orig,
    HTE_type            = NA,
    eta_type            = NA,
    K                   = K_FOLDS,
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
}

message("Fitting TV-CSL with marginal propensity ...")
res_marg <- run_tvcsl_lib("cox-linear-censored-only")

message("Fitting TV-CSL with time-varying propensity ...")
res_tvp  <- run_tvcsl_lib("cox-time-varying-prop")

# Omega_hat is the asymptotic covariance of sqrt(n) * beta; n = subjects.
n_subj <- nrow(df_orig)
res_marg_fmt <- group_effects(res_marg$beta_HTE, res_marg$Omega_hat / n_subj,
                              with_p = FALSE)
res_tvp_fmt  <- group_effects(res_tvp$beta_HTE,  res_tvp$Omega_hat / n_subj,
                              with_p = FALSE)

# ── 8. Write tables ───────────────────────────────────────────────────────────

tables_dir <- OUTPUT_DIR
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

## Table 1: Effect of ignoring treatment time

rows1 <- list(
  c("Used CC Assist (senior)",
    res_fixed$senior$ci,  res_fixed$senior$pval,
    res_tvcox$senior$ci,  res_tvcox$senior$pval),
  c("Used CC Assist (junior)",
    res_fixed$junior$ci,  res_fixed$junior$pval,
    res_tvcox$junior$ci,  res_tvcox$junior$pval),
  c("Junior $-$ senior difference",
    res_fixed$inter$ci,   res_fixed$inter$pval,
    res_tvcox$inter$ci,   res_tvcox$inter$pval)
)

tab1_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\caption{Cox Estimates That Ignore and Include CC Assist Initiation Time}",
  "\\label{tab:copilot-ignore-time}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{l|cc|cc}",
  "\\toprule",
  " & \\multicolumn{2}{c|}{Ignore treatment time}",
  " & \\multicolumn{2}{c}{Time-varying treatment} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Effect & Estimate (95\\% CI) & $p$-value",
  "       & Estimate (95\\% CI) & $p$-value \\\\",
  "\\midrule"
)
for (r in rows1) {
  tab1_lines <- c(tab1_lines,
    paste0(r[1], " & ", r[2], " & ", r[3], " & ", r[4], " & ", r[5], " \\\\"))
}
tab1_lines <- c(tab1_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\par\\smallskip",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize",
  "Estimates are on the log-hazard-ratio scale; positive values indicate a higher",
  "instantaneous rate of becoming a regular Copilot user. The junior$-$senior difference",
  "and its confidence interval use the estimated covariance of the two group effects.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab1_lines,
           file.path(tables_dir, "tab_copilot_ignore_time.tex"))

## Table 2: TV-CSL comparison

rows2 <- list(
  c("Used CC Assist (senior)",
    res_tvcox$senior$ci,
    res_marg_fmt$senior,
    res_tvp_fmt$senior),
  c("Used CC Assist (junior)",
    res_tvcox$junior$ci,
    res_marg_fmt$junior,
    res_tvp_fmt$junior),
  c("Junior $-$ senior difference",
    res_tvcox$inter$ci,
    res_marg_fmt$inter,
    res_tvp_fmt$inter)
)

tab2_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Comparison of Time-Varying Cox and TV-CSL Estimates}",
  "\\label{tab:copilot-tvcsl}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  " & Time-varying & \\multicolumn{2}{c}{TV-CSL} \\\\",
  "\\cmidrule(lr){3-4}",
  "Effect & Cox & Marginal propensity & Time-varying propensity \\\\",
  " & Estimate (95\\% CI) & Estimate (95\\% CI) & Estimate (95\\% CI) \\\\",
  "\\midrule"
)
for (r in rows2) {
  tab2_lines <- c(tab2_lines,
    paste0(r[1], " & ", r[2], " & ", r[3], " & ", r[4], " \\\\"))
}
tab2_lines <- c(tab2_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\par\\smallskip",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize",
  "Estimates are on the log-hazard-ratio scale; positive values indicate a higher",
  "instantaneous rate of becoming a regular Copilot user. TV-CSL (time-varying propensity)",
  "intervals use the cross-fitted sandwich variance estimator. Marginal-propensity",
  "intervals are working sandwich intervals.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab2_lines,
           file.path(tables_dir, "tab_copilot_tvcsl.tex"))

message("Tables written to ", tables_dir)
