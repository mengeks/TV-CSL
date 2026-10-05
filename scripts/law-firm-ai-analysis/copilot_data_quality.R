## copilot_data_quality.R
## Data-quality checks for the nuisance models of copilot_analysis_tvcsl_lib.R.
##
## 1. Baseline weekly hours: how much is missing, where, whether missingness
##    is related to the outcome or to CC Assist adoption, and how much the
##    effect estimates move when a missing-hours indicator is added.
## 2. Sparse levels and overlap (rank, office, practice group): per level,
##    lawyers, events, and events that are treated vs untreated at the event
##    time (what the time-varying propensity model is fitted on), overall and
##    in every cross-fitting training fold; fitted propensities near 0 or 1;
##    signs of separation in both propensity models.
## 3. Collapsing: office / practice-group levels that fail the thresholds are
##    merged into "Other" and the primary analysis is refitted for comparison.
##
## Writes report.txt and CSVs to OUTPUT_DIR.

# ── User configuration ────────────────────────────────────────────────────────

INPUT_FILE  <- here::here("scripts/law-firm-ai-analysis/survival_ccassist_SYNTHETIC.csv")
OUTPUT_DIR  <- here::here("scripts/law-firm-ai-analysis/data-quality")
K_FOLDS     <- 5L
PRIMARY     <- c(cohort = "new_user", ties = "lt")

MIN_LAWYERS        <- 30           # per level
MIN_EVENTS_PER_ARM <- 10           # treated and untreated events per level
PROP_BOUNDS        <- c(0.05, 0.95)

# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(survival)
  library(mgcv)
  library(here)
  source(here::here("R/TV-CSL.R"))
  source(here::here("R/data-handler.R"))
  source(here::here("scripts/law-firm-ai-analysis/copilot_common.R"))
})

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
options(width = 200, dplyr.summarise.inform = FALSE)

report_file <- file.path(OUTPUT_DIR, "report.txt")
sink(report_file, split = TRUE)

header <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
sub_header <- function(x) cat("\n--", x, "--\n")
show <- function(df, digits = 3) print(as.data.frame(df), digits = digits, row.names = FALSE)
pct <- function(x) round(100 * mean(x), 1)

df_raw <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
FACTORS <- c("rank", "office", "practice_group")

# Person-level + counting-process data for one cohort under the primary ties
build <- function(cohort, d = prepare_cohort(df_raw, cohort),
                  covariates = cohort_covariates(cohort)) {
  des <- add_design(d, covariates)
  dat <- make_tv(des$data, PRIMARY[["ties"]])
  dat$orig <- dat$orig %>% mutate(treated_event = as.integer(Delta == 1 & A_eff < U))
  c(dat, list(xvars = des$xvars, labels = des$labels))
}
cohorts <- setNames(lapply(names(COHORT_LABELS), build), names(COHORT_LABELS))

header("Sample flow")
print(t(cohort_counts(df_raw) %>% tibble::column_to_rownames("cohort")))

# ═════════════════════════════════════════════════════════════════════════════
header("1. Baseline weekly hours")
# ═════════════════════════════════════════════════════════════════════════════

hours_by_level <- list()
for (ch in names(cohorts)) {
  o <- cohorts[[ch]]$orig
  obs <- df_raw$baseline_weekly_hours[match(o$person_id, df_raw$person_id)]
  sub_header(sprintf("%s: %d of %d missing (%.1f%%)", COHORT_LABELS[[ch]],
                     sum(o$baseline_hours_missing), nrow(o), pct(o$baseline_hours_missing)))

  cat("Observed values:\n"); print(summary(obs))
  cat(sprintf("Outside 0-100 hours: %d\n", sum(obs < 0 | obs > 100, na.rm = TRUE)))

  cat("\nMissing % by level (sorted):\n")
  by_lvl <- bind_rows(lapply(c("seniority", FACTORS), function(v) {
    o %>% group_by(level = as.character(.data[[v]])) %>%
      summarise(n = n(), n_missing = sum(baseline_hours_missing),
                pct_missing = pct(baseline_hours_missing)) %>%
      mutate(variable = v, .before = 1)
  })) %>% mutate(cohort = ch, .before = 1)
  hours_by_level[[ch]] <- by_lvl
  show(by_lvl %>% select(-cohort) %>% arrange(desc(pct_missing)), 3)

  cat("\nOutcome and treatment by missing status:\n")
  show(o %>% group_by(hours_missing = baseline_hours_missing) %>% summarise(
    n = n(),
    pct_event = pct(Delta),
    median_followup_wk = median(U - entry_week),
    pct_adopt_during_fu = pct(Delta_A),
    pct_events_treated = 100 * sum(treated_event) / sum(Delta)))

  # Adjusted associations: does missingness predict the outcome hazard or
  # treatment status at the event, beyond the other covariates?
  tv <- cohorts[[ch]]$tv
  xv <- cohorts[[ch]]$xvars
  hrs_var <- names(cohorts[[ch]]$labels)[cohorts[[ch]]$labels == "baseline_weekly_hours"]
  f_out <- as.formula(paste("Surv(tstart, tstop, Delta) ~ W + baseline_hours_missing +",
                            paste(xv, collapse = " + ")))
  m_out <- coxph(f_out, data = tv, ties = "breslow")
  ev <- tv %>% filter(Delta == 1)
  f_trt <- as.formula(paste("W ~ s(tstop, k = 5) + baseline_hours_missing +",
                            paste(xv, collapse = " + ")))
  m_trt <- gam(f_trt, family = binomial, data = ev)
  pt <- summary(m_trt)$p.table
  co <- summary(m_out)$coefficients
  cat("\nAdjusted associations (log scale):\n")
  show(data.frame(
    model = c("outcome Cox: missing indicator", "outcome Cox: hours (per 10 h)",
              "treatment-at-event GAM: missing indicator", "treatment-at-event GAM: hours (per 10 h)"),
    estimate = c(co["baseline_hours_missing", "coef"], 10 * co[hrs_var, "coef"],
                 pt["baseline_hours_missing", "Estimate"], 10 * pt[hrs_var, "Estimate"]),
    se = c(co["baseline_hours_missing", "se(coef)"], 10 * co[hrs_var, "se(coef)"],
           pt["baseline_hours_missing", "Std. Error"], 10 * pt[hrs_var, "Std. Error"]),
    p = c(co["baseline_hours_missing", "Pr(>|z|)"], co[hrs_var, "Pr(>|z|)"],
          pt["baseline_hours_missing", "Pr(>|z|)"], pt[hrs_var, "Pr(>|z|)"])), 3)
}
write.csv(bind_rows(hours_by_level), file.path(OUTPUT_DIR, "hours_missing_by_level.csv"),
          row.names = FALSE)

sub_header("Impact on the primary analysis: add a missing-hours indicator")
ch <- PRIMARY[["cohort"]]
d_primary <- prepare_cohort(df_raw, ch)
res_base <- run_methods(d_primary, cohort_covariates(ch), PRIMARY[["ties"]], K_FOLDS)
res_mind <- run_methods(d_primary, c(cohort_covariates(ch), "baseline_hours_missing"),
                        PRIMARY[["ties"]], K_FOLDS)
cmp_hours <- res_base %>% select(method, effect, est_mean_impute = est, se) %>%
  left_join(res_mind %>% select(method, effect, est_with_indicator = est), by = c("method", "effect")) %>%
  mutate(change = est_with_indicator - est_mean_impute, change_in_se = change / se)
show(cmp_hours, 3)
write.csv(cmp_hours, file.path(OUTPUT_DIR, "hours_indicator_impact.csv"), row.names = FALSE)

# ═════════════════════════════════════════════════════════════════════════════
header("2. Sparse levels and overlap")
# ═════════════════════════════════════════════════════════════════════════════
cat(sprintf(paste0(
  "Thresholds: >= %d lawyers per level; >= %d treated and >= %d untreated events per level\n",
  "(treated = already adopted at the event time, the response of the time-varying\n",
  "propensity model); no level without treated or untreated events in any of the\n",
  "%d training folds; fitted propensities inside [%.2f, %.2f].\n"),
  MIN_LAWYERS, MIN_EVENTS_PER_ARM, MIN_EVENTS_PER_ARM, K_FOLDS, PROP_BOUNDS[1], PROP_BOUNDS[2]))

level_table <- function(o, var) {
  folds <- cut(seq_len(nrow(o)), breaks = K_FOLDS, labels = FALSE)  # as in TV_CSL()
  overall <- o %>% group_by(level = .data[[var]], .drop = FALSE) %>% summarise(
    n = n(), juniors = sum(is_junior), events = sum(Delta),
    treated_events = sum(treated_event), untreated_events = events - treated_events,
    pct_events_treated = round(100 * treated_events / pmax(events, 1), 1),
    adoptions_among_events = sum(Delta_A[Delta == 1]))
  by_fold <- bind_rows(lapply(seq_len(K_FOLDS), function(k) {
    o[folds != k, ] %>% group_by(level = .data[[var]], .drop = FALSE) %>%
      summarise(te = sum(treated_event), ue = sum(Delta) - sum(treated_event),
                ad = sum(Delta_A[Delta == 1]))
  })) %>% group_by(level) %>%
    summarise(fold_min_treated = min(te), fold_min_untreated = min(ue),
              fold_min_adoptions = min(ad))
  overall %>% left_join(by_fold, by = "level") %>%
    mutate(variable = var, .before = 1) %>%
    mutate(level = as.character(level),
           flag_small = n < MIN_LAWYERS,
           flag_arm   = treated_events < MIN_EVENTS_PER_ARM |
                        untreated_events < MIN_EVENTS_PER_ARM,
           flag_fold  = fold_min_treated == 0 | fold_min_untreated == 0 |
                        fold_min_adoptions == 0,
           flag_any   = flag_small | flag_arm | flag_fold)
}

level_tables <- list()
prop_overlap <- list()
for (ch in names(cohorts)) {
  cc <- cohorts[[ch]]
  sub_header(sprintf("%s (n = %d, events = %d)", COHORT_LABELS[[ch]],
                     nrow(cc$orig), sum(cc$orig$Delta)))
  lt <- bind_rows(lapply(FACTORS, function(v) level_table(cc$orig, v))) %>%
    mutate(cohort = ch, .before = 1)
  level_tables[[ch]] <- lt
  show(lt %>% select(-cohort))

  cat("\nEffect modifier: treated / untreated events by seniority\n")
  show(cc$orig %>% filter(Delta == 1) %>% group_by(is_junior) %>%
         summarise(events = n(), treated = sum(treated_event), untreated = n() - sum(treated_event)))

  # Time-varying propensity, fitted on all event rows (diagnostic, not cross-fitted)
  ev <- cc$tv %>% filter(Delta == 1)
  m_tvp <- gam(as.formula(paste("W ~ s(tstop, k = 5) +", paste(cc$xvars, collapse = " + "))),
               family = binomial, data = ev)
  ev$a_hat <- as.vector(fitted(m_tvp))
  pt <- summary(m_tvp)$p.table[-1, , drop = FALSE]
  # Marginal propensity: Cox for adoption among events; a = 1 - exp(-t exp(X alpha))
  ev_o <- cc$orig %>% filter(Delta == 1)
  m_marg <- coxph(as.formula(paste("Surv(U_A, Delta_A) ~", paste(cc$xvars, collapse = " + "))),
                  data = ev_o, ties = "breslow")
  alpha <- coef(m_marg); alpha[is.na(alpha)] <- 0
  a_marg <- 1 - exp(-ev_o$U * exp(as.vector(as.matrix(ev_o[, cc$xvars]) %*% alpha)))

  cat("\nFitted propensity at event times (all events):\n")
  ov <- bind_rows(
    tibble::tibble(model = "time-varying", treated = ev$W, a = ev$a_hat),
    tibble::tibble(model = "marginal", treated = ev_o$treated_event, a = a_marg)) %>%
    group_by(model, treated) %>%
    summarise(n = n(), min = min(a), q05 = quantile(a, .05), median = median(a),
              q95 = quantile(a, .95), max = max(a),
              pct_below = pct(a < PROP_BOUNDS[1]), pct_above = pct(a > PROP_BOUNDS[2])) %>%
    mutate(cohort = ch, .before = 1)
  prop_overlap[[ch]] <- ov
  show(ov %>% select(-cohort))

  cat("\nCovariate coefficients suggesting separation (|coef| > 5, SE > 5 or NA):\n")
  sep <- bind_rows(
    tibble::tibble(model = "time-varying GAM", term = cc$labels[rownames(pt)],
                   coef = pt[, "Estimate"], se = pt[, "Std. Error"]),
    tibble::tibble(model = "marginal Cox", term = cc$labels[names(coef(m_marg))],
                   coef = coef(m_marg), se = sqrt(diag(vcov(m_marg))))) %>%
    filter(is.na(coef) | abs(coef) > 5 | se > 5)
  if (nrow(sep)) show(sep) else cat("none\n")
}
write.csv(bind_rows(level_tables), file.path(OUTPUT_DIR, "level_overlap.csv"), row.names = FALSE)
write.csv(bind_rows(prop_overlap), file.path(OUTPUT_DIR, "propensity_overlap.csv"), row.names = FALSE)

# ═════════════════════════════════════════════════════════════════════════════
header("3. Collapsing sparse office / practice-group levels (primary analysis)")
# ═════════════════════════════════════════════════════════════════════════════
# Rank is not collapsed: it defines seniority, the effect modifier.

flagged <- level_tables[[PRIMARY[["cohort"]]]] %>%
  filter(variable %in% c("office", "practice_group"), flag_any, level != "Other")
if (nrow(flagged) == 0) {
  cat("No office or practice-group level fails the thresholds; nothing to collapse.\n")
} else {
  cat("Levels merged into \"Other\":\n")
  show(flagged %>% select(variable, level, n, treated_events, untreated_events,
                          fold_min_treated, fold_min_untreated))
  d_coll <- d_primary
  for (v in unique(flagged$variable)) {
    lv <- flagged$level[flagged$variable == v]
    d_coll[[v]] <- droplevels(factor(ifelse(as.character(d_coll[[v]]) %in% lv, "Other",
                                            as.character(d_coll[[v]]))))
  }
  res_coll <- run_methods(d_coll, cohort_covariates(PRIMARY[["cohort"]]),
                          PRIMARY[["ties"]], K_FOLDS)
  cmp_coll <- res_base %>% select(method, effect, est_original = est, se) %>%
    left_join(res_coll %>% select(method, effect, est_collapsed = est), by = c("method", "effect")) %>%
    mutate(change = est_collapsed - est_original, change_in_se = change / se)
  cat("\nEffect estimates, original vs collapsed levels:\n")
  show(cmp_coll, 3)
  write.csv(cmp_coll, file.path(OUTPUT_DIR, "collapse_impact.csv"), row.names = FALSE)
}

sink()
message("Report written to ", report_file)
