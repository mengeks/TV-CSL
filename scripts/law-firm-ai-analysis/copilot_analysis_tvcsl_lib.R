## copilot_analysis_tvcsl_lib.R
## Copilot (CC Assist) analysis for Section 6, with every estimator taken from
## R/TV-CSL.R:
##   - Cox fits (fixed and time-varying)  -> S_lasso(regressor_spec = "linear",
##                                           HTE_spec = "linear"), the unpenalized
##                                           coxph path with W * (1 + is_junior)
##   - TV-CSL (both propensity choices)   -> TV_CSL(lasso_type = "m-regression")
##       marginal propensity      = prop_score_spec "cox-linear-censored-only"
##       time-varying propensity  = prop_score_spec "cox-time-varying-prop"
##
## Analyses (cohorts, tie conventions and estimator wrappers are in copilot_common.R):
##   primary      new-user cohort,  W(t) = 1(A < t)
##   sensitivity  new-user cohort,  W(t) = 1(A <= t)
##                all lawyers,      W(t) = 1(A < t)
##                all lawyers,      W(t) = 1(A <= t)
##
## Produces in OUTPUT_DIR:
##   tab_copilot_ignore_time.tex   primary analysis, Cox ignoring vs including time
##   tab_copilot_tvcsl.tex         primary analysis, time-varying Cox vs TV-CSL
##   tab_copilot_sensitivity.tex   all four analyses side by side
##   copilot_cohort_counts.csv     sample flow for both cohorts
##   copilot_results_all.csv       every estimate, SE, CI and p-value

# ── User configuration ────────────────────────────────────────────────────────

INPUT_FILE  <- here::here("scripts/law-firm-ai-analysis/survival_ccassist_SYNTHETIC.csv")
OUTPUT_DIR  <- here::here("scripts/law-firm-ai-analysis/tables-tvcsl-lib")
K_FOLDS     <- 5L
PRIMARY     <- c(cohort = "new_user", ties = "lt")

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

# ── 1. Load data and report sample flow ───────────────────────────────────────

df_raw <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

counts <- cohort_counts(df_raw)
write.csv(counts, file.path(OUTPUT_DIR, "copilot_cohort_counts.csv"),
          row.names = FALSE)
message("Sample flow:")
print(as.data.frame(t(counts %>% tibble::column_to_rownames("cohort"))))

# ── 2. Run all four analyses ──────────────────────────────────────────────────

analyses <- expand.grid(cohort = names(COHORT_LABELS), ties = names(TIES_LABELS),
                        stringsAsFactors = FALSE)
results <- bind_rows(lapply(seq_len(nrow(analyses)), function(j) {
  message("Fitting cohort = ", analyses$cohort[j], ", ties = ", analyses$ties[j], " ...")
  ch <- analyses$cohort[j]
  run_methods(prepare_cohort(df_raw, ch), cohort_covariates(ch),
              analyses$ties[j], K = K_FOLDS) %>%
    mutate(cohort = ch, .before = 1)
})) %>%
  mutate(primary = cohort == PRIMARY[["cohort"]] & ties == PRIMARY[["ties"]])

write.csv(results, file.path(OUTPUT_DIR, "copilot_results_all.csv"),
          row.names = FALSE)

# ── 3. Formatting helpers ─────────────────────────────────────────────────────

fmt_ci <- function(r) sprintf("%.3f (%.3f, %.3f)", r$est, r$lo, r$hi)
fmt_p  <- function(p) ifelse(p < 0.001, "$<$0.001", sprintf("%.3f", p))

pick <- function(res, m, e) res[res$method == m & res$effect == e, ]

primary  <- results %>% filter(primary)
n_lab    <- function(ch) format(counts$n_analysed[counts$cohort == ch], big.mark = ",")
primary_note <- sprintf(
  "Primary analysis: %s cohort ($n = %s$), %s; same-week adoption and event counts as untreated.",
  tolower(COHORT_LABELS[[PRIMARY[["cohort"]]]]), n_lab(PRIMARY[["cohort"]]),
  TIES_LABELS[[PRIMARY[["ties"]]]])

# ── 4. Table 1: effect of ignoring treatment time (primary) ───────────────────

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
for (e in names(EFFECTS)) {
  a <- pick(primary, "fixed_cox", e); b <- pick(primary, "tv_cox", e)
  tab1_lines <- c(tab1_lines, paste0(
    EFFECTS[[e]], " & ", fmt_ci(a), " & ", fmt_p(a$p),
    " & ", fmt_ci(b), " & ", fmt_p(b$p), " \\\\"))
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
  primary_note,
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab1_lines, file.path(OUTPUT_DIR, "tab_copilot_ignore_time.tex"))

# ── 5. Table 2: time-varying Cox vs TV-CSL (primary) ──────────────────────────

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
for (e in names(EFFECTS)) {
  tab2_lines <- c(tab2_lines, paste0(
    EFFECTS[[e]], " & ", fmt_ci(pick(primary, "tv_cox", e)),
    " & ", fmt_ci(pick(primary, "tvcsl_marg", e)),
    " & ", fmt_ci(pick(primary, "tvcsl_tvp", e)), " \\\\"))
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
  primary_note,
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab2_lines, file.path(OUTPUT_DIR, "tab_copilot_tvcsl.tex"))

# ── 6. Sensitivity table: all four analyses side by side ──────────────────────
# Rows: method x effect. Columns: cohort x tie convention, primary first.

col_order <- analyses %>%
  mutate(primary = cohort == PRIMARY[["cohort"]] & ties == PRIMARY[["ties"]],
         cohort_first = cohort == PRIMARY[["cohort"]]) %>%
  arrange(desc(cohort_first), desc(ties == PRIMARY[["ties"]]))

col_head <- ifelse(col_order$primary,
                   paste0(TIES_LABELS[col_order$ties], "$^\\dagger$"),
                   TIES_LABELS[col_order$ties])
cohort_runs <- rle(col_order$cohort)
cohort_head <- paste(sprintf("\\multicolumn{%d}{c}{%s ($n = %s$)}",
                             cohort_runs$lengths,
                             COHORT_LABELS[cohort_runs$values],
                             sapply(cohort_runs$values, n_lab)),
                     collapse = " & ")
cmid <- {
  ends <- cumsum(cohort_runs$lengths) + 1
  starts <- ends - cohort_runs$lengths + 1
  paste(sprintf("\\cmidrule(lr){%d-%d}", starts, ends), collapse = "")
}

tab3_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Sensitivity of CC Assist Effect Estimates to Cohort and Tie Convention}",
  "\\label{tab:copilot-sensitivity}",
  "\\resizebox{\\linewidth}{!}{%",
  paste0("\\begin{tabular}{l", strrep("c", nrow(col_order)), "}"),
  "\\toprule",
  paste0(" & ", cohort_head, " \\\\"),
  cmid,
  paste0("Effect & ", paste(col_head, collapse = " & "), " \\\\"),
  "\\midrule"
)
for (m in names(METHODS)) {
  tab3_lines <- c(tab3_lines,
    sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", nrow(col_order) + 1, METHODS[[m]]))
  for (e in names(EFFECTS)) {
    cells <- sapply(seq_len(nrow(col_order)), function(j) {
      r <- results %>% filter(cohort == col_order$cohort[j], ties == col_order$ties[j],
                              method == m, effect == e)
      fmt_ci(r)
    })
    tab3_lines <- c(tab3_lines,
      paste0("\\quad ", EFFECTS[[e]], " & ", paste(cells, collapse = " & "), " \\\\"))
  }
  if (m != tail(names(METHODS), 1)) tab3_lines <- c(tab3_lines, "\\addlinespace")
}
tab3_lines <- c(tab3_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\par\\smallskip",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize",
  "Estimates (95\\% CI) on the log-hazard-ratio scale. $^\\dagger$Primary analysis.",
  "New users have no recorded CC Assist use before time zero; the all-lawyers cohort",
  "adds prior users and adjusts for prior use. Both cohorts exclude lawyers who became",
  "regular users in their entry week. $1(A < t)$ counts adoption in the event week as",
  "after the event; $1(A \\le t)$ counts it as before. The ignore-time Cox model uses the",
  "ever-used indicator and so does not depend on the tie convention.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab3_lines, file.path(OUTPUT_DIR, "tab_copilot_sensitivity.tex"))

message("Tables written to ", OUTPUT_DIR)
