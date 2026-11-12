## stanford_analysis.R
## Stanford Heart Transplant analysis for the appendix, following
## scripts/law-firm-ai-analysis/copilot_analysis_tvcsl_lib.R.
## Prior surgery is an adjustment covariate only; age and year are the
## effect modifiers.
##
## Produces in OUTPUT_DIR:
##   tab_stanford_ignore_time.tex   Cox ignoring vs including transplant time
##   tab_stanford_tvcsl.tex         time-varying Cox vs TV-CSL
##   stanford_fold_check.csv        prior-surgery counts per training fold
##   stanford_results.csv           every estimate, SE, CI and p-value
##   stanford_split_sensitivity.csv TV-CSL estimates over N_SPLITS fold splits

# ── User configuration ────────────────────────────────────────────────────────

OUTPUT_DIR <- here::here("scripts/heart-transplant-analysis/tables-stanford")
K_FOLDS    <- 5L
SEED       <- 2026L
N_SPLITS   <- 20L

# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(survival)
  library(mgcv)
  library(here)
  source(here::here("R/TV-CSL.R"))
  source(here::here("R/data-handler.R"))
  source(here::here("scripts/heart-transplant-analysis/stanford_common.R"))
})

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Data and fold check ────────────────────────────────────────────────────

set.seed(SEED)
d <- shuffle_rows(load_stanford())
message(sprintf("n = %d, deaths = %d, transplanted during follow-up = %d, prior surgery = %d",
                nrow(d), sum(d$Delta), sum(d$A < d$U), sum(d$X.3)))

fc <- fold_check(d, K_FOLDS)
write.csv(fc, file.path(OUTPUT_DIR, "stanford_fold_check.csv"), row.names = FALSE)
message("Prior-surgery counts per training fold:")
print(as.data.frame(fc))
if (!all(fc$ok)) warning("A training fold lacks treated or untreated prior-surgery events.")

# ── 2. Fit all four methods ───────────────────────────────────────────────────

results <- run_methods(d, K = K_FOLDS)
write.csv(results, file.path(OUTPUT_DIR, "stanford_results.csv"), row.names = FALSE)

# TV-CSL depends on the random fold split; record how much.
split_res <- bind_rows(lapply(seq_len(N_SPLITS), function(s) {
  set.seed(SEED + s)
  run_methods(shuffle_rows(d), K = K_FOLDS) %>%
    filter(method %in% c("tvcsl_marg", "tvcsl_tvp")) %>%
    mutate(split = s, .before = 1)
}))
write.csv(split_res, file.path(OUTPUT_DIR, "stanford_split_sensitivity.csv"),
          row.names = FALSE)
message("TV-CSL estimates across ", N_SPLITS, " fold splits:")
print(as.data.frame(split_res %>% group_by(method, effect) %>%
  summarise(min = min(est), median = median(est), max = max(est),
            .groups = "drop")), digits = 3)

pick <- function(m, e) results[results$method == m & results$effect == e, ]
cell_ci <- function(r) fmt_ci(r$est, r$lo, r$hi)

note_lines <- c(
  "Estimates are log hazard ratios for death; negative values indicate lower mortality",
  "after transplant. Age and year of acceptance are standardized. All models adjust",
  "for age, year of acceptance and prior surgery; prior surgery is not an effect modifier."
)

# ── 3. Table 1: effect of ignoring transplant time ────────────────────────────

tab1 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Cox Estimates That Ignore and Include Transplant Time}",
  "\\label{tab:impact-time-varying-HTE}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{l|cc|cc}",
  "\\toprule",
  " & \\multicolumn{2}{c|}{Ignore transplant time}",
  " & \\multicolumn{2}{c}{Include transplant time} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Effect & Estimate (95\\% CI) & $p$-value & Estimate (95\\% CI) & $p$-value \\\\",
  "\\midrule"
)
for (e in names(EFFECTS)) {
  a <- pick("fixed_cox", e); b <- pick("tv_cox", e)
  tab1 <- c(tab1, paste0(EFFECTS[[e]], " & ", cell_ci(a), " & ", fmt_p(a$p),
                         " & ", cell_ci(b), " & ", fmt_p(b$p), " \\\\"))
}
tab1 <- c(tab1, "\\bottomrule", "\\end{tabular}", "}", "\\par\\smallskip",
          "\\begin{minipage}{\\linewidth}", "\\footnotesize", note_lines,
          "\\end{minipage}", "\\end{table}")
writeLines(tab1, file.path(OUTPUT_DIR, "tab_stanford_ignore_time.tex"))

# ── 4. Table 2: time-varying Cox vs TV-CSL ────────────────────────────────────

tab2 <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Time-Varying Cox and TV-CSL Estimates for the Stanford Heart Transplant Data}",
  "\\label{tab:stanford-tvcsl}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  " & Time-varying & \\multicolumn{2}{c}{TV-CSL} \\\\",
  "\\cmidrule(lr){3-4}",
  "Effect & Cox & Marginal propensity & Time-varying propensity \\\\",
  "\\midrule"
)
for (e in names(EFFECTS)) {
  tab2 <- c(tab2, paste0(EFFECTS[[e]], " & ", cell_ci(pick("tv_cox", e)),
                         " & ", cell_ci(pick("tvcsl_marg", e)),
                         " & ", cell_ci(pick("tvcsl_tvp", e)), " \\\\"))
}
tab2 <- c(tab2, "\\bottomrule", "\\end{tabular}", "}", "\\par\\smallskip",
          "\\begin{minipage}{\\linewidth}", "\\footnotesize", note_lines,
          sprintf("TV-CSL uses %d-fold cross-fitting and sandwich standard errors.", K_FOLDS),
          "\\end{minipage}", "\\end{table}")
writeLines(tab2, file.path(OUTPUT_DIR, "tab_stanford_tvcsl.tex"))

message("Tables written to ", OUTPUT_DIR)
