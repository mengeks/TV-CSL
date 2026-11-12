## stanford_plasmode.R
## Semi-synthetic (plasmode) study at the Stanford sample size: the real 103
## patients' covariates and censoring dates are kept, and transplant and death
## times are redrawn from models fitted to the real data, so the true effects
## are known.
##
## Truth (fitted on the real data):
##   death:      time-varying Cox  X.1 + X.2 + X.3 + W * (1 + X.1 + X.2)
##               -> gamma (baseline risk) and beta (transplant effects)
##   transplant: Cox for time to transplant on X.1 + X.2 + X.3
## Both baselines are Breslow estimates, linearly interpolated so that times
## can be drawn by exact inversion of the cumulative hazard (including the
## hazard change at transplant). Censoring is each patient's real
## administrative censoring date (study close), capped at the last death time.
##
## Produces in OUTPUT_DIR:
##   tab_stanford_plasmode.tex   bias and coverage by method
##   plasmode_truth.csv          true coefficients
##   plasmode_reps.csv           every replicate's estimates
##   plasmode_summary.csv        bias, SD, RMSE and coverage by method
##   plasmode_fold_check.csv     prior-surgery fold counts per replicate

# ── User configuration ────────────────────────────────────────────────────────

OUTPUT_DIR <- here::here("scripts/heart-transplant-analysis/tables-stanford")
K_FOLDS    <- 5L
N_REPS     <- 500L
N_CORES    <- 12L
SEED       <- 2026L

# ─────────────────────────────────────────────────────────────────────────────

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(survival)
  library(mgcv)
  library(parallel)
  library(here)
  source(here::here("R/TV-CSL.R"))
  source(here::here("R/data-handler.R"))
  source(here::here("scripts/heart-transplant-analysis/stanford_common.R"))
})

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Truth from the real data ───────────────────────────────────────────────

# Breslow cumulative hazard at zero covariates, linearly interpolated
cumhaz <- function(fit) {
  # survfit warns about interactions; with centered = FALSE the curve is at
  # zero covariates, which is what is needed here.
  bh <- suppressWarnings(basehaz(fit, centered = FALSE))
  t  <- c(0, bh$time)
  H  <- c(0, bh$hazard)
  keep <- !duplicated(H)
  t <- t[keep]; H <- H[keep]
  list(
    H    = approxfun(t, H, rule = 2),
    Hinv = function(h) ifelse(h > max(H), Inf, approx(H, t, xout = h)$y),
    tmax = max(t)
  )
}

d_real   <- load_stanford()
dat_real <- make_tv(d_real)

fit_death <- coxph(Surv(tstart, tstop, Delta) ~ X.1 + X.2 + X.3 + W + W:X.1 + W:X.2,
                   data = dat_real$tv, ties = "breslow")
fit_tx    <- coxph(Surv(pmax(U_A, 0.5), Delta_A) ~ X.1 + X.2 + X.3,
                   data = dat_real$orig, ties = "breslow")

truth <- list(
  gamma = coef(fit_death)[XVARS],
  beta  = setNames(coef(fit_death)[c("W", "X.1:W", "X.2:W")], names(EFFECTS)),
  alpha = coef(fit_tx)[XVARS],
  H_T   = cumhaz(fit_death),
  H_A   = cumhaz(fit_tx)
)
write.csv(data.frame(effect = names(EFFECTS), truth = truth$beta),
          file.path(OUTPUT_DIR, "plasmode_truth.csv"), row.names = FALSE)
message("True effects: ", paste(names(truth$beta), round(truth$beta, 3),
                                sep = " = ", collapse = ", "))

# ── 2. Data generation ────────────────────────────────────────────────────────

simulate_once <- function(d, truth) {
  n <- nrow(d)
  X <- as.matrix(d[, XVARS])
  # Local names differ from d's columns: inside mutate(), A would refer to
  # the real transplant times.
  A_sim <- truth$H_A$Hinv(rexp(n) / exp(as.vector(X %*% truth$alpha)))

  # Cumulative hazard is H0(t) e^eta0 before A and gains e^(eta0 + tau) after.
  eta0  <- as.vector(X %*% truth$gamma)
  tau   <- as.vector(cbind(1, X[, MODIFIERS]) %*% truth$beta)
  E     <- rexp(n) / exp(eta0)
  HA    <- truth$H_T$H(pmin(A_sim, truth$H_T$tmax))
  T_sim <- ifelse(E <= HA, truth$H_T$Hinv(E),
                  truth$H_T$Hinv(HA + (E - HA) * exp(-tau)))

  C_sim <- pmin(d$C_admin, truth$H_T$tmax)
  d$A     <- A_sim
  d$U     <- pmin(T_sim, C_sim)
  d$Delta <- as.integer(T_sim <= C_sim)
  d
}

# ── 3. Replicates ─────────────────────────────────────────────────────────────

one_rep <- function(r) {
  set.seed(SEED + r)
  d <- shuffle_rows(simulate_once(d_real, truth))
  fc <- fold_check(d, K_FOLDS)
  list(
    est = run_methods(d, K = K_FOLDS) %>% mutate(rep = r, .before = 1),
    fc  = tibble::tibble(rep = r, deaths = sum(d$Delta),
                         transplanted = sum(d$A < d$U),
                         min_surgery_treated_events   = min(fc$surgery_treated_events),
                         min_surgery_untreated_events = min(fc$surgery_untreated_events),
                         all_folds_ok = all(fc$ok))
  )
}

message(sprintf("Running %d replicates on %d cores ...", N_REPS, N_CORES))
t0 <- Sys.time()
reps <- mclapply(seq_len(N_REPS), one_rep, mc.cores = N_CORES)
message(sprintf("Done in %.1f min", as.numeric(Sys.time() - t0, units = "mins")))

failed <- vapply(reps, inherits, logical(1), "try-error")
if (any(failed)) warning(sum(failed), " replicates failed entirely")
reps <- reps[!failed]

est <- bind_rows(lapply(reps, `[[`, "est"))
fcs <- bind_rows(lapply(reps, `[[`, "fc"))
write.csv(est, file.path(OUTPUT_DIR, "plasmode_reps.csv"), row.names = FALSE)
write.csv(fcs, file.path(OUTPUT_DIR, "plasmode_fold_check.csv"), row.names = FALSE)

message(sprintf(
  "Simulated data: deaths %.1f (real %d), transplanted %.1f (real %d); all folds ok in %.1f%% of replicates",
  mean(fcs$deaths), sum(d_real$Delta), mean(fcs$transplanted), sum(d_real$A < d_real$U),
  100 * mean(fcs$all_folds_ok)))

# ── 4. Summary ────────────────────────────────────────────────────────────────

summ <- est %>%
  mutate(truth = truth$beta[effect]) %>%
  group_by(method, effect) %>%
  summarise(
    truth      = first(truth),
    n_ok       = sum(!is.na(est)),
    mean_est   = mean(est, na.rm = TRUE),
    bias       = mean(est - truth, na.rm = TRUE),
    emp_sd     = sd(est, na.rm = TRUE),
    mean_se    = mean(se, na.rm = TRUE),
    rmse       = sqrt(mean((est - truth)^2, na.rm = TRUE)),
    cover      = mean(abs(est - truth) <= 1.96 * se, na.rm = TRUE),
    cover_info = mean(abs(est - truth) <= 1.96 * se_info, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(method, names(METHODS)), match(effect, names(EFFECTS)))
write.csv(summ, file.path(OUTPUT_DIR, "plasmode_summary.csv"), row.names = FALSE)
print(as.data.frame(summ), digits = 3)

# ── 5. Table ──────────────────────────────────────────────────────────────────
# Columns: bias and 95% coverage for each effect. TV-CSL coverage is
# information / sandwich, as in the main simulation tables.

cell <- function(m, e) {
  r <- summ[summ$method == m & summ$effect == e, ]
  cov <- if (m %in% c("tvcsl_marg", "tvcsl_tvp"))
    sprintf("%.2f/%.2f", r$cover_info, r$cover) else sprintf("%.2f", r$cover)
  c(sprintf("%.3f", r$bias), cov)
}

head_eff <- paste(sprintf("\\multicolumn{2}{c}{%s ($%.2f$)}", EFFECTS, truth$beta),
                  collapse = " & ")
tab <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Semi-Synthetic Bias and Coverage at the Stanford Sample Size}",
  "\\label{tab:stanford-plasmode}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  paste0(" & ", head_eff, " \\\\"),
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  paste0("Method", strrep(" & Bias & Coverage", 3), " \\\\"),
  "\\midrule"
)
for (m in names(METHODS)) {
  tab <- c(tab, paste0(METHODS[[m]], " & ",
                       paste(unlist(lapply(names(EFFECTS), cell, m = m)), collapse = " & "),
                       " \\\\"))
}
tab <- c(tab, "\\bottomrule", "\\end{tabular}", "}", "\\par\\smallskip",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize",
  sprintf("Results over %d replicates of $n = %d$. True values (in parentheses) are the",
          nrow(fcs), nrow(d_real)),
  "time-varying Cox estimates from the real data. Covariates and censoring dates are the",
  "real ones; transplant and death times are redrawn from Cox models fitted to the real data.",
  "Coverage is of nominal 95\\% intervals; TV-CSL entries are information/sandwich.",
  "\\end{minipage}",
  "\\end{table}")
writeLines(tab, file.path(OUTPUT_DIR, "tab_stanford_plasmode.tex"))

message("Plasmode outputs written to ", OUTPUT_DIR)
