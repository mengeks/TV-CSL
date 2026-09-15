## copilot_analysis.R
## Produces two LaTeX table files for Section 6 of the paper:
##   writeups/tables/tab_copilot_ignore_time.tex
##   writeups/tables/tab_copilot_tvcsl.tex
##
## Uses survival::coxph and mgcv::gam directly; sources R/TV-CSL.R for
## package loading only.

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(mgcv)
  library(here)
})

# ── 1. Load data ──────────────────────────────────────────────────────────────

df_raw <- read.csv(
  here::here("scripts/law-firm-ai-analysis/survival_ccassist_SYNTHETIC.csv"),
  stringsAsFactors = FALSE
)

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

# Dummy-code adjustment covariates; make valid R names
mm <- model.matrix(
  ~ rank + office + practice_group + baseline_weekly_hours + cca_pre_washout - 1,
  data = df_base
)
xcols <- make.names(colnames(mm))
colnames(mm) <- xcols

df_orig <- df_base
for (j in xcols) df_orig[[j]] <- mm[, j]

# ── 3. Counting-process dataset (delayed entry handled) ───────────────────────
# W_i(t) = 1(A_i < t); convention: same-week event and treatment → untreated

make_tv_data <- function(d) {
  result <- vector("list", nrow(d))
  for (i in seq_len(nrow(d))) {
    r  <- d[i, , drop = FALSE]
    t0 <- r$entry_week
    t1 <- r$U
    ev <- r$Delta
    tx <- r$A          # Inf when never treated

    if (tx >= t1) {
      r$tstart <- t0; r$tstop <- t1; r$Delta <- ev; r$W <- 0L
      result[[i]] <- r
    } else if (tx <= t0) {
      r$tstart <- t0; r$tstop <- t1; r$Delta <- ev; r$W <- 1L
      result[[i]] <- r
    } else {
      r1 <- r; r2 <- r
      r1$tstart <- t0; r1$tstop <- tx; r1$Delta <- 0L; r1$W <- 0L
      r2$tstart <- tx; r2$tstop <- t1; r2$Delta <- ev; r2$W <- 1L
      result[[i]] <- rbind(r1, r2)
    }
  }
  bind_rows(result) %>% filter(tstop > tstart)
}

df_tv <- make_tv_data(df_orig)

# ── 4. Fixed-treatment Cox (ignore treatment time) ────────────────────────────
# used_cca is the ever-treated indicator; adjusted for rank, office, pg, hours, cca_pre

adj_str <- paste(xcols, collapse = " + ")

fit_fixed <- coxph(
  as.formula(paste(
    "Surv(entry_week, U, Delta) ~ used_cca + used_cca:is_junior +", adj_str
  )),
  data = df_orig, ties = "breslow"
)

# ── 5. Time-varying Cox ───────────────────────────────────────────────────────

fit_tvcox <- coxph(
  as.formula(paste(
    "Surv(tstart, tstop, Delta) ~ W + W:is_junior +", adj_str
  )),
  data = df_tv, ties = "breslow"
)

# ── 6. Extract estimates from a coxph fit ─────────────────────────────────────
# Returns formatted list(senior, inter, junior) each with $ci and $pval

extract_cox <- function(fit, beta_nm, inter_nm) {
  cf <- coef(fit)
  vc <- vcov(fit)

  b   <- cf[beta_nm];  se  <- sqrt(vc[beta_nm,  beta_nm])
  bi  <- cf[inter_nm]; sei <- sqrt(vc[inter_nm, inter_nm])
  bs  <- b + bi
  ses <- sqrt(vc[beta_nm, beta_nm] + vc[inter_nm, inter_nm] +
              2 * vc[beta_nm, inter_nm])

  fmt <- function(est, se_est) {
    pv <- 2 * pnorm(abs(est / se_est), lower.tail = FALSE)
    list(
      ci   = sprintf("%.3f (%.3f, %.3f)",
                     est, est - 1.96 * se_est, est + 1.96 * se_est),
      pval = if (pv < 0.001) "$<$0.001" else sprintf("%.3f", pv)
    )
  }
  list(senior = fmt(b, se), inter = fmt(bi, sei), junior = fmt(bs, ses))
}

res_fixed  <- extract_cox(fit_fixed,  "used_cca", "used_cca:is_junior")
res_tvcox  <- extract_cox(fit_tvcox,  "W",         "W:is_junior")

# ── 7. TV-CSL: direct K-fold DML implementation ───────────────────────────────
# Nuisance: marginal baseline Cox (ignoring treatment).
# Propensity: either marginal GH-style or time-varying event-regression.
# Final model: coxph with offset(nu) and orthogonalized regressors.

run_tvcsl <- function(df_tv, df_orig, xcols, prop_type = "marginal",
                      K = 5L, seed = 42L) {
  set.seed(seed)
  n     <- nrow(df_orig)
  folds <- cut(seq_len(n), breaks = K, labels = FALSE)

  beta_list  <- vector("list", K)
  info_list  <- vector("list", K)
  score_list <- vector("list", K)

  f_base <- as.formula(
    paste("Surv(tstart, tstop, Delta) ~", paste(xcols, collapse = " + "))
  )

  for (k in seq_len(K)) {
    in_test  <- folds == k
    test_ids <- df_orig$id[in_test]

    train_orig <- df_orig[!in_test, ]
    train_tv   <- df_tv[!(df_tv$id %in% test_ids), ]
    test_tv    <- df_tv[df_tv$id %in% test_ids, ]

    # Nuisance: baseline Cox (η̂₀, ignoring treatment)
    m_base <- coxph(f_base, data = train_tv, ties = "breslow")
    nu_hat <- as.vector(predict(m_base, newdata = test_tv, type = "lp"))

    # Propensity â(t, X)
    if (prop_type == "marginal") {
      # GH: Cox on treatment initiation among event subjects only
      ev_orig <- train_orig %>% filter(Delta == 1)
      f_prop  <- as.formula(
        paste("Surv(U_A, Delta_A) ~", paste(xcols, collapse = " + "))
      )
      m_prop  <- coxph(f_prop, data = ev_orig, ties = "breslow")
      alpha   <- coef(m_prop)
      alpha[is.na(alpha)] <- 0  # empty category in fold → zero contribution
      lp_p    <- as.vector(as.matrix(test_tv[, xcols]) %*% alpha)
      # P(A < t | X) approximated via exponential baseline
      a_hat   <- 1 - exp(-test_tv$tstop * exp(lp_p))
    } else {
      # Time-varying: logistic/GAM on W at event time
      ev_tv  <- train_tv %>% filter(Delta == 1)
      cov_s  <- paste(xcols, collapse = " + ")
      f_gam  <- as.formula(paste("W ~ s(tstop, k = 5) +", cov_s))
      f_glm  <- as.formula(paste("W ~ tstop +", cov_s))
      m_prop <- tryCatch(
        mgcv::gam(f_gam, family = binomial, data = ev_tv),
        error = function(e) glm(f_glm, family = binomial, data = ev_tv)
      )
      a_hat  <- as.vector(
        plogis(predict(m_prop, newdata = test_tv, type = "link"))
      )
    }
    a_hat <- pmax(0.01, pmin(0.99, a_hat))

    # Orthogonalized regressors: q(t) = (W(t) - a(t,X)) * p(X),  p(X) = (1, is_junior)
    d        <- test_tv
    d$nu_X   <- nu_hat
    d$q_int  <- d$W - a_hat
    d$q_jun  <- (d$W - a_hat) * d$is_junior

    m_final <- coxph(
      Surv(tstart, tstop, Delta) ~ q_int + q_jun + offset(nu_X),
      data = d, ties = "breslow"
    )
    beta_list[[k]] <- coef(m_final)

    # Subject-level score residuals for sandwich SE
    vcov_k    <- vcov(m_final)
    info_list[[k]] <- solve(vcov_k)
    sr <- residuals(m_final, type = "score")
    if (!is.matrix(sr)) sr <- matrix(sr, ncol = 1)
    sr_df    <- as.data.frame(sr)
    sr_df$id <- d$id
    score_list[[k]] <- sr_df %>%
      group_by(id) %>%
      summarise(across(everything(), sum), .groups = "drop") %>%
      dplyr::select(-id) %>%
      as.matrix()
  }

  beta_HTE  <- rowMeans(do.call(cbind, beta_list))
  all_sc    <- do.call(rbind, score_list)
  n_subj    <- nrow(all_sc)
  Sigma_hat <- crossprod(all_sc) / n_subj
  J_hat     <- Reduce("+", info_list) / n_subj
  J_inv     <- solve(J_hat)
  Omega_hat <- J_inv %*% Sigma_hat %*% J_inv
  se_sw     <- sqrt(diag(Omega_hat) / n_subj)

  list(beta_HTE = beta_HTE, Omega_hat = Omega_hat,
       se_sw = se_sw, n_subj = n_subj)
}

message("Fitting TV-CSL with marginal propensity ...")
res_marg <- run_tvcsl(df_tv, df_orig, xcols, prop_type = "marginal")

message("Fitting TV-CSL with time-varying propensity ...")
res_tvp  <- run_tvcsl(df_tv, df_orig, xcols, prop_type = "time-varying")

# ── 8. Extract TV-CSL estimates ───────────────────────────────────────────────

extract_tvcsl <- function(res) {
  beta      <- res$beta_HTE    # c(senior, inter)
  Omega     <- res$Omega_hat
  n         <- res$n_subj

  fmt_sw <- function(b, idx, c_vec = NULL) {
    if (is.null(c_vec)) {
      se <- res$se_sw[idx]
    } else {
      se <- sqrt(drop(t(c_vec) %*% (Omega / n) %*% c_vec))
    }
    sprintf("%.3f (%.3f, %.3f)", b, b - 1.96 * se, b + 1.96 * se)
  }

  b_senior <- beta[1]
  b_inter  <- beta[2]
  b_junior <- beta[1] + beta[2]

  list(
    senior = fmt_sw(b_senior, 1L),
    inter  = fmt_sw(b_inter,  2L),
    junior = fmt_sw(b_junior, NA, c(1, 1))
  )
}

res_marg_fmt <- extract_tvcsl(res_marg)
res_tvp_fmt  <- extract_tvcsl(res_tvp)

# ── 9. Write tables ───────────────────────────────────────────────────────────

tables_dir <- here::here("writeups/tables")
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
