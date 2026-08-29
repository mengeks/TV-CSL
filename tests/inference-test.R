suppressPackageStartupMessages(library(tidyverse))
source("R/old/TV-CSL.R")
source("R/old/data-handler.R")
source("R/datagen-helper.R")
source("scripts/TV-CSL/time-varying-estimate.R")

set.seed(42)

# ---- Generate small dataset -----------------------------------------------
n        <- 500
K        <- 2
eta_type <- "linear"
HTE_type <- "linear"
p        <- 3

datagen_params <- list(
  is_time_varying       = TRUE,
  light_censoring       = FALSE,
  lambda_C              = 0.1,
  p                     = p,
  X_distribution        = "normal",
  X_cov_type            = "identity",
  tx_difficulty         = "simple",
  linear_HTE_multiplier = 1,
  linear_intercept      = 0,
  linear_slope_multiplier = 2.5
)

data_dir <- here::here("data", paste0(eta_type, "_", HTE_type))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Use iteration 999 to avoid collision with experiment runs
i <- 999
generate_and_save_data(i = i,     n = n, path_for_sim_data = data_dir,
                       params = c(datagen_params, list(eta_type = eta_type, HTE_type = HTE_type)),
                       verbose = 0)
generate_and_save_data(i = i+100, n = n, path_for_sim_data = data_dir,
                       params = c(datagen_params, list(eta_type = eta_type, HTE_type = HTE_type)),
                       verbose = 0)

single_data <- read_single_simulation_data(n=n, i=i,     eta_type=eta_type, HTE_type=HTE_type)$data
test_data   <- read_single_simulation_data(n=n, i=i+100, eta_type=eta_type, HTE_type=HTE_type)$data

# True HTE coefficients for HTE_spec="linear", p=3:
# regressor = cbind(W, W*X.1, W*X.2, W*X.3)  --> beta_true ~ c(0, 1, 1, 1)
true_beta <- c(0, 1, 1, 1)

cat("\n========== S-Cox inference ==========\n")
train_data_cox <- preprocess_data(single_data, run_time_varying = TRUE)
scox_ret <- S_cox(train_data = train_data_cox, test_data = test_data,
                  regressor_spec = "linear", HTE_spec = "linear", verbose = 0)

cat("beta_HTE:", round(scox_ret$beta_HTE, 4), "\n")
cat("se_HTE:  ", round(scox_ret$se_HTE, 4), "\n")
stopifnot("S-Cox: se_HTE not finite" = all(is.finite(scox_ret$se_HTE)))
stopifnot("S-Cox: se_HTE not positive" = all(scox_ret$se_HTE > 0))

ci_lo <- scox_ret$beta_HTE - 1.96 * scox_ret$se_HTE
ci_hi <- scox_ret$beta_HTE + 1.96 * scox_ret$se_HTE
covered <- true_beta >= ci_lo & true_beta <= ci_hi
cat("95% CI covers true_beta:", covered, "(diagnostic only -- single sample)\n")

cat("\n========== S-Lasso linear/linear inference ==========\n")
slasso_ret <- S_lasso(train_data = train_data_cox, test_data = test_data,
                      regressor_spec = "linear", HTE_spec = "linear",
                      verbose = 0)

cat("beta_HTE:", round(slasso_ret$beta_HTE, 4), "\n")
cat("se_HTE:  ", round(slasso_ret$se_HTE, 4), "\n")
stopifnot("S-Lasso: se_HTE not finite" = all(is.finite(slasso_ret$se_HTE)))
stopifnot("S-Lasso: se_HTE not positive" = all(slasso_ret$se_HTE > 0))

ci_lo <- slasso_ret$beta_HTE - 1.96 * slasso_ret$se_HTE
ci_hi <- slasso_ret$beta_HTE + 1.96 * slasso_ret$se_HTE
covered <- true_beta >= ci_lo & true_beta <= ci_hi
cat("95% CI covers true_beta:", covered, "(diagnostic only -- single sample)\n")

cat("\n========== S-Lasso complex/linear inference (expect NA SEs) ==========\n")
slasso_cx_ret <- S_lasso(train_data = train_data_cox, test_data = test_data,
                         regressor_spec = "complex", HTE_spec = "linear",
                         verbose = 0)
cat("se_HTE:", slasso_cx_ret$se_HTE, "\n")
stopifnot("S-Lasso complex: se_HTE should be NA" = all(is.na(slasso_cx_ret$se_HTE)))

cat("\n========== TV-CSL inference (naive + sandwich) ==========\n")
train_data_pseudo      <- create_pseudo_dataset(survival_data = single_data)
single_data_orig       <- single_data %>%
  mutate(U_A = pmin(A, U), Delta_A = A <= U)

tvcsl_ret <- TV_CSL(
  train_data          = train_data_pseudo,
  test_data           = test_data,
  train_data_original = single_data_orig,
  HTE_type            = HTE_type,
  eta_type            = eta_type,
  K                   = K,
  prop_score_spec     = "cox-intercept-only",
  lasso_type          = "S-lasso",
  regressor_spec      = "linear",
  final_model_method  = "lasso_coxph",
  HTE_spec            = "linear",
  i                   = i
)

cat("beta_HTE:   ", round(tvcsl_ret$beta_HTE, 4), "\n")
cat("se_naive:   ", round(tvcsl_ret$se_naive, 4), "\n")
cat("se_sandwich:", round(tvcsl_ret$se_sandwich, 4), "\n")

stopifnot("TV-CSL: se_naive not finite"    = all(is.finite(tvcsl_ret$se_naive)))
stopifnot("TV-CSL: se_naive not positive"  = all(tvcsl_ret$se_naive > 0))
stopifnot("TV-CSL: se_sandwich not finite" = all(is.finite(tvcsl_ret$se_sandwich)))
stopifnot("TV-CSL: se_sandwich not positive" = all(tvcsl_ret$se_sandwich > 0))

ci_lo_naive <- tvcsl_ret$beta_HTE - 1.96 * tvcsl_ret$se_naive
ci_hi_naive <- tvcsl_ret$beta_HTE + 1.96 * tvcsl_ret$se_naive
ci_lo_sand  <- tvcsl_ret$beta_HTE - 1.96 * tvcsl_ret$se_sandwich
ci_hi_sand  <- tvcsl_ret$beta_HTE + 1.96 * tvcsl_ret$se_sandwich

covered_naive <- true_beta >= ci_lo_naive & true_beta <= ci_hi_naive
covered_sand  <- true_beta >= ci_lo_sand  & true_beta <= ci_hi_sand
cat("Naive CI covers true_beta:", covered_naive, "(diagnostic only -- single sample)\n")
cat("Sandwich CI covers true_beta:", covered_sand, "(diagnostic only -- single sample)\n")

cat("\n========== TV-CSL inference: cox-direct-event ==========\n")
tvcsl_de <- TV_CSL(
  train_data          = train_data_pseudo,
  test_data           = test_data,
  train_data_original = single_data_orig,
  HTE_type            = HTE_type,
  eta_type            = eta_type,
  K                   = K,
  prop_score_spec     = "cox-direct-event",
  lasso_type          = "S-lasso",
  regressor_spec      = "linear",
  final_model_method  = "lasso_coxph",
  HTE_spec            = "linear",
  i                   = i
)

cat("beta_HTE:   ", round(tvcsl_de$beta_HTE, 4), "\n")
cat("se_naive:   ", round(tvcsl_de$se_naive, 4), "\n")
cat("se_sandwich:", round(tvcsl_de$se_sandwich, 4), "\n")

stopifnot("cox-direct-event: se_naive not finite"      = all(is.finite(tvcsl_de$se_naive)))
stopifnot("cox-direct-event: se_naive not positive"    = all(tvcsl_de$se_naive > 0))
stopifnot("cox-direct-event: se_sandwich not finite"   = all(is.finite(tvcsl_de$se_sandwich)))
stopifnot("cox-direct-event: se_sandwich not positive" = all(tvcsl_de$se_sandwich > 0))

cat("\n========== TV-CSL inference: cox-risk-set-adjusted-oracle ==========\n")
tvcsl_oracle <- TV_CSL(
  train_data          = train_data_pseudo,
  test_data           = test_data,
  train_data_original = single_data_orig,
  HTE_type            = HTE_type,
  eta_type            = eta_type,
  K                   = K,
  prop_score_spec     = "cox-risk-set-adjusted-oracle",
  lasso_type          = "S-lasso",
  regressor_spec      = "linear",
  final_model_method  = "lasso_coxph",
  HTE_spec            = "linear",
  i                   = i
)

cat("beta_HTE:   ", round(tvcsl_oracle$beta_HTE, 4), "\n")
cat("se_naive:   ", round(tvcsl_oracle$se_naive, 4), "\n")
cat("se_sandwich:", round(tvcsl_oracle$se_sandwich, 4), "\n")

stopifnot("oracle: se_naive not finite"      = all(is.finite(tvcsl_oracle$se_naive)))
stopifnot("oracle: se_naive not positive"    = all(tvcsl_oracle$se_naive > 0))
stopifnot("oracle: se_sandwich not finite"   = all(is.finite(tvcsl_oracle$se_sandwich)))
stopifnot("oracle: se_sandwich not positive" = all(tvcsl_oracle$se_sandwich > 0))

cat("\n========== All assertions passed ==========\n")
