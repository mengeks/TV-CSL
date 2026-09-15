# tests/test-tvcsl.R
# Basic unit tests for TV-CSL core functions.
# Run from project root: Rscript tests/test-tvcsl.R
# All tests use small n to run quickly (no SLURM needed).

suppressPackageStartupMessages({
  library(survival)
  library(tidyverse)
})
source("R/datagen-helper.R")
source("scripts/TV-CSL/time-varying-estimate.R")

pass <- function(msg) cat(sprintf("PASS: %s\n", msg))
fail <- function(msg) stop(sprintf("FAIL: %s", msg))

# ---- 1. calculate_eX: output in [0, 1] ----------------------------------------
{
  set.seed(42)
  alpha <- c(0.5, -0.3, 0.2)
  X     <- matrix(rnorm(30), nrow = 10, ncol = 3)
  t     <- runif(10, min = 0.01, max = 5)
  vals  <- calculate_eX(alpha, X, t)
  if (!all(vals >= 0 & vals <= 1)) fail("calculate_eX out of [0,1]")
  pass("calculate_eX values in [0,1]")
}

# ---- 2. generate_simulated_data: correct column names and dimensions ----------
{
  dat <- generate_simulated_data(
    n = 50, is_time_varying = TRUE, light_censoring = FALSE,
    lambda_C = 0.1, p = 3, eta_type = "linear", HTE_type = "linear",
    seed_value = 1
  )
  required_cols <- c("A", "U", "Delta", "id", "HTE", "X.1", "X.2", "X.3")
  missing <- setdiff(required_cols, names(dat))
  if (length(missing) > 0) fail(paste("Missing columns:", paste(missing, collapse = ", ")))
  if (nrow(dat) != 50) fail(paste("Expected 50 rows, got", nrow(dat)))
  pass("generate_simulated_data columns and dimensions")
}

# ---- 3. Delta is binary (0/1) and A, U are non-negative ----------------------
{
  dat <- generate_simulated_data(
    n = 100, is_time_varying = TRUE, lambda_C = 0.1, p = 3,
    eta_type = "linear", HTE_type = "linear", seed_value = 7
  )
  if (!all(dat$Delta %in% c(0, 1))) fail("Delta is not binary")
  if (any(dat$A < 0, na.rm = TRUE)) fail("Treatment time A has negative values")
  if (any(dat$U <= 0)) fail("Event time U has non-positive values")
  pass("Delta binary, A and U non-negative")
}

# ---- 4. get_X_vars helper (requires TV-CSL.R to be sourced) ------------------
{
  df <- data.frame(X.1 = 1:3, X.2 = 4:6, Z = 7:9)
  # NULL → all X. columns
  result_null <- get_X_vars(df, NULL)
  if (!setequal(names(result_null), c("X.1", "X.2"))) fail("get_X_vars(NULL) wrong cols")
  # explicit → only specified
  result_spec <- get_X_vars(df, "X.1")
  if (!identical(names(result_spec), "X.1")) fail("get_X_vars('X.1') wrong cols")
  pass("get_X_vars selects correct columns")
}

# ---- 5. TV_CSL_nuisance: propensity scores in (0, 1) -------------------------
{
  set.seed(99)
  dat <- generate_simulated_data(
    n = 80, is_time_varying = TRUE, lambda_C = 0.1, p = 3,
    eta_type = "linear", HTE_type = "linear", seed_value = 99
  )
  pseudo <- create_pseudo_dataset(dat)
  dat    <- dat |> dplyr::mutate(U_A = pmin(A, U), Delta_A = as.numeric(A <= U))

  # Use the first 40 subjects as train, last 40 as test
  train_ids <- dat$id[1:40]
  test_ids  <- dat$id[41:80]
  fold_train <- pseudo |> dplyr::filter(id %in% train_ids)
  fold_test  <- pseudo |> dplyr::filter(id %in% test_ids)
  train_orig <- dat    |> dplyr::filter(id %in% train_ids)

  ret <- TV_CSL_nuisance(
    fold_train          = fold_train,
    fold_test           = fold_test,
    train_data_original = train_orig,
    prop_score_spec     = "cox-linear-censored-only",
    lasso_type          = "S-lasso",
    regressor_spec      = "linear",
    HTE_spec            = "linear"
  )

  a_t_X <- ret$fold_test_final$a_t_X
  if (any(is.na(a_t_X))) fail("TV_CSL_nuisance produced NA propensity scores")
  if (!all(a_t_X >= 0 & a_t_X <= 1)) fail("TV_CSL_nuisance propensity outside [0,1]")
  pass("TV_CSL_nuisance propensity scores in [0,1]")
}

cat("\nAll tests passed.\n")
