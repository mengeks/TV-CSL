library(jsonlite)
suppressPackageStartupMessages(library(tidyverse))

source("R/old/data-handler.R")
source("scripts/TV-CSL/time-varying-estimate.R")


#' Run a Single Iteration of the Experiment and Save Results to CSV
#'
#' @param i          Iteration number (maps to a pre-generated dataset and a PRNG seed).
#' @param json_file  Path to the JSON configuration file.
#' @param eta_type   Baseline hazard type ("linear" or "non-linear").
#' @param HTE_type   Treatment effect type ("constant", "linear", ...).
#' @param n          Sample size.
#' @param verbose    Verbosity level (0 = silent, 1 = progress, 2 = detailed).
run_experiment_iteration <-
  function(i, json_file, eta_type, HTE_type, n, verbose = 0) {

  if (verbose >= 1)
    message("Running iteration ", i)

  config  <- fromJSON(json_file)
  methods <- config$methods
  K       <- ifelse(is.null(config$K), 5, config$K)

  input_setting <- paste0(eta_type, "_", HTE_type)
  seed_value    <- 123 + 11 * i
  set.seed(seed_value)

  if (verbose >= 2) {
    message("Configuration Parameters:")
    message("n: ", n, "\neta_type: ", eta_type, "\nHTE_type: ", HTE_type,
            "\nK: ", K, "\nSeed: ", seed_value)
  }

  # ---- Load data ------------------------------------------------------------
  start_time  <- Sys.time()
  loaded_data <- read_single_simulation_data(
    n = n, i = i, eta_type = eta_type, HTE_type = HTE_type
  )
  test_data <- read_single_simulation_data(
    n = n, i = i + 100, eta_type = eta_type, HTE_type = HTE_type
  )$data
  end_time <- Sys.time()

  if (verbose >= 1)
    message("Time to load dataset: ",
            as.numeric(difftime(end_time, start_time, units = "secs")), " seconds")

  single_data <- loaded_data$data

  # ---- Cox ------------------------------------------------------------------
  beta_estimates_cox <- mse_estimates_cox <- time_taken_cox <- list()
  is_running_cox <- !is.null(methods$cox) && methods$cox$enabled

  if (is_running_cox) {
    start_time <- Sys.time()
    cox_results <- run_cox_estimation(
      single_data, methods$cox, HTE_type = HTE_type, eta_type = eta_type
    )
    end_time <- Sys.time()

    for (config_name in names(cox_results)) {
      beta_estimates_cox[[config_name]] <- cox_results[[config_name]]$beta_estimate
      time_taken_cox[[config_name]]     <- cox_results[[config_name]]$time_taken
      mse_estimates_cox[[config_name]]  <- calculate_mse(
        beta_estimates_cox[[config_name]], n, i, HTE_type, eta_type
      )
      print(paste0("MSE of Cox config '", config_name, "': "))
      print(mse_estimates_cox[[config_name]])
    }

    if (verbose >= 1)
      message("Time to run Cox: ",
              as.numeric(difftime(end_time, start_time, units = "secs")), " seconds")
  }

  result_df_cox <- data.frame(
    Method        = rep("Cox", length(mse_estimates_cox)),
    Specification = names(mse_estimates_cox),
    MSE_Estimate  = unlist(mse_estimates_cox),
    Time_Taken    = unlist(time_taken_cox)
  )

  # ---- S-lasso --------------------------------------------------------------
  time_taken_lasso <- mse_estimates_lasso <- list()
  is_running_lasso <- !is.null(methods$lasso) && methods$lasso$enabled

  if (is_running_lasso) {
    start_time <- Sys.time()
    lasso_results <- run_lasso_estimation(
      single_data   = single_data,
      i             = i,
      methods_lasso = methods$lasso,
      HTE_type      = HTE_type,
      eta_type      = eta_type
    )
    end_time <- Sys.time()

    for (config_name in names(lasso_results)) {
      time_taken_lasso[[config_name]]    <- lasso_results[[config_name]]$time_taken
      mse_estimates_lasso[[config_name]] <- lasso_results[[config_name]]$MSE
    }

    if (verbose >= 1)
      message("Time to run lasso: ",
              as.numeric(difftime(end_time, start_time, units = "secs")), " seconds")
    if (verbose >= 2) { message("Lasso Results:"); print(mse_estimates_lasso) }
  }

  result_df_lasso <- data.frame(
    Method        = rep("Lasso", length(mse_estimates_lasso)),
    Specification = names(mse_estimates_lasso),
    MSE_Estimate  = unlist(mse_estimates_lasso),
    Time_Taken    = unlist(time_taken_lasso)
  )

  # ---- S-Cox ----------------------------------------------------------------
  time_taken_s_cox <- mse_estimates_s_cox <- list()
  is_running_s_cox <- !is.null(methods$s_cox) && methods$s_cox$enabled

  if (is_running_s_cox) {
    start_time <- Sys.time()
    s_cox_results <- run_s_cox_estimation(
      single_data   = single_data,
      i             = i,
      methods_s_cox = methods$s_cox,
      HTE_type      = HTE_type,
      eta_type      = eta_type
    )
    end_time <- Sys.time()

    for (config_name in names(s_cox_results)) {
      time_taken_s_cox[[config_name]]    <- s_cox_results[[config_name]]$time_taken
      mse_estimates_s_cox[[config_name]] <- s_cox_results[[config_name]]$MSE
    }

    if (verbose >= 1)
      message("Time to run S-Cox: ",
              as.numeric(difftime(end_time, start_time, units = "secs")), " seconds")
    if (verbose >= 2) { message("S-Cox Results:"); print(mse_estimates_s_cox) }
  }

  result_df_s_cox <- data.frame(
    Method        = rep("S-Cox", length(mse_estimates_s_cox)),
    Specification = names(mse_estimates_s_cox),
    MSE_Estimate  = unlist(mse_estimates_s_cox),
    Time_Taken    = unlist(time_taken_s_cox)
  )

  # ---- TV-CSL ---------------------------------------------------------------
  mse_estimates_TV_CSL <- time_taken_TV_CSL <- list()
  is_running_TV_CSL <- !is.null(methods$TV_CSL) && methods$TV_CSL$enabled

  if (is_running_TV_CSL) {
    start_time <- Sys.time()

    temp_result_csv_file <- paste0(
      "scripts/TV-CSL/results/temp/",
      input_setting, "-n_", n,
      "-iteration_", i, "-seed_", seed_value, ".csv"
    )

    TV_CSL_results <- run_TV_CSL_estimation(
      train_data_original  = single_data,
      test_data            = test_data,
      methods_TV_CSL       = methods$TV_CSL,
      i                    = i,
      K                    = K,
      HTE_type             = HTE_type,
      eta_type             = eta_type,
      temp_result_csv_file = temp_result_csv_file
    )
    end_time <- Sys.time()

    for (config_name in names(TV_CSL_results)) {
      time_taken_TV_CSL[[config_name]]    <- TV_CSL_results[[config_name]]$time_taken
      mse_estimates_TV_CSL[[config_name]] <- TV_CSL_results[[config_name]]$MSE
    }

    if (verbose >= 1)
      message("Time to run TV-CSL: ",
              as.numeric(difftime(end_time, start_time, units = "secs")), " seconds")
    if (verbose >= 2) { message("TV-CSL Results:"); print(mse_estimates_TV_CSL) }
  }

  result_df_TV_CSL <- data.frame(
    Method        = rep("TV_CSL", length(mse_estimates_TV_CSL)),
    Specification = names(mse_estimates_TV_CSL),
    MSE_Estimate  = unlist(mse_estimates_TV_CSL),
    Time_Taken    = unlist(time_taken_TV_CSL)
  )

  # ---- Combine and save -----------------------------------------------------
  result_df <- rbind(
    result_df_cox,
    result_df_lasso,
    result_df_s_cox,
    result_df_TV_CSL
  )

  result_csv_file <- generate_output_path(
    results_dir       = RESULTS_DIR,
    is_running_cox    = is_running_cox,
    is_running_lasso  = is_running_lasso,
    is_running_s_cox  = is_running_s_cox,
    is_running_TV_CSL = is_running_TV_CSL,
    eta_type          = eta_type,
    HTE_type          = HTE_type,
    n                 = n,
    i                 = i,
    seed_value        = seed_value
  )

  write.csv(result_df, result_csv_file, row.names = FALSE)

  if (verbose >= 1)
    message("Results for iteration ", i, " saved to ", result_csv_file)
}


calculate_mse <- function(beta_estimate, n, i, HTE_type, eta_type) {
  test_data <- read_single_simulation_data(
    n = n, i = i + 100,
    eta_type = eta_type, HTE_type = HTE_type
  )$data

  if (HTE_type == "linear") {
    HTE_est <- beta_estimate[1] * test_data$X.1 + beta_estimate[2] * test_data$X.10
  } else if (HTE_type == "constant") {
    HTE_est <- rep(beta_estimate, nrow(test_data))
  }

  mean((HTE_est - test_data$HTE)^2)
}
