source("R/old/data-handler.R")
library(jsonlite)
library(ggplot2)
library(dplyr)

#' Load Experiment Configuration
#'
#' Loads the experiment configuration from a JSON file.
#' 
#' @param json_file The path to the JSON configuration file.
#' @return A list containing the experiment configuration.
load_experiment_config <- function(json_file) {
  config <- fromJSON(json_file)
  return(config)
}


# Helper function to read individual CSV files
read_single_iteration_result <- function(csv_file) {
  df <- read.csv(csv_file)
  return(df)
}

get_is_running_flags <- function(methods) {
  list(
    is_running_cox    = !is.null(methods$cox)    && methods$cox$enabled,
    is_running_lasso  = !is.null(methods$lasso)  && methods$lasso$enabled,
    is_running_s_cox  = !is.null(methods$s_cox)  && methods$s_cox$enabled,
    is_running_TV_CSL = !is.null(methods$TV_CSL) && methods$TV_CSL$enabled
  )
}

# Function to aggregate all CSV files and calculate final metrics
process_all_iterations <- function(config,eta_type, HTE_type, results_dir, n) {
  methods <- config$methods
  # eta_type <- config$eta_type
  # HTE_type <- config$HTE_type
  R <- config$R
  
  running_flags <- get_is_running_flags(methods)

  is_running_cox    <- running_flags$is_running_cox
  is_running_lasso  <- running_flags$is_running_lasso
  is_running_s_cox  <- running_flags$is_running_s_cox
  is_running_TV_CSL <- running_flags$is_running_TV_CSL

  
  all_results <- list()
  all_times <- list()
  
  for (i in 1:R) {
    seed_value <- 123 + 11 * i
    
    result_csv_file <- generate_output_path(
      results_dir       = results_dir,
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
    
    
    if (file.exists(result_csv_file)) {
      iteration_result <- read_single_iteration_result(result_csv_file)
      all_results[[i]] <- iteration_result
    } else {
      cat("Warning: File not found for iteration", i, "\n")
    }
  }
  
  combined_results <- do.call(rbind, all_results)
  combined_results$Specification <- 
    gsub("S_lasso", "S-lasso", combined_results$Specification)
  combined_results$Specification <- 
    gsub("T_lasso", "T-lasso", combined_results$Specification)
  
  
  aggregated_metrics <- combined_results %>%
    group_by(Method, Specification) %>%
    summarise(
      MSE          = mean(MSE_Estimate),
      MCSE_MSE     = sqrt(var(MSE_Estimate) / n()),
      Mean_Time_s  = mean(Time_Taken),
      n_iterations = n()
    )
  
  return(aggregated_metrics)
}


get_method_setting <- function(config){
  methods <- config$methods
  is_running_cox    <- !is.null(methods$cox)    && methods$cox$enabled
  is_running_lasso  <- !is.null(methods$lasso)  && methods$lasso$enabled
  is_running_s_cox  <- !is.null(methods$s_cox)  && methods$s_cox$enabled
  is_running_TV_CSL <- !is.null(methods$TV_CSL) && methods$TV_CSL$enabled
  method_setting <- paste0(
    ifelse(is_running_cox,    "cox_",    ""),
    ifelse(is_running_lasso,  "lasso_",  ""),
    ifelse(is_running_s_cox,  "s-cox_",  ""),
    ifelse(is_running_TV_CSL, "TV-CSL_", "")
  )
  return(method_setting)
}

get_table_csv_file <- function(config, 
                               eta_type,
                               HTE_type, 
                               output_csv_dir = "scripts/TV-CSL/tables-and-plots",
                               file_suffix = "_est_quality.csv" ){
  eta_type_folder_name <- 
    paste0(eta_type, "_", HTE_type)
  
  
  method_setting <- get_method_setting(config)
  
  table_csv_file <- 
    file.path(output_csv_dir, 
              paste0(method_setting, 
                     eta_type_folder_name, 
                     file_suffix ) )
  return(table_csv_file)
}

#' Combine Individual MSE Results to a table by method averages
#' 
#' @description 
#' Takes DGP configuration a JSON configuration file, loop over
#' different sample sizes, gather the estimation results and saves the aggregated metrics to a CSV file.
#' 
#' @param json_file character, Path to the JSON configuration file containing 
#'   experiment parameters
#' @param n_list numeric vector, Sample sizes to process (default: c(200, 500, 1000, 2000))
#' @param output_csv_dir character, Directory path where the output CSV will be saved
#'   (default: "scripts/TV-CSL/tables-and-plots")
#' 
#' @return None (invisible). Writes results to a CSV file and prints confirmation message.
#' 
#' @details 
#' The function performs the following steps:
#' 1. Loads experiment configuration from JSON
#' 2. Processes results for each sample size in n_list
#' 3. Combines results into a single data frame
#' 4. Saves aggregated results to CSV
#' 
#' @dependencies 
#' Requires load_experiment_config(), process_all_iterations(), and get_table_csv_file()
#' 
#' @examples 
#' \dontrun{
#' process_results_to_csv(
#'   json_file = "path/to/config.json",
#'   n_list = c(100, 300, 500),
#'   output_csv_dir = "output/directory"
#' )
#' }
#' 
#' @export
process_results_to_csv <- function(json_file,
                                   eta_type,
                                   HTE_type,
                                   n_list = c(200, 500, 1000, 2000),
                                   output_csv_dir = "scripts/TV-CSL/tables-and-plots",
                                   results_dir=RESULTS_DIR) {
  
  config <- load_experiment_config(json_file)
  
  aggregated_metrics <- NULL
  
  for (n in n_list){
    aggregated_metrics_n <- 
      process_all_iterations(
        config=config, 
        eta_type = eta_type,
        HTE_type = HTE_type,
        results_dir=results_dir,
        n = n
      )
    aggregated_metrics_n$n <- n
    if (is.null(aggregated_metrics)){
      aggregated_metrics <- aggregated_metrics_n
    }else{
      aggregated_metrics <- rbind(aggregated_metrics, aggregated_metrics_n)
    }
  }
  
  
  table_csv_file <- get_table_csv_file(config = config,
                                       eta_type = eta_type,
                                       HTE_type = HTE_type,
                                       output_csv_dir = output_csv_dir)
  
  write.csv(
    aggregated_metrics, 
    table_csv_file, 
    row.names = FALSE)
  
  cat("Aggregated results saved to", table_csv_file, "\n")
}

get_plot_file <- function(config, 
                          output_csv_dir = "scripts/TV-CSL/tables-and-plots"){
  eta_type_folder_name <- 
    paste0(config$eta_type, "_", config$HTE_type)
  
  
  method_setting <- get_method_setting(config)
  
  plot_file <- 
    file.path(output_csv_dir, 
              paste0(method_setting, 
                     eta_type_folder_name, 
                     "_MSE_plots.png") )
  return(plot_file)
}

# Function to calculate sample skewness
calc_skewness <- function(T_values, mean_T, sd_T) {
  R <- length(T_values)
  sum((T_values - mean_T)^3) / (R * sd_T^3)
}

# Function to calculate sample kurtosis
calc_kurtosis <- function(T_values, mean_T, sd_T) {
  R <- length(T_values)
  sum((T_values - mean_T)^4) / (R * sd_T^4)
}

# Function to calculate MCSE of MSE
calc_mcse_mse <- function(T_values, theta) {
  R <- length(T_values)
  mean_T <- mean(T_values)
  sd_T <- sd(T_values)
  
  skewness_T <- calc_skewness(T_values, mean_T, sd_T)
  kurtosis_T <- calc_kurtosis(T_values, mean_T, sd_T)
  
  term1 <- sd_T^4 * (kurtosis_T - 1)
  term2 <- 4 * sd_T^3 * skewness_T * (mean_T - theta)
  term3 <- 4 * sd_T^2 * (mean_T - theta)^2
  
  mcse_mse <- sqrt((1 / R) * (term1 + term2 + term3))
  return(mcse_mse)
}


evaluate_HTE_metrics <- function(
    json_file, 
    results_dir, 
    eta_type, 
    HTE_type, 
    n_list, 
    HTE_spec,
    output_csv_dir = "scripts/TV-CSL/tables-and-plots") {
  # Load the experiment configuration
  config <- load_experiment_config(json_file)
  methods <- config$methods
  
  # Get running flags
  running_flags <- get_is_running_flags(methods)
  is_running_cox    <- running_flags$is_running_cox
  is_running_lasso  <- running_flags$is_running_lasso
  is_running_s_cox  <- running_flags$is_running_s_cox
  is_running_TV_CSL <- running_flags$is_running_TV_CSL

  # Initialize aggregated metrics
  aggregated_metrics <- NULL

  for (n in n_list) {
    # Generate method settings and output folder
    method_setting <- paste0(
      ifelse(is_running_cox,    "cox_",    ""),
      ifelse(is_running_lasso,  "lasso_",  ""),
      ifelse(is_running_s_cox,  "s-cox_",  ""),
      ifelse(is_running_TV_CSL, "TV-CSL_", "")
    )
    
    output_folder <- generate_output_folder(
      results_dir = results_dir,
      method_setting = method_setting, 
      eta_type = eta_type, 
      HTE_type = HTE_type, 
      n = n
    )
    
    # Path to the HTE coefficient file
    HTE_coef_file <- paste0("HTE-spec-", HTE_spec, "_beta-HTE.csv")
    result_csv_file <- file.path(output_folder, HTE_coef_file)
    
    # Check if the file exists
    if (!file.exists(result_csv_file)) {
      cat("Warning: File not found for n =", n, "\n")
      next
    }
    
    # Load HTE coefficient data
    df_HTE_coef <- read.csv(result_csv_file)
    
    
    
    # Read a single test set
    test_data <- read_single_simulation_data(
      n = 2000, 
      i = 101, 
      eta_type = eta_type,
      HTE_type = HTE_type
    )$data
    
    # Extract HTE_true and transform test data
    HTE_true <- test_data$HTE
    X_HTE_test <- transform_X(
      single_data = test_data,
      transform_spec = HTE_spec
    )
    test_regressor_HTE <- cbind(1, X_HTE_test)
    # print(paste0("dim(test_regressor_HTE)=",dim(test_regressor_HTE)))
    
    p_regressor_HTE <- ncol(test_regressor_HTE)
    n_col_df <- ncol(df_HTE_coef)
    
    
    ## Compute the first stage average
    compute_first_stage_average <- function(df_HTE_coef, p_regressor_HTE, n_col_df) {
      
      coefficient_colnames <- colnames(df_HTE_coef)[(n_col_df - p_regressor_HTE + 1):n_col_df]
      # print(coefficient_colnames)
      
      first_stage_df <- df_HTE_coef %>% 
        filter(stage == "first")
      
      first_stage_avg <- first_stage_df %>%
        group_by(iteration, lasso_type, eta_spec, HTE_spec, prop_score_spec, stage) %>%
        summarise(across(all_of(coefficient_colnames), mean, na.rm = TRUE)) %>%
        ungroup()
      
      first_stage_avg$k <- 0
      return(first_stage_avg)
    }
    
    # This is average over k
    first_stage_avg <- 
      compute_first_stage_average(df_HTE_coef=df_HTE_coef, 
                                  p_regressor_HTE=p_regressor_HTE,
                                  n_col_df = n_col_df)
    

    df_HTE_coef <- rbind(df_HTE_coef,first_stage_avg)
    
    
    # Compute MSE for each row of df_HTE_coef
    # HTE_true is a vector of length 2000
    df_HTE_coef$MSE <- apply(df_HTE_coef, 1, function(row) {
      beta_HTE <- as.numeric(row[(n_col_df - p_regressor_HTE + 1):n_col_df])
      HTE_est <- as.vector(test_regressor_HTE %*% beta_HTE)
      mean((HTE_est - HTE_true)^2)
    })
    
    # ## This is averaging over all "iterations", i.e., R
    # ## Add MCSE of MSE here
    # aggregated_metrics_n <- df_HTE_coef %>%
    #   group_by(lasso_type, eta_spec, HTE_spec, prop_score_spec, stage, k) %>%
    #   summarise(MSE = mean(MSE), n_iterations = n()) %>%
    #   ungroup()
    aggregated_metrics_n <- df_HTE_coef %>%
      group_by(lasso_type, eta_spec, HTE_spec, prop_score_spec, stage, k) %>%
      summarise(
        MSE_mean = mean(MSE),
        n_iterations = n(),
        MCSE_MSE = sqrt(var(MSE) / n())
        # MCSE_MSE = {
        #   mse_values <- MSE
        #   calc_mcse_mse(mse_values, theta = HTE_true) # Adjust theta if needed
        # }
      ) %>%
      ungroup() %>%
      rename(MSE = "MSE_mean")
    
    # Add n as a variable
    aggregated_metrics_n$n <- n
    
    # Combine with overall aggregated metrics
    if (is.null(aggregated_metrics)) {
      aggregated_metrics <- aggregated_metrics_n
    } else {
      aggregated_metrics <- rbind(aggregated_metrics, aggregated_metrics_n)
    }
  }
  
  table_csv_file <- get_table_csv_file(config = config,
                                       eta_type = eta_type,
                                       HTE_type = HTE_type,
                                       output_csv_dir = output_csv_dir)
  
  write.csv(
    aggregated_metrics, 
    table_csv_file, 
    row.names = FALSE)
  
  return(aggregated_metrics)
}


load_and_process_table_data_single_spec <- function(json_file,
                                        eta_type,
                                        HTE_type,
                                        output_csv_dir) {
  config <- load_experiment_config(json_file)
  table_csv_file <- get_table_csv_file(config,eta_type, HTE_type, output_csv_dir)
  data <- read.csv(table_csv_file)
  data <- data %>%
    mutate(Specification = as.character(Specification),
           Model = sapply(strsplit(Specification, "_"), `[`, 1),
           eta_model = sapply(strsplit(Specification, "_"), `[`, 2),
           HTE_model = sapply(strsplit(Specification, "_"), `[`, 3))
  return(data)
}

load_and_process_table_data_multi_spec <- function(json_file,
                                               eta_type,
                                               HTE_type,
                                               MSE_csv_dir) {
  config <- load_experiment_config(json_file)
  table_csv_file <- get_table_csv_file(config,eta_type, HTE_type, output_csv_dir = MSE_csv_dir)
  data <- read.csv(table_csv_file)
  return(data)
}

### TODO: 
### 0. Make the file name better
### 1. Remove "regressor-spec" from eta_model label
### 2. Remove "HTE-spec" from eta_model label
### input: the table
make_HTE_by_eta_plots <- function(json_file_lasso,
                                 json_file_TV_CSL, 
                                 eta_type,
                                 HTE_type,
                                 output_csv_dir = "scripts/TV-CSL/tables-and-plots") {
  data_lasso <- 
    load_and_process_table_data_single_spec(json_file_lasso,eta_type,HTE_type, output_csv_dir)
  
  data_TV_CSL <- 
    load_and_process_table_data_single_spec(json_file_TV_CSL,eta_type,HTE_type, output_csv_dir)
  
  data_TV_CSL <- data_TV_CSL %>% filter(!grepl("m-regression", Specification))
  
  combined_data <- bind_rows(data_lasso, data_TV_CSL, .id = "source")
  
  
  print((combined_data))
  
  hte_labels <- c(
    "HTE-spec-complex" = "HTE: Complex",
    "HTE-spec-linear" = "HTE: Linear"
  )
  
  eta_labels <- c(
    "regressor-spec-complex" = "Baseline: Linear",
    "regressor-spec-linear" = "Baseline: Complex"
  )
  
  p <- ggplot(combined_data <- combined_data %>%
                mutate(Method = ifelse(Method == "TV_CSL", "TV-CSL", 
                                       ifelse(Method == "Lasso", "S-Lasso", Method))), 
              aes(x = as.factor(n), y = MSE, color = Method)) +
    geom_point() +
    # geom_line() +
    geom_line(aes(group = interaction(Method, eta_model, HTE_model))) +
    geom_errorbar(aes(ymin = MSE - 1.96 * MCSE_MSE, 
                      ymax = MSE + 1.96 * MCSE_MSE), 
                  width = 0.2) +  # Width of the error bar caps
    facet_grid(
      HTE_model ~ eta_model,
      labeller = labeller(
        HTE_model = hte_labels,
        eta_model = eta_labels
      )
    ) + 
    labs(x = "n", y = "MSE", title = "MSE by Method and Model Specifications") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.1),
          strip.text = element_text(face = "bold"))
  
  output_plot_path <- 
    file.path(output_csv_dir, 
              paste0("both-methods_", 
                       paste0(eta_type, "_", HTE_type), 
                     "_MSE_plots.png") )
  ggsave(output_plot_path, plot = p, width = 8, height = 6)
  
  return(output_plot_path)
}

