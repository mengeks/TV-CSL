read_single_simulation_data <- 
  function(n, is_time_varying=NULL, 
           i = NULL, 
           eta_type = "linear-interaction", 
           baseline_type = "linear",
           HTE_type = "constant",
           folder_name = "data") {
  
  path_for_sim_data <- here::here(
    folder_name, 
    paste0(eta_type, "_", HTE_type)
  )
  
  
  file_name <- paste0("sim_data_", i, 
                      "_n_", n,
                      ".rds")
  
  file_path <- file.path(path_for_sim_data, file_name)
  
  if (!file.exists(file_path)) {
    stop(paste("File does not exist:", file_path))
  }
  
  data <- readRDS(file_path)
  
  cat("Loaded dataset", i, "from", file_path, "\n")
  
  return(data)
}






RESULTS_DIR <- "scripts/TV-CSL/results/"

generate_output_folder <- function(results_dir, method_setting, eta_type, HTE_type, n) {
  dgp_setting <- paste0("eta-", eta_type, "_HTE-", HTE_type)
  output_folder <- paste0(results_dir, method_setting, dgp_setting, "_n-", n, "/")
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  return(output_folder)
}


generate_output_folder_heart_transplant <- 
  function(results_dir = "scripts/heart-transplant-analysis/results/") {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  return(results_dir)
}


#' Generate file paths for experiment results and save data to CSV
#'
#' @param is_running_cox Logical, whether cox method is running
#' @param is_running_lasso Logical, whether lasso method is running
#' @param is_running_TV_CSL Logical, whether TV-CSL method is running
#' @param eta_type Character, the eta type
#' @param HTE_type Character, the HTE type
#' @param n Integer, the sample size
#' @param i Integer, the iteration number
#' @param seed_value Integer, the seed value for reproducibility
#' 
#' @return The full path of the saved CSV file
generate_output_path <- function(results_dir = "scripts/TV-CSL/results/",
                                 is_running_cox = FALSE,
                                 is_running_lasso = FALSE,
                                 is_running_s_cox = FALSE,
                                 is_running_TV_CSL = FALSE,
                                 eta_type,
                                 HTE_type,
                                 n,
                                 i,
                                 seed_value) {

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
  
  result_csv_file <- paste0(
    output_folder, 
    "result-iteration_", i, "-seed_", seed_value, ".csv"
  )
  
  
  return(result_csv_file)
}


library(readr)
save_res_to_csv<-
  function(curr_res,
           FNAME){
    curr_res_df <- as.data.frame(curr_res, stringsAsFactors = FALSE)
    
    if (file.exists(FNAME)) {
      write_csv(curr_res_df, FNAME, append=TRUE)
    } else {
      write_csv(curr_res_df, FNAME)
    }
    print(paste("Result for config_name:", curr_res_df$config_name, "saved to", FNAME))
  } 


save_lasso_beta <- function(lasso_ret, 
                            output_folder, 
                            i, 
                            lasso_type,
                            eta_spec, 
                            HTE_spec, 
                            prop_score_spec,
                            stage = "final",
                            k = 0,
                            method = "NA") {
  # fname_HTE <- paste0(output_folder, "/", "eta-spec-", eta_spec,"_HTE-spec-",HTE_spec, "_beta-HTE.csv")
  # fname_eta_0 <- paste0(output_folder,"/", "eta-spec-", eta_spec,"_HTE-spec-",HTE_spec, "_beta-HTE.csv")
  fname_HTE <- paste0(output_folder, "/", "HTE-spec-",HTE_spec, "_beta-HTE.csv")
  curr_res_beta_HTE <- c(iteration = i, 
                         method = method,
                         lasso_type = lasso_type, 
                         eta_spec = eta_spec, 
                         HTE_spec = HTE_spec, 
                         prop_score_spec = prop_score_spec, 
                         stage = stage, k = k,  
                         lasso_ret$beta_HTE)
  curr_res_beta_HTE_df <- as.data.frame(t(curr_res_beta_HTE), stringsAsFactors = FALSE)
  save_res_to_csv(curr_res_beta_HTE_df, fname_HTE)
  
  if (stage != "final"){
    dir.create(file.path(output_folder, "eta-0"), recursive = TRUE)
    fname_eta_0 <- paste0(output_folder,"/eta-0/", "eta-spec-", eta_spec, "_beta-0.csv")
    if (lasso_type == "m-regression"){
      curr_res_beta_eta_0 <- c(iteration = i, 
                               method = method,
                               lasso_type = lasso_type, 
                               eta_spec = eta_spec, 
                               HTE_spec = HTE_spec, 
                               prop_score_spec = prop_score_spec, 
                               stage = stage, 
                               k = k, 
                               lasso_ret$m_beta)
    }else{
      curr_res_beta_eta_0 <- c(iteration = i, 
                               method = method,
                               lasso_type = lasso_type, 
                               eta_spec = eta_spec, 
                               HTE_spec = HTE_spec, 
                               prop_score_spec = prop_score_spec, 
                               stage = stage, 
                               k = k, 
                               lasso_ret$beta_eta_0)
    }
    
    curr_res_beta_eta_0_df <- as.data.frame(t(curr_res_beta_eta_0), stringsAsFactors = FALSE)
    save_res_to_csv(curr_res_beta_eta_0_df, fname_eta_0)
  }
  
}

save_lasso_MSE <- function(lasso_ret, 
                           HTE_true, 
                           output_folder, 
                           i, 
                           lasso_type,
                           eta_spec, 
                           HTE_spec, 
                           prop_score_spec,
                           stage = "final",
                           k = 0,
                           method = "NA") {
  HTE_est <- lasso_ret$HTE_est
  MSE <- mean((HTE_true - HTE_est)^2)
  
  fname_MSE <- paste0(output_folder, "/", "MSE.csv")
  
  curr_res_MSE <- data.frame(
    iteration = i,
    method = method,
    lasso_type = lasso_type,
    eta_spec = eta_spec,
    HTE_spec = HTE_spec,
    prop_score_spec = prop_score_spec, 
    stage = stage,
    k = k,
    MSE = MSE,
    stringsAsFactors = FALSE
  )
  
  save_res_to_csv(curr_res_MSE, fname_MSE)
}
