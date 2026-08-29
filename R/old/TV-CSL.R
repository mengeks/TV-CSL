library(survival)
library(glmnet)
library(tidyverse)
source(here::here("R/old/cox-loglik.R"))


#' Estimate Cox model
#'
#' @param single_data The input dataset.
#' @param methods_cox A list containing model specification and run_time_varying flag.
#' @param light_censoring Optional parameter for light censoring.
#' @return The estimated coefficient for W in the Cox model.
#' @export
cox_model_estimation <- function(
    single_data, 
    methods_cox, 
    light_censoring = NULL,
    HTE_type = "constant",
    eta_type = "linear-interaction"
  ) {
  model_spec <- methods_cox$model_spec
  run_time_varying <- methods_cox$run_time_varying
  
  data_to_use <- 
    preprocess_data(single_data, 
                    run_time_varying)
  
  formula <- 
    create_cox_formula(model_spec = model_spec, 
                       run_time_varying = run_time_varying,
                       HTE_type = HTE_type,
                       eta_type = eta_type)
  
  cox_model <- 
    coxph(formula, 
          data = data_to_use, 
          ties = "breslow")
  
  # tau_est_cox <- cox_model$coefficients["W"]
  beta_est_cox <- extract_W_coefficients(
    fit=cox_model, 
    HTE_type = "constant")
  
  # return(tau_est_cox)
  return(beta_est_cox)
}



extract_W_coefficients <- function(fit, HTE_type = "constant") {
  all_coefs <- coef(fit)
  
  if (HTE_type == "constant") {
    W_part <- "W"
  } else if (HTE_type == "linear") {
    W_part <- "W:X.1|W:X.10"  # Modify as per how W interacts with X covariates in your model
  } else {
    stop("Unsupported HTE_type")
  }
  
  W_coefs <- all_coefs[grep(W_part, names(all_coefs))]
  
  return(W_coefs)
}



#' #' Estimate Lasso model
#' #'
#' #' @param single_data The input dataset.
#' #' @param methods_cox A list containing model specification and run_time_varying flag.
#' #' @param light_censoring Optional parameter for light censoring.
#' #' @return The estimated coefficient for W in the Cox model.
#' #' @export
#' lasso_model_estimation <- function(
#'     single_data, 
#'     i,
#'     methods_lasso,
#'     HTE_type = "constant",
#'     eta_type = "linear-interaction"
#' ) {
#'   regressor_spec <- methods_lasso$regressor_spec
#'   lasso_type <- methods_lasso$lasso_type
#'   HTE_spec <- methods_lasso$HTE_spec
#'   
#'   train_data <- 
#'     preprocess_data(single_data, 
#'                     run_time_varying = T)
#'   
#'   n <- nrow(single_data)
#'   test_data <- 
#'     read_single_simulation_data(
#'       n = n, 
#'       i = i + 100, 
#'       eta_type = eta_type,
#'       HTE_type = HTE_type)$data
#'   
#'   T_lasso_ret <- 
#'     T_lasso(train_data = train_data, 
#'             test_data = test_data,
#'             lasso_type)
#'   
#'   
#'   # lasso_model <- 
#'   #   glmnet(transformed_X, Surv(tstart, tstop, Delta), 
#'   #          family = "cox",
#'   #          data = data_to_use)
#'   
#'   
#'   # 
#'   # cox_model <- 
#'   #   coxph(formula, 
#'   #         data = data_to_use, 
#'   #         ties = "breslow")
#'   # beta_est_cox <- extract_W_coefficients(
#'   #   fit=cox_model, 
#'   #   HTE_type = "constant")
#'   
#'   return(beta_est_cox)
#' }


#' Cox-based HTE Estimation with Flexible Regressor Specification
#'
#' This function fits a Cox proportional hazards regression model for estimating Conditional Average Treatment Effects (HTE) 
#' using flexible regressor specifications, including natural splines and interaction terms.
#'
#' @param train_data A data frame containing the training data. Must include columns for treatment (`W`), 
#'   time-varying covariates (`X.`), and survival times (`tstart`, `tstop`, `Delta`).
#' @param test_data A data frame containing the test data for evaluating the model. Must include columns 
#'   for HTE and survival-related variables.
#' @param regressor_spec A character string specifying the type of regressor transformation:
#'   - "linear": Uses linear terms only.
#'   - "complex": Includes natural splines, square terms, and pairwise interactions for continuous variables.
#' @param HTE_spec A character string specifying the specification for HTE estimation:
#'   - "correctly-specified": Uses true variables (`X.1` and `X.10`) for HTE estimation.
#'   - "linear": Uses all variables starting with `X.` for HTE estimation with linear terms.
#'   - "flexible": Uses all transformed variables (splines, interactions, etc.) for flexible HTE estimation.
#' @param verbose A numeric flag controlling print statements for debugging (0 = no output, 1 = progress output, 2 = detailed output).
#'
#' @return A list containing the following components:
#'   - `m`: The fitted Cox model.
#'   - `m_beta`: The estimated coefficients from the model.
#'   - `beta_HTE`: Coefficients for HTE estimation.
#'   - `beta_eta_0`: Coefficients for baseline hazard estimation.
#'   - `y_0_pred`: Predicted baseline outcomes for the control group.
#'   - `y_1_pred`: Predicted outcomes for the treatment group.
#'   - `HTE_est`: Estimated HTE values.
#'   - `HTE_true`: True HTE values (from test data).
#'   - `MSE`: Mean squared error of HTE estimation.
#'
#' @examples
#' # Fit the model with "linear" regressor and "correctly-specified" HTE
#' result <- S_cox(train_data = df_train, test_data = df_test, regressor_spec = "linear", HTE_spec = "correctly-specified")
#'
#' @export
S_cox <- function(train_data, 
                  test_data, 
                  regressor_spec, 
                  HTE_spec, 
                  verbose = 0) {
  
  transformed_X <- transform_X(
    single_data = train_data, 
    transform_spec = regressor_spec
  )
  
  if (HTE_spec == "correctly-specified") {
    regressor_HTE <- train_data$W * cbind(train_data$X.1, train_data$X.10)
  } else if (HTE_spec == "linear") {
    regressor_HTE <- cbind(train_data$W, train_data$W * train_data %>% select(starts_with("X.")))
  } else if (HTE_spec == "flexible") {
    regressor_HTE <- cbind(train_data$W, train_data$W * transformed_X)
  }
  
  regressor <- cbind(transformed_X, regressor_HTE)
  regressor <- as.matrix(regressor)
  
  # Fit Cox proportional hazards model instead of Lasso
  m <- coxph(Surv(train_data$tstart, train_data$tstop, train_data$Delta) ~ regressor, data = train_data)
  
  m_beta <- coef(m)  # Get the coefficients of the Cox model
  
  n_transformed_X <- ncol(transformed_X)
  beta_HTE <- m_beta[(n_transformed_X + 1):length(m_beta)] 
  beta_eta_0 <- m_beta[1:(n_transformed_X)] 
  
  if (verbose >= 1){
    print("Finished fitting the Cox model.")
  }
  
  X_HTE_test <- transform_X(
    single_data = test_data,
    transform_spec = regressor_spec
  )
  
  if (verbose >= 1){
    print("Start prediction.")
  }
  
  if (HTE_spec == "correctly-specified") {
    test_regressor_HTE <- cbind(test_data$X.1, test_data$X.10)
  } else if (HTE_spec == "linear") {
    test_regressor_HTE <- cbind(1, as.matrix(test_data %>% select(starts_with("X."))))
  } else if (HTE_spec == "flexible") {
    test_regressor_HTE <- cbind(1, as.matrix(X_HTE_test))
  }
  
  HTE_est <- as.vector(test_regressor_HTE %*% beta_HTE)
  HTE_true <- test_data$HTE
  
  # Added y_0_pred
  y_0_pred <- as.vector(X_HTE_test %*% beta_eta_0)
  y_1_pred <- y_0_pred + HTE_est
  
  if (verbose >= 1){
    print("Finished prediction.")
  }
  
  MSE <- mean((HTE_est - HTE_true)^2)

  if (verbose == 2){
    print("MSE:")
    print(MSE)
  }

  # Inference: SE for beta_HTE from vcov of the joint coxph model
  vcov_full <- vcov(m)
  idx_HTE   <- (n_transformed_X + 1):nrow(vcov_full)
  vcov_HTE  <- vcov_full[idx_HTE, idx_HTE, drop = FALSE]
  se_HTE    <- sqrt(diag(vcov_HTE))

  ret <- list(
    m = m,
    m_beta = m_beta,
    beta_HTE = beta_HTE,
    beta_eta_0 = beta_eta_0,
    y_0_pred = y_0_pred,
    y_1_pred = y_1_pred,
    HTE_est = HTE_est,
    HTE_true = HTE_true,
    MSE = MSE,
    se_HTE   = se_HTE,
    vcov_HTE = vcov_HTE
  )

  if (verbose >= 1){
    print("m_beta:")
    print(m_beta)
    print("beta_HTE:")
    print(beta_HTE)
  }

  class(ret) <- "scox"
  ret
}



T_lasso <- function(train_data, 
                    test_data, 
                    regressor_spec = "complex") {

  
  # Get row indexes for control and treatment groups
  index_co <- which(train_data$W == 0)
  index_tx <- which(train_data$W == 1)
  
  
  transformed_X <- transform_X(
    single_data = train_data, 
    transform_spec = regressor_spec
  )
  # Separate transformed_X for control and treatment groups
  transformed_X_co <- transformed_X[index_co, ]
  transformed_X_tx <- transformed_X[index_tx, ]
  
  # Separate train_data for control and treatment groups
  data_co <- train_data[index_co, ]
  data_tx <- train_data[index_tx, ]
  
  # Fit Cox models for control and treatment groups
  eta_0 <- cv.glmnet(transformed_X_co, Surv(data_co$tstart, data_co$tstop, data_co$Delta), family = "cox")
  eta_1 <- cv.glmnet(transformed_X_tx, Surv(data_tx$tstart, data_tx$tstop, data_tx$Delta), family = "cox")
  
  print("eta_0 coefficient from T-lasso: ")
  print(coef(eta_0, s = "lambda.min"))
  print("eta_1 coefficient from T-lasso: ")
  print(coef(eta_1, s = "lambda.min"))
  
  
  
  # Transform the test data
  X_HTE_test <- transform_X(
    single_data = test_data, 
    transform_spec = regressor_spec
  )
  
  # Predict on test data using both models
  y_1_pred <- predict(eta_1, newx = X_HTE_test, s = "lambda.min")
  y_0_pred <- predict(eta_0, newx = X_HTE_test, s = "lambda.min")
  
  # Compute HTE estimate
  HTE_est <- y_1_pred - y_0_pred
  HTE_true <- test_data$HTE
  
  # print("head(HTE_est)")
  # print(head(HTE_est))
  # print("head(HTE_true)")
  # print(head(HTE_true))
  
  # Calculate Mean Squared Error
  MSE <- mean((HTE_est - HTE_true)^2)
  
  # Create and return the result list
  ret <- list(
    eta_0 = eta_0,
    eta_1 = eta_1,
    y_1_pred = y_1_pred,
    y_0_pred = y_0_pred,
    HTE_est = HTE_est,
    HTE_true = HTE_true,
    MSE = MSE
  )
  
  class(ret) <- "tlasso"
  return(ret)
}

#' Lasso-based HTE Estimation with Flexible Regressor Specification
#'
#' This function fits a Lasso Cox regression model for estimating Conditional Average Treatment Effects (HTE) 
#' using flexible regressor specifications, including natural splines and interaction terms.
#'
#' @param train_data A data frame containing the training data. Must include columns for treatment (`W`), 
#'   time-varying covariates (`X.`), and survival times (`tstart`, `tstop`, `Delta`).
#' @param test_data A data frame containing the test data for evaluating the model. Must include columns 
#'   for HTE and survival-related variables.
#' @param regressor_spec A character string specifying the type of regressor transformation:
#'   - "linear": Uses linear terms only.
#'   - "complex": Includes natural splines, square terms, and pairwise interactions for continuous variables.
#' @param HTE_spec A character string specifying the specification for HTE estimation:
#'   - "correctly-specified": Uses true variables (`X.1` and `X.10`) for HTE estimation.
#'   - "linear": Uses all variables starting with `X.` for HTE estimation with linear terms.
#'   - "flexible": Uses all transformed variables (splines, interactions, etc.) for flexible HTE estimation.
#' (TODO) change the S_lasso code to specify the names of the covariates. If covariates names is NULL, then we assumerain_data has columns named X.1, X.2, ..., X.n as covariates
#' @return A list containing the following components:
#'   - `m`: The fitted Lasso Cox model.
#'   - `m_beta`: The estimated coefficients from the model.
#'   - `beta_HTE`: Coefficients for HTE estimation.
#'   - `beta_eta_0`: Coefficients for baseline hazard estimation.
#'   - `y_0_pred`: Predicted baseline outcomes for the control group.
#'   - `y_1_pred`: Predicted outcomes for the treatment group.
#'   - `HTE_est`: Estimated HTE values.
#'   - `HTE_true`: True HTE values (from test data).
#'   - `MSE`: Mean squared error of HTE estimation.
#'
#' @examples
#' # Fit the model with "linear" regressor and "correctly-specified" HTE
#' result <- S_lasso(train_data = df_train, test_data = df_test, regressor_spec = "linear", HTE_spec = "correctly-specified")
#'
#' @export
S_lasso <- function(train_data, 
                    test_data, 
                    regressor_spec, 
                    HTE_spec, 
                    verbose = 0) {
  
  X_eta <- transform_X(
    single_data = train_data, 
    transform_spec = regressor_spec
  )
  
  # complex_X <- transform_X(
  #   single_data = train_data, 
  #   transform_spec = "complex"
  # )
  
  # if (HTE_spec == "correctly-specified") {
  #   regressor_HTE <- train_data$W * cbind(train_data$X.1, train_data$X.10)
  # } else if (HTE_spec == "linear") {
  # if (HTE_spec == "linear") {
  #   regressor_HTE <- cbind(train_data$W, train_data$W * train_data %>% select(starts_with("X.")))
  # } else if (HTE_spec == "flexible") {
  #   regressor_HTE <- cbind(train_data$W, train_data$W * transformed_X)
  # } else if (HTE_spec == "complex") {
  #   regressor_HTE <- cbind(train_data$W, train_data$W * complex_X)
  # }
  
  X_HTE <- transform_X(
    single_data = train_data, 
    transform_spec = HTE_spec # "linear" or "complex"
  )
  regressor_HTE <- cbind(train_data$W, train_data$W * X_HTE)
  
  regressor <- cbind(X_eta, regressor_HTE)
  regressor <- as.matrix(regressor)
  
  if (HTE_spec == "linear" & regressor_spec == "linear"){
    # m <- coxph(
    #   Surv(tstart, tstop, Delta) ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8 + X.9 + X.10 + 
    #     W * (1 + X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8 + X.9 + X.10),
    #   data = train_data
    # )
    # m <- coxph(
    #   Surv(tstart, tstop, Delta) ~ X.1 + X.2 + X.3 + 
    #     W * (1 + X.1 + X.2 + X.3),
    #   data = train_data
    # )
    num_covariates <- ncol(train_data %>% select(starts_with("X.")))
    
    covariate_terms <- paste("X.", 1:num_covariates, sep = "", collapse = " + ")
    interaction_terms <- paste("X.", 1:num_covariates, sep = "", collapse = " + ")
    interaction_formula <- paste("W * (1 +", interaction_terms, ")")
    
    full_formula <- paste("Surv(tstart, tstop, Delta) ~", covariate_terms, "+", interaction_formula)
    
    m <- coxph(as.formula(full_formula), data = train_data)
    
    m_beta <- coef(m)
  }else{
    m <- cv.glmnet(regressor, Surv(train_data$tstart, train_data$tstop, train_data$Delta), 
                   # lambda = 10^seq(-4, 1, length = 100),
                   family = "cox"
    )
    
    m_beta <- coef(m, s = "lambda.min")
  }
  

  n_X_eta <- ncol(X_eta)
  beta_HTE <- m_beta[(n_X_eta + 1):length(m_beta)] 
  beta_eta_0 <- m_beta[1:(n_X_eta)] 
  
  if (verbose >= 1){
    print( "Finished fitting the lasso model. ")
  }
  
  
  # X_HTE_test <- transform_X(
  #   single_data = test_data,
  #   transform_spec = regressor_spec)
  # test_complex_X <- transform_X(
  #   single_data = test_data,
  #   transform_spec = "complex")
  
  X_HTE_test <- transform_X(
    single_data = test_data,
    transform_spec = HTE_spec)
  
  # test_regressor_HTE <- cbind(1, train_data$W * X_HTE_test)
  test_regressor_HTE <- cbind(1, X_HTE_test)
  
  if (verbose >= 1){
    print( "Start prediction. ")
  }
  
  # if (HTE_spec == "correctly-specified") {
  #   test_regressor_HTE <- cbind(test_data$X.1, test_data$X.10)
  # } else if (HTE_spec == "linear") {
  #   test_regressor_HTE <- cbind(1, as.matrix(test_data %>% select(starts_with("X.")) ) )
  # } else if (HTE_spec == "flexible") {
  #   test_regressor_HTE <- cbind(1, as.matrix(X_HTE_test))
  # } else if (HTE_spec == "complex") {
  #   test_regressor_HTE <- cbind(1, as.matrix(test_complex_X))
  # }
  
  HTE_est <- as.vector(test_regressor_HTE %*% beta_HTE)
  HTE_true <- test_data$HTE
  
  # Added y_0_pred
  X_eta_test <- transform_X(
    single_data = test_data,
    transform_spec = regressor_spec)
  
  y_0_pred <- as.vector(X_eta_test %*% beta_eta_0)
  y_1_pred <- y_0_pred + HTE_est
  
  if (verbose >= 1){
    print( "Finished prediction. ")
  }
  
  
  MSE <- mean((HTE_est - HTE_true)^2)
  if (verbose == 2){
    print( "MSE: ")
    print( MSE )
  }

  # Inference: SE for beta_HTE.
  # Only available for the unpenalized coxph path (linear+linear); NA otherwise.
  if (HTE_spec == "linear" && regressor_spec == "linear") {
    vcov_full <- vcov(m)
    idx_HTE   <- (n_X_eta + 1):nrow(vcov_full)
    vcov_HTE  <- vcov_full[idx_HTE, idx_HTE, drop = FALSE]
    se_HTE    <- sqrt(diag(vcov_HTE))
  } else {
    se_HTE   <- rep(NA_real_, length(beta_HTE))
    vcov_HTE <- NULL
  }

  ret <- list(
    m = m,
    m_beta = m_beta,
    beta_HTE = beta_HTE,
    beta_eta_0 = beta_eta_0,
    y_0_pred = y_0_pred,
    y_1_pred = y_1_pred,
    HTE_est = HTE_est,
    HTE_true = HTE_true,
    MSE = MSE,
    se_HTE   = se_HTE,
    vcov_HTE = vcov_HTE
  )

  if (verbose >= 1){
    print("m_beta for S-lasso: ")
    print(m_beta)
    print("beta_HTE for S-lasso: ")
    print(beta_HTE)
  }
  
  class(ret) <- "slasso"
  ret
}

#' Lasso-based HTE Estimation with Flexible Regressor Specification
#'
#' This function fits a Lasso Cox regression model for estimating Conditional Average Treatment Effects (HTE) 
#' using flexible regressor specifications, including natural splines and interaction terms.
#'
#' @param train_data A data frame containing the training data. Must include columns for treatment (`W`), 
#'   time-varying covariates (`X.`), and survival times (`tstart`, `tstop`, `Delta`).
#' @param regressor_spec A character string specifying the type of regressor transformation:
#'   - "linear": Uses linear terms only.
#'   - "complex": Includes natural splines, square terms, and pairwise interactions for continuous variables.
#'
#'
#' @examples
#' # Fit the model with "linear" regressor and "correctly-specified" HTE
#' result <- S_lasso(train_data = df_train, test_data = df_test, regressor_spec = "linear", HTE_spec = "correctly-specified")
#'
#' @export
m_regression <- function(train_data, 
                         test_data,
                    regressor_spec,
                    lambda = NULL,
                    verbose = 0) {
  
  transformed_X <- transform_X(
    single_data = train_data, 
    transform_spec = regressor_spec
  )
  
  regressor <- as.matrix(transformed_X)
  
  if (regressor_spec == "linear"){
    num_covariates <- ncol(train_data %>% select(starts_with("X.")))
    
    covariate_terms <- paste("X.", 1:num_covariates, sep = "", collapse = " + ")
    
    full_formula <- paste("Surv(tstart, tstop, Delta) ~", covariate_terms)
    
    m <- coxph(as.formula(full_formula), data = train_data)
    # m <- coxph(
    #   Surv(tstart, tstop, Delta) ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8 + X.9 + X.10,
    #   data = train_data
    # )
    
    m_beta <- coef(m)
  }else if (regressor_spec == "complex"){
    m <- cv.glmnet(regressor, 
                   Surv(train_data$tstart, train_data$tstop, train_data$Delta), 
                   family = "cox",
                   lambda = lambda
    )
    
    m_beta <- coef(m, s = "lambda.min")
  }
  
  
  X_HTE_test <- transform_X(
    single_data = test_data,
    transform_spec = regressor_spec)
  
  
  y_pred <- as.vector(X_HTE_test %*% m_beta)
  y_0_pred <- y_1_pred <-  y_pred 
  # in m_regression, there is no need to have y_0_pred, y_1_pred
  # We have them to conform with the existing TV_CSL coding logic
  #   to construct nu easily
  
  ret <- list(
    m = m,
    m_beta = m_beta,
    y_pred = y_pred,
    y_0_pred = y_0_pred,
    y_1_pred = y_1_pred
  )
  
  print("m_beta for m_regression: ")
  print(m_beta)
  
  class(ret) <- "m-regression"
  ret
}




#' Transform X Variables for Lasso Regression
#'
#' This function performs transformations on variables starting with "X." from the input data.
#' Depending on the `transform_spec`, the function creates either simple linear terms or a more complex set of transformations including natural splines, square terms, and interaction terms.
#'
#' @param single_data A data frame that contains the regressors, with column names starting with "X.".
#' @param transform_spec A character string specifying the type of transformation. 
#'   - "linear": Generates linear terms only.
#'   - "complex": Generates linear terms, 3 natural splines for continuous variables, square terms, and pairwise interaction terms for both continuous and binary variables.
#'
#' @details
#' - If `transform_spec` is "linear", the function will return a matrix of linear terms from all variables starting with "X.".
#' - If `transform_spec` is "complex", the function applies the following transformations:
#'   - Linear terms for all variables.
#'   - 3 natural spline terms for continuous variables.
#'   - Square terms for continuous variables.
#'   - Pairwise interaction terms between all variables (both continuous and binary).
#'
#' @return A matrix of transformed regressors, including linear terms, splines, square terms, and interactions depending on the specified transformation type.
#' 
#' @examples
#' # Using linear transformation
#' transformed_X_linear <- transform_X(single_data = df_time_var, transform_spec = "linear")
#'
#' # Using complex transformation
#' transformed_X_complex <- transform_X(single_data = df_time_var, transform_spec = "complex")
#'
#' @export
transform_X <- function(single_data, transform_spec = "linear") {
  library(splines)
  library(dplyr)
  
  X_vars <- single_data %>% select(starts_with("X."))
  # X_vars <- single_data[grep("^X", names(single_data))]
  
  transformed_X <- matrix(nrow = nrow(X_vars), ncol = 0)
  
  if (transform_spec == "complex") {
    
    # transformed_X <- as.matrix(X_vars)
    continuous_vars <- names(X_vars)[sapply(X_vars, function(x) length(unique(x)) > 2)]
    
    # Apply natural splines to continuous variables
    for (var_name in continuous_vars) {
      spline_terms <- ns(X_vars[[var_name]], df = 3)
      colnames(spline_terms) <- paste0(var_name, "_spline", 1:3)
      transformed_X <- cbind(transformed_X, as.matrix(spline_terms))
    }
    
    # Add square terms for continuous variables
    for (var_name in continuous_vars) {
      square_term <- X_vars[[var_name]]^2
      transformed_X <- cbind(transformed_X, square_term)
      colnames(transformed_X)[ncol(transformed_X)] <- paste0(var_name, "_squared")
    }
    
    # Add interaction terms between all pairs of variables (both continuous and binary)
    var_names <- names(X_vars)
    for (i in 1:(length(var_names) - 1)) {
      for (j in (i + 1):length(var_names)) {
        interaction_term <- X_vars[[var_names[i]]] * X_vars[[var_names[j]]]
        transformed_X <- cbind(transformed_X, interaction_term)
        colnames(transformed_X)[ncol(transformed_X)] <- paste0(var_names[i], "_x_", var_names[j])
      }
    }
    
  } else if (transform_spec == "linear") {
    
    transformed_X <- as.matrix(X_vars)
    
  }
  
  return(transformed_X)
}



#' Run Cox Model Estimation for Different Specifications
#'
#' This function runs Cox model estimation for the provided dataset 
#' using different model specifications from the \code{methods_cox}.
#'
#' @param single_data A data frame containing the survival data.
#' @param methods_cox A list containing the model specification and time-varying flag for Cox model.
#' @return A list of results containing the tau estimates and time taken for each model specification.
#' @export
run_cox_estimation <- function(
    single_data, methods_cox, HTE_type, eta_type) {
  results <- list()
  
  if (methods_cox$enabled) {
    for (spec in methods_cox$model_specifications) {
      for (run_time_varying in methods_cox$run_time_varying) {
        
        config_name <- paste(spec, run_time_varying, sep = "_")
        
        start_time <- Sys.time()
        
        # Create a temporary methods_cox object that contains the current specification and run_time_varying
        current_methods_cox <- methods_cox
        current_methods_cox$model_spec <- spec
        current_methods_cox$run_time_varying <- run_time_varying
        
        # Run the Cox model estimation using the current configuration
        # tau_est_cox <- 
        beta_est_cox <- 
          cox_model_estimation(
            single_data = single_data, 
            methods_cox = current_methods_cox,
            HTE_type = HTE_type,
            eta_type = eta_type
          )
        
        end_time <- Sys.time()
        
        time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
        
        results[[config_name]] <- list(
          # tau_estimate = tau_est_cox,
          beta_estimate = beta_est_cox,
          time_taken = time_taken
        )
      }
    }
  }
  
  return(results)
}


#'
#' This function runs lasso model estimation for the provided dataset 
#' using different model specifications from the \code{methods_lasso}.
#'
#' @param single_data A data frame containing the survival data.
#' @param methods_cox A list containing the model specification and time-varying flag for Cox model.
#' @return A list of results containing the tau estimates and time taken for each model specification.
#' @export
run_lasso_estimation <- function(
    single_data, i, methods_lasso, HTE_type, eta_type) {
  results <- list()
  n <- nrow(single_data)
  test_data <- 
    read_single_simulation_data(
      n = n, 
      i = i + 100, 
      eta_type = eta_type,
      HTE_type = HTE_type)$data
  train_data <- 
    preprocess_data(single_data, 
                    run_time_varying = T)
  
    for (regressor_spec in methods_lasso$regressor_specs) {
      for (lasso_type in methods_lasso$lasso_types) {
        
        if (lasso_type == "T-lasso"){
          # config_name <- paste(lasso_type, regressor_spec, sep = "_")
          config_name <- paste("lasso-type-",lasso_type, "_regressor-spec_", regressor_spec)
          start_time <- Sys.time()
          lasso_ret <-
            T_lasso(train_data = train_data,
                    test_data = test_data,
                    regressor_spec = regressor_spec)
          
          end_time <- Sys.time()
          
          
          time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
          
          results[[config_name]] <- list(
            HTE_est = lasso_ret$HTE_est,
            HTE_true = lasso_ret$HTE_true,
            MSE = lasso_ret$MSE,
            time_taken = time_taken
          )
          print(paste0("config_name: ", config_name, ". MSE: ", lasso_ret$MSE, ". time_taken: ", time_taken))
          
        }else if (lasso_type == "S-lasso"){
          for (HTE_spec in methods_lasso$HTE_specs){
            
            start_time <- Sys.time()
            # config_name <- paste(lasso_type, HTE_spec, regressor_spec, sep = "_")
            config_name <- paste0("lasso-type-",lasso_type, "_regressor-spec-", regressor_spec, "_HTE-spec-", HTE_spec)
            lasso_ret <-
              S_lasso(train_data = train_data,
                      test_data = test_data,
                      regressor_spec = regressor_spec,
                      HTE_spec = HTE_spec)
            
            end_time <- Sys.time()
            
            time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
            print(paste0("config_name: ", config_name, "; MSE: ", lasso_ret$MSE, "; time_taken: ", time_taken))
            
            ## save lasso_ret$beta_HTE and lasso_ret$beta_eta_0
            output_folder <- generate_output_folder(
              results_dir = RESULTS_DIR,
              method_setting = "lasso_", 
              eta_type = eta_type, 
              HTE_type = HTE_type, 
              n = n
            )
            
            save_lasso_beta(lasso_ret = lasso_ret,
                            output_folder = output_folder,
                            i = i,
                            lasso_type = lasso_type,
                            eta_spec = regressor_spec,
                            HTE_spec = HTE_spec,
                            prop_score_spec = "NA")
            
            
            results[[config_name]] <- list(
              HTE_est  = lasso_ret$HTE_est,
              HTE_true = lasso_ret$HTE_true,
              MSE      = lasso_ret$MSE,
              time_taken = time_taken,
              beta_HTE = lasso_ret$beta_HTE,
              se_HTE   = lasso_ret$se_HTE
            )
          }
        }else{
          stop("lasso_type should be T_lasso or S_lasso")
        }
        
      } # End looping lasso_types
    } # End looping regressor_specs

  return(results)
}


#' Run S-Cox Estimation
#'
#' Runs the naive single Cox model (S-Cox) for multiple regressor/HTE-spec configs.
#' S-Cox fits one time-varying Cox model with W as a predictor (no propensity debiasing).
#'
#' @param single_data Data frame with original survival data.
#' @param i Iteration index (used to locate the held-out test dataset).
#' @param methods_s_cox List with `regressor_specs` and `HTE_specs`.
#' @param HTE_type Type of HTE ("constant", "linear", ...).
#' @param eta_type Type of baseline hazard ("linear", "non-linear", ...).
#' @return A list of results per configuration containing MSE and time_taken.
#' @export
run_s_cox_estimation <- function(single_data, i, methods_s_cox, HTE_type, eta_type) {
  results <- list()
  n <- nrow(single_data)
  test_data <-
    read_single_simulation_data(
      n = n, i = i + 100,
      eta_type = eta_type, HTE_type = HTE_type
    )$data
  train_data <- preprocess_data(single_data, run_time_varying = TRUE)

  for (regressor_spec in methods_s_cox$regressor_specs) {
    for (HTE_spec in methods_s_cox$HTE_specs) {
      config_name <- paste0("regressor-spec-", regressor_spec, "_HTE-spec-", HTE_spec)
      start_time <- Sys.time()

      s_cox_ret <- S_cox(
        train_data    = train_data,
        test_data     = test_data,
        regressor_spec = regressor_spec,
        HTE_spec      = HTE_spec
      )

      end_time <- Sys.time()
      time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

      results[[config_name]] <- list(
        MSE      = s_cox_ret$MSE,
        time_taken = time_taken,
        beta_HTE = s_cox_ret$beta_HTE,
        se_HTE   = s_cox_ret$se_HTE
      )
      print(paste0("S-Cox config: ", config_name,
                   "; MSE: ", s_cox_ret$MSE, "; time: ", time_taken, "s"))
    }
  }
  return(results)
}


#' Run TV_CSL Estimation
#'
#' This function runs the Time-Varying Conditional Survival Learning (TV-CSL) model estimation with multiple configurations. It performs k-fold cross-validation for different specifications of the model, including the type of lasso, regressor, and propensity score method.
#'
#' @param train_data_original A data frame containing the original training data.
#' @param test_data A data frame containing the test data.
#' @param methods_TV_CSL A list specifying the various configurations to use in TV-CSL, including the propensity score specification, lasso types, and final model methods.
#' @param K An integer representing the number of folds for cross-validation. Default is 5.
#'
#' @details
#' The function takes the training and test data, and for each configuration of propensity score specification, lasso type, regressor specification, and final model method, it performs k-fold cross-validation. It estimates the Conditional Average Treatment Effect (HTE) for each fold and returns the results, including the MSE and computation time for each configuration.
#'
#' @return A list where each entry corresponds to a configuration and contains the HTE estimates, true HTE values, MSE, and the time taken.
#'
#' @examples
#' \dontrun{
#' # Example of how to run the function
#' results <- run_TV_CSL_estimation(train_data_original, test_data, methods_TV_CSL, K = 5)
#' }
#'
#' @export
run_TV_CSL_estimation <- function(
    train_data_original, 
    test_data,
    methods_TV_CSL,
    i, 
    K,
    HTE_type,
    eta_type, 
    temp_result_csv_file,
    verbose = 2
) {
  
  results <- list()
  n <- nrow(train_data_original)
  
  train_data <- create_pseudo_dataset(survival_data = train_data_original)
  train_data_original <- train_data_original %>% 
    mutate(U_A = pmin(A,U),
           Delta_A = A <= U)
  
  for (prop_score_spec in methods_TV_CSL$prop_score_specs) {
    for (regressor_spec in methods_TV_CSL$regressor_specs) {
      for (lasso_type in methods_TV_CSL$lasso_types) {
        for (final_model_method in methods_TV_CSL$final_model_methods) {
          for (HTE_spec in methods_TV_CSL$HTE_specs){
          
            # config_name <- paste0("lasso-type-",lasso_type, "_regressor-spec-", regressor_spec, "_HTE-spec-", HTE_spec)
            config_name <- paste0("lasso-type-",lasso_type, "_regressor-spec-", regressor_spec, "_HTE-spec-", HTE_spec, "_prop-score-spec-", prop_score_spec)
            print(config_name)
            start_time <- Sys.time()
            
            TV_CSL_ret <- TV_CSL(train_data = train_data, 
                                 test_data = test_data, 
                                 train_data_original = train_data_original, 
                                 HTE_type = HTE_type,
                                 eta_type = eta_type,
                                 K = K, 
                                 prop_score_spec = prop_score_spec, 
                                 lasso_type = lasso_type, 
                                 regressor_spec = regressor_spec, 
                                 final_model_method = final_model_method,
                                 HTE_spec = HTE_spec,
                                 i = i)
            end_time <- Sys.time()
            
            time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
          
            if (verbose == 2){
              print("Saving K = 0 results to beta and MSE csv.")
              output_folder <- generate_output_folder(
                results_dir = RESULTS_DIR,
                method_setting = "TV-CSL_", 
                eta_type = eta_type, 
                HTE_type = HTE_type, 
                n = n
              )
              save_lasso_beta(lasso_ret = TV_CSL_ret, 
                              output_folder = output_folder, 
                              i = i, 
                              k = 0,
                              lasso_type = lasso_type, 
                              eta_spec = regressor_spec, 
                              HTE_spec = HTE_spec,
                              prop_score_spec = prop_score_spec,
                              stage = "final") 
              
              save_lasso_MSE(lasso_ret = TV_CSL_ret, 
                             HTE_true = test_data$HTE , 
                             output_folder = output_folder, 
                             i = i,
                             k = 0,
                             lasso_type = lasso_type, 
                             eta_spec = regressor_spec, 
                             HTE_spec = HTE_spec,
                             prop_score_spec = prop_score_spec,
                             stage = "final")
              print("Finished K = 0 results to beta and MSE csvs.")
            }
            
            
            results[[config_name]] <- TV_CSL_ret 
            MSE <- TV_CSL_ret$MSE ## this is 5-fold result
            print(paste0("config_name: ", config_name, ". MSE: ", MSE, ". time_taken: ", time_taken))
          
          
          } # End looping over HTE_specs
        } # End looping over final_model_methods
      } # End looping over lasso_types
    } # End looping over regressor_specs
  } # End looping over prop_score_specs
  
  return(results)
}


calculate_eX <- function(alpha_estimate, X, t) {
  hazard_adapt_date <- exp(X %*% alpha_estimate)
  trt_prob <- 1 - exp(-t * hazard_adapt_date)
  return(trt_prob)
}

#' TV_CSL Nuisance Estimation
#'
#' This function estimates the nuisance components for the time-varying conditional survival model (TV-CSL). It calculates the propensity score and performs the lasso-based nuisance estimation.
#'
#' @param fold_train A data frame containing the training set.
#' @param fold_test A data frame containing the test set.
#' @param prop_score_spec A character string specifying the method for propensity score estimation (e.g., "cox-linear-censored-only").
#' @param lasso_type A character string specifying the type of lasso model ("T-lasso" or "S-lasso").
#' @param regressor_spec A character string specifying the regressor used in the lasso model.
#'
#' @details
#' The function first estimates the propensity score based on the specified propensity score model. If the model is based on the Cox proportional hazards model (e.g., "cox-linear-censored-only"), it estimates the treatment probability using a Cox regression. It then splits the test set based on time intervals, computes the treatment probabilities within each interval, and fits a lasso model (either "T-lasso" or "S-lasso") to estimate the nuisance parameters for each fold.
#'
#' @return A data frame containing the augmented test set, including the nuisance estimates.
#'
#' @examples
#' # Example usage of TV_CSL_nuisance
#' fold_train <- your_train_data
#' fold_test <- your_test_data
#' result <- TV_CSL_nuisance(fold_train, fold_test, "cox-linear-censored-only", "T-lasso", "linear")
#'
#' @importFrom dplyr mutate left_join filter
#' @importFrom survival coxph Surv
#' @importFrom glmnet cv.glmnet
#' @export
TV_CSL_nuisance <- function(fold_train, 
                            fold_test,
                            train_data_original, # for estimating the propensity score
                            prop_score_spec,
                            lasso_type,
                            regressor_spec,
                            HTE_spec, 
                            id_var = "id") {
  
  # 1. Estimate the propensity score
  if (prop_score_spec == "cox-intercept-only") {
    # Intercept-only misspecification: no covariates used.
    alpha_estimate <- 0
  } else if (prop_score_spec %in% c("cox-risk-set-adjusted", "cox-time-varying-oracle",
                                    "cox-time-varying-prop")) {
    # Propensity specs estimated after fold_test_final is built.
    alpha_estimate <- NULL
  } else if (grepl("^cox", prop_score_spec)) {
    if (prop_score_spec == "cox-linear-all-data") {
      df_prop_score <- train_data_original
    } else {
      df_prop_score <- train_data_original %>% filter(Delta == 1)
    }

    if (prop_score_spec == "cox-linear-mis-specification") {
      formula <- as.formula("Surv(U_A, Delta_A) ~ X.1")
    } else {
      num_covariates <- ncol(df_prop_score %>% select(starts_with("X.")))
      covariate_terms <- paste("X.", 1:num_covariates, sep = "", collapse = " + ")
      formula <- as.formula(paste("Surv(U_A, Delta_A) ~", covariate_terms))
    }

    treatment_model <- coxph(formula, data = df_prop_score, ties = "breslow")
    alpha_estimate <- treatment_model$coefficients

  } else {
    stop(paste0("Unknown prop_score_spec: ", prop_score_spec))
  }
  
  # 2. Obtain the propensity score at each interval

  
  split_within_intervals <- function(row, cut_points) {
    original_tstart <- row[["tstart"]]
    original_tstop <- row[["tstop"]]
    
    valid_cuts <- cut_points[cut_points > original_tstart & cut_points < original_tstop]
    final_cuts <- sort(c(original_tstart, valid_cuts, original_tstop))
    
    new_intervals <- data.frame(
      tstart = head(final_cuts, -1),
      tstop = tail(final_cuts, -1),
      id = row[["id"]],
      W = row[["W"]]
    )
    return(new_intervals)
  }
  
  # Save the true eta_0 from DGP before the lasso step overwrites it.
  # The oracle needs the population eta_0(x), not the first-stage estimate.
  if ("eta_0" %in% names(fold_test)) {
    fold_test$true_eta_0 <- fold_test$eta_0
  }

  granular_cut_points <- unique(fold_test$tstop)
  fold_test_split <- map_dfr(seq_len(nrow(fold_test)), function(i) {
    split_within_intervals(fold_test[i, ], granular_cut_points)
  })



  # 3. Obtain nu
  if (lasso_type == "T-lasso") {
    lasso_ret <- T_lasso(
      train_data = fold_train,  
      test_data = fold_test,  
      regressor_spec = regressor_spec
    )
  } else if (lasso_type == "S-lasso") {
    lasso_ret <- S_lasso(
      train_data = fold_train,
      test_data = fold_test,
      regressor_spec = regressor_spec,
      # HTE_spec = "linear" # we hard code this
      HTE_spec = HTE_spec
    )
  }else if (lasso_type == "m-regression"){
    lasso_ret <- m_regression(
      train_data = fold_train,
      test_data = fold_test,
      regressor_spec = regressor_spec
    )
  }
  fold_test$eta_1 <- lasso_ret$y_1_pred
  fold_test$eta_0 <- lasso_ret$y_0_pred
  
  
  # 3.2 Join the data
  fold_test_final <- 
    left_join(fold_test_split, 
              fold_test %>% select(-tstart, -tstop), 
              by = c("id", "W") )
  
  fold_test_final <- fold_test_final %>%
    arrange(id, tstop)
  
  fold_test_final <- fold_test_final %>%
    group_by(id) %>%
    mutate(Delta = if_else(tstop < max(tstop), 0, Delta)) %>%
    ungroup()

  if (prop_score_spec == "cox-risk-set-adjusted") {
    # Risk-set-adjusted propensity: a_0(t,x) = expit(logit(π_R(t,x)) + η₁(x) - η₀(x))
    # π_R estimated by logistic regression on the training pseudo-rows (all at-risk by construction).
    num_cov    <- ncol(fold_train %>% select(starts_with("X.")))
    cov_str    <- paste0("X.", 1:num_cov, collapse = " + ")
    pi_R_mod   <- glm(as.formula(paste("W ~ tstart +", cov_str)),
                      family = binomial, data = fold_train)
    logit_pi_R <- predict(pi_R_mod, newdata = fold_test_final, type = "link")
    # Plug-in: use estimated τ̂(x) = η̂₁(x) - η̂₀(x) from first-stage lasso.
    eta_diff   <- fold_test_final$eta_1 - fold_test_final$eta_0
    a0_hat     <- plogis(logit_pi_R + eta_diff)
    fold_test_final <- fold_test_final %>% mutate(a_t_X = a0_hat)

  } else if (prop_score_spec == "cox-time-varying-oracle") {
    # True analytical oracle for a_t(x) = P{W(t)=1 | event at t, X=x}.
    # DGP: A|X ~ Exp(r(x)) with r(x)=exp(x2+x3), h0(t)=t, tau(x)=x1+x2+x3.
    # Formula (same for linear and non-linear eta_0; only eta_0(x) value differs):
    #   logit a_t(x) = x1+2*x2+2*x3 + r*t - d*t^2 + log I(t,x)
    # where r = exp(x2+x3), d(x) = 0.5*exp(eta_0(x))*(exp(tau(x))-1),
    #       I(t,x) = integral_0^t exp(d*s^2 - r*s) ds.
    # Uses true_eta_0 (DGP value) rather than the first-stage lasso estimate.
    # Numerical stabilisation: subtract M = max(0, d*t^2 - r*t) inside the exponent.
    compute_true_a_t <- function(t, x1, x2, x3, eta0_true, tau) {
      if (is.na(t) || t <= 0) return(0.5)
      r <- exp(x2 + x3)
      d <- 0.5 * exp(eta0_true) * (exp(tau) - 1)
      M <- max(0, d * t^2 - r * t)
      stable_int <- tryCatch(
        integrate(function(s) exp(d * s^2 - r * s - M),
                  lower = 0, upper = t, rel.tol = 1e-5)$value,
        error = function(e) NA_real_
      )
      if (is.na(stable_int) || stable_int <= 0) return(0.5)
      plogis(x1 + 2*x2 + 2*x3 + r * t - d * t^2 + M + log(stable_int))
    }
    a0_hat <- mapply(compute_true_a_t,
                     t         = fold_test_final$tstop,
                     x1        = fold_test_final$X.1,
                     x2        = fold_test_final$X.2,
                     x3        = fold_test_final$X.3,
                     eta0_true = fold_test_final$true_eta_0,
                     tau       = fold_test_final$HTE)
    fold_test_final <- fold_test_final %>% mutate(a_t_X = as.vector(a0_hat))

  } else if (prop_score_spec == "cox-time-varying-prop") {
    # Direct event-regression estimator.
    # a_t(x) = P{W(t)=1 | event at t, X=x} — no η₀ or η₁ needed.
    # Training data: event subjects only (Delta==1 rows of the training pseudo-dataset).
    # Response: W(U_i) = treatment status at the time of the event = W column on event row.
    # Time covariate: tstop = U_i for event rows.
    # Model: additive GAM with smooth for time plus linear covariate terms.
    # Time-covariate interactions exist in the truth but are a refinement; an additive
    # GAM already improves on logistic regression by allowing a flexible time trend.
    library(mgcv)
    events  <- fold_train %>% filter(Delta == 1)
    num_cov <- ncol(fold_train %>% select(starts_with("X.")))
    cov_str <- paste0("X.", 1:num_cov, collapse = " + ")
    event_mod <- tryCatch(
      gam(as.formula(paste("W ~ s(tstop, k=5) +", cov_str)),
          family = binomial, data = events),
      error = function(e) {
        # Fallback to logistic regression if GAM fails (too few unique event times).
        glm(as.formula(paste("W ~ tstop +", cov_str)),
            family = binomial, data = events)
      }
    )
    # Evaluate â_t(X_j) for each pseudo-row (t = tstop, X = X_j) in the test fold.
    # as.vector() strips the dim attribute that predict.gam adds to its output;
    # without it, a_t_X becomes a 1D array and subsequent matrix arithmetic fails.
    a_hat <- as.vector(plogis(predict(event_mod, newdata = fold_test_final, type = "link")))
    fold_test_final <- fold_test_final %>% mutate(a_t_X = a_hat)

  } else {
    print(paste("alpha_estimate: ", alpha_estimate))
    if (prop_score_spec == "cox-intercept-only") {
      test_X         <- matrix(0, nrow = nrow(fold_test_final), ncol = 1)
      alpha_estimate <- 0
    } else if (prop_score_spec == "cox-linear-mis-specification") {
      test_X <- cbind(fold_test_final$X.1)
    } else {
      test_X <- as.matrix(fold_test_final %>% select(starts_with("X.")))
    }
    prop_scores <- calculate_eX(
      alpha_estimate = alpha_estimate,
      X = test_X,
      t = fold_test_split$tstop
    )
    fold_test_final <- fold_test_final %>% mutate(a_t_X = prop_scores)
  }

  fold_test_final <- fold_test_final %>%
    mutate(nu_X = a_t_X * eta_1 + (1 - a_t_X) * eta_0)
  
  return( list(fold_test_final = fold_test_final,
               lasso_ret = lasso_ret) )
}


fit_TV_CSL <- function(fold_causal_fitted, 
                       test_data, 
                       beta_HTE_first_stage = NULL,
                       HTE_spec = "linear",
                       final_model_method = "lasso_coxph") {
  
  beta_HTE <- NULL
  if (is.null(beta_HTE_first_stage)){
    n_regressors <- ncol(fold_causal_fitted %>% select(starts_with("X."))) ## we hand assume HTE type is linear
    beta_HTE_first_stage <- rep(0, n_regressors+1)
  }
  
  
  if (HTE_spec == "linear") {
    regressors <- cbind(1, fold_causal_fitted %>% select(starts_with("X.")) )
    n_regressors <- ncol(regressors) - 1
    
    interaction_terms <- (fold_causal_fitted$W - fold_causal_fitted$a_t_X) * regressors
    colnames(interaction_terms) <- c("intercept", paste0("interaction_X", 1:n_regressors)) 
    
    
    fold_causal_fitted <- cbind(fold_causal_fitted, interaction_terms)
    
    interaction_formula <- paste(colnames(interaction_terms), collapse = " + ")
    final_formula <- as.formula(paste("Surv(tstart, tstop, Delta) ~", interaction_formula, "+ offset(nu_X)"))
    
    if (final_model_method == "lasso_coxph"){
      final_model <- coxph(
        final_formula,
        data = fold_causal_fitted,
        ties = "breslow",
        init = beta_HTE_first_stage
      )
      beta_HTE <- coef(final_model)

      # Inference quantities: subject-level score residuals and vcov.
      # residuals(type="score") is one row per pseudo-row; sum within subject.
      vcov_k    <- vcov(final_model)
      info_k    <- solve(vcov_k)
      score_raw <- residuals(final_model, type = "score")
      if (!is.matrix(score_raw)) score_raw <- matrix(score_raw, ncol = 1)
      score_df  <- as.data.frame(score_raw)
      score_df$id <- fold_causal_fitted$id
      score_subj <- score_df %>%
        group_by(id) %>%
        summarise(across(everything(), sum), .groups = "drop") %>%
        select(-id) %>%
        as.matrix()
    }
    
    
    
    
    # devtools::install_github("CYGUBICKO/pcoxtime-pkg")
    # The paper: https://arxiv.org/pdf/2102.02297
    # library(pcoxtime)
    # final_model <- pcoxtimecv(
    #   final_formula,
    #   data = fold_causal_fitted,
    #   alphas = 0.001,  # Very close to 0 for ridge-like behavior
    #   lambdas = NULL,
    #   lamfract = 0.6,  # Use 60% of lambda sequence for speed
    #   nclusters = 4    # Use parallel processing
    # )
    # beta_HTE <- coef(final_model, s = "lambda.min")
    
    if (final_model_method == "pcoxtime"){
      library(pcoxtime)
      # Skip cross-validation by using pcoxtime directly with a single lambda
      final_model <- pcoxtime(
        final_formula,
        data = fold_causal_fitted,
        alpha = 0.001,  # Very close to 0 for ridge-like behavior
        lambda = 0.1    # Single moderate lambda value - adjust based on your needs
      )
      
      beta_HTE <- coef(final_model)
    }
    
    
    
  } else if (HTE_spec == "complex") {
    
    complex_X <- transform_X(
      single_data = fold_causal_fitted, 
      transform_spec = "complex"
    )
    
    regressor_TV_CSL <- (fold_causal_fitted$W - as.vector(fold_causal_fitted$a_t_X)) *
      cbind(1, complex_X)
    
    final_model <- cv.glmnet(
      regressor_TV_CSL, 
      Surv(fold_causal_fitted$tstart, fold_causal_fitted$tstop, fold_causal_fitted$Delta), 
      offset = fold_causal_fitted$nu_X,
      family = "cox"
    )

    beta_HTE <- as.vector(coef(final_model, s = "lambda.min"))  
  }
  
  print("beta_HTE")
  print(beta_HTE)
  
  if (HTE_spec == "linear") {
    test_regressor_HTE <- cbind(1, as.matrix(test_data %>% select(starts_with("X.")) ) )
  } else if (HTE_spec == "complex") {
    test_complex_X <- transform_X(
      single_data = test_data,
      transform_spec = "complex")
    test_regressor_HTE <- cbind(1, as.matrix(test_complex_X))
  }
  
  HTE_est <- as.vector(test_regressor_HTE %*% beta_HTE)

  # Collect inference objects if available (only for linear coxph path)
  if (!exists("vcov_k"))    vcov_k    <- NULL
  if (!exists("info_k"))    info_k    <- NULL
  if (!exists("score_subj")) score_subj <- NULL

  ret <- list(
    HTE_est    = HTE_est,
    beta_HTE   = beta_HTE,
    final_model_method = final_model_method,
    vcov_k     = vcov_k,
    info_k     = info_k,
    score_subj = score_subj
  )

  return(ret)
}

#' Time-Varying Conditional Survival Learning (TV-CSL) with K-Fold Cross-Fitting
#'
#' This function performs K-fold cross-fitting to estimate the Conditional Average Treatment Effect (HTE) using 
#' Time-Varying Conditional Survival Learning (TV-CSL). The procedure includes fitting nuisance parameters 
#' and the final model on different folds of the data.
#'
#' @param train_data A data frame containing the training data. Must include variables needed for model training and cross-fitting, including an `id` column.
#' @param test_data A data frame containing the test data. Must include the true HTE values (`HTE`).
#' @param train_data_original A data frame with the original training data used for model fitting.
#' @param folds A vector or factor specifying the fold assignment for each observation in `train_data`.
#' @param K An integer specifying the number of folds for cross-fitting.
#' @param prop_score_spec A character string specifying the type of propensity score specification to use in the nuisance model.
#' @param lasso_type A character string specifying the type of Lasso to use in the model (e.g., "linear", "complex").
#' @param regressor_spec A character string specifying the regressor transformation ("linear" or "complex").
#'
#' @return A list containing the following components:
#'   - `HTE_est`: Estimated HTE values for the test data.
#'   - `HTE_true`: True HTE values from the test data.
#'   - `MSE`: Mean squared error between the estimated and true HTE.
#'   - `time_taken`: Time taken to run the K-fold cross-fitting procedure.
#'
#' @examples
#' # Example usage of TV_CSL function
#' result <- TV_CSL(
#'   train_data = df_train, 
#'   test_data = df_test, 
#'   train_data_original = df_train_original,
#'   folds = folds_vector,
#'   K = 5, 
#'   prop_score_spec = "logit",
#'   lasso_type = "complex",
#'   regressor_spec = "linear",
#'   final_model_method = "cox"
#' )
#'
#' @export
TV_CSL <- function(train_data, 
                   test_data, 
                   train_data_original, 
                   HTE_type,
                   eta_type,
                   K, 
                   prop_score_spec, 
                   lasso_type, 
                   regressor_spec, 
                   final_model_method,
                   HTE_spec,
                   i = 0,
                   id_var = "id",
                   lasso_warmstart = 1,
                   verbose = 2) {
  
  n <- nrow(test_data)
  folds <- cut(seq(1, nrow(train_data_original)), breaks = K, labels = FALSE)
  
  HTE_ests <- matrix(NA, nrow = n, ncol = K)
  
  # Measure time taken for the procedure
  start_time <- Sys.time()
  
  {
    train_data <- train_data %>% mutate(id = !!sym(id_var))
    test_data <- test_data %>% mutate(id = !!sym(id_var))
    train_data_original <- train_data_original %>% mutate(id = !!sym(id_var))
  }
  
  
  
  # Perform K-fold cross-fitting
  first_stage_lassos <- list()
  fit_TV_CSL_rets    <- list()   # store per-fold fit for inference
  for (k in 1:K) {
    
    # Get IDs for nuisance and causal splits
    nuisance_ids <- train_data_original[folds != k, id_var]
    causal_ids <- train_data_original[folds == k, id_var]
    
    # Split the training data
    fold_nuisance <- train_data[train_data[[id_var]] %in% nuisance_ids, ]
    fold_causal <- train_data[train_data[[id_var]] %in% causal_ids, ]
    train_data_original_nuisance <- 
      train_data_original[train_data_original[[id_var]] %in% nuisance_ids, ]
    
    
    object_causal_fitted <- TV_CSL_nuisance(
      fold_train = fold_nuisance, 
      fold_test = fold_causal, 
      train_data_original = train_data_original_nuisance,
      prop_score_spec = prop_score_spec,
      lasso_type = lasso_type,
      regressor_spec = regressor_spec,
      HTE_spec = HTE_spec
    )
    fold_causal_fitted <- object_causal_fitted$fold_test_final
    first_stage_lassos[[k]] <- first_stage_lasso <- 
      object_causal_fitted$lasso_ret
    
    
    if (verbose == 2){
      output_folder <- generate_output_folder(
        results_dir = RESULTS_DIR,
        method_setting = "TV-CSL_", 
        eta_type = eta_type, 
        HTE_type = HTE_type, 
        n = n
      )
      
      if (lasso_type == "T-lasso" | lasso_type == "S-lasso"){ # this rules out m-regression
        
        save_lasso_beta(lasso_ret = first_stage_lasso, 
                        output_folder = output_folder, 
                        i = i, 
                        k = k,
                        lasso_type = lasso_type,
                        # eta_type = regressor_spec, 
                        # HTE_type = HTE_spec, 
                        eta_spec = regressor_spec, 
                        HTE_spec = HTE_spec,
                        prop_score_spec = prop_score_spec,
                        stage = "first") 
        
        save_lasso_MSE(lasso_ret = first_stage_lasso, 
                       HTE_true = test_data$HTE , 
                       output_folder = output_folder, 
                       i = i,
                       k = k,
                       lasso_type = lasso_type,
                       # eta_type = regressor_spec, 
                       # HTE_type = HTE_spec, 
                       eta_spec = regressor_spec, 
                       HTE_spec = HTE_spec,
                       prop_score_spec = prop_score_spec,
                       stage = "first")
      }
    }
    
    
    if (lasso_warmstart){
      beta_HTE_first_stage <- 
        if ( lasso_type == "T-lasso" | lasso_type == "S-lasso") 
          first_stage_lasso$beta_HTE else NULL
      if (k == 1){
        beta_HTE_first_stages <- matrix(NA, nrow = length(beta_HTE_first_stage), ncol = K)
      }
      
      beta_HTE_first_stages[, k] <- beta_HTE_first_stage
    }else{
      beta_HTE_first_stage <- NULL
      beta_HTE_first_stages <- NULL
    }
    
    fit_TV_CSL_ret <- fit_TV_CSL(
      fold_causal_fitted = fold_causal_fitted,
      test_data = test_data,
      HTE_spec = HTE_spec,
      beta_HTE_first_stage = beta_HTE_first_stage,
      final_model_method = final_model_method
    )
    fit_TV_CSL_rets[[k]] <- fit_TV_CSL_ret
    
    
    fit_TV_CSL_ret$beta_eta_0 <- 
      if ( lasso_type == "T-lasso" | lasso_type == "S-lasso") 
        first_stage_lasso$beta_eta_0 else NA
    
    if (k == 1){
      beta_HTEs <- matrix(NA, nrow = length(fit_TV_CSL_ret$beta_HTE), ncol = K)
    }
    
    beta_HTEs[, k] <- fit_TV_CSL_ret$beta_HTE
    HTE_ests[, k] <- fit_TV_CSL_ret$HTE_est
    HTE_true <- test_data$HTE 
    MSE_this_round <- mean((HTE_true - HTE_ests[, k])^2)
    
    
    if (verbose == 2){
      # Save beta_HTE from second stage to the file 
      save_lasso_beta(lasso_ret = fit_TV_CSL_ret, 
                      output_folder = output_folder, 
                      i = i, 
                      k = k,
                      lasso_type = lasso_type, 
                      # eta_type = regressor_spec, 
                      # HTE_type = HTE_spec, 
                      eta_spec = regressor_spec, 
                      HTE_spec = HTE_spec,
                      prop_score_spec = prop_score_spec,
                      stage = "final") 
      
      save_lasso_MSE(lasso_ret = fit_TV_CSL_ret, 
                     HTE_true = test_data$HTE , 
                     output_folder = output_folder, 
                     i = i,
                     k = k,
                     lasso_type = lasso_type, 
                     # eta_type = regressor_spec, 
                     # HTE_type = HTE_spec, 
                     eta_spec = regressor_spec, 
                     HTE_spec = HTE_spec,
                     prop_score_spec = prop_score_spec,
                     stage = "final")
      
    }
    if (verbose >= 1){
      print("MSE_this_round: ")
      print(MSE_this_round)
    }
    
  
  }
  
  # Average HTE estimates across folds
  beta_HTE <- rowMeans(beta_HTEs, na.rm = TRUE)
  HTE_est  <- rowMeans(HTE_ests,  na.rm = TRUE)
  HTE_true <- test_data$HTE
  MSE      <- mean((HTE_true - HTE_est)^2)
  time_taken <- as.numeric(Sys.time() - start_time, units = "secs")

  # ---- Inference -------------------------------------------------------
  # Only computed when all folds used the linear coxph path.
  has_inference <- !is.null(fit_TV_CSL_rets[[1]]$vcov_k)

  if (has_inference) {
    vcov_list  <- lapply(fit_TV_CSL_rets, `[[`, "vcov_k")
    info_list  <- lapply(fit_TV_CSL_rets, `[[`, "info_k")
    score_list <- lapply(fit_TV_CSL_rets, `[[`, "score_subj")

    # Naive SE: Var(beta_hat) = Var((1/K) sum beta_k) = (1/K^2) sum Vcov_k
    vcov_naive <- Reduce("+", vcov_list) / K^2
    se_naive   <- sqrt(diag(vcov_naive))

    # Sandwich SE: stack subject-level scores across folds
    all_scores <- do.call(rbind, score_list)   # n_subjects x d
    n_subj     <- nrow(all_scores)
    Sigma_hat  <- t(all_scores) %*% all_scores / n_subj
    # J_hat = (1/n) * sum of per-fold observed information
    # solve(vcov_k) = observed info for n/K subjects, so sum / n gives J
    J_hat      <- Reduce("+", info_list) / n_subj
    J_inv      <- solve(J_hat)
    Omega_hat  <- J_inv %*% Sigma_hat %*% J_inv
    se_sandwich <- sqrt(diag(Omega_hat) / n_subj)
  } else {
    d           <- length(beta_HTE)
    se_naive    <- rep(NA_real_, d)
    se_sandwich <- rep(NA_real_, d)
    vcov_naive  <- NULL
    Omega_hat   <- NULL
  }

  return(list(
    first_stage_lassos    = first_stage_lassos,
    beta_HTE_first_stages = beta_HTE_first_stages,
    beta_HTE    = beta_HTE,
    beta_HTEs   = beta_HTEs,
    HTE_est     = HTE_est,
    HTE_ests    = HTE_ests,
    HTE_true    = HTE_true,
    MSE         = MSE,
    time_taken  = time_taken,
    se_naive    = se_naive,
    se_sandwich = se_sandwich,
    vcov_naive  = vcov_naive,
    Omega_hat   = Omega_hat
  ))
}

