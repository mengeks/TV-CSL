#' Estimate the Counterfactual Outcome E[Y(W)|X]
#'
#' This function estimates the expected outcome E[Y(W)|X] given a trained model, 
#' a binary treatment indicator W, and a set of confounders X. 
#' The function can be used to obtain counterfactual estimates for either the 
#' treated (W = 1) or control (W = 0) group.
#'
#' @param model A trained model object that predicts E[Y|X, W]. 
#' This could be any model that has a `predict` method.
#' @param W A binary value (0 or 1) indicating whether to return the 
#' counterfactual estimate for the control group (W = 0) or the treated group (W = 1).
#' This is not a vector, but a single binary number.
#' @param X A data frame or matrix of confounders (covariates) for which the 
#' counterfactual estimate is required. The number of rows should correspond 
#' to the number of observations.
#'
#' @return A numeric vector containing the predicted counterfactual outcomes.
#' @export
#'
#' @examples 
#' # Assuming 'trained_model' is a model object and 'confounders' is a data frame of covariates
#' get_counterfactual_estimation(trained_model, W = 1, X = confounders)
get_counterfactual_estimation <- function(model, W, X) {
  counterfactual_est <- predict(
    model, 
    newdata = data.frame(X = X, W = W)
  )
  return(counterfactual_est)
}


get_outcome_predictions <- function(model, X) {
  pred_control <- get_counterfactual_estimation(model, W = 0, X = X)
  pred_treated <- get_counterfactual_estimation(model, W = 1, X = X)
  return(list(pred_control = pred_control, pred_treated = pred_treated))
}

get_bias <- function(model){
  model_summary <- 
    summary(model)
  model_bias <- (model_summary$coefficients[2, 1] - beta_t)
  return(model_bias)
}