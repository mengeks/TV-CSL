library(ggplot2)

create_plot <- function(data, naive_pred, adjusted_pred) {
  p <- ggplot(data, aes(x = X, y = Y, color = as.factor(W))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("blue", "red"), labels = c("Control", "Treated")) +
    labs(color = "Group", x = expression(x[1]), y = "Y", title = "Additive Error Model") +
    geom_line(aes(y = naive_pred$pred_control), color = "green", size = 1) +
    geom_line(aes(y = naive_pred$pred_treated), color = "green", size = 1) +
    geom_line(aes(y = adjusted_pred$pred_control), color = "yellow", size = 1) +
    geom_line(aes(y = adjusted_pred$pred_treated), color = "yellow", size = 1) +
    theme_minimal()
  
  p <- p + 
    scale_color_manual(values = c("blue", "red"),
                       labels = c("Control", "Treated")) +
    guides(color = guide_legend(override.aes = list(
      linetype = c(0, 0),  # No lines for points
      shape = c(16, 16),   # Points for Control and Treated
      size = c(3, 3)       # Size of points
    ))) +
    geom_line(aes(x = X, y = naive_pred$pred_control), color = "green", linetype = 1, size = 1) +
    geom_line(aes(x = X, y = naive_pred$pred_treated), color = "green", linetype = 1, size = 1) +
    geom_line(aes(x = X, y = adjusted_pred$pred_control), color = "yellow", linetype = 1, size = 1) +
    geom_line(aes(x = X, y = adjusted_pred$pred_treated), color = "yellow", linetype = 1, size = 1)
  
  p
}


#' Title
#'
#' @param X confounder
#' @param Y outcome
#' @param W treatment
#'
#' @return A plot showing how confounding will add bias to ATE estimation
#' @export
#'
#' @examples
illustrate_confounding_ATE <- function(X,Y,W) {
  
  data <- data.frame(X = X, Y = Y, W = W)
  
  naive_model <- lm(Y ~ W, data = data)
  adjusted_model <- lm(Y ~ W + X, data = data)
  
  cat("Bias in Naive Model (ignoring confounder):")
  cat(naive_bias <- get_bias(naive_model))
  
  cat("Bias in Adjusted Model (ignoring confounder):")
  cat(adjusted_bias <- get_bias(adjusted_model))
  
  naive_pred <- 
    get_outcome_predictions(naive_model, X)
  adjusted_pred <- 
    get_outcome_predictions(adjusted_model, X)
  p <- create_plot(data, naive_pred, adjusted_pred)
  
  print(p)
  return(p)
}
