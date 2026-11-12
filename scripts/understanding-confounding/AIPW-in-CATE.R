# Load necessary libraries
library(MASS)
library(dplyr)
library(estimatr)

# Simulation parameters
set.seed(123)
n <- 1000  # Sample size
beta <- 1  # True treatment effect
gamma <- 0.5  # Coefficient for propensity score model
n_sim <- 1000  # Number of simulations

# Storage for results
gcomp_results <- matrix(NA, nrow = n_sim, ncol = 100)
aipw_results <- matrix(NA, nrow = n_sim, ncol = 100)
X_vals <- seq(-3, 3, length.out = 100)  # Range of X values

# Run the simulation n_sim times
for (i in 1:n_sim) {
  print(i)
  # Generate covariate and treatment assignment
  X <- rnorm(n)  # Covariate
  p <- 1 / (1 + exp(-gamma * X))  # Propensity score
  W <- rbinom(n, 1, p)  # Treatment assignment
  
  # Generate potential outcomes
  Y0 <- X + rnorm(n)  # Outcome under control
  Y1 <- X + beta + rnorm(n)  # Outcome under treatment
  
  # Observed outcome
  Y <- ifelse(W == 1, Y1, Y0)
  
  # Misspecified outcome model (ignores interaction with X)
  outcome_model <- lm(Y ~ W)  # Misspecified model
  
  # Predicted outcomes under treatment and control using the misspecified model
  Y_hat_0_mis <- predict(outcome_model, newdata = data.frame(W = 0, X = X))
  Y_hat_1_mis <- predict(outcome_model, newdata = data.frame(W = 1, X = X))
  
  # Correct propensity score model
  ps_model <- glm(W ~ X, family = binomial)
  ps <- predict(ps_model, type = "response")
  
  # Estimate HTE for a range of X values
  cate_gcomp_estimates <- sapply(X_vals, function(x) {
    mean(Y_hat_1_mis[abs(X - x) < 0.1]) - mean(Y_hat_0_mis[abs(X - x) < 0.1])
  })
  
  cate_aipw_estimates <- sapply(X_vals, function(x) {
    ind <- which(abs(X - x) < 0.1)
    term1 <- (W[ind] / ps[ind]) * (Y[ind] - Y_hat_1_mis[ind]) + Y_hat_1_mis[ind]
    term2 <- ((1 - W[ind]) / (1 - ps[ind])) * (Y[ind] - Y_hat_0_mis[ind]) + Y_hat_0_mis[ind]
    mean(term1 - term2)
  })
  
  # Store results
  gcomp_results[i, ] <- cate_gcomp_estimates
  aipw_results[i, ] <- cate_aipw_estimates
}

# Compute means across simulations
mean_gcomp <- colMeans(gcomp_results)
mean_aipw <- colMeans(aipw_results)

# Plot the average HTE estimates
plot(X_vals, mean_gcomp, type = "l", col = "red", lwd = 2,
     xlab = "X", ylab = "HTE", main = "Comparison of AIPW and g-computation for HTE Estimation")
lines(X_vals, mean_aipw, col = "blue", lwd = 2)
abline(h = beta, col = "green", lty = 2)  # True HTE line

legend("topright", legend = c("g-computation", "AIPW", "True HTE"),
       col = c("red", "blue", "green"), lty = 1, lwd = 2)

# Summary of the results
cat("Mean HTE from g-computation across X values: \n")
print(mean(mean_gcomp))
cat("Mean HTE from AIPW across X values: \n")
print(mean(mean_aipw))
