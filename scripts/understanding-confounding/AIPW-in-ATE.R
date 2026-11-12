# Load necessary libraries
library(MASS)       # For mvrnorm (generating multivariate normal data)
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting
library(estimatr)   # For AIPW estimator

# Simulation parameters
n <- 1000              # Sample size
beta <- 1              # True treatment effect
gamma <- 0.5           # True coefficient for covariate in propensity score model

# Generate covariate and treatment assignment
set.seed(123)
X <- rnorm(n)          # Covariate
p <- 1 / (1 + exp(-gamma * X))  # Propensity score
W <- rbinom(n, 1, p)   # Treatment assignment

# Generate potential outcomes
Y0 <- X + rnorm(n)     # Outcome under control
Y1 <- X + beta + rnorm(n)  # Outcome under treatment

# Observed outcome
Y <- ifelse(W == 1, Y1, Y0)

# Misspecified outcome model (ignores X)
Y_hat_0_mis <- mean(Y[W == 0])  # Constant prediction
Y_hat_1_mis <- mean(Y[W == 1])  # Constant prediction

# Correct propensity score model
ps_model <- glm(W ~ X, family = binomial)
ps <- predict(ps_model, type = "response")

# g-computation estimator (using misspecified model)
ate_gcomp <- 
  mean(Y_hat_1_mis) - mean(Y_hat_0_mis)

# AIPW estimator
Y_hat_0 <- Y_hat_0_mis
Y_hat_1 <- Y_hat_1_mis

aipw_estimator <- function(Y, W, ps, Y_hat_0, Y_hat_1) {
  term1 <- (W / ps) * (Y - Y_hat_1) + Y_hat_1
  term2 <- ((1 - W) / (1 - ps)) * (Y - Y_hat_0) + Y_hat_0
  mean(term1 - term2)
}

ate_aipw <- 
  aipw_estimator(Y, W, ps, Y_hat_0, Y_hat_1)

# Print results
cat("True ATE: ", beta, "\n")
cat("ATE from g-computation (misspecified outcome model): ", ate_gcomp, "\n")
cat("ATE from AIPW (correct PS, misspecified outcome model): ", ate_aipw, "\n")

# Repeat the simulation many times to evaluate bias
n_sim <- 1000
gcomp_results <- numeric(n_sim)
aipw_results <- numeric(n_sim)

for (i in 1:n_sim) {
  print(i)
  # Generate new data
  X <- rnorm(n)
  p <- 1 / (1 + exp(-gamma * X))
  W <- rbinom(n, 1, p)
  Y0 <- X + rnorm(n)
  Y1 <- X + beta + rnorm(n)
  Y <- ifelse(W == 1, Y1, Y0)
  
  # Misspecified outcome model
  Y_hat_0_mis <- mean(Y[W == 0])
  Y_hat_1_mis <- mean(Y[W == 1])
  
  # Correct propensity score model
  ps_model <- glm(W ~ X, family = binomial)
  ps <- predict(ps_model, type = "response")
  
  # g-computation estimator
  gcomp_results[i] <- mean(Y_hat_1_mis) - mean(Y_hat_0_mis)
  
  # AIPW estimator
  Y_hat_0 <- Y_hat_0_mis
  Y_hat_1 <- Y_hat_1_mis
  aipw_results[i] <- aipw_estimator(Y, W, ps, Y_hat_0, Y_hat_1)
}

# Plot results
results_df <- data.frame(
  Method = rep(c("g-computation", "AIPW"), each = n_sim),
  Estimate = c(gcomp_results, aipw_results)
)

ggplot(results_df, aes(x = Method, y = Estimate, fill = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = beta, linetype = "dashed", color = "red") +
  labs(title = "Comparison of AIPW and g-computation Estimators",
       y = "ATE Estimate",
       x = "Estimator") +
  theme_minimal()

# Summary of results
summary_gcomp <- mean(gcomp_results)
summary_aipw <- mean(aipw_results)

cat("Mean ATE from g-computation across simulations: ", summary_gcomp, "\n")
cat("Mean ATE from AIPW across simulations: ", summary_aipw, "\n")
