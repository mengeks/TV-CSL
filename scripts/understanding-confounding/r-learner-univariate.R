set.seed(42)

# Generate data
n <- 2000
X <- runif(n, -2, 2)
epsilon <- rnorm(n, 0, 1)

tau <- 1

fn_eta0 <- function(x){
  1.5 *(1 + x+3*x^2)
}

eta0 <- fn_eta0(X) 
eta1 <- eta0 +tau 

# Generate Y(0) and Y(1)
Y0 <- eta0 + epsilon
Y1 <- eta1 + epsilon

# Generate treatment indicator
propensity_score <- exp(X) / (1 + exp(X))
treatment <- rbinom(n, 1, propensity_score)
W <- treatment

# Observed outcomes
Y <- ifelse(treatment == 1, Y1, Y0)

png("./figures/understanding-confounding/separate-nuisance.png",
    width=480,
    height=480)
# Plot the data
plot(X[treatment == 0], Y[treatment == 0], col = "blue", pch = 16, xlab = expression(x[1]), ylab = "Y", main = "Additive Error Model with Interaction")
points(X[treatment == 1], Y[treatment == 1], col = "red", pch = 16)
# curve(x^2, add = TRUE, col = "black", lwd = 2) # true function curve

# Fit regression model with common slope and different intercepts
model_lin <- lm(Y ~ 1 + X + W)
cat("Bias of the Naive estimate: \n")
print(summary(model_lin)$coefficients[3,1] - tau)

model_nuisance <- lm(Y ~ 1 + X)
m_x <- predict(
  model_nuisance,
  newdata = data.frame(X=X,W=W))

X_seq <- seq(min(X), max(X), length.out = 100)

# Predict values for control and treated groups
pred_control <- 
  predict(
    model_lin, 
    newdata = data.frame(X = X_seq, W = 0)
    )
pred_treated <- 
  predict(
    model_lin, 
    newdata = data.frame(X = X_seq, W = 1)
  )

# Add regression lines to the plot
lines(X_seq, pred_control, col = "green", lwd = 2)
lines(X_seq, pred_treated, col = "green", lwd = 2)

### 
W_minus_eX <- treatment - propensity_score
# Fit a linear model of Y ~ W_minus_eX
model_r_learner <- lm(
  Y ~ W_minus_eX, 
  offset = m_x)
cat("Bias of the R-learner estimate: \n")
print(summary(model_r_learner)$coefficients[2,1] - tau)

library(rlearner)
rlasso_fit = rlasso(x, w, y)

legend("topleft", legend = c("Control", "Treated", "True Function", "Fitted Regression with Interaction"),
       col = c("blue", "red", "black", "green"), pch = c(16, 16, NA, NA), lwd = c(NA, NA, 2, 2))
dev.off()


#########
#### R-learner illustration
#########

# Plot the data
# Save plot to a PNG file
png("./figures/understanding-confounding/R-learner-illustration.png", width = 800, height = 600)

# Plot the data
# plot(W_minus_eX, Y-m_x,
plot(W_minus_eX, fn_eta0(X) + tau * propensity_score,
     col = ifelse(treatment == 1, 
                  "red", "blue"), pch = 16, 
     xlab = expression(W - e(X)), ylab = "Y", 
     main = expression(paste("Plot of ", Y, " vs ", W - e(X))),
     # ylim=c(-5,5)
     )
legend("topleft", legend = c("Control", "Treated"),
       col = c("blue", "red"), pch = 16)

# Fit a linear model of Y ~ W_minus_eX
model <- lm(Y ~ 0+W_minus_eX, offset = m_x)
summary(model)$coefficient[1,1] - tau 
abline(model, col = "green", lwd = 2)
# So the problem is that the quality of the offset makes a difference...

model_correct_m <- lm(Y ~ 0+W_minus_eX, offset = fn_eta0(X))
summary(model_correct_m)$coefficient[1,1] - tau 

# Close the device
dev.off()

# png("./figures/confound-and-de-confound/plot-w-resid-mx.png", width = 800, height = 600)
# plot(W_minus_eX, X^2 + 1 * propensity_score, 
#      col = ifelse(treatment == 1, 
#                   "red", "blue"), pch = 16, 
#      xlab = expression(W - e(X)), ylab = "Y", 
#      main = expression(paste("Plot of ", Y, " vs ", W - e(X))))
# dev.off()

# The slope for this regression should be 
#   zero because 

# set.seed(42)
# 
# # Generate data
# n <- 2000
# X <- runif(n, -2, 2)
# epsilon <- rnorm(n, 0, 0.05)
# 
# # Define eta functions
# eta0 <- X^2
# eta1 <- X^2 + 1
# 
# # Generate treatment indicator
# propensity_score <- exp(X) / (1 + exp(X))
# treatment <- rbinom(n, 1, propensity_score)

# # Observed outcomes
# Y0 <- eta0 + epsilon
# Y1 <- eta1 + epsilon
# Y <- ifelse(treatment == 1, Y1, Y0)

