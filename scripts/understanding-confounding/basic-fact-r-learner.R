set.seed(123)

# Generate data
n <- 500
X <- rnorm(n)  # Generate X from a normal distribution
g_X <- X^2 -X + 3     # Define Y as a function of X, e.g., Y = X^2
epsilon <- rnorm(n)  # Noise term independent of X

# Define tilde(X) such that E[tilde(X) | X] = 0
tilde_X <- epsilon  # epsilon is independent noise, hence E[tilde(X) | X] = 0

# Define Y as some function of X
Y <- g_X  # For instance, Y = X^2

# Plot tilde(X) vs Y
png("./figures/confound-and-de-confound/basic-fact.png", width = 800, height = 600)
plot(Y, tilde_X, 
     xlab = "Y = g(X)", 
     ylab = expression(tilde(X)), 
     main = "Scatter Plot of tilde(X) vs Y",
     pch = 16, col = "blue")
abline(h = 0, col = "red", lwd = 2)  # Horizontal line at y = 0

# Add a legend
legend("topright", legend = c("tilde(X) vs Y", "E[tilde(X) | X] = 0"),
       col = c("blue", "red"), pch = c(16, NA), lwd = c(NA, 2))
dev.off()


#########
########
#######
set.seed(123)

# Generate X from a normal distribution N(1, 1)
n <- 1000
X <- rnorm(n, mean = 0, sd = 1)

# Generate Z from a mixture of N(1, 1), N(2, 1), N(-3, 1)
# Define the mixture components
mixture_components <- c(10, 20, -30)
mixture_sd <- 1
mixture_prob <- c(1/3, 1/3, 1/3)

# Sample which component to use for each observation
component <- sample(mixture_components, size = n, replace = TRUE, prob = mixture_prob)

# Generate Z based on the sampled components
Z <- rnorm(n, mean = component, sd = mixture_sd)


# Calculate tilde_X = X * Z
tilde_X <- X * Z

# Y = X
Y <- X

# Plot tilde_X against Y
plot(tilde_X, Y, 
     xlab = expression(tilde(X)), ylab = "Y = X", 
     main = expression(paste("Scatter Plot of ", tilde(X), " vs Y")),
     col = "blue", pch = 16,
     xlim=c(-12,12))

# Add a line to show the mean of Y
abline(h = mean(Y), col = "red", lwd = 2)
abline(lm(Y~tilde_X), col="green")
