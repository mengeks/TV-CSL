library(MASS)
library(here)
utils_folder <- "utils"
source(here(utils_folder,"plotter.R"))
source(here(utils_folder,"model.R"))
set.seed(123)

n <- 1000  
beta_c <- 2 
beta_t <- 3 
tau <- 1.5

X <- rnorm(n, mean = 0, sd = 1)
W <- rbinom(n, 1, prob = plogis(beta_c * X))
Y <- beta_t * W + tau * X + rnorm(n)

illustrate_confounding_ATE(X=X, Y=Y, W=W)
ggsave(
  here("figures", "understanding-confounding-ATE.png"),
    width = 4.8,
    height = 3.6
  )

