# =============================================================
# 01_utils_dgp.R
# DGP for CSF Simulation Study (log-HR treatment effect,
# RMST-CATE target, Weibull event and censoring times)
# =============================================================

expit <- function(z) 1 / (1 + exp(-z))

# Propensity score
efun <- function(X) {
  expit(-0.2 + 1.2 * (X[, 1] - 0.5) -
          1.0 * (X[, 2] - 0.5) +
          0.8 * sin(2 * pi * X[, 3]))
}

# Baseline log-hazard
bfun <- function(X) {
  1.2 * sin(2 * pi * X[, 2]) +
    0.8 * (X[, 3] - 0.5)^2 +
    0.7 * X[, 4] * X[, 5] -
    0.6 * as.numeric(X[, 6] > 0.5)
}

# Treatment effect on log-HR scale: larger X1 → stronger benefit
deltafun <- function(X) {
  0.4 + 0.8 * X[, 1]
}

# Censoring prognostic score
cfun <- function(X) {
  -0.4 + 0.8 * X[, 2] - 0.5 * X[, 3] +
    0.5 * as.numeric(X[, 7] > 0.5)
}

DGP_PARAMS <- list(
  p        = 15,    # covariate dimension
  kappa    = 1.5,   # event Weibull shape
  lambda   = 0.25,  # event Weibull scale
  kappa_c  = 1.2,   # censoring Weibull shape
  lambda_c = 0.15,  # censoring Weibull scale
  h        = 2.0    # RMST horizon
)

# Simulate n observations
# Returns list: X, A, U (obs time), Delta (event indicator),
#               plus true e, b, delta for debugging
simulate_data <- function(n, params = DGP_PARAMS) {
  p        <- params$p
  kappa    <- params$kappa
  lambda   <- params$lambda
  kappa_c  <- params$kappa_c
  lambda_c <- params$lambda_c

  X <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  e_x     <- efun(X)
  A       <- rbinom(n, 1, e_x)

  b_x     <- bfun(X)
  delta_x <- deltafun(X)

  # T ~ Weibull:  T = (E / (lambda * exp(b - A*delta)))^(1/kappa), E ~ Exp(1)
  eta_T   <- b_x - A * delta_x
  T_event <- (rexp(n) / (lambda * exp(eta_T)))^(1 / kappa)

  # C ~ Weibull with arm-dependent hazard
  eta_C   <- cfun(X) + 0.3 * A
  T_cens  <- (rexp(n) / (lambda_c * exp(eta_C)))^(1 / kappa_c)

  list(
    X       = X,
    A       = A,
    U       = pmin(T_event, T_cens),
    Delta   = as.integer(T_event <= T_cens),
    e_true  = e_x,
    b_true  = b_x,
    delta_true = delta_x
  )
}
