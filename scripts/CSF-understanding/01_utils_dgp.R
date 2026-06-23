# =============================================================
# 01_utils_dgp.R
# DGP1: log-HR treatment effect, Weibull events (RMST-CATE is
#        nonlinear in X1 on RMST scale)
# DGP2: RMST-CATE is exactly linear in X1 by construction
#        (exponential events, solve arm-1 rate numerically)
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

# =============================================================
# DGP2: Exactly linear RMST-CATE
#
# Event model: T(a)|X=x ~ Exp(r_a(x))
# Baseline rate r_0(x) is complex (drives complex baseline RMST).
# r_1(x) is solved so that RMST_1(x;h) - RMST_0(x;h) = gamma0 + gamma1*X1
# exactly.  Censoring is the same Weibull model as DGP1.
#
# True tau (analytical): tau_RMST(x;h) = gamma0 + gamma1*X1
# =============================================================

DGP2_PARAMS <- list(
  p        = 15,
  h        = 2.0,
  gamma0   = 0.10,   # RMST-CATE intercept
  gamma1   = 0.30,   # RMST-CATE slope in X1; tau in [0.10, 0.40]
  kappa_c  = 1.2,    # censoring (same as DGP1)
  lambda_c = 0.15
)

# Complex baseline rate: drives RMST_0 variation across X2-X6.
# Constructed so RMST_0(x) <= h - max(tau) = 1.60 always,
# guaranteeing RMST_1 = RMST_0 + tau < h.
r0fun_dgp2 <- function(X) {
  b2 <- -0.3 +
    0.5 * sin(2 * pi * X[, 2]) +
    0.5 * (X[, 3] - 0.5)^2 +
    0.5 * X[, 4] * X[, 5] -
    0.4 * as.numeric(X[, 6] > 0.5)
  exp(b2)   # range ~[0.30, 2.28]; RMST_0 in ~[0.43, 1.50]
}

# Invert RMST_exp(r, h) = (1-exp(-r*h))/r to find r given target RMST.
# Function is strictly decreasing in r: h at r->0, 0 at r->Inf.
solve_r1_dgp2 <- function(r0, tau_rmst, h) {
  target <- (1 - exp(-r0 * h)) / r0 + tau_rmst
  target <- min(target, h - 1e-6)   # safety: RMST < h
  uniroot(
    function(r) (1 - exp(-r * h)) / r - target,
    lower = 1e-8, upper = 1e4, tol = 1e-10
  )$root
}

simulate_data_dgp2 <- function(n, params = DGP2_PARAMS) {
  p <- params$p; h <- params$h
  kappa_c <- params$kappa_c; lambda_c <- params$lambda_c

  X <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  e_x  <- efun(X)          # same propensity as DGP1
  A    <- rbinom(n, 1, e_x)

  r0   <- r0fun_dgp2(X)
  tau_x <- params$gamma0 + params$gamma1 * X[, 1]
  r1   <- mapply(solve_r1_dgp2, r0, tau_x, MoreArgs = list(h = h))

  r_a     <- ifelse(A == 1L, r1, r0)
  T_event <- rexp(n) / r_a          # Exp(r_a) event times

  eta_C  <- cfun(X) + 0.3 * A       # same censoring as DGP1
  T_cens <- (rexp(n) / (lambda_c * exp(eta_C)))^(1 / kappa_c)

  list(
    X       = X,
    A       = A,
    U       = pmin(T_event, T_cens),
    Delta   = as.integer(T_event <= T_cens),
    e_true  = e_x,
    r0      = r0,
    r1      = r1
  )
}

# True RMST-CATE for DGP2: analytical
compute_true_tau_dgp2 <- function(X_test, params = DGP2_PARAMS) {
  params$gamma0 + params$gamma1 * X_test[, 1]
}
