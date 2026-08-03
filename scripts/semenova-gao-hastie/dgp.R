#!/usr/bin/env Rscript
# DGP: Cox PH with linear hazard-scale CATE
#
#   X1, X2 ~ N(0,1)
#   W | X  ~ Bernoulli(expit(0.3*X1 - 0.3*X2))
#   T | W,X ~ Exp(0.1 * exp(0.5*X2 + W*(-0.5 - 0.5*X1)))
#   C | W,X ~ Exp(0.05 * exp(0.2*W + 0.2*X2))
#   Y = min(T,C),  Delta = 1(T <= C)
#
#   tau(X) = -0.5 - 0.5*X1  =  p(X)' beta,  p(X) = (1, X1),  beta = (-0.5, -0.5)

DGP_PARAMS <- list(
  beta_true  = c(-0.5, -0.5),   # true CATE coefficients
  lam0       = 0.1,              # baseline hazard
  lam_c      = 0.05,             # baseline censoring hazard
  eta0_coef  = 0.5,              # eta0(X) = eta0_coef * X2
  e_coef     = c(0.3, -0.3),    # propensity logit = e_coef[1]*X1 + e_coef[2]*X2
  cens_coef  = c(0.2, 0.2)      # censoring log-rate = cens_coef[1]*W + cens_coef[2]*X2
)

# Generate one dataset of size n
simulate_data <- function(n, params = DGP_PARAMS, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params

  X   <- matrix(rnorm(n * 2), n, 2)
  e   <- plogis(p$e_coef[1] * X[, 1] + p$e_coef[2] * X[, 2])
  W   <- rbinom(n, 1, e)

  tau  <- p$beta_true[1] + p$beta_true[2] * X[, 1]
  eta0 <- p$eta0_coef * X[, 2]

  Tevt <- rexp(n, rate = p$lam0 * exp(eta0 + W * tau))
  Cens <- rexp(n, rate = p$lam_c * exp(p$cens_coef[1] * W + p$cens_coef[2] * X[, 2]))

  Y   <- pmin(Tevt, Cens)
  Del <- as.integer(Tevt <= Cens)

  list(X = X, W = W, Y = Y, Del = Del, e_true = e, tau_true = tau,
       eta0_true = eta0, event_rate = mean(Del))
}

# Basis p(X) = (1, X1)
make_basis <- function(X) cbind(1, X[, 1])

# True nuisances (oracle)
oracle_nuisances <- function(dat, params = DGP_PARAMS) {
  p   <- params
  X   <- dat$X
  W   <- dat$W

  e    <- plogis(p$e_coef[1] * X[, 1] + p$e_coef[2] * X[, 2])
  eta0 <- p$eta0_coef * X[, 2]
  tau  <- p$beta_true[1] + p$beta_true[2] * X[, 1]
  eta1 <- eta0 + tau

  # R_fn_vec(t_vec, w): n x length(t_vec) matrix of R_w(t_k, X_i)
  # S_w(t|X) = exp(-lam0 * t * exp(eta_w(X)))
  # G_w(t|X) = exp(-lam_c * t * exp(cens_coef[1]*w + cens_coef[2]*X2))
  R_fn_vec <- function(t_vec, w) {
    eta_w  <- eta0 + w * tau
    lam_cw <- p$lam_c * exp(p$cens_coef[1] * w + p$cens_coef[2] * X[, 2])
    S_w    <- exp(-outer(p$lam0 * exp(eta_w), t_vec))   # n x K
    G_w    <- exp(-outer(lam_cw, t_vec))                 # n x K
    S_w * G_w
  }

  list(e = e, eta0 = eta0, eta1 = eta1, R_fn_vec = R_fn_vec)
}
