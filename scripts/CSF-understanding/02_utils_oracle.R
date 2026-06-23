# =============================================================
# 02_utils_oracle.R
# Oracle (true-DGP) nuisance functions and ground-truth tau
# Used for: (1) computing true RMST-CATE on test set,
#           (2) oracle linear CSF estimating equation.
# =============================================================

# Assumes 01_utils_dgp.R has been sourced already.

# ---- True event and censoring survival functions ----

# S_T(t | x, w) = exp(-lambda * t^kappa * exp(b(x) - w*delta(x)))
S_T_fun <- function(t, xrow, w, params = DGP_PARAMS) {
  xmat <- matrix(xrow, nrow = 1)
  phi  <- exp(bfun(xmat) - w * deltafun(xmat))
  exp(-params$lambda * t^params$kappa * phi)
}

# S_C(t | x, w) = exp(-lambda_c * t^kappa_c * exp(c(x) + 0.3*w))
S_C_fun <- function(t, xrow, w, params = DGP_PARAMS) {
  xmat  <- matrix(xrow, nrow = 1)
  eta_c <- cfun(xmat) + 0.3 * w
  exp(-params$lambda_c * t^params$kappa_c * exp(eta_c))
}

# Cumulative censoring hazard Λ_C(t|x,w) = -log S_C(t|x,w)  [analytical]
Lambda_C_fun <- function(t, xrow, w, params = DGP_PARAMS) {
  xmat  <- matrix(xrow, nrow = 1)
  eta_c <- cfun(xmat) + 0.3 * w
  params$lambda_c * t^params$kappa_c * exp(eta_c)
}

# ---- True RMST for arm w ----

rmst_true <- function(xrow, w, h = DGP_PARAMS$h, params = DGP_PARAMS) {
  xmat  <- matrix(xrow, nrow = 1)
  phi   <- exp(bfun(xmat) - w * deltafun(xmat))
  alpha <- params$lambda * phi
  kappa <- params$kappa
  integrate(
    function(t) exp(-alpha * t^kappa),
    lower = 1e-10, upper = h, rel.tol = 1e-5
  )$value
}

# True RMST-CATE: tau(x; h) = rmst1(x) - rmst0(x)
tau_rmst_true <- function(xrow, h = DGP_PARAMS$h, params = DGP_PARAMS) {
  rmst_true(xrow, 1, h, params) - rmst_true(xrow, 0, h, params)
}

# Vectorized true RMST-CATE over test matrix (slow, used once per MC)
compute_true_tau <- function(X_test, h = DGP_PARAMS$h, params = DGP_PARAMS) {
  apply(X_test, 1, function(xrow) tau_rmst_true(xrow, h, params))
}

# ---- Q_w(s | x, h): conditional expected T∧h given T∧h > s ----
# Q_w(s) = s + ∫_s^h exp(-alpha*(t^kappa - s^kappa)) dt
# Build an approxfun over s ∈ [0, h] to avoid nested integration

make_Q_approx <- function(xrow, w, h = DGP_PARAMS$h,
                          n_grid = 80, params = DGP_PARAMS) {
  xmat  <- matrix(xrow, nrow = 1)
  phi   <- exp(bfun(xmat) - w * deltafun(xmat))
  alpha <- params$lambda * phi
  kappa <- params$kappa

  s_grid <- seq(0, h, length.out = n_grid + 1)

  Q_grid <- sapply(s_grid, function(s) {
    if (s >= h - 1e-10) return(h)
    val <- tryCatch(
      integrate(
        function(t) exp(-alpha * (t^kappa - s^kappa)),
        lower = s, upper = h,
        rel.tol = 1e-4, abs.tol = 1e-8
      )$value,
      error = function(e) {
        # Fallback: rough trapezoidal rule
        tt <- seq(s, h, length.out = 30)
        sum(diff(tt) * exp(-alpha * (head(tt, -1)^kappa - s^kappa)))
      }
    )
    s + val
  })

  approxfun(s_grid, Q_grid, rule = 2)
}

# ---- Marginal mean E[T∧h | X=x] = e(x)*rmst1 + (1-e(x))*rmst0 ----

m_true <- function(xrow, h = DGP_PARAMS$h, params = DGP_PARAMS) {
  xmat <- matrix(xrow, nrow = 1)
  e    <- efun(xmat)
  e * rmst_true(xrow, 1, h, params) + (1 - e) * rmst_true(xrow, 0, h, params)
}

# ---- Oracle B_i and H_i for one training observation ----
#
# H_i = 1/S_C(ũ) - Λ_C(ũ)                [fully analytical]
# B_i = (Q_w(ũ) + Δ^h*(y^h - Q_w(ũ)) - m) / S_C(ũ) - ∫_0^ũ (dΛ_C/ds) * (Q(s)-m) ds
#
# where ũ = min(U_i, h), Δ^h = Δ_i OR (U_i >= h)

compute_BH_one_oracle <- function(i, X, U, Delta, W,
                                  h = DGP_PARAMS$h, params = DGP_PARAMS) {
  xrow    <- X[i, ]
  w       <- W[i]
  u_h     <- min(U[i], h)
  delta_h <- as.numeric(Delta[i] == 1 || U[i] >= h)
  y_h     <- u_h   # y(T) = T ∧ h

  # Pre-compute censoring quantities
  xmat      <- matrix(xrow, nrow = 1)
  eta_c     <- as.numeric(cfun(xmat)) + 0.3 * w
  exp_eta_c <- exp(eta_c)
  kappa_c   <- params$kappa_c
  lambda_c  <- params$lambda_c

  sc_u  <- S_C_fun(u_h, xrow, w, params)   # S_C(ũ)
  lam_u <- Lambda_C_fun(u_h, xrow, w, params)  # Λ_C(ũ) = -log S_C(ũ)

  # H_i (analytical)
  H <- 1 / sc_u - lam_u

  # Q approxfun and m
  Q_approx <- make_Q_approx(xrow, w, h, params = params)
  m_x      <- m_true(xrow, h, params)
  q_u      <- Q_approx(u_h)

  # ∫_0^ũ (dΛ_C/ds) * (Q(s) - m) ds
  # dΛ_C/ds = kappa_c * lambda_c * s^(kappa_c-1) * exp(eta_c)
  int_B <- 0
  if (u_h > 1e-10) {
    int_B <- tryCatch(
      integrate(
        function(s) {
          dL <- kappa_c * lambda_c * s^(kappa_c - 1) * exp_eta_c
          dL * (Q_approx(s) - m_x)
        },
        lower = 1e-10, upper = u_h,
        rel.tol = 1e-4, abs.tol = 1e-8
      )$value,
      error = function(e) 0
    )
  }

  B <- (q_u + delta_h * (y_h - q_u) - m_x) / sc_u - int_B

  c(B = B, H = H)
}

# Compute BH for all training observations (sequential)
compute_BH_oracle <- function(X, U, Delta, W,
                              h = DGP_PARAMS$h, params = DGP_PARAMS,
                              verbose = FALSE) {
  n      <- nrow(X)
  result <- matrix(NA_real_, n, 2, dimnames = list(NULL, c("B", "H")))
  for (i in seq_len(n)) {
    if (verbose && i %% 200 == 0)
      cat("  BH:", i, "/", n, "\n")
    result[i, ] <- tryCatch(
      compute_BH_one_oracle(i, X, U, Delta, W, h, params),
      error = function(e) c(B = NA_real_, H = NA_real_)
    )
  }
  result
}
