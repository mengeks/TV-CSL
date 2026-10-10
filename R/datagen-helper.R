library(simsurv)
suppressPackageStartupMessages(library(tidyverse))
library(parallel)  # Load parallel here
library(survival)
library(here)

generate_X <- 
  function(n, p = 5, 
           distribution = "uniform", 
           cov_type = "identity") {
  if (distribution == "uniform") {
    X <- matrix(runif(n * p, min = -1, max = 1), n, p)
  } else if (distribution == "normal") {
    if (cov_type == "identity") {
      Sigma <- diag(p)
    } else if (cov_type == "toeplitz") {
      Sigma <- toeplitz((p:1)/p)
    } else {
      stop("Unsupported covariance type")
    }
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else {
    stop("Unsupported distribution")
  }
  return(X)
}

calculate_eta <- function(x, 
                          eta_type = "linear-interaction",
                          linear_intercept = 0,
                          linear_slope_multiplier = 2.5) {
  # Hardcoded betas and delta
  
  
  p <- length(x)
  
  if (eta_type == "non-linear") {
    sigma <- function(x) 2 / (1 + exp(-12 * (x - 0.5)))
    eta_0 <- (-3/4) * sigma(x[1]) * sigma(x[2])
  } else if (eta_type == "linear") {
    eta_0 <- linear_intercept + linear_slope_multiplier * sum(sapply(1:p, function(j) x[j] / j))
  }
  
  
  if (p == 5 | p == 10){
    betas <- rep(1, 5)
    delta <- 0.5
    X1 <- x[1]
    X2 <- x[2]
    X3 <- x[3]
    X4 <- x[4]
    X5 <- x[5]
    X6 <- x[6]
    X10 <- x[10]
    
    if (eta_type == "linear-interaction") {
      eta_0 <- X1 * betas[1] +
        X2 * betas[2] + 
        X3 * betas[3] + 
        X4 * betas[4] + 
        X5 * betas[5] + 
        delta * X1 * X2
      
    } else if (eta_type == "10-dim-non-linear") {
      eta_0 <- -2 + (2/3) * (
        - 0.5 * sqrt(abs(X1 * X2)) +
          0.5 * sqrt(abs(X10)) -
          0.5 * cos(X5) +
          0.5 * cos(X5) * cos(X6)
      )
      
    } else if (eta_type == "log") {
      eta_0 <- 2 * log(1 + exp(X1 + X2 + X3))
      
    } else {
      stop("Unsupported eta_type")
    }
    
  }# end p == 5 | p == 10
  
  
  return(eta_0)
}

HTE_function <- function(x, 
                          HTE_type = "zero",
                          linear_HTE_multiplier = 1) {
  
  # Handle different cases based on 'HTE_type'
  if (HTE_type == "zero") {
    return(0)
    
  } else if (HTE_type == "constant") {
    return(2)
    
  } else if (HTE_type == "ReLU") {
    # Ensure that x has at least 5 dimensions
    if (length(x) < 5) stop("ReLU HTE_type requires at least 5 dimensions")
    return(max(x[1] + x[2] + x[3], 0) - max(x[4] + x[5], 0))
    
  } else if (HTE_type == "linear") {
    p <- length(x)
    
    beta <- rep(1, p)

    # beta <- rep(0, p)
    # beta[1] <- beta[p] <- 1
    
    return(sum(x * beta * linear_HTE_multiplier))
    
  } else if (HTE_type == "non-linear") {
    sigma <- function(x) 2 / (1 + exp(-12 * (x - 0.5)))
    
    if (length(x) < 2) stop("Non-linear HTE_type requires at least 2 dimensions")
    
    # return(sigma(x[1]) * sigma(x[10]))
    return(sigma(x[1]) * sigma(x[2]))
    
  } else {
    stop("Unsupported HTE_type")
  }
}

calculate_baseline_hazard <- function(t, baseline_type = "cosine") {
  if (baseline_type == "linear") {
    baseline_hazard <- t
  } else if (baseline_type == "constant") {
    baseline_hazard <- 1
  } else {  # Default is "cosine"
    baseline_hazard <- (cos(t * 3) + 1) / 2
  }
  
  return(baseline_hazard)
}

calculate_hazard <- function(t, 
                             x, 
                             is_time_varying, 
                             baseline_type = "linear"
) {
  baseline_hazard <- 
    calculate_baseline_hazard(
      t = t, 
      baseline_type = baseline_type
    )
  
  if (is_time_varying) {
    hazard <- exp(x[, "eta_0"] + (t >= x[, "A"]) * x[, "HTE"]) * baseline_hazard
  } else {
    # hazard <- exp(eta_0 + x[, "W"] * x[, "HTE"]) * baseline_hazard
    hazard <- exp(x[, "eta_0"] + x[, "W"] * x[, "HTE"]) * baseline_hazard
  }
  
  return(hazard)
}

# Cumulative baseline hazard B(t) = int_0^t h0(s) ds, matching calculate_baseline_hazard().
calculate_cumulative_baseline_hazard <- function(t, baseline_type = "cosine") {
  if (baseline_type == "linear") {
    t^2 / 2
  } else if (baseline_type == "constant") {
    t
  } else {  # Default is "cosine"
    (t + sin(3 * t) / 3) / 2
  }
}

# Inverse of B(t). The cosine baseline has no closed form, so it is solved numerically;
# B(t) >= (t - 1/3) / 2 gives the upper bracket 2y + 1/3.
invert_cumulative_baseline_hazard <- function(y, baseline_type = "cosine") {
  if (baseline_type == "linear") {
    sqrt(2 * y)
  } else if (baseline_type == "constant") {
    y
  } else {
    vapply(y, function(yi) {
      if (!is.finite(yi)) return(Inf)
      uniroot(function(t) calculate_cumulative_baseline_hazard(t, baseline_type) - yi,
              lower = 0, upper = 2 * yi + 1/3, tol = 1e-12)$root
    }, numeric(1))
  }
}

# Draw event times exactly by inverting the cumulative hazard of calculate_hazard().
# (simsurv's 15-node quadrature is inaccurate across the hazard jump at A.)
# Time-varying: H(t) = exp(eta_0) * [B(min(t, A)) + exp(HTE) * (B(t) - B(min(t, A)))].
# Returns the same columns as simsurv: id, eventtime (truncated at maxt), status.
simulate_event_times <- function(covariates, is_time_varying, baseline_type = "linear", maxt = 20) {
  n <- nrow(covariates)
  y <- rexp(n) * exp(-covariates$eta_0)          # = B(T) if treatment had no effect
  if (is_time_varying) {
    B_A <- calculate_cumulative_baseline_hazard(covariates$A, baseline_type)
    after <- y > B_A
    y[after] <- B_A[after] + (y[after] - B_A[after]) * exp(-covariates$HTE[after])
  } else {
    y <- y * exp(-covariates$W * covariates$HTE)
  }
  eventtime <- invert_cumulative_baseline_hazard(y, baseline_type)
  data.frame(id        = covariates$id,
             eventtime = pmin(eventtime, maxt),
             status    = as.numeric(eventtime < maxt))
}




generate_treatment <-
  function(n, X, is_time_varying, difficulty = "simple") {

  simple_fn <- function(X) X[, 2] + X[, 3]
  
  hard_fn <- function(X) sin(pi * X[, 1] * X[, 2])
  
  super_hard_fn <- function(X) {
    L1 <- X[, 1]
    L2 <- X[, 2]
    L5 <- X[, 5]
    L6 <- X[, 6]
    L10 <- ifelse(ncol(X) >= 10, X[, 10], 0) # Handle if X has only 5 columns
    
    exp(-2 * sqrt(abs(L1 * L2)) + 
          2 * sqrt(abs(L10)) - 
          2 * cos(L5) + 
          2 * cos(L5) * cos(L6))
  }
  
  complex_fn <- function(X)   2 * ( sqrt(abs(X[,1] * X[,2])) -  cos (X[,3]) )
  
  # Treatment generation logic
  if (difficulty == "simple") {
    fn <- simple_fn
  } else if (difficulty == "hard") {
    fn <- hard_fn
  } else if (difficulty == "super-hard") {
    fn <- super_hard_fn
  } else if (difficulty == "complex"){
    fn <- complex_fn
  } else if (difficulty == "constant"){
    if (is_time_varying) {
      fn <- function(X) log(1/ 3)
    } else {
      fn <- function(X) 0.5
    }
  }
  
  
  
  if (difficulty == "uniform") {
    if (is_time_varying) {
      A <- runif(n, min = 0, max = 20)
    } else {
      W <- rbinom(n, 1, prob = 0.5)
    }
  }  else {
    if (is_time_varying) {
      A <- rexp(n, rate = exp(fn(X)))
    } else {
      expit <- function(x) 1 / (1 + exp(-x))  # Logistic function
      W <- rbinom(n, 1, prob = expit(fn(X)))
    }
  }
  
  # Return the generated treatment
  if (is_time_varying) {
    return(A)
  } else {
    return(W)
  }
}

generate_covariates <- function(n, 
                                p = 5, 
                                X_distribution = "uniform", 
                                X_cov_type = "identity", 
                                is_time_varying = TRUE, 
                                tx_difficulty = "simple", 
                                HTE_type = "constant",
                                linear_HTE_multiplier = 1,
                                eta_type = "linear",
                                linear_intercept = 0,
                                linear_slope_multiplier = 2.5) {
  
  # Step 1: Generate X using the generate_X function
  X <- generate_X(n = n, 
                  p = p, 
                  distribution = X_distribution, 
                  cov_type = X_cov_type)
  
  
  # Step 2: Generate HTE using the HTE_function
  HTE <- apply(X, 1, function(row) HTE_function(x = row, HTE_type = HTE_type, linear_HTE_multiplier=linear_HTE_multiplier))
  
  # Step 3: Generate eta_0
  eta_0 <- apply(X, 1, function(row) calculate_eta(x = row, 
                         eta_type = eta_type,
                         linear_intercept = linear_intercept,
                         linear_slope_multiplier = linear_slope_multiplier) )
  
  
  # Step 3: Generate treatment based on is_time_varying
  if (is_time_varying) {
    A <- generate_treatment(n = n, 
                            X = X, 
                            is_time_varying = TRUE, 
                            difficulty = tx_difficulty)
    return(data.frame(id=1:n,X = X, A = A, HTE = HTE, eta_0 = eta_0))  # Return X, A, and HTE for time-varying case
  } else {
    W <- generate_treatment(n = n, 
                            X = X, 
                            is_time_varying = FALSE, 
                            difficulty = tx_difficulty)
    return(data.frame(id=1:n, X = X, W = W, HTE = HTE, eta_0 = eta_0))  # Return X, W, and HTE for non-time-varying case
  }
}


generate_censoring_times <- 
  function(n, light_censoring, lambda_C) {
  if (light_censoring) {
    return(rep(20, n))
  } else {
    return(rexp(n, rate = lambda_C))
  }
}




generate_simulated_data <- 
  function(n, 
           is_time_varying = TRUE, 
           light_censoring = FALSE,
           lambda_C = 0.1,
           p = 3,  
           baseline_type = "linear",
           eta_type = "linear-interaction",
           X_distribution = "normal", 
           X_cov_type = "identity",
           tx_difficulty = "simple",
           HTE_type = "constant",
           linear_intercept = 0,
           linear_slope_multiplier = 2.5,
           linear_HTE_multiplier = 1, 
           max_censoring_time = 20,
           seed_value = 123,
           verbose = 0) {
  
    verbose_print <- function(message, level = 1) {
      if (verbose >= level) {
        cat(message, "\n")
      }
    }
  verbose_print("Step 0: Setting seed", 2)
  set.seed(seed_value)
    
  verbose_print("Step 1: Generating covariates...", 2)
  covariates <- generate_covariates(
    n = n, 
    p = p, 
    X_distribution = X_distribution, 
    X_cov_type = X_cov_type, 
    tx_difficulty = tx_difficulty,
    is_time_varying = is_time_varying,
    HTE_type = HTE_type,
    eta_type = eta_type,
    linear_HTE_multiplier = linear_HTE_multiplier,
    linear_intercept = linear_intercept,
    linear_slope_multiplier = linear_slope_multiplier
  )
  
  verbose_print("Step 2: Generating censoring times...", 2)
  start_time <- Sys.time()
  C_gen <- 
    generate_censoring_times(n, light_censoring, lambda_C)
  # verbose_print( paste0("C_gen is: ", ) , 2)
  C <- pmin(max_censoring_time, C_gen)
  end_time <- Sys.time()
  verbose_print(sprintf("Censoring times generated in %.2f seconds.", as.numeric(end_time - start_time)), 1)
  
  
  verbose_print(paste0("baseline_type passed to generate_simulated_data is:    ", baseline_type), 2)
  verbose_print(head(covariates), 1)

  # Step 3: Simulate survival data by exact inversion of the cumulative hazard
  verbose_print("Step 3: Simulating survival data...", 2)
  start_time <- Sys.time()
  simulated_data <- simulate_event_times(
    covariates      = covariates,
    is_time_varying = is_time_varying,
    baseline_type   = baseline_type,
    maxt            = max_censoring_time
  )
  end_time <- Sys.time()
  verbose_print(sprintf("Survival data simulated in %.2f seconds.", as.numeric(end_time - start_time)), 1)
  
  processed_data <- 
    post_process(simulated_data = simulated_data, 
                covariates = covariates, 
                C = C, 
                C_gen = C_gen,
                verbose = verbose)
  
  verbose_print("Data generation completed.", 1)
  
  return(processed_data)
}

post_process <- function(simulated_data, 
                         covariates,
                         C,
                         C_gen,
                         verbose = 0) {
  
  # Helper function for printing based on verbosity level
  verbose_print <- function(message, level = 1) {
    if (verbose >= level) {
      cat(message, "\n")
    }
  }
  
  # Step 1: Merge covariates with simulated data
  verbose_print("Step 1: Merging covariates with simulated data...", 2)
  processed_data <- 
    merge(simulated_data, 
          covariates, by = "id")
  
  # Step 2: Calculate observed times and event indicators
  verbose_print("Step 2: Calculating observed times and event indicators...", 2)
  processed_data <- processed_data %>%
    mutate(U = pmin(eventtime, C),  # Observed time is the minimum of event time and censoring time
           # Event before censoring; status == 0 means still alive at maxt (eventtime was truncated
           # to maxt), which must stay censored even when C == maxt.
           Delta = as.numeric(eventtime <= C & status == 1))
  
  # Step 3: Rename eventtime to T
  verbose_print("Step 3: Renaming eventtime to T...", 2)
  processed_data <- processed_data %>%
    rename(T = "eventtime")
  
  # Step 4: Add censoring times to the final dataset
  verbose_print("Step 4: Adding censoring times...", 2)
  processed_data <- 
    cbind(processed_data, C = C, C_gen = C_gen)
  
  return(processed_data)
}


generate_and_save_data <- 
  function(i, 
           n, 
           path_for_sim_data, 
           params,
           verbose = 0) {
    
    start_time <- Sys.time()  # Start timing
    
    seed_value <- 123 + 11 * i
    params$seed <- seed_value
    
    # Print iteration details
    cat("Iteration:", i, 
        "Data size (n):", n, 
        "eta_type:", params$eta_type, 
        "HTE_type:", params$HTE_type, "\n")
    
    # Optional params fall back to generate_simulated_data()'s defaults when absent.
    param_or <- function(name, default) if (is.null(params[[name]])) default else params[[name]]

    simulated_data_i_n <-
      generate_simulated_data(
        n,
        is_time_varying = params$is_time_varying,
        light_censoring = params$light_censoring,
        lambda_C = params$lambda_C,
        p = params$p,
        baseline_type = param_or("baseline_type", "linear"),
        eta_type = params$eta_type,
        X_distribution = params$X_distribution,
        X_cov_type = params$X_cov_type,
        tx_difficulty = params$tx_difficulty,
        HTE_type = params$HTE_type,
        linear_intercept = param_or("linear_intercept", 0),
        linear_slope_multiplier = param_or("linear_slope_multiplier", 2.5),
        linear_HTE_multiplier = param_or("linear_HTE_multiplier", 1),
        seed_value = params$seed,
        verbose = verbose
      )
    
    dataset_with_params <- list(
      data = simulated_data_i_n,
      params = params
    )
    
    file_name <- paste0("sim_data_", i, "_n_", n, ".rds")
    file_path <- file.path(path_for_sim_data, file_name)
    saveRDS(dataset_with_params, file_path)
    
    end_time <- Sys.time()  # End timing
    time_used <- end_time - start_time  # Calculate time used
    
    # Print time taken
    cat("Saved dataset", i, "to", file_path, "\n")
    cat("Time used:", time_used, "\n")
  }

