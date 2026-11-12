library(survival)
library(glmnet)
library(tidyverse)
library(here)

source(here::here("R/cox-loglik.R"))
source(here::here("R/data-handler.R"))

# ---------------------------------------------------------------------------
# Data preprocessing helpers
# ---------------------------------------------------------------------------

#' Transform a survival dataset into the (tstart, tstop, Delta, W) pseudo format
#' required for time-varying Cox models.
create_pseudo_dataset <- function(survival_data) {
  pseudo_dataset <- tibble(
    tstart = numeric(),
    tstop  = numeric(),
    Delta  = numeric(),
    W      = numeric()
  )

  covariates <- setdiff(colnames(survival_data), c("U", "Delta", "A", "id"))

  # coxph rejects (near-)zero-length intervals. Never drop an event for this: if adoption
  # falls within eps of U, move the split to U - eps (shifts A by at most eps).
  eps <- 1e-6

  for (i in 1:nrow(survival_data)) {
    U_i     <- survival_data$U[i]
    Delta_i <- survival_data$Delta[i]
    A_i     <- survival_data$A[i]
    id_i    <- survival_data$id[i]
    covariate_values <- survival_data[i, covariates, drop = FALSE]

    if (U_i <= A_i && A_i <= Inf) {
      new_rows <- tibble(tstart = 0, tstop = U_i, Delta = Delta_i, W = 0)
    } else if (A_i < U_i) {
      A_split <- min(A_i, U_i - eps)
      new_rows <- tibble(
        tstart = c(0,       A_split),
        tstop  = c(A_split, U_i),
        Delta  = c(0,       Delta_i),
        W      = c(0,       1)
      )
    }

    new_rows <- new_rows %>%
      mutate(id = id_i) %>%
      bind_cols(covariate_values)

    pseudo_dataset <- bind_rows(pseudo_dataset, new_rows)
  }

  # Drop only event-free intervals shorter than eps (negligible at-risk time).
  pseudo_dataset %>%
    filter(tstop - tstart > eps | Delta == 1)
}

#' Convert raw survival data to the format expected by Cox/lasso estimators.
preprocess_data <- function(single_data, run_time_varying) {
  if (run_time_varying) {
    create_pseudo_dataset(survival_data = single_data)
  } else {
    single_data %>% mutate(W = as.numeric(A <= U))
  }
}

#' Generate the covariate part of a Cox formula string.
generate_regressor_part <- function(model_spec  = "correctly-specified",
                                    HTE_type    = "constant",
                                    eta_type    = "linear-interaction") {
  if (HTE_type == "constant") {
    W_part <- "W"
  } else if (HTE_type == "linear") {
    W_part <- "W:X.1 + W:X.10"
  }

  if (model_spec == "correctly-specified") {
    if (eta_type == "10-dim-non-linear") {
      return(paste(W_part, "+ sqrt(abs(X.1 * X.2)) + sqrt(abs(X.10)) + cos(X.5) + cos(X.5) * cos(X.6)"))
    } else {
      return(paste(W_part, "+ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8"))
    }
  } else if (model_spec == "mildly-mis-specified") {
    return(paste(W_part, "+ X.1 + X.2 + X.3 + X.4 + X.5"))
  } else if (model_spec == "quite-mis-specified") {
    return(paste(W_part, "+ X.1 + X.2"))
  }
}

#' Build a coxph formula for a Cox model estimation.
create_cox_formula <- function(model_spec,
                               run_time_varying,
                               HTE_type = "constant",
                               eta_type = "linear-interaction") {
  regressor_part <- generate_regressor_part(
    model_spec = model_spec,
    HTE_type   = HTE_type,
    eta_type   = eta_type
  )

  if (run_time_varying) {
    as.formula(paste("Surv(tstart, tstop, Delta) ~", regressor_part))
  } else {
    as.formula(paste("Surv(U, Delta) ~", regressor_part))
  }
}

# ---------------------------------------------------------------------------
# Main TV-CSL estimation code (TV_CSL, S_cox, T_lasso, S_lasso, run_* fns)
# ---------------------------------------------------------------------------
source(here::here("R/TV-CSL.R"))
