cox_loglik_old <- function(
    beta, time, status, covar, strata, 
    offset = 0, weights = NULL, method = 0) {
  
  n <- length(time)  # Number of observations
  
  # Check if covar is a vector (1D) and handle accordingly
  if (is.vector(covar)) {
    covar <- matrix(covar, ncol = 1)  # Convert to a matrix with one column
  }
  
  nvar <- ncol(covar)  # Number of covariates
  
  if (is.null(weights)) weights <- rep(1, n)
  
  loglik <- 0  # Initialize log-likelihood
  
  # Check if offset is scalar or vector
  is_offset_scalar <- length(offset) == 1
  
  for (i in 1:n) {
    if (strata[i] == 1) {  # Only process the first event in each strata
      risk_set <- which(time >= time[i])  # Define risk set
      event_set <- which(time == time[i] & status == 1)  # Define event set
      
      for (event in event_set) {
        
        # Handle case where offset is scalar or vector
        if (is_offset_scalar) {
          offset_event <- offset  # Use the scalar offset
          offset_risk_set <- offset
        } else {
          offset_event <- offset[event]  # Use the vector offset for event
          offset_risk_set <- offset[risk_set]  # Use the vector offset for the risk set
        }
        
        if (length(beta) == 1) {
          zbeta <- offset_event + covar[event, , drop = FALSE] * beta
          risk_scores <- exp(offset_risk_set + covar[risk_set, , drop = FALSE] * beta)
        } else {
          zbeta <- offset_event + covar[event, , drop = FALSE] %*% beta
          risk_scores <- exp(offset_risk_set + covar[risk_set, , drop = FALSE] %*% beta)
        }
        
        # If there is an issue with zbeta or risk_scores, print values
        if (any(!is.finite(zbeta))) {
          cat("Non-finite zbeta encountered\n")
          print(zbeta)
        }
        
        if (any(!is.finite(risk_scores))) {
          cat("Non-finite risk_scores encountered\n")
          print(risk_scores)
        }
        
        risk_sum <- sum(risk_scores)
        
        # Check for non-finite risk_sum before taking the log
        if (!is.finite(risk_sum) || risk_sum <= 0) {
          cat("Non-finite or non-positive risk_sum encountered\n")
          print(risk_sum)
        } else {
          loglik <- loglik + (zbeta) - log(risk_sum)
        }
        
        # Check for non-finite log-likelihood
        if (!is.finite(loglik)) {
          cat("Non-finite log-likelihood encountered\n")
          print(loglik)
          print(beta)
        }
      }
    }
  }
  
  return(-loglik)  # Return negative log-likelihood (for minimization in optim)
}

cox_loglik <- function(
    beta, time, status, covar, strata, 
    offset = 0, weights = NULL, method = 0, verbose = 0) {
  n <- length(time)
  
  if (is.vector(covar)) {
    covar <- matrix(covar, ncol = 1)
  }
  
  nvar <- ncol(covar)
  
  if (is.null(weights)) weights <- rep(1, n)
  
  loglik <- 0
  is_offset_scalar <- length(offset) == 1
  
  for (i in 1:n) {
    if (strata[i] == 1) {
      risk_set <- which(time >= time[i])
      event_set <- which(time == time[i] & status == 1)
      
      for (event in event_set) {
        if (is_offset_scalar) {
          offset_event <- offset
          offset_risk_set <- offset
        } else {
          offset_event <- offset[event]
          offset_risk_set <- offset[risk_set]
        }
        
        if (length(beta) == 1) {
          zbeta <- offset_event + covar[event, , drop = FALSE] * beta
          risk_scores <- exp(offset_risk_set + covar[risk_set, , drop = FALSE] * beta)
        } else {
          zbeta <- offset_event + covar[event, , drop = FALSE] %*% beta
          risk_scores <- exp(offset_risk_set + covar[risk_set, , drop = FALSE] %*% beta)
        }
        
        if (verbose > 0) {
          if (abs(zbeta) > 1000) {
            cat("Large zbeta detected for event:", event, "\n")
            cat("Beta:", beta, "\n")
            cat("Covariate values:", covar[event, , drop = FALSE], "\n")
          }
          
          if (!is.finite(zbeta)) cat("Non-finite zbeta encountered\n")
        }
        
        risk_sum <- sum(risk_scores)
        
        if (!is.finite(risk_sum) || risk_sum <= 0) {
          next
        } else {
          loglik <- loglik + (zbeta) - log(risk_sum)
        }
        
        if (verbose > 0) {
          if (!is.finite(loglik)) {
            cat("Non-finite log-likelihood encountered\n")
          }
        }
      }
    }
  } # End loop i
  
  return(-loglik)
}

cox_loglik_new <- function(beta, time, status, covar, strata, offset = 0) {
  loglik <- 0
  
  for (i in 1:n) {
    if (strata[i] == 1) {
      risk_set <- which(time >= time[i])
      event_set <- which(time == time[i] & status == 1)
      
      for (event in event_set) {
        # Calculate linear predictors
        lin_pred <- covar[risk_set, ] %*% beta + offset
        
        # Find maximum for stability
        max_pred <- max(lin_pred)
        
        # Stabilized computation
        risk_scores <- exp(lin_pred - max_pred)
        log_risk_sum <- log(sum(risk_scores)) + max_pred
        
        # Update loglik
        loglik <- loglik + (covar[event, ] %*% beta + offset) - log_risk_sum
      }
    }
  }
  return(-loglik)
}


cox_loglik_tv <- function(
    beta, start_time, stop_time, status, covar, strata, 
    offset = 0, 
    weights = NULL, 
    method = 0) {

  unique_stop_time <- unique(stop_time)
  n <- length(unique_stop_time)
  
  if (is.vector(covar)) {
    covar <- matrix(covar, ncol = 1)
  }
  
  nvar <- ncol(covar)
  
  if (is.null(weights)) weights <- rep(1, n)
  
  loglik <- 0
  is_offset_scalar <- length(offset) == 1
  
  for (i in 1:n) {
    if (strata[i] == 1) {
      stop_time_i <- unique_stop_time[i]
      event_set <- which(stop_time == stop_time_i & status == 1)
      risk_set <- which(stop_time == stop_time_i)
      
      
      for (event in event_set) {
        if (is_offset_scalar) {
          offset_event <- offset
          offset_risk_set <- offset
        } else {
          offset_event <- offset[event]
          offset_risk_set <- offset[risk_set]
        }
        
        if (length(beta) == 1) {
          zbeta <- offset_event + covar[event, , drop = FALSE] * beta
          risk_scores <- exp(offset_risk_set + covar[risk_set, , drop = FALSE] * beta)
        } else {
          zbeta <- offset_event + covar[event, , drop = FALSE] %*% beta
          risk_scores <- exp(offset_risk_set + covar[risk_set, , drop = FALSE] %*% beta)
        }
        
        if (abs(zbeta) > 1000) {
          cat("Large zbeta detected for event:", event, "\n")
          cat("Beta:", beta, "\n")
          cat("Covariate values:", covar[event, , drop = FALSE], "\n")
        }
        
        if (!is.finite(zbeta)) cat("Non-finite zbeta encountered\n")
        
        risk_sum <- sum(risk_scores)
        
        if (!is.finite(risk_sum) || risk_sum <= 0) {
          next
        } else {
          loglik <- loglik + (zbeta) - log(risk_sum)
        }
        
        if (!is.finite(loglik)) {
          cat("Non-finite log-likelihood encountered\n")
        }
      }
    }
  }
  
  return(-loglik)
}


fit_custom_cox_model <- 
  function(init_beta, time, status, covar, 
           strata, offset = 0, verbose = 0) {
    
    track_optim <- function(par, time, status, covar, strata, offset) {
      print("par value: ")
      print(signif(par, 3) )
      
      coxlik_value <- cox_loglik(par, time, status, covar, strata, offset)
      
      print(paste0("coxlik_value: ", coxlik_value))
      
      return(coxlik_value)
    }
    
    if (verbose == 0){
      result <- optim(par = init_beta, 
                      fn = cox_loglik,
                      # fn = track_optim,
                      time = time, 
                      status = status, 
                      covar = covar, 
                      strata = strata, 
                      offset = offset,  # Include offset parameter
                      method = "BFGS")
    }else if (c == 1){
      result <- optim(par = init_beta, 
                      # fn = cox_loglik, 
                      fn = track_optim,
                      time = time, 
                      status = status, 
                      covar = covar, 
                      strata = strata, 
                      offset = offset,  # Include offset parameter
                      method = "BFGS")
    }else {
      stop("verbose must be 0 or 1. Detected other values")
    }
    
    
    
    # Return the result
    return(result)
  }




