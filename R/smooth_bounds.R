tmle_smooth <- function(A, Y, mu0, mu1, pi, threshold, smoothness, parameter = "trimmed", maxiter = 1e2) {
  fluctuation <- \(epsilon, mu0, mu1, pi) {
    cleverA <- (mu1 * s_gt_dot(pi, threshold, smoothness) - mu0 * s_lt_dot(pi, 1 - threshold, smoothness))
    if(parameter == "lower") cleverA <- cleverA + s_lt_dot(pi, 1 - threshold, smoothness)
    if(parameter == "upper") cleverA <- cleverA - s_gt_dot(pi, threshold, smoothness)
    
    list(
      mu0 = plogis(qlogis(mu0) - epsilon / (1 - pi) * s_lt(pi, 1 - threshold, smoothness)),
      mu1 = plogis(qlogis(mu1) + epsilon / pi * s_gt(pi, threshold, smoothness)),
      pi  = plogis(qlogis(pi)  + epsilon * cleverA)
    )
  }
  
  loss <- \(params, mu0, mu1, pi) {
    f <- fluctuation(params, mu0, mu1, pi) 
    mean(
      -A * log(f$pi) - (1 - A) * log(1 - f$pi) + ifelse(A == 1, -Y * log(f$mu1) - (1 - Y) * log(1 - f$mu1), -Y * log(f$mu0) - (1 - Y) * log(1 - f$mu0))
    )
  }
  
  # Start at initial estimators
  mu0_star <- mu0
  mu1_star <- mu1
  pi_star  <- pi
  for(iter in 1:maxiter) {
    epsilon_star <- optimize(loss, interval = c(-1, 1), mu0 = mu0_star, mu1 = mu1_star, pi = pi_star)$minimum
    f <- fluctuation(epsilon_star, mu0_star, mu1_star, pi_star)
    
    mu0_star <- f$mu0
    mu1_star <- f$mu1
    pi_star  <- f$pi
    
    if(abs(epsilon_star) < 1e-2) break
  }
  psi_trimmed <- mean(mu1_star * s_gt(pi_star, threshold, smoothness) - mu0_star * s_lt(pi_star, 1 - threshold, smoothness))
  
  if(parameter == "trimmed") {
    psi <- psi_trimmed
    eif <- eif_trimmed(A, Y, mu0_star, mu1_star, pi_star, threshold, smoothness)
  }
  else if(parameter == "upper") {
    psi <- psi_trimmed + 1 - mean(s_gt(pi_star, threshold, smoothness))
    eif <- eif_upper(A, Y, mu0_star, mu1_star, pi_star, threshold, smoothness)
  }
  else if(parameter == "lower") {
    psi <- psi_trimmed - 1 + mean(s_lt(pi_star, 1 - threshold, smoothness))
    eif <- eif_lower(A, Y, mu0_star, mu1_star, pi_star, threshold, smoothness)
  }
  
  list(
    psi = psi,
    eif = eif,
    mu0 = mu0,
    mu1 = mu1,
    pi = pi,
    ci = psi + qnorm(c(0.025, 0.975)) * sd(eif) / sqrt(length(Y))
  )
}

smooth_bounds <- function(data, X, A, Y, learners_trt = c("glm"), learners_outcome = c("glm"), method = "tmle", folds = 5, thresholds = c(0, 10^seq(-4, 0, 0.05)), smoothness = 0, bootstrap = TRUE, bootstrap_draws = 1e3, nuisance = NULL) {
  method <- str_to_lower(method)
  if(!(method %in% c("tmle", "onestep"))) stop("Method must be 'tmle' or 'onestep'")
  if(!all(X %in% names(data))) stop("All X columns must be in input data")
  if(!(A %in% names(data))) stop(glue::glue("treatment column '{A}' must be in input data"))
  if(!(Y %in% names(data))) stop(glue::glue("outcome column '{Y}' must be in input data"))
  if(any(sort(unique(data[[Y]])) != c(0, 1))) stop("Outcome must be binary")
  
  K <- length(thresholds)
  N <- nrow(data)
  
  # Cross-fitted nuisance models
  if(!is.null(nuisance)) {
    pi_hat <- nuisance$pi_hat
    mu0_hat <- nuisance$mu0_hat
    mu1_hat <- nuisance$mu1_hat
    mu_hat <- nuisance$mu_hat
  }
  else {
    data0 <- data1 <- data
    data0[[A]] <- 0
    data1[[A]] <- 1
    pi_hat <- mu0_hat <- mu1_hat <- numeric(N)
    
    cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = folds)
    if(folds > 1) {
      for(fold in seq_along(cv)) {
        training   <- cv[[fold]]$training_set
        validation <- cv[[fold]]$validation_set
        pi_model   <- mlr3superlearner::mlr3superlearner(data[training, c(X, A)], target = A, library = learners_trt, outcome_type = "binomial", folds = 5)
        mu_model   <- mlr3superlearner::mlr3superlearner(data[training, c(X, A, Y)], target = Y, library = learners_outcome, outcome_type = "binomial", folds = 5)
        
        pi_hat[validation]  <- predict(pi_model, newdata = data[validation, ])
        mu0_hat[validation] <- predict(mu_model, newdata = data0[validation, ])
        mu1_hat[validation] <- predict(mu_model, newdata = data1[validation, ])
      }
    }
    else {
      pi_model   <- mlr3superlearner::mlr3superlearner(data[, c(X, A)], target = A, library = learners_trt, outcome_type = "binomial", folds = 5)
      mu_model   <- mlr3superlearner::mlr3superlearner(data[, c(X, A, Y)], target = Y, library = learners_outcome, outcome_type = "binomial", folds = 5)
      
      pi_hat  <- predict(pi_model, newdata = data)
      mu0_hat <- predict(mu_model, newdata = data0)
      mu1_hat <- predict(mu_model, newdata = data1)
    }
    mu_hat <- ifelse(data[[A]] == 1, mu1_hat, mu0_hat)
  }

  results <- lapply(smoothness, \(smoothness) {
    # Set up output
    trimmed <- lower <- upper <- numeric(K)
    trimmed_ci <- lower_ci <- upper_ci <- matrix(ncol = 2, nrow = K)
    trimmed_eif <- lower_eif <- upper_eif <- matrix(nrow = N, ncol = K)
    
    #
    # Plugin and one-step estimators
    #
    if(method == "onestep") {
      for(index in seq_along(thresholds)) {
        threshold <- thresholds[index]
        
        trimmed_plugin <- mean(mu1_hat * s_gt(pi_hat, threshold, smoothness) - mu0_hat * s_lt(pi_hat, 1 - threshold, smoothness))
        lower_plugin <- trimmed_plugin - mean(1 - s_lt(pi_hat, 1 - threshold, smoothness))
        upper_plugin <- trimmed_plugin + mean(1 - s_gt(pi_hat, threshold, smoothness))
        
        trimmed_eif[, index] <- eif_trimmed(data[[A]], data[[Y]], mu0_hat, mu1_hat, pi_hat, threshold, smoothness)
        lower_eif[, index]   <- eif_lower(data[[A]], data[[Y]], mu0_hat, mu1_hat, pi_hat, threshold, smoothness)
        upper_eif[, index]   <- eif_upper(data[[A]], data[[Y]], mu0_hat, mu1_hat, pi_hat, threshold, smoothness)
        
        trimmed[index] <- onestep(trimmed_plugin, trimmed_eif[, index])
        lower[index]   <- onestep(lower_plugin, lower_eif[, index])
        upper[index]   <- onestep(upper_plugin, upper_eif[, index])
        
        trimmed_ci[index, ]  <- ci(trimmed[index], trimmed_eif[, index], N)
        lower_ci[index, ]    <- ci(lower[index], lower_eif[, index], N)
        upper_ci[index, ]    <- ci(upper[index], upper_eif[, index], N)
      }
    }
   
    #
    # TMLE
    #
    else if(method == "tmle") {
      for(index in seq_along(thresholds)) {
        threshold <- thresholds[index]
        
        tmle_lower   <- tmle_smooth(data[[A]], data[[Y]], mu0_hat, mu1_hat, pi_hat, threshold, smoothness, parameter = "lower")
        tmle_upper   <- tmle_smooth(data[[A]], data[[Y]], mu0_hat, mu1_hat, pi_hat, threshold, smoothness, parameter = "upper")
        tmle_trimmed <- tmle_smooth(data[[A]], data[[Y]], mu0_hat, mu1_hat, pi_hat, threshold, smoothness, parameter = "trimmed")
        
        trimmed[index] <- tmle_trimmed$psi
        lower[index]   <- tmle_lower$psi
        upper[index]   <- tmle_upper$psi
        
        trimmed_eif[, index] <- tmle_trimmed$eif
        lower_eif[, index]   <- tmle_lower$eif
        upper_eif[, index]   <- tmle_upper$eif
        
        trimmed_ci[index, 1] <- tmle_trimmed$ci[1]
        trimmed_ci[index, 2] <- tmle_trimmed$ci[2]
        lower_ci[index, 1]   <- tmle_lower$ci[1]
        lower_ci[index, 2]   <- tmle_lower$ci[2]
        upper_ci[index, 1]   <- tmle_upper$ci[1]
        upper_ci[index, 2]   <- tmle_upper$ci[2]
      } 
    }
    
    # Multiplier bootstrap
    uniform_ci <- matrix(NA, K, 2)
    if(bootstrap == TRUE) {
      uniform_ci <- multiplier_bootstrap(lower, upper, lower_eif, upper_eif, draws = bootstrap_draws)
    }
    
    list(
      lower = lower_ci[, 1],
      upper = upper_ci[, 2],
      lower_uniform = uniform_ci[, 1],
      upper_uniform = uniform_ci[, 2],
      thresholds = thresholds
    )
  })

  list(
    bounds = results,
    smoothness = smoothness,
    thresholds = thresholds,
    nuisance = list(
      pi_hat = pi_hat,
      mu_hat = mu_hat,
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat
    )
  )
}
