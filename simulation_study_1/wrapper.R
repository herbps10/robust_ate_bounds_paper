library(effectbounds)

wrapper <- function(index, N, alpha, beta, effect_size, pscore_threshold, data, path) {
  # Check cache
  if(!file.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  if(file.exists(path)) {
    cached_res <- read_rds(path)
    if(sum(cached_res$index == index & cached_res$N == N & cached_res$alpha == alpha & cached_res$beta == beta & cached_res$effect_size == effect_size & cached_res$pscore_threshold == pscore_threshold) > 0) return()
  }
  else {
    cached_res <- NULL
  }

  print(glue::glue("Starting: {index} {N} {alpha} {beta} {effect_size} {pscore_threshold}"))

  smoothness <- c(0.001, 0.01, 0.1)

  set.seed(10016)
  # Run statistical analysis
  thresholds <- c(10^seq(-4, log10(0.05), 0.05))
  fit <- ate_bounds(
    data = data, X = c("X1", "X2"), A = "A", Y = "Y", 
    thresholds = thresholds,
    smoothness = smoothness, 
    bootstrap_draws = 1e3, 
  )

  critical_value <- fit$uniform_critical_value
  true_ate <- calculate_ate(effect_size)
  max_pscore <- max(c(1 / fit$nuisance$pi_hat[data$A == 1], 1 / (1 - fit$nuisance$pi_hat[data$A == 0])))

  res <- map_df(seq_along(smoothness), \(sx) {
    bounds <- fit$bounds[[sx]]

    best_lower_uniform_bound <- which.min(bounds$upper_uniform)
    best_upper_uniform_bound <- which.max(bounds$lower_uniform)

    uniform_test <- any(bounds$upper_uniform < 0) | any(bounds$lower_uniform > 0)
    test <- any(bounds$upper < 0) | any(bounds$lower_uniform > 0)

    tibble(
      # Save simulation parameters
      index = index,
      N = N,
      alpha = alpha,
      beta = beta,
      effect_size = effect_size,
      smoothness = smoothness[sx],
      pscore_threshold = pscore_threshold,

      # Save estimator results
      true_ate = true_ate,
      max_pscore = max_pscore,
      test = test,
      uniform_test = uniform_test,
      critical_value = critical_value,

      # One-step 
      lower_onestep = fit$onestep$lower,
      upper_onestep = fit$onestep$upper,
      width_onestep = fit$onestep$upper - fit$onestep$lower,
      test_onestep  = fit$onestep$test,

      # Uniform bounds
      lower_uniform = bounds$lower_uniform[best_lower_uniform_bound],
      upper_uniform = bounds$upper_uniform[best_upper_uniform_bound],
      width_uniform = bounds$upper_uniform[best_upper_uniform_bound] - bounds$lower_uniform[best_lower_uniform_bound]
    )
  })


  if(is.null(cached_res)) {
    write_rds(res, path)
  }
  else {
    write_rds(bind_rows(cached_res, res), path)
  }

  return()
}
