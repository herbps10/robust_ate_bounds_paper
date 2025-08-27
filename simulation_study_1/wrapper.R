wrapper <- function(index, N, alpha, beta, effect_size, smoothness, data, path) {
  # Check cache
  if(!file.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  if(file.exists(path)) {
    return(read_rds(path))
  }

  # Run statistical analysis
  thresholds <- c(0, 10^seq(-4, -1, 0.01))
  bounds <- smooth_bounds(
    data = data, X = c("X1", "X2"), A = "A", Y = "Y", 
    smoothness = smoothness, 
    method = "tmle", 
    bootstrap_draws = 5e3, 
    thresholds = thresholds
  )
  
  shortest_bound <- which.min(bounds$upper - bounds$lower)
  shortest_uniform_bound <- which.min(bounds$upper_uniform - bounds$lower_uniform)
  uniform_test <- any(bounds$upper_uniform < 0) | any(bounds$lower_uniform > 0)
  dr_test <- bounds$upper[1] < 0 | bounds$lower[1] > 0
  test <- any(bounds$upper < 0) | any(bounds$lower_uniform > 0)
  
  max_pscore <- max(c(1 / bounds$nuisance$pi_hat[data$A == 1], 1 / (1 - bounds$nuisance$pi_hat[data$A == 0])))

  true_ate <- calculate_ate(effect_size)
  
  res <- tibble(
    # Save simulation parameters
    index = index,
    N = N,
    alpha = alpha,
    beta = beta,
    effect_size = effect_size,
    smoothness = smoothness,

    # Save estimator results
    true_ate = true_ate,
    max_pscore = max_pscore,
    test = test,
    uniform_test = uniform_test,
    dr_test = dr_test,

    # CI with c=0
    lower_dr = bounds$lower[1],
    upper_dr = bounds$upper[1],
    width_dr = bounds$upper[1] - bounds$lower[1],

    # Pointwise bounds
    lower = bounds$lower[shortest_bound],
    upper = bounds$upper[shortest_bound],
    width = bounds$upper[shortest_bound] - bounds$lower[shortest_bound],
    shortest_threshold = thresholds[shortest_bound],

    # Uniform bounds
    lower_uniform = bounds$lower_uniform[shortest_uniform_bound],
    upper_uniform = bounds$upper_uniform[shortest_uniform_bound],
    width_uniform = bounds$upper_uniform[shortest_uniform_bound] - bounds$lower_uniform[shortest_uniform_bound],
    shortest_uniform_threshold = thresholds[shortest_uniform_bound]
  )

  write_rds(res, path, compress = "gz")

  return(res)
}
