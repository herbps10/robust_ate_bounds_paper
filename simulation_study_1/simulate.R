#' Simulate data for simulation study
#'
#' @param seed random number seed
#' @param N number of observations
#' @param alpha strength of confounding
#' @param effect_size treatment effect size
#'
#' @return data frame of simulated data. 
#' 
simulate_data <- function(seed = 1, N = 1e3, alpha = 1, beta = 0.05, effect_size = 1, pscore_threshold = Inf) {
  set.seed(seed)
  attempts <- 0
  repeat { 
    X1 <- runif(N, -1, 1)
    X2 <- extraDistr::rcat(N, c(beta / 2, 1 - beta, beta / 2)) - 2
    pscore <- plogis(X1 + alpha * X2)
    A  <- rbinom(N, 1, pscore)

    Y  <- rbinom(N, 1, plogis(X1 + (A - 0.5) * effect_size))

    min_pscore <- min(c(pscore[A == 1], (1 - pscore)[A == 0]))

    if(min_pscore < pscore_threshold) break
    attempts <- attempts + 1
    if(attempts > 1e3) stop("Failed to simulate dataset")
  }
  
  data.frame(X1 = X1, X2 = X2, A = A, Y = Y, pscore = pscore)
}

#' Calculate true ATE under simulation DGP
#'
#' @param effect_size treatment effect size
calculate_ate <- function(effect_size) {
  integrate(\(x) plogis(x + 0.5 * effect_size) * 1/2, -1, 1)$value - integrate(\(x) plogis(x - 0.5 * effect_size) * 1/2, -1, 1)$value
}
