#' Simulate data for simulation study
#'
#' @param seed random number seed
#' @param N number of observations
#' @param alpha strength of confounding
#' @param effect_size treatment effect size
#'
#' @return data frame of simulated data. 
#' 
simulate_data <- function(seed = 1, N = 1e3, gamma0 = 3, effect_size = 1, pscore_threshold = Inf) {
  set.seed(seed)
  attempts <- 0
  a <- b <- gamma0 - 1

  repeat { 
    U <- rnorm(N)

    pscore <- qbeta(pnorm(U), a, b)
    A  <- rbinom(N, 1, pscore)

    mu <- plogis(0.5 * U + effect_size * (A - 0.5))
    Y  <- rbinom(N, 1, mu)

    X1 <- U + rnorm(N, 0, 0.25)
    X2 <- U^2 + rnorm(N, 0, 0.25)
    X3 <- rnorm(N, 0, 1)

    min_pscore <- min(c(pscore[A == 1], (1 - pscore)[A == 0]))

    if(min_pscore < pscore_threshold) break
    attempts <- attempts + 1
    if(attempts > 1e3) stop("Failed to simulate dataset")
  }
  
  data.frame(X1 = X1, X2 = X2, X3, A = A, Y = Y, pscore = pscore)
}

#' Calculate true ATE under simulation DGP
#'
#' @param effect_size treatment effect size
calculate_ate <- function(effect_size) {
  EY0 <- integrate(\(x) plogis(0.5 * x - 0.5 * effect_size) * dnorm(x), -Inf, Inf)$value
  EY1 <- integrate(\(x) plogis(0.5 * x + 0.5 * effect_size) * dnorm(x), -Inf, Inf)$value 

  EY1 - EY0
}
