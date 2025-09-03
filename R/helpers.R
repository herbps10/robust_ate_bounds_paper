# Approximate I[x > t] as gamma -> 0
# Note that s(x) <= I[x > t]
s_gt <- function(x, t, gamma) ifelse(x - t <= 0, 0, ifelse(x - t >= 0 + gamma, 1, 1 - exp(1 / (((x - t) / gamma)^2 - 1))/ exp(-1)))
s_gt_dot <- function(x, t, gamma) {
  ifelse(x - t <= 0, 0, ifelse(x - t >= 0 + gamma, 0, 
    2 * gamma^2 * exp(1 / ((t - x)^2 / gamma^2 - 1) + 1) * (x - t) / (gamma^2 - (t - x)^2)^2))
}

# Approximate I[x < t] as gamma -> 0
# Note that s(x) <= I[x < t]
s_lt <- function(x, t, gamma) ifelse(x <= t - gamma, 1, ifelse(x >= t, 0,  1 - exp(1 / (((x - t) / gamma)^2 - 1))/ exp(-1)))
s_lt_dot <- function(x, t, gamma) {
  ifelse(x <= t - gamma, 0, ifelse(x >= t, 0, 
                               2 * gamma^2 * exp(1 / ((t - x)^2 / gamma^2 - 1) + 1) * (x - t) / (gamma^2 - (t - x)^2)^2))
}

#' Compute one-step estimtaor
#' 
#' @param plugin plugin point estimate
#' @param eif EIF vector
onestep <- function(plugin, eif) plugin + mean(eif)

ci <- function(estimate, eif, N, alpha = 0.05) {
  estimate + qnorm(c(alpha / 2, 1 - alpha / 2)) * sd(eif) / sqrt(N)
}

multiplier_bootstrap <- function(lower, upper, eif_lower, eif_upper, draws = 1e3, alpha = 0.05) {
  N <- nrow(eif_lower)
  K <- ncol(eif_lower)

  se_lower <- apply(eif_lower, 2, sd) |> matrix(ncol = K, nrow = N, byrow = TRUE)
  se_upper <- apply(eif_upper, 2, sd) |> matrix(ncol = K, nrow = N, byrow = TRUE)
  eif_lower_scaled <- eif_lower / se_lower
  eif_upper_scaled <- eif_upper / se_upper

  zs <- matrix(2 * rbinom(draws * N, 1, 0.5) - 1, nrow = N, ncol = draws)
  
  T_lower <- map(1:draws, \(draw) colSums(zs[, draw] * eif_lower_scaled) / sqrt(N))
  T_upper <- map(1:draws, \(draw) colSums(zs[, draw] * eif_upper_scaled) / sqrt(N))
  
  T_max <- map(1:draws, \(draw) pmax(T_lower[[draw]], -T_upper[[draw]]))
  
  max_max <- map_dbl(T_max, max, na.rm = TRUE)
  
  calpha <- quantile(max_max, 1 - alpha)
  print(calpha)

  matrix(c(
    lower - calpha * apply(eif_lower, 2, sd) / sqrt(N), 
    upper + calpha * apply(eif_upper, 2, sd) / sqrt(N)
  ), ncol = 2, nrow = K, byrow = FALSE)
}
