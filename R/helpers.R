# Approximate I[x > t] as gamma -> 0
# Note that s(x) <= I[x > t]
s_gt <- function(x, t, gamma) ifelse(x - t <= 0, 0, ifelse(x - t >= 0 + gamma, 1, 1 - exp(1 / (((x - t) / gamma)^2 - 1))/ exp(-1)))
s_gt_dot <- function(x, t, gamma) {
  ifelse(x - t <= 0, 0, ifelse(x - t >= 0 + gamma, 0, 
    2 * gamma^2 * exp(1 / ((t - x)^2 / gamma^2 - 1) + 1) * (x - t) / (gamma^2 - (t - x)^2)^2))
}

s_gte <- function(x, t, gamma) ifelse(x - t <= -gamma, 0, ifelse(x - t >= 0, 1, 1 - exp(1 / (((x - t + gamma) / gamma)^2 - 1))/ exp(-1)))
s_gte_dot <- function(x, t, gamma) {
  ifelse(x - t <= 0, 0, ifelse(x - t >= 0 + gamma, 0, 
    2 * gamma^2 * exp(1 / ((t - x + gamma)^2 / gamma^2 - 1) + 1) * (x - t) / (gamma^2 - (t - x)^2)^2))
}

# Approximate I[x < t] as gamma -> 0
# Note that s(x) <= I[x < t]
s_lt <- function(x, t, gamma) 1 - s_gte(x, t, gamma) 
s_lt_dot <- function(x, t, gamma) -s_gte_dot(x, t, gamma)

#' Compute one-step estimtaor
#' 
#' @param plugin plugin point estimate
#' @param eif EIF vector
onestep <- function(plugin, eif) plugin + mean(eif)

ci <- function(estimate, eif, N, alpha = 0.05) {
  estimate + qnorm(c(alpha / 2, 1 - alpha / 2)) * sd(eif) / sqrt(N)
}

multiplier_bootstrap <- function(onestep, eif, draws = 1e3, alpha = 0.05) {
  N <- nrow(eif)
  se <- apply(eif, 2, sd) |> matrix(ncol = ncol(eif), nrow = N, byrow = TRUE)
  mu <- apply(eif, 2, mean) |> matrix(ncol = ncol(eif), nrow = nrow(eif), byrow = TRUE)
  eif2 <- (eif - mu) / se
  zs <- matrix(2 * rbinom(draws * N, 1, 0.5) - 1, nrow = N, ncol = draws)
  
  x <- map_dbl(1:draws, \(draw) max(apply(zs[, draw] * eif2, 2, \(x) abs(sum(x) / sqrt(N))), na.rm = T))
  calpha <- quantile(x, 1 - alpha)
  
  matrix(c(onestep - calpha * apply(eif, 2, sd) / sqrt(N), onestep + calpha * apply(eif, 2, sd) / sqrt(N)), ncol = 2, nrow = ncol(eif), byrow = FALSE)
}
