eif_ate <- function(X, A, Y, mu0, mu1, pi) {
  mu <- ifelse(A == 1, mu1, mu0)
  mu1 - mu0 - mean(mu1 - mu0) + (A / pi - (1 - A) / (1 - pi)) * (Y - mu)
}

eif_trimmed <- function(A, Y, mu0, mu1, pi, threshold, smoothness) {
  w1 <- s_gt(pi, threshold, smoothness)
  w1_dot <- s_gt_dot(pi, threshold, smoothness)
  w0 <- s_lt(pi, 1 - threshold, smoothness)
  w0_dot <- s_lt_dot(pi, 1 - threshold, smoothness)
  
  A / pi * w1 * (Y - mu1) - (1 - A) / (1 - pi) * w0 * (Y - mu0) +
    (mu1 * w1_dot - mu0 * w0_dot) * (A - pi) +
    mu1 * w1 - mu0 * w0 - mean(mu1 * w1 - mu0 * w0)
}

eif_lower <- function(A, Y, mu0, mu1, pi, threshold, smoothness) {
  eif_trimmed(A, Y, mu0, mu1, pi, threshold, smoothness) + 
    s_lt_dot(pi, 1 - threshold, smoothness) * (A - pi) + s_lt(pi, 1 - threshold, smoothness) - mean(s_lt(pi, 1 - threshold, smoothness))
}

eif_upper <- function(A, Y, mu0, mu1, pi, threshold, smoothness) {
  eif_trimmed(A, Y, mu0, mu1, pi, threshold, smoothness) - 
    s_gt_dot(pi, threshold, smoothness) * (A - pi) - s_gt(pi, threshold, smoothness) + mean(s_gt(pi, threshold, smoothness))
}
