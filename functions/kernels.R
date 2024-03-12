logit <- function(t) {
  log(t / (1 - t))
}

invLogit <- function(v) {
  exp(v) / (1 + exp(v))
}

k_t <- function(t, mu, theta) {
  dnorm(logit(t), mu, sqrt(exp(theta))) / (t * (1 - t))
}

k_s <- function(x1, x2, mu1, mu2, alpha, beta) {
  cov_matrix <- diag(c(exp(alpha), exp(beta)))
  dmvnorm(c(logit(x1), logit(x2)), c(mu1, mu2), cov_matrix) / (x1 * (1 - x1) * x2 * (1 - x2))
}

k <- function(t, x1, x2, mu, theta, mu1, mu2, alpha, beta) {
  k_t(t, mu, theta) * k_s(x1, x2, mu1, mu2, alpha, beta)
}