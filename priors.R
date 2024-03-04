# mu, mu1, and mu2 (normal)
logPriorMu <- function(mu, mu0 = 0, sig0 = 1) {
  sum(dnorm(mu, mu0, sig0, log = TRUE))
}

# sigma.sq, tau.sq1, and tau.sq2 (inverse gamma)
logPriorSigma2 <- function(sigma2, a = 1, b = 1) {
  sum(a * log(b) - lgamma(a) - (a + 1) * log(sigma2) - b / sigma2)
}

# pi (beta)
logPriorPi <- function(p, a = 0.5, b = 0.5) {
  sum(dbeta(p, a, b, log = TRUE))
}

# mu1


# mu2


# tau.sq1


# tau.sq2


