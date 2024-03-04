source("kernels.R")
source("gradients.R")
source("priors.R")
source("likelihood.R")
source("langevin.R")
library(mvtnorm)

# Dimensions
N <- 500
n <- 100

# Number of basis functions
K <- 5

# Times
t <- seq(0, 1, length = (n+2))[-c(1, n+2)]

# Spatial locations
x1 <- runif(n)
x2 <- runif(n)
x <- matrix(c(x1, x2), ncol = 2)

# Temporal parameters
mu <- rnorm(K)
sig.sq <- rgamma(K, 1, 2)
theta <- log(sig.sq)

# Spatial parameters
mu1 <- rnorm(K)
mu2 <- rnorm(K)
tau.sq1 <- rgamma(K, 1, 2)
tau.sq2 <- rgamma(K, 1, 2)
alpha <- log(tau.sq1)
beta <- log(tau.sq2)

# Weights and transformed weights
p <- c(0.1, 0.2, 0.25, 0.15, 0.3)
v <- log(p / (1 - p))

# Test out likelihood functions
log(lik(t, x, v, mu, theta, mu1, mu2, alpha, beta))
loglik(t, x, v, mu, theta, mu1, mu2, alpha, beta)

# Test out gradient functions
gradient_v(v)
gradient_mu(mu)
gradient_theta(theta)
gradient_mu1(mu1)
gradient_mu2(mu2)
gradient_alpha(alpha)
gradient_beta(beta)

p
mu
sig.sq

starting <- list(v = v, mu = mu, theta = theta, mu1 = mu1, mu2 = mu2, alpha = alpha, beta = beta)
step_sizes <- c(0.0005, 0.0001, 1e-4, 1e-5, 1e-16, 1e-4, 1e-5)
langevin_pp(x, t, N, K, starting, step = step_sizes, nIter = 100, nBurn = 10, nThin = 2)
