source("kernels.R")
source("gradients.R")
source("priors.R")
source("likelihood.R")
source("langevin.R")
library(mvtnorm)
library(spatstat)
str(lansing)
summary(lansing$x)
summary(lansing$y)

N <- lansing$n
n <- 100
subset <- sample(1:N, n)
x <- cbind(x1 = lansing$x[subset], x2 = lansing$y[subset])
K <- 10
t <- rep(0.5, n)

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
p <- rep(1/K, K)
v <- log(p / (1 - p))

starting <- list(v = v, mu = mu, theta = theta, mu1 = mu1, mu2 = mu2, alpha = alpha, beta = beta)
step_sizes <- c(5e-4, 1e-4, 3e-6, 1e-6, 5e-16, 1e-8, 5e-8)
model <- langevin_pp(x = x, t, N, K, starting, step = step_sizes, nIter = 50, nBurn = 10, nThin = 2)

f_true <- sapply(1:n, function(i) {
  sum(sapply(1:K, function(j) {
    k(t[i], x[i,1], x[i,2],
      mu[j], log(sig.sq[j]),
      mu1[j], mu2[j],
      log(tau.sq1[j]), log(tau.sq2[j])) * p[j]
  }))
})

f_est <- sapply(1:n, function(i) {
  sum(sapply(1:K, function(j) {
    k(t[i], x[i,1], x[i,2],
      model$posteriorMeans$mu[j], log(model$posteriorMeans$sigma.sq[j]),
      model$posteriorMeans$mu1[j], model$posteriorMeans$mu2[j],
      log(model$posteriorMeans$tau.sq1[j]), log(model$posteriorMeans$tau.sq2[j])) * 
      model$posteriorMeans$pi[j]
  }))
})

df_true <- data.frame(x=x[,1], y=x[,2], f=f_true)
df_est <- data.frame(x=x[,1], y=x[,2], f=f_est)
lims <- c(0, 52)

ggplot(data=df_true, aes(x, y, height=0.05, width=0.05)) +
  geom_tile(aes(fill = f)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = lims) +
  theme_classic()

ggplot(data=df_est, aes(x, y, height=0.05, width=0.05)) +
  geom_tile(aes(fill = f)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = lims) +
  theme_classic()

library(interp)
df_true <- interp(x[,1], x[,2], f_true, nx = 500, ny = 500) |> 
  interp2xyz() |> 
  as.data.frame()
ggplot(data = df_true, aes(x, y)) +
  geom_raster(aes(fill = z)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA) + 
  theme_classic() +
  ggtitle("Linear interpolation")