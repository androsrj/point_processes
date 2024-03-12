source("functions/kernels.R")
source("functions/gradients.R")
source("functions/priors.R")
source("functions/likelihood.R")
source("functions/langevin.R")
source("functions/parallel.R")

library(mvtnorm)
library(ggplot2)
library(spatstat)
library(parallel)
library(doParallel)
library(foreach)

index_rm <- lansing$x != 1 & lansing$y != 0
lansing$x <- lansing$x[index_rm]
lansing$y <- lansing$y[index_rm]
mySeed <- 123
nCores <- 20
N <- lansing$n - 4
n <- 100

subsets <- vector("list", nCores)
for (i in 1:nCores) {
  set.seed(mySeed)
  subsets[[i]] <- ((i-1)*n+1):(i*n)
}

#x <- cbind(x1 = lansing$x[subset], x2 = lansing$y[subset])
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

#starting <- list(v = v, mu = mu, theta = theta, mu1 = mu1, mu2 = mu2, alpha = alpha, beta = beta)
starting <- list(v = v, mu = mu + rnorm(K, 0, 0.25), 
                 theta = theta + rnorm(K, 0, 0.25), 
                 mu1 = mu1 + rnorm(K, 0, 0.25), mu2 = mu2 + rnorm(K, 0, 0.25), 
                 alpha = alpha + rnorm(K, 0, 0.25), beta = beta + rnorm(K, 0, 0.25))
step_sizes <- c(1e-3, 1e-4, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6)

# Parallel
cl <- makeCluster(nCores)
registerDoParallel(cl)
strt<-Sys.time()
set.seed(mySeed)
obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "spatstat")) %dopar% pp_parallel(i)  
final.time <- Sys.time() - strt 
stopCluster(cl)

str(obj)

wass_pi <- rowMeans(sapply(1:nCores, function(z) obj[[z]]$posteriorMedians$pi))
wass_mu <- rowMeans(sapply(1:nCores, function(z) obj[[z]]$posteriorMedians$mu))
wass_sigma.sq <- rowMeans(sapply(1:nCores, function(z) obj[[z]]$posteriorMedians$sigma.sq))
wass_mu1 <- rowMeans(sapply(1:nCores, function(z) obj[[z]]$posteriorMedians$mu1))
wass_mu2 <- rowMeans(sapply(1:nCores, function(z) obj[[z]]$posteriorMedians$mu2))
wass_tau.sq1 <- rowMeans(sapply(1:nCores, function(z) obj[[z]]$posteriorMedians$tau.sq1))
wass_tau.sq2 <- rowMeans(sapply(1:nCores, function(z) obj[[z]]$posteriorMedians$tau.sq2))

#model <- langevin_pp(x = x, t, N, K, starting, step = step_sizes, nIter = 1000, nBurn = 10, nThin = 2)

x <- cbind(x1 = lansing$x[unlist(subsets)], x2 = lansing$y[unlist(subsets)])
t <- rep(0.5, nrow(x))

f_true <- sapply(1:nrow(x), function(i) {
  sum(sapply(1:K, function(j) {
    k(t[i], x[i,1], x[i,2],
      mu[j], log(sig.sq[j]),
      mu1[j], mu2[j],
      log(tau.sq1[j]), log(tau.sq2[j])) * p[j]
  }))
})

f_est <- sapply(1:nrow(x), function(i) {
  sum(sapply(1:K, function(j) {
    k(t[i], x[i,1], x[i,2],
      wass_mu[j], log(wass_sigma.sq[j]),
      wass_mu1[j], wass_mu2[j],
      log(wass_tau.sq1[j]), log(wass_tau.sq2[j])) * 
      wass_pi[j]
  }))
})



df_true <- data.frame(x=x[,1], y=x[,2], f=f_true)
df_est <- data.frame(x=x[,1], y=x[,2], f=f_est)
lims <- c(0, max(f_est))

ggplot(data=df_true, aes(x, y, height=0.05, width=0.05)) +
  geom_tile(aes(fill = f)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = lims) +
  theme_classic()

ggplot(data=df_est, aes(x, y, height=0.05, width=0.05)) +
  geom_tile(aes(fill = f)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = lims) +
  theme_classic()

f_true

library(interp)
df_true <- interp(x[,1], x[,2], f_true, nx = 50, ny = 50) |> 
  interp2xyz() |> 
  as.data.frame()
ggplot(data = df_true, aes(x, y)) +
  geom_raster(aes(fill = z)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = lims) + 
  theme_classic() +
  ggtitle("True f")
ggsave(filename = "figures/surf_truth.pdf", height = 5)

df_est <- interp(x[,1], x[,2], f_est, nx = 50, ny = 50) |> 
  interp2xyz() |> 
  as.data.frame()
ggplot(data = df_est, aes(x, y)) +
  geom_raster(aes(fill = z)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = lims) + 
  theme_classic() +
  ggtitle("Estimated f")
ggsave(filename = "figures/surf_est.pdf", height = 5)

model$acceptance
starting
mu[1]
sig.sq[1]
mu1[1]
mu2[1]
tau.sq1[1]
tau.sq2[1]
