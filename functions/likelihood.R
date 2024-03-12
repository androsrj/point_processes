# Likelihood
lik <- function(t, x, v, mu, theta, mu1, mu2, alpha, beta) {
  prod(sapply(1:n, function(i) {
    sum(sapply(1:K, function(j) {
      k(t[i], x[i,1], x[i,2], mu[j], theta[j], mu1[j], mu2[j], alpha[j], beta[j]) * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 
    }))
  })^(N / n))
}

# Log likelihood
loglik <- function(t, x, v, mu, theta, mu1, mu2, alpha, beta) {
  sum(log(sapply(1:n, function(i) {
    sum(sapply(1:K, function(j) {
      k(t[i], x[i,1], x[i,2], mu[j], theta[j], mu1[j], mu2[j], alpha[j], beta[j]) * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 
    }))
  }))) * (N / n) 
}
