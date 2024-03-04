# Denominator used in most gradient calculations
den <- function(i) {
  sum(sapply(1:K, function(j) {
    k(t[i], x[i,1], x[i,2], mu[j], theta[j], mu1[j], mu2[j], alpha[j], beta[j]) * 
      exp(2*v[j]) / (1 + exp(v[j]))^3
  }))
}

gradient_v <- function(v) {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      k(t[i], x[i,1], x[i,2], mu[j], theta[j], mu1[j], mu2[j], alpha[j], beta[j]) * 
        (2*exp(2*v[j]) - exp(3*v[j])) / (1 + exp(v[j])^4) / den(i)
        
    }))
  })
}

gradient_mu <- function(mu) {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
     dkt_dmu(i, j) * k_s(x[i,1], x[i,2], mu1[j], mu2[j], alpha[j], beta[j])  * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 / den(i)
    }))
  })
}

gradient_theta <- function(theta) {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dkt_dtheta(i, j) * k_s(x[i,1], x[i,2], mu1[j], mu2[j], alpha[j], beta[j])  * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 / den(i)
    }))
  })
}

gradient_mu1 <- function(mu1) {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_dmu1(i, j) * k_t(x[i,1], mu[j], theta[j])  * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 / den(i)
    }))
  })
}

gradient_mu2 <- function(mu2) {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_dmu2(i, j) * k_t(x[i,1], mu[j], theta[j])  * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 / den(i)
    }))
  })
}

gradient_alpha <- function(alpha) {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_da(i, j) * k_t(x[i,1], mu[j], theta[j])  * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 / den(i)
    }))
  })
}

gradient_beta <- function(beta) {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_db(i, j) * k_t(x[i,1], mu[j], theta[j])  * 
        exp(2*v[j]) / (1 + exp(v[j]))^3 / den(i)
    }))
  })
}

dkt_dmu <- function(i, j) {
  (sqrt(2 * pi * exp(theta[j])) * t[i] * (1-t[i]))^(-1) * mu[j] * 
    exp(-0.5 * (logit(t[i]) - mu[j])^2 / exp(theta[j]))
}

dkt_dtheta <- function(i, j) {
  (2 * sqrt(2 * pi) * t[i] * (1-t[i]))^(-1) * (1 - exp(-theta[j]) * (logit(t[i]) - mu[j])^2) * 
    exp(0.5 * theta[j] - 0.5 * exp(-theta[j]) * (logit(t[i]) - mu[j])^2)
}

# Jacobian for x1, x2
g <- function(x_i) {
  (x_i[1] * (1 - x_i[1]) * x_i[2] * (1 - x_i[2]))^(-1)
}

dks_dmu1 <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) / (2 * pi) * sqrt(exp(beta[j]) / exp(alpha[j])) * mu1[j] * 
    exp(-0.5 * ( (y_i[1] - mu1[j])^2 / exp(alpha[j]) + (y_i[2] - mu2[j])^2 / exp(beta[j]) ))
}

dks_dmu2 <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) / (2 * pi) * sqrt(exp(alpha[j]) / exp(beta[j])) * mu2[j] * 
    exp(-0.5 * ( (y_i[1] - mu1[j])^2 / exp(alpha[j]) + (y_i[2] - mu2[j])^2 / exp(beta[j]) ))
}

dks_da <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) * sqrt(exp(beta[j])) / (2 * pi) * (0.5 + 0.5 * (y_i[1] - mu1[j])^2 / exp(alpha[j])) * 
    exp(0.5 * alpha[j] - 0.5 * ( (y_i[1] - mu1[j])^2 / exp(alpha[j]) + (y_i[2] - mu2[j])^2 / exp(beta[j]) ))
}

dks_db <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) * sqrt(exp(alpha[j])) / (2 * pi) * (0.5 + 0.5 * (y_i[2] - mu2[j])^2 / exp(beta[j])) * 
    exp(0.5 * beta[j] - 0.5 * ( (y_i[1] - mu1[j])^2 / exp(alpha[j]) + (y_i[2] - mu2[j])^2 / exp(beta[j]) ))
}

