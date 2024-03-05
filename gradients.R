# General gradient function
gradient <- function(param, v, mu, theta, mu1, mu2, alpha, beta) {
  pars <<- list(v=v, mu=mu, theta=theta, mu1=mu1, mu2=mu2, alpha=alpha, beta=beta)
  if (param == "v") { return(gradient_v())
  } else if (param == "mu")    { return(gradient_mu())
  } else if (param == "theta") { return(gradient_theta())
  } else if (param == "mu1")   { return(gradient_mu1())
  } else if (param == "mu2")   { return(gradient_mu2())
  } else if (param == "alpha") { return(gradient_alpha())
  } else if (param == "beta")  { return(gradient_beta())
  } else { stop("Wrong parameter specified for gradient.")
  }
}

# Denominator used in most gradient calculations
den <- function(i) {
  sum(sapply(1:K, function(j) {
    k(t[i], x[i,1], x[i,2], pars$mu[j], pars$theta[j], pars$mu1[j], pars$mu2[j], pars$alpha[j], pars$beta[j]) * 
      exp(2*pars$v[j]) / (1 + exp(pars$v[j]))^3
  }))
}

gradient_v <- function() {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      k(t[i], x[i,1], x[i,2], pars$mu[j], pars$theta[j], pars$mu1[j], pars$mu2[j], pars$alpha[j], pars$beta[j]) * 
        (2*exp(2*pars$v[j]) - exp(3*pars$v[j])) / (1 + exp(pars$v[j])^4) / den(i)
        
    }))
  })
}

gradient_mu <- function() {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
     dkt_dmu(i, j) * k_s(x[i,1], x[i,2], pars$mu1[j], pars$mu2[j], pars$alpha[j], pars$beta[j])  * 
        exp(2*pars$v[j]) / (1 + exp(pars$v[j]))^3 / den(i)
    }))
  })
}

gradient_theta <- function() {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dkt_dtheta(i, j) * k_s(x[i,1], x[i,2], pars$mu1[j], pars$mu2[j], pars$alpha[j], pars$beta[j])  * 
        exp(2*pars$v[j]) / (1 + exp(pars$v[j]))^3 / den(i)
    }))
  })
}

gradient_mu1 <- function() {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_dmu1(i, j) * k_t(x[i,1], pars$mu[j], pars$theta[j])  * 
        exp(2*pars$v[j]) / (1 + exp(pars$v[j]))^3 / den(i)
    }))
  })
}

gradient_mu2 <- function() {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_dmu2(i, j) * k_t(x[i,1], pars$mu[j], pars$theta[j])  * 
        exp(2*pars$v[j]) / (1 + exp(pars$v[j]))^3 / den(i)
    }))
  })
}

gradient_alpha <- function() {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_da(i, j) * k_t(x[i,1], pars$mu[j], pars$theta[j])  * 
        exp(2*pars$v[j]) / (1 + exp(pars$v[j]))^3 / den(i)
    }))
  })
}

gradient_beta <- function() {
  sapply(1:K, function(j) {
    (N / n) * sum(sapply(1:n, function(i) {
      dks_db(i, j) * k_t(x[i,1], pars$mu[j], pars$theta[j])  * 
        exp(2*pars$v[j]) / (1 + exp(pars$v[j]))^3 / den(i)
    }))
  })
}

dkt_dmu <- function(i, j) {
  (sqrt(2 * pi * exp(pars$theta[j])) * t[i] * (1-t[i]))^(-1) * pars$mu[j] * 
    exp(-0.5 * (logit(t[i]) - pars$mu[j])^2 / exp(pars$theta[j]))
}

dkt_dtheta <- function(i, j) {
  (2 * sqrt(2 * pi) * t[i] * (1-t[i]))^(-1) * (1 - exp(-pars$theta[j]) * (logit(t[i]) - pars$mu[j])^2) * 
    exp(0.5 * pars$theta[j] - 0.5 * exp(-pars$theta[j]) * (logit(t[i]) - pars$mu[j])^2)
}

# Jacobian for x1, x2
g <- function(x_i) {
  (x_i[1] * (1 - x_i[1]) * x_i[2] * (1 - x_i[2]))^(-1)
}

dks_dmu1 <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) / (2 * pi) * sqrt(exp(pars$beta[j]) / exp(pars$alpha[j])) * pars$mu1[j] * 
    exp(-0.5 * ( (y_i[1] - pars$mu1[j])^2 / exp(pars$alpha[j]) + (y_i[2] - pars$mu2[j])^2 / exp(pars$beta[j]) ))
}

dks_dmu2 <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) / (2 * pi) * sqrt(exp(pars$alpha[j]) / exp(pars$beta[j])) * pars$mu2[j] * 
    exp(-0.5 * ( (y_i[1] - pars$mu1[j])^2 / exp(pars$alpha[j]) + (y_i[2] - pars$mu2[j])^2 / exp(pars$beta[j]) ))
}

dks_da <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) * sqrt(exp(pars$beta[j])) / (2 * pi) * (0.5 + 0.5 * (y_i[1] - pars$mu1[j])^2 / exp(pars$alpha[j])) * 
    exp(0.5 * pars$alpha[j] - 0.5 * ( (y_i[1] - pars$mu1[j])^2 / exp(pars$alpha[j]) + (y_i[2] - pars$mu2[j])^2 / exp(pars$beta[j]) ))
}

dks_db <- function(i, j) {
  y_i <- logit(x[i,])
  g(x[i,]) * sqrt(exp(pars$alpha[j])) / (2 * pi) * (0.5 + 0.5 * (y_i[2] - pars$mu2[j])^2 / exp(pars$beta[j])) * 
    exp(0.5 * pars$beta[j] - 0.5 * ( (y_i[1] - pars$mu1[j])^2 / exp(pars$alpha[j]) + (y_i[2] - pars$mu2[j])^2 / exp(pars$beta[j]) ))
}

