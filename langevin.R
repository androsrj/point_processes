### LANGEVIN-ADJUSTED METROPOLIS FUNCTION FOR POINT PROCESS DATA ###
# Priors, Jacobians, likelihood, etc. must already be sourced

langevin_pp <- function(x, t, N, K,
                        starting, step,
                        nIter = 1000, nBurn = 100, nThin = 2) {
  
  # Dimensions
  n <<- nrow(x)
  d <<- ncol(x)
  
  # MCMC chain properties
  nIter <- nBurn + nIter
  
  # Initialize parameter vectors and starting values
  v <- mu <- theta <- mu1 <- mu2 <- alpha <- beta <- matrix(0, nrow = nIter, ncol = K) 
  v[1,] <- starting[[1]]
  mu[1,] <- starting[[2]]
  theta[1,] <- starting[[3]]
  mu1[1,] <- starting[[4]]
  mu2[1,] <- starting[[5]]
  alpha[1,] <- starting[[6]]
  beta[1,] <- starting[[7]]
  
  # Track acceptance rates
  acc <- rep(list(rep(0, K)), 7)
  names(acc) <- c("v", "mu", "theta", "mu1", "mu2", "alpha", "beta")
  
  # Run Langevin for one chain
  for (i in 2:nIter) {
    
    cat(paste0("Beginning iteration ", i, ".\n"))
    
    ### Update v ###
    grad <- gradient(param = "v", v[i-1,], mu[i-1,], theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1,])
    v_prop <- v[i-1, ] + step[1] * grad + sqrt(2*step[1]) * rmvnorm(1, sigma = diag(K))
    MHratio <- loglik(t, x, v_prop, mu[i-1,], theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1]) - 
      loglik(t, x, v[i-1,], mu[i-1,], theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1]) + 
      logPriorPi(invLogit(v_prop)) - logPriorPi(invLogit(v[i-1,])) + 
      sum(log(invLogit(v_prop))) - sum(log(invLogit(v[i-1,]))) # Jacobians
    if(runif(1) < exp(MHratio)) {
      v[i,] <- v_prop
      acc$v <- acc$v + rep(1, K)
    } else {
      v[i,] <- v[i -1,]
    }
    
    ### Update mu ###
    grad <- gradient(param = "mu", v[i,], mu[i-1,], theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1,])
    mu_prop <- mu[i-1, ] + step[2] * grad + sqrt(2*step[2]) * rmvnorm(1, sigma = diag(K))
    MHratio <- loglik(t, x, v[i,], mu_prop, theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1]) - 
      loglik(t, x, v[i,], mu[i-1,], theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1]) + 
      logPriorMu(mu_prop) - logPriorMu(mu[i-1,])
    if(runif(1) < exp(MHratio)) {
      mu[i,] <- mu_prop
      acc$mu <- acc$mu + rep(1, K)
    } else {
      mu[i,] <- mu[i -1,]
    }
    
    ### Update theta ###
    grad <- gradient(param = "theta", v[i,], mu[i,], theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1,])
    theta_prop <- theta[i-1, ] + step[3] * grad + sqrt(2*step[3]) * rmvnorm(1, sigma = diag(K))
    MHratio <- loglik(t, x, v[i,], mu[i,], theta_prop, mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1]) - 
      loglik(t, x, v[i,], mu[i,], theta[i-1,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1]) + 
      logPriorSigma2(exp(theta_prop)) - logPriorSigma2(exp(theta[i-1,])) +
      sum(theta_prop) - sum(theta[i-1,]) # Log of Jacobian cancels out
    if(runif(1) < exp(MHratio)) {
      theta[i,] <- theta_prop
      acc$theta <- acc$theta + rep(1, K)
    } else {
      theta[i,] <- theta[i -1,]
    }
    
    ### Update mu1 ###
    grad <- gradient(param = "mu1", v[i,], mu[i,], theta[i,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1,])
    mu1_prop <- mu1[i-1, ] + step[4] * grad + sqrt(2*step[4]) * rmvnorm(1, sigma = diag(K))
    MHratio <- loglik(t, x, v[i,], mu[i,], theta[i,], mu1_prop, mu2[i-1,], alpha[i-1,], beta[i-1]) - 
      loglik(t, x, v[i,], mu[i,], theta[i,], mu1[i-1,], mu2[i-1,], alpha[i-1,], beta[i-1]) + 
      logPriorMu(mu1_prop) - logPriorMu(mu1[i-1,])
    if(runif(1) < exp(MHratio)) {
      mu1[i,] <- mu1_prop
      acc$mu1 <- acc$mu1 + rep(1, K)
    } else {
      mu1[i,] <- mu1[i -1,]
    }
    
    ### Update mu2 ###
    grad <- gradient(param = "mu2", v[i,], mu[i,], theta[i,], mu1[i,], mu2[i-1,], alpha[i-1,], beta[i-1,])
    mu2_prop <- mu2[i-1, ] + step[5] * grad + sqrt(2*step[5]) * rmvnorm(1, sigma = diag(K))
    MHratio <- loglik(t, x, v[i,], mu[i,], theta[i,], mu1[i,], mu2_prop, alpha[i-1,], beta[i-1]) - 
      loglik(t, x, v[i,], mu[i,], theta[i,], mu[i,], mu2[i-1,], alpha[i-1,], beta[i-1]) + 
      logPriorMu(mu2_prop) - logPriorMu(mu2[i-1,])
    if(runif(1) < exp(MHratio)) {
      mu2[i,] <- mu2_prop
      acc$mu2 <- acc$mu2 + rep(1, K)
    } else {
      mu2[i,] <- mu2[i-1,]
    }
    
    ### Update alpha ###
    grad <- gradient(param = "alpha", v[i,], mu[i,], theta[i,], mu1[i,], mu2[i,], alpha[i-1,], beta[i-1,])
    alpha_prop <- alpha[i-1, ] + step[6] * grad + sqrt(2*step[6]) * rmvnorm(1, sigma = diag(K))
    MHratio <- loglik(t, x, v[i,], mu[i,], theta[i,], mu1[i,], mu2[i,], alpha_prop, beta[i-1]) - 
      loglik(t, x, v[i,], mu[i,], theta[i,], mu1[i,], mu2[i,], alpha[i-1,], beta[i-1]) + 
      logPriorSigma2(exp(alpha_prop)) - logPriorSigma2(exp(alpha[i-1,])) + 
      sum(alpha_prop) - sum(alpha[i-1,]) # Log of Jacobian cancels out
    if(runif(1) < exp(MHratio)) {
      alpha[i,] <- alpha_prop
      acc$alpha <- acc$alpha + rep(1, K)
    } else {
      alpha[i,] <- alpha[i -1,]
    }
    
    ### Update beta ###
    grad <- gradient(param = "beta", v[i,], mu[i,], theta[i,], mu1[i,], mu2[i,], alpha[i,], beta[i-1,])
    beta_prop <- beta[i-1, ] + step[7] * grad + sqrt(2*step[7]) * rmvnorm(1, sigma = diag(K))
    MHratio <- loglik(t, x, v[i,], mu[i,], theta[i,], mu1[i,], mu2[i,], alpha[i,], beta_prop) - 
      loglik(t, x, v[i,], mu[i,], theta[i,], mu1[i,], mu2[i,], alpha[i,], beta[i-1]) + 
      logPriorSigma2(exp(beta_prop)) - logPriorSigma2(exp(beta[i-1,])) + 
      sum(beta_prop) - sum(beta[i-1,]) # Log of Jacobian cancels out
    if(runif(1) < exp(MHratio)) {
      beta[i,] <- beta_prop
      acc$beta <- acc$beta + rep(1, K)
    } else {
      beta[i,] <- beta[i -1,]
    }

  }
  
  # Acceptance rates (for Metropolis-sampled parameters)
  acceptance <- sapply(1:7, \(i) acc[[i]] / nIter)
  colnames(acceptance) <- names(acc)
  
  # Remove burn-in and perform thinning
  index <- seq(nBurn + 1, nIter, by = nThin)
  nSamples <- length(index)
  v <- v[index, ]
  mu <- mu[index, ]
  theta <- theta[index, ]
  mu1 <- mu1[index, ]
  mu2 <- mu2[index, ]
  alpha <- alpha[index, ]
  beta <- beta[index, ]
  
  # Back-transform
  pi <- invLogit(v)
  sigma.sq <- exp(theta)
  tau.sq1 <- exp(alpha)
  tau.sq2 <- exp(beta)
  
  # Trace plots
  pdf("figures/trace_plots.pdf")
  plot(1:nSamples, sigma.sq[,1], type = 'l', ylab = "Sigma2", main = "")
  plot(1:nSamples, mu[,1], type = 'l', ylab = "Mu", main = "")
  plot(1:nSamples, mu1[,1], type = 'l', ylab = "Mu1")
  plot(1:nSamples, mu2[,1], type = 'l', ylab = "Mu2", main = "")
  plot(1:nSamples, tau.sq1[,1], type = 'l', ylab = "Tau2", main = "")
  plot(1:nSamples, tau.sq2[,1], type = 'l', ylab = "Tau2", main = "")
  dev.off()  
  
  # Posterior mean estimates (can be somewhat skewed because of back-transformations)
  posteriorMeans <- list(pi = apply(pi, 2, mean),
                         mu = apply(mu, 2, mean),
                         sigma.sq = apply(sigma.sq, 2, mean),
                         mu1 = apply(mu1, 2, mean),
                         mu2 = apply(mu2, 2, mean),
                         tau.sq1 = apply(tau.sq1, 2, mean),
                         tau.sq2 = apply(tau.sq2, 2, mean))
  
  # Posterior median estimates (more accurate)
  posteriorMedians <- list(pi = apply(pi, 2, median),
                           mu = apply(mu, 2, median),
                           sigma.sq = apply(sigma.sq, 2, median),
                           mu1 = apply(mu1, 2, median),
                           mu2 = apply(mu2, 2, median),
                           tau.sq1 = apply(tau.sq1, 2, median),
                           tau.sq2 = apply(tau.sq2, 2, median))
  
  # 95% credible interval bounds
  #credLower <- list(sigma2 = quantile(sigma2, 0.025), 
  #                  tau2 = quantile(tau2, 0.025),
  #                  beta = apply(beta, 1, quantile, 0.025))
  #credUpper <- list(sigma2 = quantile(sigma2, 0.975), 
  #                  tau2 = quantile(tau2, 0.975),
  #                  beta = apply(beta, 1, quantile, 0.975))
  
  # Posterior predictive results for test data
  #preds <- lapply(1:nTestSubj, function(j) {
  #  apply(YPreds[[j]], 1, quantile, c(0.025, 0.5, 0.975))
  #})
  
  # Return results
  return(list(acceptance = acceptance, 
              posteriorMeans = posteriorMeans,
              posteriorMedians = posteriorMedians))
}


