#ifndef LINK_FUNCTIONS_HPP
#define LINK_FUNCTIONS_HPP


arma::vec F(
    std::string link_function,
    const arma::vec& vector) {
  // 
  // 
  // Parameters:
  //   d:
  
  int n = vector.n_elem;
  arma::vec cdf(n);
  
  if (link_function == "probit") {
    
    for (int i = 0; i < n; ++i) {
      // Compute the CDF for each element of BKst
      cdf(i) = R::pnorm(vector(i), 0.0, 1.0, true, false);
    }
    
  } else if (link_function == "sigmoid") {
    
    cdf = 1.0 / (1.0 + exp(- vector));
    
  }
  
  return cdf;
  
}


double F_scalar(
    std::string link_function,
    const double& scalar) {
  // 
  // 
  // Parameters:
  //   d:
  
  double cdf;
  
  if (link_function == "probit") {
    
    cdf = R::pnorm(scalar, 0.0, 1.0, true, false);
    
  } else if (link_function == "sigmoid") {
    
    cdf = 1.0 / (1.0 + exp(- scalar));
    
  }
  
  return cdf;
  
}


arma::vec log_F(
    std::string link_function,
    const arma::vec& vector){
  // 
  // 
  // Parameters:
  //   d:
  
  int n = vector.n_elem;
  arma::vec cdf(n);
  
  if (link_function == "probit") {
    
    for (int i = 0; i < n; ++i) {
      // Compute the CDF for each element of BKst
      cdf(i) = R::pnorm(vector(i), 0.0, 1.0, true, true);
    }
    
  } else if (link_function == "sigmoid") {
    
    cdf = - log(1.0 + exp(- vector));
    
  }
  
  return cdf;
  
}


#endif // LINK_FUNCTIONS_HPP
