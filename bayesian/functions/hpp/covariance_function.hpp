#ifndef COVARIANCE_FUNCTION_HPP
#define COVARIANCE_FUNCTION_HPP


arma::mat cov_function(
    double s20,
    double t20,
    double exponent,
    arma::mat dist){
  
  arma::mat S = s20 * exp(-(1 / (2 * t20)) * pow(dist, exponent));
  
  return S;
  
}


#endif // COVARIANCE_FUNCTION_HPP
