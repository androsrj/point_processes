#ifndef UPDATE_BETA_RETROSPECTIVE_HPP
#define UPDATE_BETA_RETROSPECTIVE_HPP


arma::vec update_beta_retrospective(
    int n1,
    arma::mat x1_locations,
    int n2,
    arma::vec x2,
    arma::mat x2_locations,
    double mu0,
    double s20,
    double t20,
    double gaussian_processes_exponent,
    arma::mat S22i){
  
  arma::vec mu, x1;
  arma::mat term1, term2, dist, dist2, S11, S12, S;

  term1 = pow(repmat(x2_locations.col(0), 1, n1) -
    repmat(x1_locations.col(0).t(), n2, 1), 2);
  term2 = pow(repmat(x2_locations.col(1), 1, n1) -
    repmat(x1_locations.col(1).t(), n2, 1), 2);
  
  dist = arma::sqrt(term1 + term2);
  S12 = cov_function(s20, t20, gaussian_processes_exponent, dist);
  
  dist2 = dmatrix(x1_locations);
  S11 = cov_function(s20, t20, gaussian_processes_exponent, dist2);
  
  mu = arma::ones(n1) * mu0 + S12.t() * S22i * x2;
  S = S11 - S12.t() * S22i * S12;
  S = (S + S.t()) / 2.0;
  
  x1 = rmvnorm(n1, mu, S);
  
  return x1;
  
}


#endif // UPDATE_BETA_RETROSPECTIVE_HPP
