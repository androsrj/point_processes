#ifndef UPDATE_LOG_POSTERIOR_HPP
#define UPDATE_LOG_POSTERIOR_HPP


double update_log_posterior_1(
    int n,
    int m,
    arma::mat BNM0,
    double lambdastar,
    double alpha,
    double beta,
    double range_log_prior,
    double mu0,
    arma::mat C,
    double gp_log_dens,
    double vol,
    bool lambdastar_estimated,
    bool gaussian_processes_estimated,
    bool ranges_estimated,
    std::string link_function){
  // 
  // 
  // Parameters:
  //   d:
  
  arma::vec mu;
  arma::mat BKtilda_l, CKl;
  arma::mat Y = BNM0.submat(0, 0, n - 1, 2);
  arma::mat Ytilda = BNM0.submat(n, 0, BNM0.n_rows - 1, 2);
  
  double augmented_likelihood =
    sum(log_F(link_function, Y.col(2))) + 
    sum(log_F(link_function, - Ytilda.col(2))) + 
    (- vol * lambdastar) + ((n + m) * log(lambdastar));
  
  double res = augmented_likelihood;
  
  if (lambdastar_estimated) {
    // lambdastar prior
    res += (beta * lambdastar) + ((alpha - 1) * log(lambdastar));
  }
  
  if (ranges_estimated && gaussian_processes_estimated) {
    // range prior and gaussian process prior
    res += range_log_prior + gp_log_dens;
    
  } else if (ranges_estimated){
    // range prior
    res += range_log_prior;
    
  } else if (gaussian_processes_estimated){
    // gaussian process prior
    mu = arma::ones(n + m) * mu0;
    res += dmvnorm_to_double(BNM0.col(2), mu, C, true);
    
  }
  
  if (std::isnan(res)) {
    Rcpp::stop("Log posterior is NA!");
  }
  
  return res;
  
}


#endif // UPDATE_LOG_POSTERIOR_HPP
