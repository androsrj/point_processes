#ifndef GET_LINEAR_PREDICTOR_HPP
#define GET_LINEAR_PREDICTOR_HPP


arma::mat get_linear_predictor(
    arma::mat B_aux,
    arma::vec lambdastar,
    arma::mat xy_vals,
    std::string link_function) {
  
  arma::mat B = B_aux.rows(find(B_aux.col(4) == 1));
  
  int n = xy_vals.n_rows;
  arma::mat betas(n, 3);
  int count = 0;
  arma::vec regions(n);
  
  for (unsigned int i = 0; i < n; ++i) {
    double x = xy_vals(i, 0);
    double y = xy_vals(i, 1);
    
    for (unsigned int j = 0; j < B.n_rows; ++j) {
      if (B(j, 0) == x && B(j, 1) == y) {
        betas(i, 0) = B(j, 2);
        regions(i) = B(j, 3);
        count += 1;
      }
    }
  }
  
  betas.col(1) = subsetVector(lambdastar, regions) % F(link_function, betas.col(0));
  betas.col(2) = regions;
  
  if (count != n) {
    Rcpp::stop("Some additional location is not at BKtilda!");
  }
  
  return betas;
  
}

arma::mat get_linear_predictor_1(
    arma::mat B,
    double lambdastar,
    arma::mat xy_vals,
    std::string link_function) {
  
  int n = xy_vals.n_rows;
  arma::mat betas(n, 2);
  int count = 0;
  
  for (unsigned int i = 0; i < n; ++i) {
    double x = xy_vals(i, 0);
    double y = xy_vals(i, 1);
    
    for (unsigned int j = 0; j < B.n_rows; ++j) {
      if (B(j, 0) == x && B(j, 1) == y) {
        betas(i, 0) = B(j, 2);
        count += 1;
      }
    }
  }
  
  betas.col(1) = lambdastar * F(link_function, betas.col(0));
  
  if (count != n) {
    Rcpp::stop("Some additional location is not in Y!");
  }
  
  return betas;
  
}


#endif // GET_LINEAR_PREDICTOR_HPP
