#ifndef UPDATE_BETA_L_HPP
#define UPDATE_BETA_L_HPP


arma::vec update_beta_l__probit(
    int n,
    int m,
    arma::mat cvm0,
    double mu0) {
  // 
  // 
  // Parameters:
  //   d:
  
  arma::mat BNM0;
  
  if (m == 0) {
    
    arma::mat WN = arma::ones(n, 1);
    arma::mat WW = WN;
    WW = Wmatrix(WW);
    
    int K = n;
    
    arma::mat G = arma::eye(K, K) + WW * cvm0 * WW.t();
    G = (G + G.t()) / 2.0;
    arma::mat Gi = inv_sympd_2(G);
    arma::mat Gc = arma::chol(G, "lower");
    
    BNM0 = rsun(arma::zeros(K, 1) + mu0, cvm0, WW, Gi, Gc);
    
  } else {
    
    arma::mat WW = join_cols(arma::ones(n, 1), -arma::ones(m, 1));
    WW = Wmatrix(WW);
    int K = n + m;
    
    arma::mat G = arma::eye(K, K) + WW * cvm0 * WW;
    G = (G + G.t()) / 2.0;
    arma::mat Gi = inv_sympd_2(G);
    arma::mat Gc = arma::chol(G, "lower");
    
    BNM0 = rsun(arma::zeros(K, 1) + mu0, cvm0, WW, Gi, Gc);
    
  }
  
  return BNM0;
  
}


arma::vec update_beta_l__sigmoid(
    int n,
    int m,
    arma::mat B_1,
    double b,
    arma::vec beta_l) {
  // 
  // 
  // Parameters:
  //   d:
  // ,
  // 
  
  arma::vec BNM0(n + m);
  
  arma::vec alpha(n + m);
  alpha = rpg_vector(1, beta_l);
  arma::mat A = arma::diagmat(alpha);
  
  arma::vec y = arma::ones(n);
  arma::vec ytilda = arma::zeros(m);
  arma::vec prob = 0.5 * arma::ones(n + m);
  
  arma::vec tau = (arma::join_vert(y, ytilda) - prob);
  
  arma::mat V_alpha = 
    arma::inv(A + B_1);
  
  arma::vec mu_alpha = 
    V_alpha * (tau + B_1 * (arma::ones(n + m) * b) );
  
  BNM0 = rmvnorm(n + m, mu_alpha, V_alpha);
  
  return BNM0;
  
}


arma::vec update_beta_l(
    int n,
    int m,
    arma::mat CK,
    double mu0,
    arma::vec beta_l,
    arma::mat CKi,
    std::string link_function) {
  // 
  // 
  // Parameters:
  //   d:
  
  arma::vec BNM0;
  
  if (link_function == "probit") {
    
    BNM0 = update_beta_l__probit(n, m, CK, mu0);
    
  } else if (link_function == "sigmoid") {
    
    BNM0 = update_beta_l__sigmoid(
      n,
      m,
      CKi,
      mu0,
      beta_l);
    
  }
  
  return BNM0;
  
}


#endif // UPDATE_BETA_L_HPP
