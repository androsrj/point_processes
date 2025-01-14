#ifndef UPDATE_LATENT_PROCESSES_HPP
#define UPDATE_LATENT_PROCESSES_HPP


List update_latent_processes_1(
    int n,
    arma::mat BNM0,
    double lambdastr,
    arma::mat C11i,
    arma::mat C11,
    arma::mat Ci,
    double mu0,
    double s20,
    double t20,
    double gaussian_processes_exponent,
    double H,
    double W,
    arma::mat dist_Y,
    std::string link_function) {
  
  arma::mat term1, term2, dist, dist2, C12, C22, C;
  
  int K = BNM0.n_rows;
  int Kst = R::rpois(lambdastr * H * W);
  
  // sampling the process X
  arma::mat sm(Kst, 2);
  sm.col(0) = runif_cpp(Kst, 0, H);
  sm.col(1) = runif_cpp(Kst, 0, W);
  
  arma::vec BKst = update_beta_retrospective(
    Kst, sm, K, BNM0.col(2), BNM0, mu0, s20, t20,
    gaussian_processes_exponent, Ci);
  
  // thinned process
  arma::vec unretained = rbinom_cpp(Kst, F(link_function, - BKst));
  arma::uvec ind = arma::find(unretained);
  arma::mat Ytilda = sm.rows(ind);
  int M = Ytilda.n_rows;
  Ytilda = arma::join_rows(Ytilda, BKst(ind));
  
  arma::mat Y = BNM0.submat(0, 0, n - 1, BNM0.n_cols - 1);
  term1 = pow(repmat(Y.col(0), 1, M) - repmat(Ytilda.col(0).t(), n, 1), 2);
  term2 = pow(repmat(Y.col(1), 1, M) - repmat(Ytilda.col(1).t(), n, 1), 2);
  dist = arma::sqrt(term1 + term2);
  C12 = cov_function(s20, t20, gaussian_processes_exponent, dist);
  
  dist2 = dmatrix(Ytilda.cols(0, 1));
  C22 = cov_function(s20, t20, gaussian_processes_exponent, dist2);;
  
  C = join_vert(join_horiz(C11, C12),
                join_horiz(C12.t(), C22));
  C = (C + C.t()) / 2.0;
  
  arma::mat SS_aux = C22 - C12.t() * C11i * C12;
  SS_aux = (SS_aux + SS_aux.t()) / 2.0;
  arma::mat SS = inv_sympd_2(SS_aux);
  arma::mat bd = C11i * C12;
  arma::mat CC12 = -bd * SS;
  arma::mat CCi = join_vert(join_horiz(C11i + bd * SS * bd.t(), CC12),
                            join_horiz(CC12.t(), SS));
  CCi = (CCi + CCi.t()) / 2.0;
  
  arma::mat dist_Y_Ytilda = join_vert(join_horiz(dist_Y, dist),
                                      join_horiz(dist.t(), dist2));
  
  // updating Ytilda with the generated Ytilda_l
  BNM0 = arma::join_cols(Y, Ytilda);
  
  return List::create(M,
                      Ytilda.cols(0, 1),
                      C,
                      CCi,
                      dist_Y_Ytilda,
                      BNM0);
  
}


#endif // UPDATE_LATENT_PROCESSES_HPP
