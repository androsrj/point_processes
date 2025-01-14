#ifndef GENERATE_DISTRIBUTIONS_HPP
#define GENERATE_DISTRIBUTIONS_HPP


NumericVector generate_truncated_norm(double a, double b) {
  // Use Rcpp to call the rtruncnorm function in R
  Environment pkg = Environment::namespace_env("truncnorm");
  Function rtruncnorm = pkg["rtruncnorm"];
  
  NumericVector result = rtruncnorm(1, a, b, 0, 1);
  return result;
}


arma::vec runif_cpp(
    int n,
    double min,
    double max){
  // 
  // 
  // Parameters:
  //   d:
  
  arma::vec u(n);
  for (int i = 0; i < n; i++) {
    u(i) = R::runif(0, 1) * (max - min) + min;
  }
  
  return u;
  
}


arma::vec rbinom_cpp(
    int n,
    arma::vec p){
  // 
  // 
  // Parameters:
  //   d:
  
  arma::vec u(n);
  
  for (int i = 0; i < n; i++) {
    u(i) = R::runif(0, 1);
  }
  
  arma::vec sm = (u - p) / abs(u - p); // 1 if u > p; and -1 if u < p
  sm = (1 - sm) / 2; // 0 if u > p; and 1 if u < p
  
  return sm;
  
}


arma::vec rnorm_cpp(
    int n,
    double mean,
    double sd){
  // 
  // 
  // Parameters:
  //   d:
  
  arma::vec u(n);
  for (int i = 0; i < n; i++) {
    u(i) = R::rnorm(mean, sd);
  }
  return u;
}


arma::vec rmvnorm(
    int d,
    const arma::vec& mu,
    const arma::mat& sigma2){
  // 
  // 
  // Parameters:
  //   d:
  arma::mat A;
  
  try{
    A = arma::chol(sigma2, "lower");
  } catch(...) {
    arma::mat sigma2_modified = mpd(sigma2, 1e-5);
    A = arma::chol(sigma2_modified, "lower");
  }
  
  arma::vec uncorrelated_samples = arma::randn(d, 1);
  arma::vec samples = A * uncorrelated_samples + mu;
  return samples;
  
}


double rtnorm_rejection_inf_truncated(
    double a){
  
  double z = generate_truncated_norm(a, INFINITY)[0];
  
  if (z == INFINITY) {
    Rcpp::stop("rtnorm_rejection_inf_truncated generated an infinity value");
  }
  
  return z;
  
}


double rtnorm_rejection(
    double a,
    const double b){
  
  double z = generate_truncated_norm(a, b)[0];
  
  if (z == INFINITY) {
    Rcpp::stop("rtnorm_rejection generated an infinity value");
  }
  
  return z;
  
}


double rtnorm_inversion(
    double a,
    const double b){
  
  double u = generate_truncated_norm(a, b)[0];
  
  if (u == INFINITY) {
    Rcpp::stop("rtnorm_inversion generated an infinity value");
  }
  
  return u;
  
}


arma::vec rsun(
    const arma::vec mu,
    arma::mat Sigma,
    arma::mat W,
    arma::mat iGam,
    arma::mat A){
  // 
  // 
  // Parameters:
  //   d:
  
  int m = W.n_rows;
  int d = W.n_cols;
  
  arma::vec gam = W * mu; // 0, if mu = 0;
  
  // gerando um valor inicial dentro da região
  arma::vec u0str = arma::zeros(m);
  if (A(0, 0) > 0) {
    u0str(0) = rtnorm_inversion(-gam(0) / A(0, 0), INFINITY);
  } else {
    u0str(0) = rtnorm_inversion(-INFINITY, -gam(0) / A(0, 0));
  }
  
  arma::mat A_signal = (A / abs(A) + 1.0) / 2.0; // matrix of 0 or 1
  
  for (int i = 1; i < m; ++i) {
    
    double limm0 = (- gam(i) - arma::as_scalar(
      A.row(i).subvec(0, i - 1) * u0str.subvec(0, i - 1)) ) / A(i, i);
    double ind = A_signal(i, i);
    
    arma::vec mia;
    if(ind == 1){
      mia = {limm0, INFINITY};
    } else{
      mia = {-INFINITY, limm0};
    }
    
    u0str(i) = rtnorm_inversion(mia(0), mia(1));
    
  }
  
  // embebed gibbs sampling
  for (int j = 1; j <= std::min(m, 10); ++j) {
    
    for (int i = 0; i < m; ++i) {
      
      u0str(i) = 0;
      arma::vec A_signal_i = A_signal.submat(i, i, m - 1, i);
      
      arma::vec limm = ( - gam.subvec(i, m - 1) -
        A.submat(i, 0, m - 1, m - 1) * u0str ) /
          A.submat(i, i, m - 1, i);
      
      // A_signal_i == 0.0 => min(limm)
      // A_signal_i == 1.0 => max(limm)
      
      arma::uvec A_pos = find(A_signal_i == 1.0);  
      double mi = max(limm(A_pos));
      
      if (all(A_signal_i == 1)) {
        
        double ma = INFINITY;
        
        if (mi > 8) {
          u0str(i) = rtnorm_rejection_inf_truncated(mi);
        } else {
          u0str(i) = rtnorm_inversion(mi, ma);
        }
        
      } else {
        
        arma::uvec A_neg = find(A_signal_i == 0.0);  
        double ma = min(limm(A_neg));
        
        if (mi > 7) {
          u0str(i) = rtnorm_rejection(mi, ma);
        } else if (ma < -7) {
          u0str(i) = -rtnorm_rejection(-ma, -mi);
        } else {
          u0str(i) = rtnorm_inversion(mi, ma);
        }
        
      }
      
    }
  }
  
  // Algorithm 4.1
  arma::mat tDel = W * Sigma;
  arma::mat Del = tDel.t();
  
  arma::vec u0 = A * u0str;
  arma::vec mu1 = Del * iGam * u0;
  arma::mat Sigma1 = Sigma - Del * iGam * tDel;
  Sigma1 = (Sigma1 + Sigma1.t()) / 2.0;
  
  arma::vec zstr = rmvnorm(d, mu1, Sigma1);
  
  arma::vec z = zstr + mu;
  
  return z;
  
}


arma::vec rpg_vector(const double h, const arma::vec z) {
  
  int n = z.size();
  arma::vec omega(n);
  
  for (int i = 0; i < n; ++i) {
    omega(i) = pg::rpg_scalar_hybrid(h, z(i));
    
    if (omega(i) == 0) {
      omega(i) = pow(10, -250);
    }
    
  }
  
  return omega;
  
}


double repulsive_prior_log_density(
    arma::mat U,
    int L) {
  // 
  // 
  // Parameters:
  //   d:
  
  arma::mat distMatrix = arma::zeros(L, L);
  arma::vec point1, point2;
  double result = 0;
  double value, term;
  
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < i; j++) {
      arma::vec point1 = U.row(i).t();
      arma::vec point2 = U.row(j).t();
      value = euclideanDistance(point1, point2);
      term = log(1.0 - exp(-1.5 * pow(value, 4)));
      result += term;
    }
  }
  
  return result;
  
}


static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}


arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::vec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat chol_sigma;
  try{
    chol_sigma = arma::chol(sigma);
  } catch(...) {
    arma::mat sigma_aux = mpd(sigma, 1e-5);
    chol_sigma = arma::chol(sigma_aux);
  }
  arma::mat const rooti = arma::inv(trimatu(chol_sigma));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean.t());
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }
  
  if (logd)
    return out;
  return exp(out);
}


double dmvnorm_to_double(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& S,
    const bool log_p = false){
  //
  //
  // Parameters:
  //   d:
  
  arma::mat x_matrix = x.t();
  return arma::as_scalar(dmvnrm_arma_fast(x_matrix, mu, S, log_p));
  
}


#endif // GENERATE_DISTRIBUTIONS_HPP
