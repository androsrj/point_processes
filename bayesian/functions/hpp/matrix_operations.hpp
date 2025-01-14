#ifndef MATRIX_OPERATIONS_HPP
#define MATRIX_OPERATIONS_HPP


arma::mat mpd(const arma::mat& m, const double tol) {
  // Make Positive Definite (MPD)
  int n = m.n_rows;
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, m);
  
  double delta = 2 * tol;
  arma::vec tau = arma::zeros(n);
  for (int i = 0; i < n; ++i) {
    if (delta - eigval(i) >= 0) {
      tau(i) = delta - eigval(i);
    }
  }
  
  arma::mat dm = eigvec * diagmat(tau) * eigvec.t();
  
  arma::mat M = m + dm;
  M = (M + M.t()) / 2.0;
  
  return M;
  
}


arma::mat inv_sympd_2(arma::mat SS_aux){

  arma::mat SS;
  double tol = 1e-5;
  
  try{
    SS = inv_sympd(SS_aux);
  } catch(...) {
    SS_aux = mpd(SS_aux, tol);
    SS_aux = (SS_aux + SS_aux.t()) / 2.0;
    SS = inv_sympd(SS_aux);
  }
  
  SS = (SS + SS.t()) / 2.0;
  
  return SS;
  
}


arma::mat Wmatrix(
    const arma::mat& W){
  // 
  // 
  // Parameters:
  //   d:
  
  int p = W.n_cols - 1;
  
  // Initialize WW with the first column of W as a diagonal matrix
  arma::mat WW = arma::diagmat(W.col(0));
  
  // If there are additional columns, concatenate diagonal blocks
  if (p >= 1) {
    for (int i = 1; i <= p; ++i) {
      WW = arma::join_cols(WW, arma::diagmat(W.col(i)));
    }
  }
  
  return WW;
}


arma::mat dmatrix(
    const arma::mat& s){
  // 
  // 
  // Parameters:
  //   d:
  
  arma::mat tt1 = repmat(s.col(0), 1, s.n_rows);
  arma::mat tt2 = repmat(s.col(1), 1, s.n_rows);
  
  arma::mat diff1 = tt1 - tt1.t();
  arma::mat diff2 = tt2 - tt2.t();
  
  arma::mat WW = sqrt(pow(diff1, 2) + pow(diff2, 2));
  
  return WW;
}


#endif // MATRIX_OPERATIONS_HPP
