#ifndef LOCATION_FUNCTIONS_HPP
#define LOCATION_FUNCTIONS_HPP


arma::mat expand_grid(
    int intensity_function_grid_size,
    int W,
    int H){
  // 
  // 
  // Parameters:
  //   d:
  
  int numRows = intensity_function_grid_size * intensity_function_grid_size;
  
  // Create an arma::mat to store combinations
  arma::mat result(numRows, 2);
  
  arma::vec vec1 = arma::linspace<arma::vec>(0, W, intensity_function_grid_size);
  arma::vec vec2 = arma::linspace<arma::vec>(0, H, intensity_function_grid_size);
  
  // Generate all combinations
  int idx = 0;
  for (int i = 0; i < intensity_function_grid_size; ++i) {
    for (int j = 0; j < intensity_function_grid_size; ++j) {
      result(idx, 0) = vec1(i);
      result(idx, 1) = vec2(j);
      ++idx;
    }
  }
  
  return(result);
}


arma::uvec closer_location (
  arma::mat M1,
  arma::mat M2
){
  
  arma::mat term1, term2, dist;
  
  int n1 = M1.n_rows;
  int n2 = M2.n_rows;
  
  term1 = pow(repmat(M1.col(0), 1, n2) -
                repmat(M2.col(0).t(), n1, 1), 2);
  term2 = pow(repmat(M1.col(1), 1, n2) -
                repmat(M2.col(1).t(), n1, 1), 2);
  dist = arma::sqrt(term1 + term2);
  arma::uvec indexes = arma::index_min(dist, 1);
  
  return indexes;
  
}


NumericVector whichTileCpp(
    arma::mat centers,
    arma::mat points){
  // 
  // 
  // Parameters:
  //   d:
  
  if(points.n_cols > 2){
    arma::uvec two_first = {0, 1};
    points = points.cols(two_first);
  }
  
  arma::mat pos = arma::join_vert(centers, points);
  int L = centers.n_rows;
  int n = points.n_rows;
  double dx,dy;
  NumericMatrix distm(n, L);
  NumericVector mins(n);
  for(int i = L; i < n + L; i++){
    for(int j = 0; j < L; j++){
      dx = pos(i, 0) - pos(j, 0);
      dy = pos(i, 1) - pos(j, 1);
      distm(i - L, j) = sqrt(dx*dx + dy*dy);
    }
    mins(i - L) = which_min(distm.row(i - L));
  }
  return(mins + 1);
  
}


arma::vec nTileCpp(
    arma::vec tiles,
    int L){
  // 
  // 
  // Parameters:
  //   d:
  
  arma::vec kl(L);
  
  for(int l = 1; l < L; l++){
    arma::uvec ids = find(tiles == l);
    arma::vec aux = tiles(ids);
    kl(l - 1) = aux.n_elem;
  }
  kl(L - 1) = tiles.n_elem - sum(kl);
  
  return(kl);
  
}


arma::uvec sample_neigh(
    arma::mat U,
    int B,
    int L){
  // 
  // 
  // Parameters:
  //   d:
  
  // choose a random row index
  int rand_idx = floor(R::runif(0, L));
  
  // calculate distances from the chosen row to all other rows
  NumericVector dist(L);
  for (int i = 0; i < L; i++) {
    double sum = 0;
    for (int j = 0; j < 2; j++) {
      sum += pow(U(i, j) - U(rand_idx, j), 2);
    }
    dist[i] = sqrt(sum);
  }
  
  // find the indices of the b-1 closest rows
  NumericVector sorted_idx = clone(dist).sort();
  arma::uvec closest_idx =  Rcpp::as<arma::uvec>(match(sorted_idx, dist));
  
  // return the sampled index and the indices of the b-1 closest rows
  arma::uvec result =  closest_idx.subvec(0, B - 1) - 1;
  
  return result;
  
}


#endif // LOCATION_FUNCTIONS_HPP
