#ifndef VECTOR_OPERATIONS_HPP
#define VECTOR_OPERATIONS_HPP


arma::vec subsetVector(
    const arma::vec& lambdastar,
    const arma::vec& ind){
  // 
  // 
  // Parameters:
  //   d:
  
  // Initialize a new arma::vec to store the subset
  arma::vec subset(ind.n_elem);
  int aux = ind.n_elem;  
  
  // Extract elements at specified indices
  for (int i = 0; i < aux; ++i) {
    // Subtract 1 from ind(i) because C++ uses 0-based indexing
    subset(i) = lambdastar(ind(i) - 1);
  }
  
  return subset;
}


IntegerVector order(arma::vec x) {
  int n = x.size();
  IntegerVector indices = seq(0, n - 1);
  
  std::sort(indices.begin(), indices.end(), [&](int i, int j) {
    return x[i] < x[j];
  });
  
  return indices;
}


bool isNotOrdered(const arma::vec& lambdastar_prop, int L){
  
  for (int i = 0; i < L - 1; ++i) {
    if (lambdastar_prop(i) > lambdastar_prop(i + 1)) {
      return true;
    }
  }
  
  return false;
}


double euclideanDistance(
    arma::vec point1,
    arma::vec point2){
  // 
  // 
  // Parameters:
  //   d:
  
  int n = point1.size();
  double sum = 0;
  
  for (int i = 0; i < n; i++) {
    double diff = point1(i) - point2(i);
    sum += diff * diff;
  }
  
  return sqrt(sum);
  
}


#endif // VECTOR_OPERATIONS_HPP
