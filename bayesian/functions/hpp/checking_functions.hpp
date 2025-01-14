#ifndef CHECKING_FUNCTIONS_HPP
#define CHECKING_FUNCTIONS_HPP


bool check__link_function(std::string link_function) {
  
  return (link_function == "probit" || link_function == "sigmoid");
  
}


bool check_BKtilda(const arma::mat& BKtilda) {
  // Check if the matrix has at least 5 columns
  
  // Extract the fifth column
  arma::vec column = BKtilda.col(4); // Column indices are 0-based in Armadillo
  
  // Check if all elements in the column are either 1 or 2
  bool valid = arma::all(column == 1 || column == 2);
  
  return valid;
}


#endif // CHECKING_FUNCTIONS_HPP
