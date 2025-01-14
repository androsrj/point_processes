// Include necessary headers
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <RcppParallel.h>
#include <iostream>
#include <mvnorm.h>
#include <pg.h>

// Rcpp annotations for dependencies
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppParallel, pg)]]
// [[Rcpp::plugins(openmp)]]

// Using namespaces for convenience
using namespace Rcpp;
using namespace RcppParallel;

#include <hpp/location_functions.hpp>
#include <hpp/checking_functions.hpp>
#include <hpp/vector_operations.hpp>
#include <hpp/matrix_operations.hpp>
#include <hpp/link_functions.hpp>
#include <hpp/generate_distributions.hpp>
#include <hpp/covariance_function.hpp>
#include <hpp/update_retrospective_gaussian_process.hpp>
#include <hpp/update_latent_processes__L1.hpp>
#include <hpp/update_gaussian_process.hpp>
#include <hpp/get_linear_predictor.hpp>
#include <hpp/update_log_posterior__L1.hpp>


inline void set_seed(unsigned int seed){
  // Routine for setting seed within cpp.

  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
  
}


std::string addLeadingZeros(int value, int width) {
  // String representation of a integer until it reaches the specified width.
  
  std::string result = std::to_string(value);
  int aux = result.size();
  while (aux < width) {
    result = "0" + result;
    aux = result.size();
  }
  
  return result;
  
}


// [[Rcpp::export]]
double mcmc(
    arma::mat Y,
    double H,
    double W,
    bool lambdastar_estimated,
    double lambdastar_init,
    double lambdastar_prior_a,
    double lambdastar_prior_b,
    std::string link_function,
    bool gaussian_processes_estimated,
    double gaussian_processes_means,
    double gaussian_processes_variances,
    double gaussian_processes_ranges,
    bool additional_locations__lp_if = false,
    NumericMatrix additional_locations = NumericMatrix(),
    double gaussian_processes_exponent = 1.9,
    bool intensity_function_estimated = false,
    int intensity_function_grid_size = 2500,
    std::string path_to_save = "",
    double niter = 1000,
    double niter_to_save = 1000,
    int seed = 1
){
  
  Timer timer;
  timer.step("");
  double tot_time;
  int counter = 1;
  arma::vec tt;
  
  // ---------------------------------------------------------------------------
  // Basic checks --------------------------------------------------------------
  
  if ((static_cast<int>(niter) % static_cast<int>(niter_to_save)) != 0){
    stop("niter_to_save must be a multiple of niter.");
  }
  
  
  // ---------------------------------------------------------------------------
  // Defining components -------------------------------------------------------
  
  set_seed(seed);
  
  int n = Y.n_rows;
  int niter_int = niter;
  
  double range_taxa = 0;
  double niter_inicial = niter_to_save;
  
  List list__latent_processes, list__range;
  
  arma::uvec two_first = {0, 1};
  arma::uvec three_first = {0, 1, 2};
  arma::mat term1, term2, dist, sm;
  
  double m;
  
  
  // ---------------------------------------------------------------------------
  // MCMC Quantities -----------------------------------------------------------
  
  double vol_S = H * W;
  
  // lambdastar
  arma::vec chain__lambdastar(niter_to_save);
  double lambdastar_prop = lambdastar_init;
  
  // additional locations
  arma::mat additional_locations_mat = Rcpp::as<arma::mat>(additional_locations);
  int n_add = additional_locations_mat.n_rows;
  arma::mat additional_locs_i(n_add, 2);
  arma::mat chain__additional_locations__lp(n_add, niter_to_save);
  arma::mat chain__additional_locations__if(n_add, niter_to_save);
  
  // ranges
  bool ranges_estimated = false;
  double range_log_prior = 0;
  
  // intensity function
  arma::mat grid = expand_grid(ceil(sqrt(intensity_function_grid_size)), W, H);
  arma::vec lambda = arma::zeros(grid.n_rows);
  arma::vec squared_lambda = arma::zeros(grid.n_rows);
  int grid_rows = grid.n_rows;
  arma::vec lambda_i(grid.n_rows);
  double total_lambda = 0;
  arma::mat intensity_function_mean;
  arma::mat squared_intensity_function_mean;
  
  // log likelihood
  double one_cell_grid_area = (vol_S / grid.n_rows);
  arma::vec chain__log_likelihood = arma::ones(niter_to_save) * arma::datum::nan;
  arma::vec chain__integrated_if(niter_to_save);
  arma::vec chain__sum_log_lambda_y(niter_to_save);
  arma::mat Y_aux;
  double sum_log_lambda_y;
  
  // log posterior
  arma::vec chain__log_posterior(niter_to_save);
  double gp_log_dens = 0;
  double log_posterior_i;
  
  
  // ---------------------------------------------------------------------------
  // Inicial configuration -----------------------------------------------------
  
  arma::mat dist_Y, dist_Y_Ytilda, BNM0, C11, C11i, C, Ci;
  
  BNM0 = Y.cols(three_first);
  Y = Y.cols(two_first);
  dist_Y = dmatrix(Y);
  
  // autocorrelation function
  C11 = cov_function(
    gaussian_processes_variances,
    gaussian_processes_ranges,
    gaussian_processes_exponent,
    dist_Y);
  C11 = (C11 + C11.t()) / 2.0;
  C11i = arma::inv(C11);
  C = C11;
  Ci = C11i;
  
  list__latent_processes = update_latent_processes_1(
    n, BNM0, lambdastar_init, C11i, C11, Ci,
    gaussian_processes_means,
    gaussian_processes_variances,
    gaussian_processes_ranges,
    gaussian_processes_exponent,
    H, W, dist_Y, link_function);
  
  m = list__latent_processes[0];
  sm = Rcpp::as<arma::mat>(list__latent_processes[1]);
  C = Rcpp::as<arma::mat>(list__latent_processes[2]);
  Ci = Rcpp::as<arma::mat>(list__latent_processes[3]);
  dist_Y_Ytilda = Rcpp::as<arma::mat>(list__latent_processes[4]);
  BNM0 = Rcpp::as<arma::mat>(list__latent_processes[5]);
  
  
  // ---------------------------------------------------------------------------
  // MCMC ----------------------------------------------------------------------
  
  int i_aux = 0;
  std::string i_str;
  arma::vec width_aux = {
    5, static_cast<double>(
        std::to_string(static_cast<int>(niter / niter_to_save)).size())};
  int width = max(width_aux);
  
  for (int i = 1; i <= niter; ++i) {
    
    // -------------------------------------------------------------------------
    // Updating the latent processes -------------------------------------------
    
    list__latent_processes = update_latent_processes_1(
      n, BNM0, lambdastar_prop, C11i, C11, Ci,
      gaussian_processes_means,
      gaussian_processes_variances,
      gaussian_processes_ranges,
      gaussian_processes_exponent,
      H, W, dist_Y, link_function);
    
    m = list__latent_processes[0];
    sm = Rcpp::as<arma::mat>(list__latent_processes[1]);
    C = Rcpp::as<arma::mat>(list__latent_processes[2]);
    Ci = Rcpp::as<arma::mat>(list__latent_processes[3]);
    dist_Y_Ytilda = Rcpp::as<arma::mat>(list__latent_processes[4]);
    BNM0 = Rcpp::as<arma::mat>(list__latent_processes[5]);
    
    
    // -------------------------------------------------------------------------
    // Updating lambdastar -----------------------------------------------------
    
    if (lambdastar_estimated) {
      
      lambdastar_prop = R::rgamma(lambdastar_prior_a + (n + m),
                                  1/(lambdastar_prior_b + (H * W)));
      chain__lambdastar(i - i_aux - 1) = lambdastar_prop;
      
    }
    
    
    // -------------------------------------------------------------------------
    // Updating the Gaussian processes -----------------------------------------
    
    if (gaussian_processes_estimated){
      
      BNM0 = update_beta_l(
        n, m, C, gaussian_processes_means, BNM0.col(2), Ci, link_function);
      BNM0 = join_horiz(join_vert(Y, sm), BNM0);
      
    }
    
    
    // -------------------------------------------------------------------------
    // Updating lambda ---------------------------------------------------------
    
    Y_aux = BNM0.submat(0, 0, n - 1, BNM0.n_cols - 1);
    
    sum_log_lambda_y = (n * log(lambdastar_prop)) +
      sum(log_F(link_function, Y_aux.col(2)));
    
    chain__sum_log_lambda_y(i - i_aux - 1) = sum_log_lambda_y;
    
    if (intensity_function_estimated){
      
      term1 = pow(repmat(grid.col(0), 1, n + m) -
        repmat(BNM0.col(0).t(), grid_rows, 1), 2);
      term2 = pow(repmat(grid.col(1), 1, n + m) -
        repmat(BNM0.col(1).t(), grid_rows, 1), 2);
      dist = arma::sqrt(term1 + term2); // grid_rows x BNM0.nrows
      arma::uvec mm = arma::index_min(dist, 1); // locations with min distance
      
      arma::vec F_link(grid_rows);
      for (int ii = 0; ii < grid_rows; ii++) {
        F_link(ii) = F_scalar(link_function, BNM0(mm(ii), 2)); //
      }
      
      lambda_i = lambdastar_prop * F_link;
      lambda += lambda_i;
      squared_lambda += pow(lambda_i, 2);
      total_lambda += 1;
      
      // ---------------------------------------------------------------------
      // Log-likelihood ------------------------------------------------------
      
      chain__log_likelihood(i - i_aux - 1) = (
        - one_cell_grid_area * sum(lambda_i)) + sum_log_lambda_y;
      
      chain__integrated_if(i - i_aux - 1) = (
        one_cell_grid_area * sum(lambda_i));
      
    }
    
    
    // -------------------------------------------------------------------------
    // Generating the GP at additional locations -------------------------------
    
    if (additional_locations__lp_if){
      
      additional_locs_i = get_linear_predictor_1(
        BNM0, lambdastar_prop, additional_locations_mat, link_function);
      
      chain__additional_locations__lp.col(i - i_aux - 1) = 
        additional_locs_i.col(0);
      
      chain__additional_locations__if.col(i - i_aux - 1) = 
        additional_locs_i.col(1);
      
    }
    
    
    // -------------------------------------------------------------------------
    // Log-posterior -----------------------------------------------------------
    
    log_posterior_i = update_log_posterior_1(
      n, m,
      BNM0,
      lambdastar_prop,
      lambdastar_prior_a,
      lambdastar_prior_b,
      range_log_prior,
      gaussian_processes_means,
      C,
      gp_log_dens,
      vol_S,
      lambdastar_estimated,
      gaussian_processes_estimated,
      ranges_estimated,
      link_function);
    
    chain__log_posterior(i - i_aux - 1) = log_posterior_i;
    
    
    // -------------------------------------------------------------------------
    // Saving files ------------------------------------------------------------
    
    if ((i - i_aux) == niter_to_save) {
      
      i_str = addLeadingZeros(static_cast<int>(i / niter_to_save), width);
      
      // Chains
      std::ofstream outfile_lambdastar(
          path_to_save + "chain__lambdastar__" + i_str + ".txt");
      outfile_lambdastar << chain__lambdastar;
      outfile_lambdastar.close();
      
      std::ofstream outfile_add_gp(
          path_to_save + "chain__additional_locations__lp__" + i_str + ".txt");
      outfile_add_gp << chain__additional_locations__lp;
      outfile_add_gp.close();
      
      std::ofstream outfile_add_if(
          path_to_save + "chain__additional_locations__if__" + i_str + ".txt");
      outfile_add_if << chain__additional_locations__if;
      outfile_add_if.close();
      
      std::ofstream outfile_log_posterior(
          path_to_save + "chain__log_posterior__" + i_str + ".txt");
      outfile_log_posterior << chain__log_posterior;
      outfile_log_posterior.close();
      
      std::ofstream outfile_log_likelihood(
          path_to_save + "chain__log_likelihood__" + i_str + ".txt");
      outfile_log_likelihood << chain__log_likelihood;
      outfile_log_likelihood.close();
      
      std::ofstream outfile_integrated_if(
          path_to_save + "chain__integrated_if__" + i_str + ".txt");
      outfile_integrated_if << chain__integrated_if;
      outfile_integrated_if.close();    
      
      intensity_function_mean = arma::join_rows(grid, lambda/total_lambda);
      std::ofstream outfile_if(
          path_to_save + "intensity_function_mean__" + i_str + ".txt");
      outfile_if << intensity_function_mean;
      outfile_if.close();
      lambda = arma::zeros(grid.n_rows);
      
      squared_intensity_function_mean = arma::join_rows(
        grid, squared_lambda/total_lambda);
      std::ofstream outfile_if2(
          path_to_save + "squared_intensity_function_mean__" + i_str + ".txt");
      outfile_if2 << squared_intensity_function_mean;
      outfile_if2.close();
      squared_lambda = arma::zeros(grid.n_rows);
      
      total_lambda = 0;
      
      // Elapsed time
      timer.step("");
      tt = Rcpp::as<arma::vec>(timer);
      tot_time = (tt(counter) - tt(counter - 1))/1e+9;
      std::ofstream outfile_time(
          path_to_save + "elapsed_time__" + i_str + ".txt");
      outfile_time << tot_time;
      outfile_time.close();
      counter += 1;
      
    }
    
    if ((i - i_aux) == niter_to_save) {
      i_aux += niter_to_save;
    }
    
    Rprintf("\rIteration: %d/%d", i, niter_int);
    
  }
  
  return 0;
  
}