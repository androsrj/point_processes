# exact-bayesian-spatial-cox-process--rcpp

**Rcpp Implementation of ‘Exact Bayesian Inference in Spatiotemporal Cox Processes Driven by Multivariate Gaussian Processes’**  

Flávio B. Gonçalves, Dani Gamerman, Exact Bayesian Inference in Spatiotemporal Cox Processes Driven by Multivariate Gaussian Processes, Journal of the Royal Statistical Society Series B: Statistical Methodology, Volume 80, Issue 1, January 2018, Pages 157–175, https://doi.org/10.1111/rssb.12237

Code implemented by: Izabel Nolau.

The code enables exact Bayesian inference for spatial Cox processes, with the models driven by multivariate Gaussian processes. The core algorithms are implemented using **Rcpp**, which integrates C++ code into R for optimized computational performance, particularly useful for large datasets and complex models.


## Folder Structure

The directory structure is organized as follows:

- `datasets/` — Contains the dataset files used for analysis.
- `function/` — Contains the Rcpp functions and auxiliary scripts necessary for the computation.
- `results/` — Directory where the results of the analysis (MCMC samples) are saved.

## How to Run

To perform the MCMC analysis, follow these steps:

1. **Prepare your environment**  
   Ensure that you have R and the necessary packages installed, including **Rcpp**.

2. **Configure the script**  
   Open the `run_mcmc.R` script and adjust the configuration settings according to your needs, such as specifying the dataset and adjusting MCMC parameters.

3. **Run the script**  
   In R, execute the script. This will run the MCMC algorithm, perform the Bayesian inference, and save the results to the `results/` directory.
