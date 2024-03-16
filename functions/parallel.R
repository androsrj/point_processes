# Helper function to run subsets in parallel for D-and-C
pp_parallel <- function(i) {
  x <<- cbind(x1 = lansing$x[subsets[[i]]], x2 = lansing$y[subsets[[i]]])
  langevin_pp(x = x, t, N, K, starting, step = step_sizes, nIter = 5000, nBurn = 500, nThin = 2)
}


