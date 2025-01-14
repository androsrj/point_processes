library(spatstat.sparse)
library(spatstat.geom)
library(spatstat.utils)

fit_homog <- function(obj, method, window, type = "trend", nCores = 20) {
  
  if(missing(window)) {
    pred_window <- obj$window
  } else {
    pred_window <- intersect.owin(obj$window, window)
  }
  
  if (method == "full") {
    fit1 <- ppm(unmark(obj) ~ 1)
    preds <- predict(fit1, type = type, window = pred_window)
    ints <- predict(fit1, interval="confidence", type = type, window = pred_window)
  } else if (method == "divided") {
    indices <- sample(1:nCores, size = obj$n, replace = TRUE)
    models <- vector("list", nCores)
    for (i in 1:nCores) {
      disk <- subset(unmark(obj), indices == i)
      models[[i]] <- ppm(disk ~ 1)
    }
    preds <- Reduce("+", lapply(1:nCores, \(i) predict.ppm(models[[i]], 
                                                           type = type,
                                                           window = pred_window)))
    se_agg <- Reduce("+", lapply(1:nCores, \(i) predict.ppm_homog(models[[i]], 
                                                           type = type,
                                                           window = pred_window,
                                                           nCores = 20,
                                                           interval='confidence'))) 
    se <- (mean(se_agg[[2]]) - mean(se_agg[[1]])) / (2 * sqrt(nCores))
    ints <- list(preds - se, preds + se)
  } else {
    stop("method must either be 'full' or 'divided'.")
  }
  return(list(preds = preds, interval = ints))
}