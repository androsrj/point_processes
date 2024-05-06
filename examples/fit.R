library(spatstat.sparse)
library(spatstat.geom)
library(spatstat.utils)

fit <- function(obj, method, window, type = "trend", nCores = 20) {
  
  if(missing(window)) {
    pred_window <- obj$window
  } else {
    pred_window <- intersect.owin(obj$window, window)
  }
  
  if (method == "full") {
    fit1 <- ppm(unmark(obj) ~ polynom(x, y, 3))
    preds <- predict(fit1, type = type, window = pred_window)
    ints <- predict(fit1, interval="confidence", type = type, window = pred_window)
  } else if (method == "divided") {
    indices <- sample(1:nCores, size = obj$n, replace = TRUE)
    models <- vector("list", nCores)
    for (i in 1:nCores) {
      disk <- subset(unmark(obj), indices == i)
      models[[i]] <- ppm(disk ~ polynom(x, y, 3))
    }
    preds <- Reduce("+", lapply(1:nCores, \(i) predict.ppm2(models[[i]], 
                                                            type = type,
                                                            window = pred_window,
                                                            N = obj$n, 
                                                            m = table(indices)[i]))) / nCores
    ints <- Reduce("+", lapply(1:nCores, \(i) predict.ppm2(models[[i]], 
                                                           type = type,
                                                           window = pred_window,
                                                           N = obj$n, 
                                                           m = table(indices)[i], 
                                                           interval='confidence'))) / nCores
  } else {
    stop("method must either be 'full' or 'divided'.")
  }
  return(list(preds = preds, interval = ints))
}
