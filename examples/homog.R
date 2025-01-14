library(spatstat)
library(ggplot2)
source("fit_homog.R")
source("predict.ppm_homog.R")
source("predconfPois_homog.R")

fire_centr <- fit_homog(rescale(clmfires, 100), method = "full", type = "trend")
centr_preds <- fire_centr$preds
centr_ints <- fire_centr$interval

fires_div <- fit_homog(rescale(clmfires, 100), method = "divided", type = "trend")
agg_preds <- fires_div$preds
agg_ints <- fires_div$interval

# Point prediction for entire window is the same for both
mean(centr_preds)
mean(agg_preds)

# Intervals
mean(centr_ints[[1]])
mean(centr_ints[[2]])
mean(agg_ints[[1]])
mean(agg_ints[[2]])

# Average interval width
mean(centr_ints[[2]] - centr_ints[[1]], na.rm=TRUE)
mean(agg_ints[[2]] - agg_ints[[1]], na.rm=TRUE)

# Count predictions
n_regions <- 5
region_radius <- 0.75
set.seed(75193)
centroids_x <- runif(n_regions, min(clmfires$x)/100, 
                     max(clmfires$x)/100)
centroids_y <- runif(n_regions, min(clmfires$y)/100,
                     max(clmfires$y)/100)

fire_centr <- fit_homog(rescale(clmfires, 100), method = "full", type = "trend")
centr_preds <- fire_centr$preds
for (i in 1:n_regions) {
  temp_window <- disc(radius = region_radius, centre = c(centroids_x[i], centroids_y[i]))
  
  fire_centr <- fit_homog(rescale(clmfires, 100), window = temp_window, 
                    method = "full", type = "count")
  centr_preds <- fire_centr$preds
  centr_ints <- fire_centr$interval
  
  fires_div <- fit_homog(rescale(clmfires, 100), window = temp_window,
                   method = "divided", type = "count")
  agg_preds <- fires_div$preds
  agg_ints <- fires_div$interval
  
  true_n <- rescale(clmfires, 100)[temp_window]$n
  cat(paste0("Region ", i, ":\n"))
  cat(paste0("True count in region ", i, ": ", true_n, "\n"))
  cat(paste0("Estimated count in region ", i, " from full model: ", 
             round(centr_preds, 2), " (", 
             round(centr_ints[1], 2), ", ", 
             round(centr_ints[2], 2), ") \n"))
  cat(paste0("Estimated count in region ", i, " from partitioned model: ", 
             round(agg_preds, 2), " (", 
             round(agg_ints[1], 2), ", ", 
             round(agg_ints[2], 2), ") \n"))
  cat("\n")
}

