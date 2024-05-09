library(spatstat)
library(ggplot2)
source("fit.R")
source("predict.ppm2.R")
source("predconfPois2.R")

fire_centr <- fit(rescale(clmfires, 100), method = "full", type = "trend")
centr_preds <- fire_centr$preds
centr_ints <- fire_centr$interval

fires_div <- fit(rescale(clmfires, 100), method = "divided", type = "trend")
agg_preds <- fires_div$preds
agg_ints <- fires_div$interval

# Predictive plots
#pdf("fires.pdf", width = 10, height = 5)
lim_lwr <- min(c(min(centr_preds), min(agg_preds)))
lim_upr <- max(c(max(centr_preds), max(agg_preds)))
bounds <- c(lim_lwr, lim_upr)
par(mfrow=c(1,2))
image(centr_preds, zlim = bounds,
      main = "Centralized Model \n (Wildfires Data)") # Whole model
image(agg_preds, zlim = bounds,
      main = "Partitioned Model \n (Wildfires Data)") # Divided model
plot(clmfires$x, clmfires$y)
#dev.off()

lim_lwr <- min(c(min(centr_ints), min(agg_ints[[1]]), min(agg_ints[[2]])))
lim_upr <- max(c(max(centr_ints), max(agg_ints[[1]]), max(agg_ints[[2]])))
bounds <- c(lim_lwr, lim_upr)
#pdf("fire_intervals.pdf")
image(centr_ints, zlim = bounds,
      main = "Centralized Model \n (Wildfires Data)")
image(agg_ints, zlim = bounds,
      main = "Partitioned Model \n (Wildfires Data)")
#dev.off()

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

fire_centr <- fit(rescale(clmfires, 100), method = "full", type = "trend")
centr_preds <- fire_centr$preds
image(centr_preds, zlim = bounds,
      main = "Centralized Model \n (Wildfires Data)") # Whole model
for (i in 1:n_regions) {
  temp_window <- disc(radius = region_radius, centre = c(centroids_x[i], centroids_y[i]))
  
  fire_centr <- fit(rescale(clmfires, 100), window = temp_window, 
                    method = "full", type = "count")
  centr_preds <- fire_centr$preds
  centr_ints <- fire_centr$interval
  
  fires_div <- fit(rescale(clmfires, 100), window = temp_window,
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
  lines(temp_window$bdry[[1]]$x, temp_window$bdry[[1]]$y)
  text(centroids_x[i], centroids_y[i], true_n)
}

