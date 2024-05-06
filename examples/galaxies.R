library(spatstat)
library(usmap)
library(ggplot2)
source("fit_example.R")
source("predict.ppm2.R")

galaxies <- fit(rescale(shapley, 100))
centr_preds <- galaxies$centr_preds
centr_ints <- galaxies$centr_ints
agg_preds <- galaxies$agg_preds
agg_ints <- galaxies$agg_ints

# Predictive plots
pdf("galaxies.pdf", width = 10, height = 5)
lim_lwr <- min(c(min(centr_preds), min(agg_preds)))
lim_upr <- max(c(max(centr_preds), max(agg_preds)))
bounds <- c(lim_lwr, lim_upr)
par(mfrow=c(1,2))
image(centr_preds, zlim = bounds, main = "") # Whole model
title("Centralized Model \n (Galaxies Data)", line = -3)
image(agg_preds, zlim = bounds, main = "") # Divided model
title("Partitioned Model \n (Galaxies Data)", line = -3)
dev.off()

lim_lwr <- min(c(min(centr_ints), min(agg_ints[[1]]), min(agg_ints[[2]])))
lim_upr <- max(c(max(centr_ints), max(agg_ints[[1]]), max(agg_ints[[2]])))
bounds <- c(lim_lwr, lim_upr)
pdf("galaxies_intervals.pdf")
image(centr_ints, zlim = bounds,
      main = "Centralized Model (Galaxies Data)")
image(agg_ints, zlim = bounds,
      main = "Partitioned Model (Galaxies Data)")
dev.off()

# Average interval width
mean(centr_ints[[2]] - centr_ints[[1]], na.rm=TRUE)
mean(agg_ints[[2]] - agg_ints[[1]], na.rm=TRUE)


