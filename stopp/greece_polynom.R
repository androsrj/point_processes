library(stopp)
library(spatstat)
library(dplyr)
library(plotrix)

nCores <- 8
df <- 4
N <- nrow(greececatalog$df)

# Fit full model
lgcp_full <- stlgcppm(greececatalog, ~ 1 + poly(x, df) + poly(y, df) + poly(t, df))

# Divide up data
indices <- sample(c(rep(1:nCores, each = N / nCores), 1:(N %% nCores)))

# Calculate volume for entire window
full_window <- owin(xrange = range(greececatalog$df$x), yrange = range(greececatalog$df$y))
vol_full <- area.owin(full_window) * diff(range(greececatalog$df$t))

# Dividing time treshold for smaller windows
div_time <- 40000
window_time <- c(39000, 41000)

# Calculate volume for smaller window (window 1)
r1 <- 1.25
window1 <- intersect.owin(disc(radius = r1, centre = c(22, 39)),
                          full_window)
index1 <- ((greececatalog$df$x - 22)^2 + (greececatalog$df$y - 39)^2 < r1^2) & (greececatalog$df$t > div_time)
vol_window1 <- area.owin(window1) * (max(greececatalog$df$t) - div_time)

# Calculate volume for medium window (window 2)
r2 <- 1.5
window2 <- intersect.owin(disc(radius = r2, centre = c(22, 36)),
                          full_window)
index2 <- (greececatalog$df$x - 22)^2 + (greececatalog$df$y - 36)^2 < r2^2 & (greececatalog$df$t < div_time)
vol_window2 <- area.owin(window2) * (div_time - min(greececatalog$df$t))

# Calculate volume for larger window (window 3)
r3 <- 2
window3 <- intersect.owin(disc(radius = r3, centre = c(26, 36)),
                          full_window)
index3 <- (greececatalog$df$x - 26)^2 + (greececatalog$df$y - 36)^2 < r3^2 & 
  (greececatalog$df$t > window_time[1]) & (greececatalog$df$t < window_time[2])
vol_window3 <- area.owin(window3) * (window_time[2] - window_time[1])

# For loop (will later change to parallel)
runtimes <- numeric(nCores)
agg_intens <- numeric(N)
coefs <- matrix(0, nrow = nCores, ncol = 1 + 3*df)
for (i in 1:nCores) {
  index <- which(indices == i)
  subset <- stp(data.frame(x = greececatalog$df$x[index], 
                           y = greececatalog$df$y[index], 
                           t = greececatalog$df$t[index]))
  lgcp_subset <- stlgcppm(subset, ~ 1 + poly(x, df) + poly(y, df) + poly(t, df))
  coefs[i , ] <- lgcp_subset$IntCoefs
  agg_intens[indices == i] <- lgcp_subset$l
  runtimes[i] <- as.numeric(substr(lgcp_subset$time, 1, nchar(lgcp_subset$time) - 8))
  print(mean(lgcp_subset$w[1:length(index)]) / mean(lgcp_full$w[index]))
}

# Predictions (agg model)
pred_agg <- c(mean(agg_intens[index1]) * vol_window1 * nCores, 
              mean(agg_intens[index2]) * vol_window2 * nCores, 
              mean(agg_intens[index3]) * vol_window3 * nCores, 
              mean(agg_intens) * vol_full * nCores)

# Predictions (full model)
pred_full <- c(mean(lgcp_full$l[index1]) * vol_window1, 
               mean(lgcp_full$l[index2]) * vol_window2, 
               mean(lgcp_full$l[index3]) * vol_window3, 
               mean(lgcp_full$l) * vol_full)

sum(index1)
sum(index2)
sum(index3)
sum(runtimes)
max(runtimes)


results <- data.frame(
  Region = 1:4,
  Truth = c(sum(index1), sum(index2), sum(index3), N),
  Pred_Aggregated = pred_agg,
  Pred_Full = pred_full,
  Time_Aggregated = paste0(max(runtimes), " minutes"),
  Time_Full = lgcp_full$time
)
results


apply(coefs, 2, mean)
lgcp_full$IntCoefs
