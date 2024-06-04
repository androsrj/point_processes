library(stopp)
library(spatstat)
library(dplyr)

N <- clmfires$n
time <- seq(0, 1, length = N)
fires <- stp(data.frame(x = clmfires$x, y = clmfires$y, t = time))
nCores <- 8
indices <- sample(rep(1:nCores, each = N / nCores))

# Calculate volume for entire window
vol_full <- diff(range(fires$df[, 1])) * 
  diff(range(fires$df[, 2])) * 
  diff(range(fires$df[, 3]))

# Calculate volume for smaller window
new_index <- fires$df$x > 250 & fires$df$y < 150 & fires$df$t > 0.25
sum(new_index)
newdf <- fires$df[new_index, ]
vol_window <- diff(range(newdf[, 1])) * 
  diff(range(newdf[, 2])) * 
  diff(range(newdf[, 3]))

# For loop (will later change to parallel)
pred_full <- pred_window <- runtimes <- numeric(nCores)
for (i in 1:nCores) {
  index <- which(indices == i)
  subset <- stp(data.frame(x = clmfires$x[index], y = clmfires$y[index], t = time[index]))
  lgcp_subset <- stlgcppm(subset, ~ 1 + poly(x, 3) + y)
  pred_full[i] <- mean(lgcp_subset$l) * vol_full
  pred_window[i] <- mean(lgcp_subset$l) * vol_window
  runtimes[i] <- as.numeric(substr(lgcp_subset$time, 1, nchar(lgcp_subset$time) - 8))
}
sum(pred_full)
sum(pred_window)
sum(runtimes)
max(runtimes)