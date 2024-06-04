library(stopp)
library(spatstat)
library(dplyr)

N <- nrow(greececatalog$df)
nCores <- 8
indices <- sample(rep(1:nCores, each = N / nCores))

# Calculate volume for entire window
vol_full <- diff(range(greececatalog$df[, 1])) * 
  diff(range(greececatalog$df[, 2])) * 
  diff(range(greececatalog$df[, 3]))

# Calculate volume for smaller window
new_index <- greececatalog$df$x > 25 & greececatalog$df$y < 36 & greececatalog$df$t > 40000
sum(new_index)
newdf <- greececatalog$df[new_index, ]
vol_window <- diff(range(newdf[, 1])) * 
  diff(range(newdf[, 2])) * 
  diff(range(newdf[, 3]))

# For loop (will later change to parallel)
pred_full <- pred_window <- runtimes <- numeric(nCores)
for (i in 1:nCores) {
  index <- which(indices == i)
  subset <- stp(data.frame(x = greececatalog$df$x[index], 
                           y = greececatalog$df$y[index], 
                           t = greececatalog$df$t[index]))
  lgcp_subset <- stlgcppm(subset, ~ 1 + poly(x, 3) + y + t)
  pred_full[i] <- mean(lgcp_subset$l) * vol_full
  pred_window[i] <- mean(lgcp_subset$l) * vol_window
  runtimes[i] <- as.numeric(substr(lgcp_subset$time, 1, nchar(lgcp_subset$time) - 8)) * 60
}
sum(pred_full)
sum(pred_window)
sum(runtimes)
max(runtimes)
