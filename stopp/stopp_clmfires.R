library(stopp)
library(spatstat)
library(dplyr)
library(plotrix)

nCores <- 8

N <- clmfires$n
time <- seq(0, 1, length = N)
fires <- stp(data.frame(x = clmfires$x, y = clmfires$y, t = time))

# Fit full model
lgcp_full <- stlgcppm(fires, ~ 1 + log(x) + poly(x, 5) + log(y) + poly(y, 5))

# Divide up data
indices <- sample(rep(1:nCores, each = N / nCores))

# Calculate volume for entire window
vol_full <- area.owin(clmfires$window)

# Calculate volume for smaller window (window 1)
r1 <- 30
window1 <- intersect.owin(disc(radius = r1, centre = c(300, 200)),
                          clmfires$window)
index1 <- (fires$df$x - 300)^2 + (fires$df$y - 200)^2 < r1^2
vol_window1 <- area.owin(window1)

# Calculate volume for larger window (window 2)
r2 <- 60
window2 <- intersect.owin(disc(radius = r2, centre = c(200, 250)),
                          clmfires$window)
index2 <- (fires$df$x - 200)^2 + (fires$df$y - 250)^2 < r2^2
vol_window2 <- area.owin(window2)

# For loop (will later change to parallel)
runtimes <- numeric(nCores)
agg_intens <- numeric(N)
for (i in 1:nCores) {
  index <- which(indices == i)
  subset <- stp(data.frame(x = clmfires$x[index], y = clmfires$y[index], t = time[index]))
  lgcp_subset <- stlgcppm(subset, ~ 1 + log(x) + poly(x, 5) + log(y) + poly(y, 5))
  agg_intens[indices == i] <- lgcp_subset$l
  runtimes[i] <- as.numeric(substr(lgcp_subset$time, 1, nchar(lgcp_subset$time) - 8))
}

# Predictions (agg model)
pred_agg <- c(mean(agg_intens[index1]) * vol_window1 * nCores, 
              mean(agg_intens[index2]) * vol_window2 * nCores, 
              mean(agg_intens) * vol_full * nCores)

# Predictions (full model)
pred_full <- c(mean(lgcp_full$l[index1]) * vol_window1, 
               mean(lgcp_full$l[index2]) * vol_window2, 
               mean(lgcp_full$l) * vol_full)

sum(index1)
sum(index2)
sum(runtimes)
max(runtimes)

plot(unmark(clmfires), main = "Castilla-La Mancha Forest Fires")
draw.circle(300, 200, r1, border = 'red', lwd = 4)
draw.circle(200, 250, r2, border = 'blue', lwd = 4)
mtext(paste0("Region 1: Red circle with radius ", r1),
      side = 1, line = 0)
mtext(paste0("Region 2: Blue circle with radius ", r2),
      side = 1, line = 1)
mtext(paste0("Region 3: Entire domain"),
      side = 1, line = 2)


results <- data.frame(
  Region = 1:3,
  Truth = c(sum(index1), sum(index2), N),
  Pred_Aggregated = pred_agg,
  Pred_Full = pred_full,
  Time_Aggregated = paste0(max(runtimes), " minutes"),
  Time_Full = lgcp_full$time
)
results
