library(stopp)
library(spatstat)
library(dplyr)
library(plotrix)

nCores <- 8

N <- nrow(greececatalog$df)

# Fit full model
lgcp_full <- stlgcppm(greececatalog, ~ 1 + log(x) + log(y) + log(t))

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
coefs <- matrix(0, nrow = nCores, ncol = 4)
for (i in 1:nCores) {
  index <- which(indices == i)
  subset <- stp(data.frame(x = greececatalog$df$x[index], 
                           y = greececatalog$df$y[index], 
                           t = greececatalog$df$t[index]))
  lgcp_subset <- stlgcppm(subset, ~ 1 + log(x) + log(y) + log(t))
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

pdf("figures/greece.pdf", width = 9, height = 6)
plot(greececatalog$df$x, greececatalog$df$y, asp=1,
     xlab = "", ylab = "", 
     main = "Greece Earthquakes")
draw.circle(22, 39, r1, border = 'red', lwd = 4)
draw.circle(22, 36, r2, border = 'blue', lwd = 4)
draw.circle(26, 36, r3, border = 'forestgreen', lwd = 4)
legend("topright", 
       fill = c("red", "blue", "forestgreen"), 
       legend = c(paste0(" Region 1 \n r = ", r1),
                  paste0(" Region 2 \n r = ", r2),
                  paste0(" Region 3 \n r = ", r3)),
       y.intersp = 2)
dev.off()

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
