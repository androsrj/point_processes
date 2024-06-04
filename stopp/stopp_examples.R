library(stopp)
library(spatstat)
library(dplyr)

catsub <- stp(greececatalog$df[1:200,])
lgcp1 <- stlgcppm(catsub, ~ 1 + poly(x, 3) + y + t)
lgcp1

vol <- diff(range(greececatalog$df[1:200, 1])) * 
  diff(range(greececatalog$df[1:200, 2])) * 
  diff(range(greececatalog$df[1:200, 3]))
pred <- mean(lgcp1$l) * vol
pred


lgcp2 <- stlgcppm(greececatalog, ~ 1 + poly(x, 3) + y + t)
lgcp2

vol <- diff(range(greececatalog$df[, 1])) * 
  diff(range(greececatalog$df[, 2])) * 
  diff(range(greececatalog$df[, 3]))
pred <- mean(lgcp2$l) * vol
pred

fires <- stp(data.frame(x = clmfires$x, y = clmfires$y, t = seq(0, 1, length = clmfires$n)))
lgcp3 <- stlgcppm(fires, ~ 1 + poly(x, 3) + y)
lgcp3

vol <- diff(range(fires$df[, 1])) * 
  diff(range(fires$df[, 2])) * 
  diff(range(fires$df[, 3]))
pred <- mean(lgcp3$l) * vol
pred


index <- fires$df$x > 250 & fires$df$y < 150 & fires$df$t > 0.25
sum(index)
newdf <- fires$df[index, ]
vol <- diff(range(newdf[, 1])) * 
  diff(range(newdf[, 2])) * 
  diff(range(newdf[, 3]))
pred <- mean(lgcp3$l[index]) * vol
pred


### Preliminary simulation: Randomly divide clmfires into 2 datasets
### Fit 2 separate models and predict occurrences for each
### Then add predictions together
n <- clmfires$n
time <- seq(0, 1, length = n)
index <- sample(1:n, n/2)
fires1 <- stp(data.frame(x = clmfires$x[index], y = clmfires$y[index], t = time[index]))
fires2 <- stp(data.frame(x = clmfires$x[-index], y = clmfires$y[-index], t = time[-index]))
lgcp_subset1 <- stlgcppm(fires1, ~ 1 + poly(x, 3) + y)
lgcp_subset2 <- stlgcppm(fires2, ~ 1 + poly(x, 3) + y)
lgcp_subset1$IntCoefs
lgcp_subset2$IntCoefs
lgcp3$IntCoefs

# Predict for each model (entire window)
vol <- diff(range(fires$df[, 1])) * 
  diff(range(fires$df[, 2])) * 
  diff(range(fires$df[, 3]))
pred1 <- mean(lgcp_subset1$l) * vol
pred2 <- mean(lgcp_subset2$l) * vol
pred1
pred2
pred1 + pred2

# Predict for each model (smaller window)
new_index <- fires$df$x > 250 & fires$df$y < 150 & fires$df$t > 0.25
sum(new_index)
newdf <- fires$df[new_index, ]
vol <- diff(range(newdf[, 1])) * 
  diff(range(newdf[, 2])) * 
  diff(range(newdf[, 3]))
pred1 <- mean(lgcp_subset1$l) * vol
pred2 <- mean(lgcp_subset2$l) * vol
pred1
pred2
pred1 + pred2


