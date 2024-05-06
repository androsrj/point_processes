library(spatstat)
library(usmap)
library(ggplot2)

# Full dataset fit
fit1 <- ppm(lansing~x+y)
fit1 <- ppm(lansing ~ marks * polynom(x,y,3))
fit1

# Partition data
nCores <- 5
indices <- sample(1:nCores, size = lansing$n, replace = TRUE)
table(indices)
weights <- lansing$n / table(indices)

models <- vector("list", nCores)
for (i in 1:nCores) {
  disk <- subset(lansing, indices == i)
  models[[i]] <- ppm(disk~marks * polynom(x,y,3))
}
agg_preds <- Reduce("+", lapply(1:nCores, \(i) predict(models[[i]]) * weights[i] )) / nCores

# Point estimates
fit1$coef
apply(sapply(1:nCores, \(i) models[[i]]$coef), 1, mean)

# Confidence intervals
summary(fit1)$coefs.SE.CI[,3:4]
Reduce("+", lapply(1:nCores, \(i) summary(models[[i]])$coefs.SE.CI[,3:4])) / nCores

# Predictive plots
plot(predict(fit1), main = "Centralized Model (Lansing Woods Data)") # Whole model
plot(agg_preds, main = "Partitioned Model (Lansing Woods Data)") # Divided model

# Lightning data
nldn <- read.csv("../nldn/daily_data/2023/05/02.csv", stringsAsFactors = TRUE)
dim(nldn)
head(nldn)
bounds <- c(min(nldn$LON), max(nldn$LON), min(nldn$LAT), max(nldn$LAT))
nldn_reduced <- nldn[,c(11,10)]
nldn_pp <- as.ppp(nldn_reduced, bounds)
nldn_fit <- ppm(nldn_pp ~ x+y, Poisson())
nldn_fit

# Partition data
nCores <- 20
indices <- sample(1:nCores, size = nldn_pp$n, replace = TRUE)
table(indices)
weights <- nldn_pp$n / table(indices)

models <- vector("list", nCores)
for (i in 1:nCores) {
  disk <- subset(nldn_pp, indices == i)
  models[[i]] <- ppm(disk~x+y, Poisson())
}
preds_nldn <- Reduce("+", lapply(1:nCores, \(i) predict(models[[i]], type = "intensity") * weights[i])) / nCores

# Coefficients
nldn_fit$coef
apply(sapply(1:nCores, \(i) models[[i]]$coef), 1, mean)

# Confidence intervals
summary(nldn_fit)$coefs.SE.CI[,3:4]
Reduce("+", lapply(1:nCores, \(i) summary(models[[i]])$coefs.SE.CI[,3:4])) / nCores

# Add prediction surface plots here 
plot(predict(nldn_fit, type="intensity"), main = "Centralized Model (NLDN Data)") # whole model
plot(preds_nldn, main = "Partitioned Model (NLDN Data)") # Divided model

# Actual lightning strike data overlayed on US map
nldn_df <- data.frame(lon = nldn_pp$x, lat = nldn_pp$y)
nldn_map <- usmap_transform(nldn_df)
plot_usmap(exclude = c("HI", "AK")) + 
  geom_sf(data = nldn_map, alpha = 0.5) + 
  labs(title = "Actual Lightning Strikes (NLDN Data)") + 
  theme(plot.title = element_text(size = 20, hjust = 0.5))

