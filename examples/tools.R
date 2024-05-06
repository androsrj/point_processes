library(spatstat)
library(usmap)
library(ggplot2)

# Full dataset fit
fit1 <- ppm(unmark(stonetools)~x+y)
fit1

# Partition data
nCores <- 20
indices <- sample(1:nCores, size = stonetools$n, replace = TRUE)
table(indices)
weights <- stonetools$n / table(indices)

models <- vector("list", nCores)
for (i in 1:nCores) {
  disk <- subset(unmark(stonetools), indices == i)
  models[[i]] <- ppm(disk~x+y)
}
agg_preds <- Reduce("+", lapply(1:nCores, \(i) predict(models[[i]]) * weights[i] )) / nCores

# Point estimates
fit1$coef
apply(sapply(1:nCores, \(i) models[[i]]$coef), 1, mean)

# Confidence intervals
summary(fit1)$coefs.SE.CI[,3:4]
Reduce("+", lapply(1:nCores, \(i) summary(models[[i]])$coefs.SE.CI[,3:4])) / nCores

# Predictive plots
pdf("tools.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
plot(predict(fit1), main = "Centralized Model (Stone Tools Data)") # Whole model
plot(agg_preds, main = "Partitioned Model (Stone Tools Data)") # Divided model
#plot(stonetools$x, stonetools$y)
dev.off()
