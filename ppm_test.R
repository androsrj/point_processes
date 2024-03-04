library(spatstat.model)

# fit the stationary Poisson process to point pattern 'nztrees'
ppm(nztrees ~ 1)

# fit the nonstationary Poisson process with intensity function lambda(x,y) = exp(a + bx)
# where x,y are the Cartesian coordinates and a,b are parameters to be estimated
fit1 <- ppm(nztrees ~ x)
fit1
coef(fit1)
coef(summary(fit1))

# Marked example with toy data
n <- 100
xy <- matrix(runif(2*n), ncol=2)
z <- as.factor(sample(c("A", "B", "C"), n, replace=T))
data <- data.frame(xy, z)
pp <- as.ppp(data, c(0,1,0,1))
fit2 <- ppm(pp ~ marks, Poisson())
fit2

# Try it on lightning data for just one date. Use polarity as the mark
nldn <- read.csv("../nldn/daily_data/2023/01/01.csv", stringsAsFactors = TRUE)
dim(nldn)
head(nldn)
bounds <- c(min(nldn$LON), max(nldn$LON), min(nldn$LAT), max(nldn$LAT))
nldn_reduced <- nldn[,c(11,10,3)]
nldn_pp <- as.ppp(nldn_reduced, bounds)
nldn_fit <- ppm(nldn_pp ~ marks, Poisson())
nldn_fit
