# libraries
library(lgcp)
library(caret)
library(spatstat)
library(spatstat.utils)
library(spatstat.geom)
library(sp)

# load data
load("btb/BTBppp.RData")
load("btb/farmspdf.RData")

# find highly correlated variables
d <- farmspdf@data[, 30:40]
d <- d[!is.na(d[, 1]), ]
findCorrelation(cor(d))

# minimum contrast estimation
W <- pppdata$window
simpW <- simplify.owin(W, 1000)
ttx <- pppdata$x[pppdata$marks == "x1" | pppdata$marks == "x4" | pppdata$marks == 
                   "x5" | pppdata$marks == "x7"]
tty <- pppdata$y[pppdata$marks == "x1" | pppdata$marks == "x4" | pppdata$marks == 
                   "x5" | pppdata$marks == "x7"]
ttm <- as.character(pppdata$marks[pppdata$marks == "x1" | pppdata$marks == "x4" | 
                                    pppdata$marks == "x5" | pppdata$marks == "x7"])
tempppp <- ppp(x = ttx, y = tty, window = simpW, marks = as.factor(ttm))
denls <- lapply(as.character(levels(tempppp$marks)), function(x) {
  density.ppp(tempppp[tempppp$marks == x])
})
mn <- minimum.contrast(tempppp, model = "exponential", method = "g", intens = denls, 
                       transform = log)

# select cell width
chooseCellwidth(pppdata, cwinit = 1800)
CELLWIDTH <- 3600

# select extension
EXT <- 2

# perform polygon overlay operations and compute computational grid
polyolay <- getpolyol(data = pppdata, pixelcovariates = farmspdf, cellwidth = CELLWIDTH, 
                      ext = EXT)

# the statistical model for the main effets, $\beta$
fl <- formulaList(list(x1 ~ K207 + K208 + K209 + K210, x4 ~ K207 + K208 + K209 + 
                         K210, x5 ~ K207 + K208 + K209 + K210, x7 ~ K207 + K208 + K209 + K210))

# specify formula for purposes of interpolation
FORM <- X ~ K207 + K208 + K209 + K210

# interpolation has already been set in farmspdf.RData

# perform interpolation of the covariates onto the
Zmat <- getZmat(formula = FORM, data = pppdata, pixelcovariates = farmspdf, cellwidth = CELLWIDTH, 
                ext = EXT, overl = polyolay)

Zmat[, "K207"] <- log(Zmat[, "K207"])
Zmat[, "K207"][is.infinite(Zmat[, "K207"])] <- min(Zmat[, "K207"][!is.infinite(Zmat[, 
                                                                                    "K207"])])

Zmat[, "K208"] <- log(Zmat[, "K208"])
Zmat[, "K208"][is.infinite(Zmat[, "K208"])] <- min(Zmat[, "K208"][!is.infinite(Zmat[, 
                                                                                    "K208"])])

Zmat[, "K209"] <- log(Zmat[, "K209"])
Zmat[, "K209"][is.infinite(Zmat[, "K209"])] <- min(Zmat[, "K209"][!is.infinite(Zmat[, 
                                                                                    "K209"])])

Zmat[, "K210"] <- log(Zmat[, "K210"])
Zmat[, "K210"][is.infinite(Zmat[, "K210"])] <- min(Zmat[, "K210"][!is.infinite(Zmat[, 
                                                                                    "K210"])])

# plot the interpolated covariates
plot(Zmat, ask = FALSE)

# define the priors
pr.mn <- log(c(1, 1500))
pr.vr <- c(0.2, 0.05)
priors <- list()
for (i in 1:4) {
  priors[[i]] <- lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = pr.mn, 
                                                                 variance = diag(pr.vr))), betaprior = PriorSpec(GaussianPrior(mean = rep(0, 
                                                                                                                                          5), variance = diag(10^6, 5))))
}
priors[[5]] <- lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = pr.mn, variance = diag(pr.vr))), 
                         betaprior = NULL)

# set initial values were not set for this example

# choose the covariance function
cfs <- list()
for (i in 1:4) {
  cfs[[i]] <- CovFunction(exponentialCovFct)
}
cfs[[5]] <- CovFunction(RandomFieldsCovFct(model = "Matern32", additionalparameters = 1))

# run the MCMC algorithm, note that as this takes a long time, we present a
# shorter run here in the article ,we used mala.length=1000000, burnin=100000,
# thin=900
BASEDR <- getwd()
lg <- lgcpPredictMultitypeSpatialPlusPars(formulaList = fl, sd = pppdata, Zmat = Zmat, 
                                          model.priorsList = priors, spatial.covmodelList = cfs, cellwidth = CELLWIDTH, 
                                          mcmc.control = mcmcpars(mala.length = 10000, burnin = 1000, retain = 10, adaptivescheme = andrieuthomsh(inith = 1, 
                                                                                                                                               alpha = 0.5, C = 1, targetacceptance = 0.574)), output.control = setoutput(gridfunction = dump2dir(dirname = file.path(BASEDR, 
                                                                                                                                                                                                                                                                      "BTB"), forceSave = TRUE)), ext = EXT)
save(list = ls(), file = file.path(BASEDR, "BTB", "BTB.RData"))

plot(lg)
pppdata$n
trueCounts <- getCounts(pppdata, M=32, N=32, ext=1)
predCounts <- lg$EY.mean$grid[[1]]
sum(trueCounts)
sum(predCounts)
sum(abs(trueCounts - predCounts))

truePos <- trueCounts >= 1
predPos <- predCounts >= 1
mean(predPos == truePos)

sum(lg$cellInside)
table(lg$cellInside)
image(lg$mcens, lg$ncens, lg$cellInside)

