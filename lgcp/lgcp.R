library(lgcp)
library(spatstat)
library(spatstat.geom)
library(spatstat.utils)
library(sp)

load("sd_liver.RData")
load("popshape_liver.RData")
source("polygons.R")

set.seed(12)
#grp1 <- sample(1:sd$n, 100)
#sd <- sd[grp1]

#minimum.contrast(sd, model = "exponential", method = "g", intens = density(sd), transform = log)
#chooseCellwidth(sd, cwinit = 300)
CELLWIDTH <- 3000
EXT <- 2
polyolay <- getpolyol(data = sd, regionalcovariates = popshape,
                      cellwidth = CELLWIDTH, ext = EXT)

popshape@data <- guessinterp(popshape@data)
popshape@data <- assigninterp(df = popshape@data,
                              vars = c("pop", "males", "females"), 
                              value = "ArealWeightedSum")

FORM <- X ~ pop + propmale + Income + Employment + Education +
  Barriers + Crime + Environment

Zmat <- getZmat(formula = FORM, data = sd, regionalcovariates = popshape,
                cellwidth = CELLWIDTH, ext = EXT, overl = polyolay)
Zmat[, "pop"] <- log(Zmat[, "pop"])
Zmat[, "pop"][is.infinite(Zmat[, "pop"])] <- min(
  Zmat[, "pop"][!is.infinite(Zmat[, "pop"])])

priors <- lgcpPrior(etaprior = PriorSpec(
  LogGaussianPrior(mean = log(c(1, 500)), variance = diag(0.15, 2))),
  betaprior = PriorSpec(
    GaussianPrior(mean = rep(0, 9), variance = diag(10^6, 9))))

INITS <- lgcpInits(etainit = log(c(sqrt(1.5), 275)), betainit = NULL)
cf <- CovFunction(exponentialCovFct)

mcmc.control <- mcmcpars(
  mala.length = 25000, burnin = 5000, retain = 10,
  adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1,
                                 targetacceptance = 0.574))

lg <- lgcpPredictSpatialPlusPars(formula = FORM, sd = sd, Zmat = Zmat,
                                 model.priors = priors, model.inits = INITS, spatial.covmodel = cf,
                                 cellwidth = CELLWIDTH, poisson.offset = NULL, mcmc.control = mcmc.control,
                                 output.control = setoutput(gridfunction = dump2dir(
                                   dirname = "liver/", forceSave = TRUE)),
                                 ext = EXT)
saveRDS(lg, file = "liver/lg.RDS")
plot(lg, ask=FALSE)
plot(sd)
lg$EY.mean$grid[[1]]
image(lg$mcens, lg$ncens, lg$EY.mean$grid[[1]])
image(lg$mcens, lg$ncens, getCounts(sd, M=16, N=8, ext=1))

sum(lg$EY.mean$grid[[1]])
sum(lg$EY.mean$grid[[1]] >= 1)
