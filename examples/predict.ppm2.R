predict.ppm2 <- function (object, nCores, window = NULL, ngrid = NULL, locations = NULL, 
          covariates = NULL, type = c("trend", "cif", "intensity", 
                                      "count"), se = FALSE, interval = c("none", "confidence", 
                                                                         "prediction"), level = 0.95, X = data.ppm(object), correction, 
          ignore.hardcore = FALSE, ..., dimyx = NULL, eps = NULL, rule.eps = c("adjust.eps", 
                                                                               "grow.frame", "shrink.frame"), new.coef = NULL, check = TRUE, 
          repair = TRUE) 
{
  interval <- match.arg(interval)
  rule.eps <- match.arg(rule.eps)
  #xarg <- xtract(...)
  #sumobj <- xarg$sumobj
  #E <- xarg$E
  #total <- xarg$total
  #getoutofjail <- xarg$getoutofjail
  seonly <- FALSE
  if (missing(type)) 
    type <- type[1]
  #else {
  #  if (length(type) > 1) 
  #    stop("Argument 'type' should be a single value")
  #  mt <- pmatch(type, typeaccept)
  #  if (is.na(mt)) 
  #    stop("Argument 'type' should be one of", commasep(sQuote(typepublic), 
  #                                                      " or "))
  #  type <- typeuse[mt]
  #  if (type == "se") {
  #    if (!getoutofjail) 
  #      message(paste("Outdated syntax:", "type='se' should be replaced by se=TRUE;", 
  #                    "then the standard error is predict(...)$se"))
  #    type <- "trend"
  #    se <- TRUE
  #    seonly <- TRUE
  #  }
  #}
  
  #if (!is.null(total)) {
  #  message("Outdated argument 'total': use 'window' and set type='count'")
  #  type <- "count"
  #  if (!is.logical(total)) 
  #    window <- if (is.tess(total)) 
  #      total
  #  else as.owin(total)
  #}
  model <- object
  verifyclass(model, "ppm")
  if (check && damaged.ppm(object)) {
    if (!repair) 
      stop("object format corrupted; try update(object, use.internal=TRUE)")
    message("object format corrupted; repairing it.")
    object <- update(object, use.internal = TRUE)
  }
  if (missing(correction) || is.null(correction)) 
    correction <- object$correction
  fitcoef <- coef(object)
  if (!is.null(new.coef)) {
    if (length(new.coef) != length(fitcoef)) 
      stop(paste("Argument new.coef has wrong length", 
                 length(new.coef), ": should be", length(fitcoef)))
    coeffs <- new.coef
  }
  else {
    coeffs <- fitcoef
  }
  sumobj <- summary(model, quick = "entries")
  poisson <- sumobj$poisson
  marked <- sumobj$marked
  multitype <- sumobj$multitype
  notrend <- sumobj$no.trend
  changedcoef <- sumobj$changedcoef || !is.null(new.coef)
  trivial <- poisson && notrend
  need.covariates <- sumobj$uses.covars
  covnames.needed <- sumobj$covars.used
  if (sumobj$antiquated) 
    warning("The model was fitted by an out-of-date version of spatstat")
  if (marked) {
    if (!multitype) 
      stop("Prediction not yet implemented for general marked point processes")
    else types <- levels(marks(sumobj$entries$data))
  }
  if (poisson && type %in% c("cif", "intensity")) 
    type <- "trend"
  if (se && type %in% c("cif", "intensity")) 
    stop(paste("Standard error for", type, "is not yet implemented"), 
         call. = FALSE)
  if (type == "count" && interval != "none" && (marked || !poisson)) {
    stop(paste0(interval, " intervals for counts are only implemented for", 
                if (marked) 
                  " unmarked"
                else "", if (!poisson) 
                  " Poisson", " models"), call. = FALSE)
  }
  if (interval == "prediction" && type != "count") 
    stop("Prediction intervals are only available for type='count'", 
         call. = FALSE)
  if (interval == "confidence" && type %in% c("intensity", 
                                              "cif")) 
    stop(paste("Confidence intervals are not yet available for", 
               type), call. = FALSE)
  estimatename <- if (interval == "none") 
    "estimate"
  else interval
  if (type == "count") {
    if (is.null(window)) {
      if (!seonly) 
        est <- predconfPois2(NULL, model, level, estimatename, nCores = nCores,
                            new.coef = new.coef)
      if (se) 
        sem <- predconfPois2(NULL, model, level, "se", nCores = nCores,
                            new.coef = new.coef)
    }
    else if (is.tess(window)) {
      tilz <- tiles(window)
      if (!seonly) {
        est <- lapply(tilz, predconfPois2, object = model, nCores = nCores,
                      level = level, what = estimatename, new.coef = new.coef)
        est <- switch(interval, none = unlist(est), confidence = , 
                      prediction = t(simplify2array(est)))
      }
      if (se) 
        sem <- sapply(tilz, predconfPois2, object = model, nCores = nCores,
                      level = level, what = "se", new.coef = new.coef)
    }
    else {
      if (!seonly) 
        est <- predconfPois2(window, model, level, estimatename, nCores = nCores,
                            new.coef = new.coef)
      if (se) 
        sem <- predconfPois2(window, model, level, "se", nCores = nCores,
                            new.coef = new.coef)
    }
    if (!se) 
      return(est)
    if (seonly) 
      return(sem)
    result <- list(est, sem)
    names(result) <- c(estimatename, "se")
    return(result)
  }
  if (interval != "none") {
    alpha2 <- (1 - level)/2
    pp <- sort(c(alpha2, 1 - alpha2))
    ci.names <- paste0(signif(100 * pp, 3), "%")
    ci.q <- qnorm(pp)
  }
  if (is.im(locations)) 
    locations <- as.owin(locations)
  if (is.null(window) && is.owin(locations) && !is.mask(locations)) {
    window <- locations
    locations <- NULL
  }
  if (!is.null(locations)) {
    offending <- c(!is.null(ngrid), !is.null(dimyx), !is.null(eps))
    if (any(offending)) {
      offenders <- c("grid", "dimyx", "eps")[offending]
      nbad <- sum(offending)
      stop(paste(ngettext(nbad, "The argument", "The arguments"), 
                 commasep(sQuote(offenders)), ngettext(nbad, "is", 
                                                       "are"), "incompatible with", sQuote("locations")), 
           call. = FALSE)
    }
  }
  if (!is.null(ngrid) && !is.null(dimyx)) 
    warning(paste("The arguments", sQuote("ngrid"), "and", 
                  sQuote("dimyx"), "are equivalent: only one should be given"), 
            call. = FALSE)
  ngrid <- ngrid %orifnull% dimyx
  if (is.null(ngrid) && is.null(locations)) 
    ngrid <- rev(spatstat.options("npixel"))
  want.image <- is.null(locations) || is.mask(locations)
  make.grid <- !is.null(ngrid)
  if (!want.image) {
    xpredict <- locations$x
    ypredict <- locations$y
    if (is.null(xpredict) || is.null(ypredict)) {
      xy <- xy.coords(locations)
      xpredict <- xy$x
      xpredict <- xy$y
    }
    if (is.null(xpredict) || is.null(ypredict)) 
      stop(paste("Don't know how to extract x,y coordinates from", 
                 sQuote("locations")))
    if (marked) {
      mpredict <- locations$marks
      if (is.null(mpredict)) 
        stop(paste("The argument", sQuote("locations"), 
                   "does not contain a column of marks", "(required since the fitted model", 
                   "is a marked point process)"))
      if (is.factor(mpredict)) {
        if (!isTRUE(all.equal(levels(mpredict), types))) {
          if (all(levels(mpredict) %in% types)) 
            mpredict <- factor(mpredict, levels = types)
          else stop(paste("The marks in", sQuote("locations"), 
                          "do not have the same levels as", "the marks in the model"))
        }
      }
      else {
        if (all(mpredict %in% types)) 
          mpredict <- factor(mpredict, levels = types)
        else stop(paste("The marks in", sQuote("locations"), 
                        "do not have the same values as the marks in the model"))
      }
    }
  }
  else {
    if (!make.grid) 
      masque <- locations
    else {
      if (!is.null(ngrid)) {
        if (!is.numeric(ngrid)) 
          stop("ngrid should be a numeric vector")
        ngrid <- ensure2vector(ngrid)
      }
      if (is.null(window)) 
        window <- sumobj$entries$data$window
      masque <- as.mask(window, dimyx = ngrid, eps = eps, 
                        rule.eps = rule.eps)
    }
    tums <- termsinformula(model$trend)
    if (any(tums == "lo(x)" | tums == "lo(y)" | tums == "lo(x,y)" | 
            tums == "lo(y,x)")) {
      gg <- model$internal$glmdata
      gxr <- range(gg$x[gg$SUBSET])
      gyr <- range(gg$y[gg$SUBSET])
      masque <- intersect.owin(masque, owin(gxr, gyr))
    }
    rxy <- rasterxy.mask(masque, drop = TRUE)
    xpredict <- rxy$x
    ypredict <- rxy$y
  }
  if (!marked) 
    newdata <- data.frame(x = xpredict, y = ypredict)
  else if (!want.image) 
    newdata <- data.frame(x = xpredict, y = ypredict, marks = mpredict)
  else {
    nt <- length(types)
    np <- length(xpredict)
    xpredict <- rep.int(xpredict, nt)
    ypredict <- rep.int(ypredict, nt)
    mpredict <- rep.int(types, rep.int(np, nt))
    mpredict <- factor(mpredict, levels = types)
    newdata <- data.frame(x = xpredict, y = ypredict, marks = mpredict)
  }
  if (need.covariates) {
    if (is.null(covariates)) {
      oldcov <- model$covariates
      if (is.null(oldcov)) 
        stop("External covariates are required, and are not available")
      if (is.data.frame(oldcov)) 
        stop(paste("External covariates are required.", 
                   "Prediction is not possible at new locations"))
      covariates <- oldcov
    }
    covariates <- if (is.data.frame(covariates)) {
      covariates[, covnames.needed, drop = FALSE]
    }
    else covariates[covnames.needed]
    covfunargs <- model$covfunargs
    covariates.df <- mpl.get.covariates(covariates, list(x = xpredict, 
                                                         y = ypredict), "prediction points", covfunargs)
    newdata <- cbind(newdata, covariates.df)
  }
  if (is.null(newdata$SUBSET)) 
    newdata$SUBSET <- rep.int(TRUE, nrow(newdata))
  if (!trivial) {
    Vnames <- model$internal$Vnames
    vnameprefix <- model$internal$vnameprefix
    glmdata <- getglmdata(model)
    glmfit <- getglmfit(model)
    if (object$method == "logi") 
      newdata$.logi.B <- rep(glmdata$.logi.B[1], nrow(newdata))
  }
  if (type == "covariates") 
    return(list(newdata = newdata, mask = if (want.image) masque else NULL))
  needSE <- se || (interval != "none")
  attribeauts <- list()
  if (trivial) {
    lambda <- exp(coeffs[[1]])
    if (needSE) {
      npts <- nobs(model)
      se.lambda <- lambda/sqrt(npts)
    }
    switch(interval, none = {
      z <- rep.int(lambda, nrow(newdata))
    }, confidence = {
      z <- matrix(lambda + se.lambda * ci.q, byrow = TRUE, 
                  nrow = nrow(newdata), ncol = 2, dimnames = list(NULL, 
                                                                  ci.names))
    }, stop("Internal error: unreached"))
    if (se) 
      zse <- rep.int(se.lambda, nrow(newdata))
  }
  else if ((type %in% c("trend", "intensity")) || poisson) {
    zeroes <- numeric(nrow(newdata))
    for (vn in Vnames) newdata[[vn]] <- zeroes
    z <- lambda <- GLMpredict(glmfit, newdata, coeffs, changecoef = changedcoef) * nCores
    if (type == "intensity") 
      z <- PoisSaddle(z, fitin(model))
    if (needSE) {
      if (inherits(glmfit, "gam")) {
        if (!is.null(new.coef)) 
          warning("new.coef ignored in standard error calculation")
        SE <- predict(glmfit, newdata = newdata, type = "response", 
                      se.fit = TRUE)[[2]]
      }
      else {
        vc <- vcov(model, new.coef = new.coef)
        fmla <- rhs.of.formula(formula(glmfit))
        mf <- model.frame(fmla, newdata, na.action = na.pass)
        mm <- model.matrix(fmla, mf, na.action = na.pass)
        if (nrow(mm) != nrow(newdata)) 
          stop("Internal error: row mismatch in SE calculation")
        if (ncol(mm) != ncol(vc)) 
          stop("Internal error: column mismatch in SE calculation")
        vv <- quadform(mm, vc)
        SE <- lambda * sqrt(vv) / sqrt(nCores)
      }
      if (se) 
        zse <- SE
      if (interval == "confidence") {
        z <- lambda + outer(SE, ci.q, "*")
        colnames(z) <- ci.names
      }
    }
  }
  else if (type == "cif" || type == "lambda") {
    inter <- model$interaction
    if (!missing(X)) 
      stopifnot(is.ppp(X))
    W <- as.owin(data.ppm(model))
    U <- ppp(newdata$x, y = newdata$y, window = W, check = FALSE)
    if (marked) 
      marks(U) <- newdata$marks
    if (is.null(E)) 
      E <- equalpairs(U, X, marked)
    Vnew <- evalInteraction(X, U, E, inter, correction = correction, 
                            splitInf = ignore.hardcore, check = check)
    if (!ignore.hardcore) {
      cif.equals.zero <- matrowany(Vnew == -Inf)
    }
    else {
      cif.equals.zero <- attr(Vnew, "-Inf") %orifnull% 
        logical(nrow(Vnew))
    }
    attribeauts <- c(attribeauts, list(isZero = cif.equals.zero))
    if (ncol(Vnew) == 1) {
      newdata[[Vnames]] <- as.vector(Vnew)
    }
    else if (is.null(avail <- colnames(Vnew))) {
      for (i in seq_along(Vnames)) newdata[[Vnames[i]]] <- Vnew[, i]
    }
    else {
      if (all(Vnames %in% avail)) {
        for (vn in Vnames) newdata[[vn]] <- Vnew[, vn]
      }
      else if (all(Vnames %in% (Pavail <- paste0(vnameprefix, 
                                                 avail)))) {
        for (vn in Vnames) newdata[[vn]] <- Vnew[, match(vn, 
                                                         Pavail)]
      }
      else stop(paste("Internal error: unable to match names", 
                      "of available interaction terms", commasep(sQuote(avail)), 
                      "to required interaction terms", commasep(sQuote(Vnames))), 
                call. = FALSE)
    }
    z <- GLMpredict(glmfit, newdata, coeffs, changecoef = changedcoef)
    if (!ignore.hardcore && any(cif.equals.zero)) 
      z[cif.equals.zero] <- 0
  }
  else stop(paste("Unrecognised type", sQuote(type)))
  if (!want.image) {
    if (!se) {
      z <- as.vector(z)
      attributes(z) <- c(attributes(z), attribeauts)
      out <- z
    }
    else if (seonly) {
      out <- as.vector(zse)
    }
    else {
      z <- as.vector(z)
      attributes(z) <- c(attributes(z), attribeauts)
      out <- list(z, as.vector(zse))
      names(out) <- c(estimatename, "se")
    }
  }
  else {
    imago <- as.im(masque, value = 1)
    if (!marked && interval == "none") {
      if (!se) {
        out <- imago
        out[] <- z
      }
      else if (seonly) {
        out <- imago
        out[] <- zse
      }
      else {
        est <- std <- imago
        est[] <- z
        std[] <- zse
        out <- list(est, std)
        names(out) <- c(estimatename, "se")
      }
    }
    else if (interval != "none") {
      if (!seonly) {
        hi <- lo <- imago
        hi[] <- z[, 1]
        lo[] <- z[, 2]
        est <- solist(hi, lo)
        names(est) <- ci.names
      }
      if (se) {
        std <- imago
        std[] <- zse
      }
      if (!se) {
        out <- est
      }
      else if (seonly) {
        out <- std
      }
      else {
        out <- list(est, std)
        names(out) <- c(estimatename, "se")
      }
    }
    else {
      out <- list()
      for (i in seq_along(types)) {
        outi <- imago
        outi[] <- z[newdata$marks == types[i]]
        out[[i]] <- outi
      }
      out <- as.solist(out)
      names(out) <- as.character(types)
    }
  }
  return(out)
}
