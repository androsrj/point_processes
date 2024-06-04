# methods for coercion to Spatial Polygons by Adrian Baddeley

owin2Polygons <- function(x, id="1") {
	#check_spatstat("spatstat.geom")
  stopifnot(spatstat.geom::is.owin(x))
  x <- spatstat.geom::as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  #check_spatstat("spatstat.utils")
  pieces <- lapply(x$bdry,
      function(p) {Polygon(coords=closering(cbind(p$x,p$y)),
                   hole=spatstat.utils::is.hole.xypolygon(p))  })
  z <- Polygons(pieces, id)
  return(z)
}

as.SpatialPolygons.tess <- function(x) {
  #check_spatstat("spatstat.geom")
  stopifnot(spatstat.geom::is.tess(x))
  y <- spatstat.geom::tiles(x)
  nam <- names(y)
  z <- list()
  for(i in seq(y)) {
    zi <- try(owin2Polygons(y[[i]], nam[i]), silent=TRUE)
    if (inherits(zi, "try-error")) {
      warning(paste("tile", i, "defective\n", as.character(zi)))
    } else {
      z[[i]] <- zi
    }
  }
  return(SpatialPolygons(z))
}

#setAs("tess", "SpatialPolygons", function(from) as.SpatialPolygons.tess(from))


as.SpatialPolygons.owin <- function(x) {
  #check_spatstat("spatstat.geom")
  stopifnot(spatstat.geom::is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

setAs("owin", "SpatialPolygons", function(from) as.SpatialPolygons.owin(from))

getZmat <- function (formula, data, regionalcovariates = NULL, pixelcovariates = NULL, 
                     cellwidth, ext = 2, inclusion = "touching", overl = NULL) 
{
  if (!is.null(overl)) {
    cat("Using 'cellwidth' and 'ext' from overl\n")
    cellwidth <- overl$cellwidth
    ext <- overl$ext
  }
  data_sf = try(st_as_sf(data), silent = TRUE)
  if (inherits(data, "SpatialPolygonsDataFrame")) {
    spatstat.options(checkpolygons = FALSE)
    W <- as(as(st_union(data_sf), "Spatial"), "owin")
    spatstat.options(checkpolygons = TRUE)
    sd <- ppp(window = W)
  }
  else {
    sd <- data
  }
  ow <- selectObsWindow(sd, cellwidth)
  sd <- ow$xyt
  M <- ow$M
  N <- ow$N
  study.region <- sd$window
  if (!is.null(overl)) {
    gridobj <- overl$gridobj
  }
  else {
    gridobj <- genFFTgrid(study.region = study.region, M = M, 
                          N = N, ext = ext, inclusion = inclusion)
  }
  del1 <- gridobj$del1
  del2 <- gridobj$del2
  Mext <- gridobj$Mext
  Next <- gridobj$Next
  mcens <- gridobj$mcens
  ncens <- gridobj$ncens
  cellarea <- gridobj$cellarea
  cellInside <- gridobj$cellInside
  ans <- cov.interp.fft2(formula = formula, W = study.region, 
                         regionalcovariates = regionalcovariates, pixelcovariates = pixelcovariates, 
                         mcens = mcens[1:M], ncens = ncens[1:N], cellInside = cellInside[1:M, 1:N], overl = overl)
  attr(ans, "gridobj") <- gridobj
  attr(ans, "inclusion") <- inclusion
  attr(ans, "ext") <- ext
  attr(ans, "cellwidth") <- cellwidth
  class(ans) <- c("lgcpZmat", "matrix")
  return(ans)
}

cov.interp.fft2 <- function (formula, W, regionalcovariates = NULL, pixelcovariates = NULL, 
                             mcens, ncens, cellInside, overl = NULL) 
{
  if (!is.null(overl)) {
    if (!(isTRUE(all.equal(mcens, overl$mcens)) & isTRUE(all.equal(ncens, 
                                                                   overl$ncens)))) {
      stop("Differing FFT grids ... check this is the correct overlay")
    }
  }
  varn <- variablesinformula(formula)[-1]
  if (!is.null(regionalcovariates)) {
    if (any(is.na(getinterp(regionalcovariates@data)))) {
      stop("No interpolation method specified for one or more regional covariates, see ?assigninterp and ?guessinterp. Note that this should also be supplied for pixelcovariates, if applicable.")
    }
  }
  if (!is.null(pixelcovariates)) {
    if (any(is.na(getinterp(pixelcovariates@data)))) {
      stop("No interpolation method specified for one or more pixel covariates, see ?assigninterp and ?guessinterp")
    }
  }
  M <- length(mcens)
  N <- length(ncens)
  if (is.null(regionalcovariates) & is.null(pixelcovariates)) {
    testf <- as.formula(X ~ 1)
    attr(testf, ".Environment") <- .GlobalEnv
    if (identical(formula, testf)) {
      Zmat <- matrix(cellInside, M * N, 1)
      colnames(Zmat) <- "(Intercept)"
      attr(Zmat, "data.frame") <- data.frame(X = rep(1, 
                                                     sum(cellInside)))
      attr(Zmat, "cellInside") <- cellInside
      attr(Zmat, "anymiss") <- NULL
      attr(Zmat, "M") <- M
      attr(Zmat, "N") <- N
      attr(Zmat, "mcens") <- mcens
      attr(Zmat, "ncens") <- ncens
      attr(Zmat, "polygonOverlay") <- NULL
      attr(Zmat, "pixelOverlay") <- NULL
      attr(Zmat, "FORM") <- formula
      attr(Zmat, "fftpoly") <- NULL
      return(Zmat)
    }
    stop("Must have either regional or pixel covariates.")
  }
  spoly <- as.SpatialPolygons.owin(W)
  if (!is.null(overl)) {
    fftpoly <- overl$fftpoly
  }
  else {
    fftpoly <- grid2spoly(mcens, ncens)
  }
  Zmat <- 0
  dmat <- as.data.frame(matrix(NA, length(fftpoly), length(varn)))
  clsvec <- c()
  polygonoverlay <- NULL
  pixeloverlay <- NULL
  if (!is.null(regionalcovariates)) {
    cat("aggregating regional covariate information ...\n")
    s <- Sys.time()
    polyareas <- sapply(1:length(regionalcovariates), function(ii) {
      sum(sapply(slot(regionalcovariates[ii, ], "polygons"), 
                 slot, "area"))
    })
    if (!is.null(overl)) {
      cat("loading polygon overlay ...\n")
      polyol <- overl$polyol
    }
    else {
      cat("performing polygon overlay operations ...\n")
      polyol <- gOverlay(fftpoly, regionalcovariates)
    }
    cidx <- match(varn, names(regionalcovariates))
    if (any(is.na(cidx))) {
      stop("There is a mismatch between the names of variables in the supplied formula and the variables in the data frame.")
    }
    cidx <- cidx[!is.na(cidx)]
    gidx <- sort(unique(polyol$info$grididx))
    dtemp <- as.data.frame(matrix(NA, length(gidx), length(cidx)))
    regionalcovariates <- regionalcovariates@data[, cidx, 
                                                  drop = FALSE]
    classes <- sapply(regionalcovariates, class)
    cat("interpolating ...\n")
    sapply(1:length(gidx), function(i) {
      dtemp[i, 1:length(cidx)] <<- aggregateCovariateInfo(cellidx = i, 
                                                          cidx = cidx, gidx = gidx, df = regionalcovariates, 
                                                          fftovl = polyol, classes = classes, polyareas = polyareas)
    })
    dmat[gidx, 1:length(cidx)] <- dtemp
    e <- Sys.time()
    cat("Time Taken: ", difftime(e, s, units = "secs"), "\n")
    clsvec <- c(clsvec, sapply(clearinterp(regionalcovariates), 
                               class))
    polygonoverlay <- polyol
  }
  if (!is.null(pixelcovariates)) {
    cat("aggregating pixel covariate information ...\n")
    s <- Sys.time()
    cat("converting SpatialPixelsDataFrame to SpatialPolygonsDataFrame ...\n")
    oldclasses <- sapply(pixelcovariates@data, class)
    polyareas <- rep(prod(pixelcovariates@grid@cellsize), 
                     length(pixelcovariates))
    pixelcovariates <- as(pixelcovariates, "SpatialPolygonsDataFrame")
    sapply(1:ncol(pixelcovariates), function(i) {
      class(pixelcovariates[[i]]) <<- oldclasses[, i]
    })
    if (!is.null(overl)) {
      cat("loading polygon overlay ...\n")
      polyol <- overl$pixol
    }
    else {
      cat("performing polygon overlay operations ...\n")
      polyol <- gOverlay(fftpoly, pixelcovariates)
    }
    cidx <- match(varn, names(pixelcovariates))
    if (any(is.na(cidx))) {
      stop("There is a mismatch between the names of variables in the supplied formula and the variables in the data frame.")
    }
    cidx <- cidx[!is.na(cidx)]
    gidx <- sort(unique(polyol$info$grididx))
    dtemp <- as.data.frame(matrix(NA, length(gidx), length(cidx)))
    pixelcovariates <- pixelcovariates@data[, cidx, drop = FALSE]
    classes <- sapply(pixelcovariates, class)
    cat("interpolating ...\n")
    sapply(1:length(gidx), function(i) {
      dtemp[i, 1:length(cidx)] <<- aggregateCovariateInfo(cellidx = i, 
                                                          cidx = cidx, gidx = gidx, df = pixelcovariates, 
                                                          fftovl = polyol, classes = classes, polyareas = polyareas)
    })
    dmat[gidx, (length(varn) - length(cidx) + 1):length(varn)] <- dtemp
    e <- Sys.time()
    cat("Time Taken: ", difftime(e, s, units = "secs"), "\n")
    clsvec <- c(clsvec, sapply(clearinterp(pixelcovariates), 
                               class))
    pixeloverlay <- polyol
  }
  cellin <- which(as.logical(cellInside))
  dmat <- dmat[cellin, , drop = FALSE]
  anymiss <- apply(dmat, 1, function(x) {
    any(is.na(x))
  })
  anymisscopy <- anymiss
  if (any(anymiss)) {
    warning(paste("There is(are)", sum(anymiss), "missing value(s) in the covariate data. Imputing with median of non-missing values."), 
            immediate. = TRUE)
    for (iii in 1:ncol(dmat)) {
      if (any(is.na(dmat[, iii]))) {
        dmat[, iii][is.na(dmat[, iii])] <- median(dmat[, 
                                                       iii], na.rm = TRUE)
      }
    }
    anymiss <- apply(dmat, 1, function(x) {
      any(is.na(x))
    })
  }
  sapply(1:length(clsvec), function(i) {
    dmat[[i]] <<- eval(call(paste("as.", clsvec[i], sep = ""), 
                            dmat[[i]]))
  })
  names(dmat) <- varn
  dmat$X <- rep(1, sum(cellInside))
  dmat <- as.data.frame(dmat)
  DM <- dmat
  dmat <- model.matrix(formula, data = dmat)
  missingind <- rep(0, dim(dmat)[2])
  Zmat <- matrix(0, M * N, dim(dmat)[2])
  colnames(Zmat) <- colnames(dmat)
  rownames(dmat) <- NULL
  idx <- as.logical(cellInside)
  idx[as.logical(cellInside)] <- !anymiss
  Zmat[idx, ] <- dmat
  if (any(anymisscopy)) {
    missingind[idx] <- as.numeric(anymisscopy)
  }
  attr(Zmat, "data.frame") <- DM
  attr(Zmat, "cellInside") <- cellInside
  attr(Zmat, "anymiss") <- anymiss
  attr(Zmat, "mcens") <- mcens
  attr(Zmat, "ncens") <- ncens
  attr(Zmat, "M") <- M
  attr(Zmat, "N") <- N
  attr(Zmat, "polygonOverlay") <- polygonoverlay
  attr(Zmat, "pixelOverlay") <- pixeloverlay
  attr(Zmat, "FORM") <- formula
  attr(Zmat, "fftpoly") <- fftpoly
  attr(Zmat, "missingind") <- missingind
  return(Zmat)
}

