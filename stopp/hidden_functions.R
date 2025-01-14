.counting.weights <- function(id, volumes, lsr = FALSE) {
  id <- as.integer(id)
  if(lsr){
    w <- volumes[id] 
  } else {
    fid <- factor(id, levels = seq_along(volumes))
    counts <<- table(fid)
    w <- volumes[id] / counts[id]
  }
  w <- as.vector(w)
  names(w) <- NULL
  return(w)
}

.default.ncube <- function(X){
  guess.ngrid <- floor((splancs::npts(X) / 2) ^ (1 / 3))
  max(5, guess.ngrid)
}

density.new.var  <-function(x,np,gx){
  nx  <-  length(x)
  n <- length(gx)
  h <- numeric(nx)
  mat.d <- as.matrix(dist((x), upper = TRUE, diag = TRUE))
  diag(mat.d) <- Inf 
  if(missing(np))  np <- trunc(nx / 2)
  
  
  for(i in 1:ncol(mat.d)){
    o <- order(mat.d[,i])
    pp <- o[np]
    h[i] <- mat.d[pp, i]
  }
  h <- h / 4
  
  ax <- outer(gx, x, "-") / h
  z <- rowSums(matrix(dnorm(ax), n, nx) / (nx * h))
  return(list(x = gx, z = z, h = h))
}

g_st <- function(param, useq, vseq, ghat, transform, power){
  uv <- expand.grid(useq, vseq)
  g.parametric <- exp(
    param[1] / ((uv[, 2] / param[3]) ^ param[5] + 1) ^ (param[6] / param[5]) *
      exp( - (uv[, 1] / param[2]) ^ param[4] / ((uv[, 2] / param[3]) ^ param[5] + 1) ^ (param[6] / param[5] * param[4] / 2))
  )
  ghat2 <- as.vector(ghat)
  return(sum((((ghat2)) ^ power - (g.parametric) ^ power) ^ 2))
}

g_st_iaco <- function (param, useq, vseq,ghat, transform, power){
  uv<-expand.grid(useq,vseq)
  g.parametric <- exp((as.numeric(param[1])*
                         (1+(uv[,1]/as.numeric(param[2]))^as.numeric(param[4])+
                            (uv[,2]/as.numeric(param[3]))^as.numeric(param[5]))^(-as.numeric(param[6]))
  ))
  ghat2 <- as.vector(ghat)
  return(sum((((ghat2)) ^ power - (g.parametric) ^ power) ^ 2))
}

.grid1.index <- function(x, xrange, nx) {
  i <- ceiling(nx * (x - xrange[1]) / diff(xrange))
  i <- pmax.int(1, i)
  i <- pmin.int(i, nx)
  i
}

.grid.index <- function(x, y, t, xrange, yrange, trange, nx, ny, nt) {
  
  ix <- .grid1.index(x, xrange, nx)
  iy <- .grid1.index(y, yrange, ny)
  it <- .grid1.index(t, trange, nt)
  
  return(list(ix = ix, iy = iy, it = it, index = as.integer((iy - 1) * nx + ix + (it - 1) * nx * ny)))
}

g.sep_st_exp_exp2 <- function(param, useq, vseq, ghat, transform, power){
  uv <- expand.grid(useq, vseq)
  g.parametric <- exp(param[1] * exp( - uv[, 2] / param[3]) * exp( - uv[, 1] / param[2]))
  ghat2 <- as.vector(ghat)
  return(sum((ghat2 ^ power - g.parametric ^ power) ^ 2))
}

kde2d.new.var <- function(x, y, gx = x, gy = y, np, var.bin = TRUE, y.var = TRUE, hx, hy){
  
  nx  <-  length(x)
  if(length(y) != nx)   stop("Data vectors must be the same length")
  if(var.bin){
    h <- numeric()
    if(y.var){
      mat.d <- as.matrix(dist(cbind(x, y), upper = TRUE, diag = TRUE))
    }
    if(!y.var){
      mat.d <- as.matrix(dist(cbind(x), upper = TRUE, diag = TRUE))
    }
    diag(mat.d) <- Inf 
    if(missing(np))  np <- trunc(nx / 2)
    np <- trunc(np)
    if(np <= 1) np <- 1
    for(i in 1:ncol(mat.d)){
      o <- order(mat.d[, i])
      pp <- o[np]
      h[i] <- mat.d[pp, i]
    }
    h <- h / 4
    hx <- h
    hy <- h
    rm(mat.d)
  } else{ np <- 0
  if(missing(hx) || missing(hy)){
    h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
    h <- h / 4 / 1.06 * 0.9
    hx <- rep(h[1], nx)
    hy <- rep(h[2], nx)
  }
  }
  n <- length(gx)
  ax <- (outer(gx, x, "-") / hx)
  ay <- (outer(gy, y, "-") / hy)
  z <- rowSums(t(t(matrix(dnorm(ax), n, nx) * matrix(dnorm(ay), n, nx)) / (nx * hx * hy)))
  return(list(x = gx, y = gy, z = z, hx = hx, hy = hy, np = np))
  
}

norm2_etas <- function(n = 1, d = 1, q = 2, x0, y0){
  theta <- runif(n, max = 2 * pi)
  U <- runif(n)
  R <- sqrt(d * (U ^ (1 / (1 - q)) - 1))
  return(cbind(x0 + R * cos(theta), y0 + R * sin(theta)))
}

permutest.stlp <- function(perm, Q_n, Xi, L0, method){
  n_perm <- length(perm)
  n_points <- sum(perm == 2) + 1
  A <- sample(perm, size = n_perm, replace = FALSE, prob = NULL)
  cod <- which(A == 2)
  Q_perm <- Q_n[cod, 1:3]
  Q <- stlnpp::as.stlpp(c(Xi$x, Q_perm$x),
                        c(Xi$y, Q_perm$y),
                        c(Xi$t, Q_perm$t),
                        L = L0)
  
  Q_dens <- rep(n_points / (spatstat.geom::volume(L0) * (Q$time[2] - Q$time[1])), n_points)
  lista_0 <- if (method == "K"){
    STLKinhom_i(Q, lambda = Q_dens)
  } else {
    STLginhom_i(Q, lambda = Q_dens)
  }
  lista_0
}

permutest.stp <- function(perm, Q_n, Xi, method){
  
  n_perm <- length(perm)
  n_points <- sum(perm == 2) + 1
  A <- sample(perm, size = n_perm, replace = FALSE, prob = NULL)
  cod <- which(A == 2)
  Q_perm <- Q_n[cod, 1:3]
  
  Q <- cbind(c(Xi$x, Q_perm$x),
             c(Xi$y, Q_perm$y),
             c(Xi$t, Q_perm$t))
  class(Q) <- "stpp"
  
  lista_0 <- if (method == "K"){
    KLISTAhat(Q)$list.KLISTA
  } else {
    LISTAhat(Q)$list.LISTA
  }
  lista_0
}

STLginhom_i <- function(X, lambda, normalize = FALSE, r = NULL, t = NULL, nxy = 10)
{
  if (!inherits(X, "stlpp"))
    stop("X should be from class stlpp")
  
  Y <- stlnpp::as.lpp.stlpp(X)
  l <- spatstat.geom::domain(Y)
  tleng <- summary(l)$totlength
  n <- spatstat.geom::npoints(Y)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- spatstat.geom::pairdist(Y)
  tdist <- as.matrix(dist(timev))
  toler <- spatstat.linnet::default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for (j in 1:n) {
    ml[-j, j] <- spatstat.linnet::countends(l, Y[-j], sdist[-j, j], toler = toler)
  }
  mtplus <- matrix(timev, n, n, byrow = T) + tdist
  mtminus <- matrix(timev, n, n, byrow = T) - tdist
  mtedge <- (mtplus <= b) + (mtminus >= a)
  diag(mtedge) <- 1
  lamden <- outer(lambda, lambda, FUN = "*")
  edgetl <- mtedge * ml * lamden
  maxs <- 0.7 * max(sdist[!is.infinite(sdist)])
  maxt <- 0.7 * (trange/2)
  if (is.null(r))
    r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
  if (is.null(t))
    t <- seq((maxt/nxy),maxt,by=(maxt-(maxt/nxy))/(nxy-1))
  
  g_local <- list()
  for(k in 1:n){
    
    g <- matrix(NA, nrow = nxy, ncol = nxy)
    no <- sdist == 0&tdist == 0|sdist == Inf|sdist > maxs|tdist > maxt
    bwl <- stats::bw.nrd0(as.numeric(sdist[!no]))
    bwt <- stats::bw.nrd0(as.numeric(tdist[!no]))
    
    for (i in 1:length(r)) {
      for (j in 1:length(t)) {
        
        outl <- spatstat.explore::dkernel(as.numeric(sdist[!no] - r[i]), sd = bwl)
        outt <- spatstat.explore::dkernel(as.numeric(tdist[!no] - t[j]), sd = bwt)
        g1 <- outl * outt / (edgetl[!no])
        no2 <- no
        no2[no] <- NA
        no2[!no] <- g1
        no2 <- no2[, k]
        g[i, j] <- sum(no2[!is.na(no2) & !is.infinite(no2)])
        
      }
    }
    g_local[[k]] <- g
  }
  
  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n-1) * (tleng * trange) / (sum(revrho[lower.tri(revrho,diag = FALSE)]) * 2)
    gval <- lapply(g_local, "*" , appx)
  }
  else {
    g_local <- lapply(g_local, FUN= function(K) (n-1) * K / (trange * tleng) )
  }
  
  arr <- array( unlist(g_local) , c(nxy, nxy, n) )
  
  return(arr)
}

STLKinhom_i <- function(X, lambda = lambda, normalize = FALSE, r = NULL, t = NULL,
                        nxy = 10)
{
  if (!inherits(X, "stlpp"))
    stop("X should be from class stlpp")
  
  Y <- as.lpp.stlpp(X)
  l <- domain(Y)
  tleng <- summary(l)$totlength
  n <- npoints(X)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- pairdist.lpp(Y)
  tdist <- as.matrix(dist(timev))
  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for (j in 1:n) {
    ml[-j, j] <- countends(l, Y[-j], sdist[-j, j], toler = toler)
  }
  mtplus <- matrix(timev, n, n, byrow = T) + tdist
  mtminus <- matrix(timev, n, n, byrow = T) - tdist
  mtedge <- (mtplus <= b) + (mtminus >= a)
  diag(mtedge) <- 1
  lamden <- outer(lambda, lambda, FUN = "*")
  diag(lamden) <- 1
  edgetl <- mtedge* ml * lamden
  
  maxs <- 0.7 * max(sdist[!is.infinite(sdist)])
  maxt <- 0.7 * (trange/2)
  if (is.null(r))
    r <- seq((maxs/nxy), maxs, by = (maxs - (maxs/nxy))/(nxy - 1))
  if (is.null(t))
    t <- seq((maxt/nxy), maxt, by = (maxt - (maxt/nxy))/(nxy - 1))
  
  K_local <- list()
  for(k in 1:n){
    K <- matrix(NA, nrow = nxy, ncol = nxy)
    for (i in 1:length(r)) {
      for (j in 1:length(t)) {
        out <- (sdist <= r[i]) * (tdist <= t[j])
        diag(out) <- 0
        kout <- out / edgetl
        kout <- kout[, k]
        K[i, j] <- sum(kout[is.finite(kout)])
      }
    }
    K_local[[k]] <- K
  }
  
  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n - 1) * (tleng * trange) / (sum(revrho[lower.tri(revrho, diag = FALSE)]) * 2)
    K_local <- lapply(K_local, "*" , appx)#
    
  }
  else {
    K_local <- lapply(K_local, FUN= function(K) (n - 1) * K / (trange * tleng) )
  }
  
  arr <- array(unlist(K_local), c(nxy, nxy, n))
  
  return(arr)
}

# aggiunte nella seconda versione
scale_to_range <- function(x, new_min = 0, new_max = 1) {
  ((x - min(x)) / (max(x) - min(x))) * (new_max - new_min) + new_min
}

distD <- function(A, B, d){
  a <- switch(as.character(d), 
              "1" = sqrt((B$x - A$xx) ^ 2),
              "2" = sqrt((B$x - A$xx) ^ 2+ (B$y - A$xy) ^ 2),
              "3" = sqrt((B$x - A$xx) ^ 2+ (B$y - A$xy) ^ 2 + (B$t - A$xt) ^ 2))
  a
}

idw3D <- function(points, covs, p, iid, d, parallel, cl){
  nV <- dim(covs)[1]
  wi <- vector(length = nV)
  dd <- if(parallel) {
    parSapply(cl, 1:nV, function(i) distD(points[iid, ], covs[i, 1:d, 
                                                              drop = FALSE], 
                                          d))
  } else {
    sapply(1:nV, function(i) distD(points[iid, ], covs[i, 1:d, 
                                                       drop = FALSE], d))
  }
  wi <- 1 / (dd) ^ p
  wi
}

interp3D <- function(points, covs, p, d, verbose = FALSE, parallel, cl){
  nU <- dim(points)[1]
  gu <- vector(length = nU)
  if(verbose) cat("Interpolating covariate values at point locations: \n")
  wi <- sapply(1:nU, function(j) idw3D(points, covs, p, iid = j, d, 
                                       parallel = parallel, cl = cl))
  gu <- apply(wi, 2 , function(p) sum(p * covs[, 4]) / sum(p))
  gu 
}

interpMin <- function(dati, covariate, parallel = parallel, cl = cl){
  # dati e covariate sono due matrici con 3 colonne l'uno (x,y,z) e righe diverse
  # colnames(dati)[1:3] <- c("xx", "xy", "xt")
  # colnames(covariate)[1:3] <- c("x", "y", "z")
  nV <- nrow(covariate)
  dd <- if(parallel) {
    parSapply(cl, 1:nV, function(i) distD(dati, covariate[i, ], d = 3))
  } else {
    sapply(1:nV, function(i) distD(dati, covariate[i, ], d = 3))
  }
  id <- apply(dd, 1, which.min)
  gu <- covariate[, 4][id]
  gu
}

cartesian_3d_N_Mark <- function(pp3_obj, box_3d, list_levels_marks) {
  
  grid_marks_coord_x <- list_levels_marks
  grid_marks_coord_x$x <- pp3_obj$data$x
  coord_x <- do.call(expand.grid, grid_marks_coord_x)
  
  grid_marks_coord_y <- list_levels_marks
  grid_marks_coord_y$y <- pp3_obj$data$y
  coord_y <- do.call(expand.grid, grid_marks_coord_y)
  
  grid_marks_coord_t <- list_levels_marks
  grid_marks_coord_t$t <- pp3_obj$data$z
  coord_t <- do.call(expand.grid, grid_marks_coord_t)
  
  result <- pp3(x = coord_x$x, y = coord_y$y, z = coord_t$t, box_3d, 
                marks = coord_x[, c(1:length(list_levels_marks)), drop = FALSE])
  return(result)
}

dummy.marked.result <- function(X, formula, dummy_points, Wdum, Wdat, ndata, ndummy) {
  
  marks_name <- names(X[,-c(1,2,3), drop=F])
  marks_in_formula <- marks_name %in% all.vars(formula)
  marks.name.formula <- marks_name[marks_in_formula]
  
  df_marks <- as.data.frame(matrix(NA, nrow = nrow(X), ncol = length(marks.name.formula)))
  colnames(df_marks) <- marks.name.formula
  j <- 1
  for (i in (3 + which(marks_in_formula))) {
    df_marks[,j] <- X[,i]
    j <- j+1
  }
  
  for(i in 1:ncol(df_marks)) {
    df_marks[,i] <- as.factor(df_marks[,i])
  }
  
  list_levels_marks0 <- lapply(df_marks, function(col) levels(as.factor(col)))
  
  n_levels_marks <- NULL
  for(i in 1:ncol(df_marks)) {
    n_levels_marks[i] <- length(levels(as.factor(df_marks[,i])))
  }
  
  n_comb_levels <- prod(n_levels_marks)
  
  b_dummy_3d <- box3(range(dummy_points$x),range(dummy_points$y),range(dummy_points$t)) 
  dummy_3d <- pp3(dummy_points$x, dummy_points$y, dummy_points$t, b_dummy_3d)
  
  dumdum <- cartesian_3d_N_Mark(dummy_3d, b_dummy_3d, list_levels_marks0)
  Wdumdum <- rep.int(Wdum, n_comb_levels)
  Idumdum <- rep.int(ndata + seq_len(ndummy), n_comb_levels)
  
  b_data_3d <- box3(range(X$x),range(X$y),range(X$t)) 
  data_3d <- pp3(X$x, X$y, X$t, b_data_3d)
  
  dumdat <- cartesian_3d_N_Mark(data_3d, b_data_3d, list_levels_marks0)
  Wdumdat <- rep.int(Wdat, n_comb_levels)
  Mdumdat <- marks(dumdat)
  Mdumdat <- as.data.frame(Mdumdat)
  Idumdat <- rep.int(1:ndata, n_comb_levels)
  
  Mrepdat <- do.call(rbind, replicate(n_comb_levels, df_marks, simplify = FALSE))
  
  casi_uguali <- apply(Mdumdat == Mrepdat, 1, all)
  
  dumdat <- dumdat[-c(which(casi_uguali == T)),]
  Wdumdat <- Wdumdat[-c(which(casi_uguali == T))]
  Idumdat <- Idumdat[-c(which(casi_uguali == T))]
  
  total_dummy_marks <- as.data.frame(matrix(NA, nrow = nrow(dumdum$data)+nrow(dumdat$data), 
                                            ncol = length(4:ncol(dumdum$data))))
  colnames(total_dummy_marks) <- colnames(dumdum$data[,4:ncol(dumdum$data)])
  
  df_dumdum <- as.data.frame(dumdum$data)
  df_dumdat <- as.data.frame(dumdat$data)
  k <- 1 
  for(i in 4:ncol(dumdum$data)) {
    total_dummy_marks[,k] <- c(df_dumdum[,i], df_dumdat[,i])
    k <- k+1
  }
  
  for(i in 1:ncol(total_dummy_marks)) {
    total_dummy_marks[,i] <- as.factor(total_dummy_marks[,i])
  }
  
  dumb <- pp3(c(dumdum$data$x, dumdat$data$x), 
              c(dumdum$data$y, dumdat$data$y),
              c(dumdum$data$z, dumdat$data$z),
              b_dummy_3d, marks = total_dummy_marks)
  
  Wdumb <- c(Wdumdum, Wdumdat)
  Idumb <- c(Idumdum, Idumdat)
  
  result.dummy.marked <- list(dumb = dumb, Wdumb = Wdumb, df_marks = df_marks, 
                              total_dummy_marks = total_dummy_marks, 
                              n_comb_levels = n_comb_levels)
  return(result.dummy.marked)
}

interp.covariate <- function(X, dummy_points, covs, formula, parallel, interp,
                             xx, xy, xt, ncores, verbose) {
  dati.interpolati <- rbind(X[,1:3], dummy_points)
  colnames(dati.interpolati) <- c("x", "y", "t")
  
  cc <- vector(length = length(covs))
  for(ki in 1:length(covs)){
    cc[ki] <- c(names(covs[[ki]]$df)[4])
  } 
  names(covs) <- cc
  ff <- which(names(covs) %in% all.vars(formula)) 
  covs <- covs[ff]                               
  
  for(k in names(covs)){  
    kk <- which(names(covs)  == k)
    if(verbose) {cat("Covariate", k, "\n")}
    if(interp) {
      covs0 <- covs[[k]]$df
      colnames(covs0) <- c("x", "y", "t", names(covs[[k]]$df)[4])
      if(parallel) {
        on.exit({stopCluster(cl)}, add = TRUE, after = TRUE)
        cl <- makeCluster(getOption("cl.cores", ncores))
        clusterExport(cl = cl, c('distD'))
        gu <- interpMin(data.frame(xx, xy, xt), covs0, parallel = TRUE, cl = cl)
      } else {
        gu <- interpMin(data.frame(xx, xy, xt), covs0, parallel = FALSE)}
    } else {
      df0 <- data.frame(covs0$x, covs0$y, covs0$t)
      colnames(df0) <- c("xx", "xy", "xt")
      covs0 <- covs[[k]]$df
      colnames(covs0) <- c("x", "y", "t", names(covs[[k]])[4])
      if(parallel) {
        on.exit({stopCluster(cl)}, add = TRUE, after = TRUE)
        cl <- makeCluster(getOption("cl.cores", ncores))
        clusterExport(cl = cl, c('distD'))
        gu <- interpMin(df0, covs0, parallel = TRUE, cl = cl)
      } else {
        gu <- interpMin(df0, covs0, parallel = FALSE)
      }
      gu <- c(gu, covs[[k]]$df[, 4])
    }
    dati.interpolati <- cbind(dati.interpolati, gu)
    colnames(dati.interpolati)[3 + kk] <- k
  }
  return(dati.interpolati)
}

as.stpp <- function(x){
  if(!inherits(x,"stp")) stop("class(x) must be stp")
  
  out <- cbind(x$df$x, x$df$y, x$df$t)
  colnames(out) <- c("x", "y", "t")
  class(out) <- "stpp"
  return(out)
}

as.stp <- function(x){
  if(!inherits(x,"stpp")) stop("class(x) must be stpp")
  
  stp(cbind(x[, 1], x[, 2], x[, 3]))
}

as.stlpp <- function(x){
  if(!inherits(x,"stlp")) stop("class(x) must be stlp")
  
  stlnpp::as.stlpp(x$df$x, x$df$y, x$df$t, x$L)
}

as.stlp <- function(x){
  if(!inherits(x,"stlpp")) stop("class(x) must be stlpp")
  
  stp(cbind(x$data$x, x$data$y, x$data$t), x$domain)
}

is.stp <- function(x){
  inherits(x,"stp")
}

is.stlp <- function(x){
  inherits(x,"stlp")
}

localplot.stlgcppm <- function(x, par = TRUE){
  
  if(inherits(x$IntCoefs, "numeric") & inherits(x$CovCoefs, "numeric")){
    stop("No local parameters to plot")
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if(!inherits(x$IntCoefs, "numeric")){
    nn <- length(names(x$IntCoefs))
    
    if(par == T){
      par(mfrow = c(1, nn))
    } else{
      par(mfrow = c(1, 1))
      par(ask = TRUE)
    }
    
    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
    plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                      theta = - 45, phi = 20,
                      col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                          range = range(x$IntCoefs$`(Intercept)`)),
                                 "stuff")$outputs,
                      ticktype = "detailed", pch = 19,
                      colvar = x$IntCoefs$`(Intercept)`,
                      xlab="x",ylab="y",zlab="t",
                      main = c("Intercept"))
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != TRUE){
      par(ask = FALSE)
    }
    
    if(nn > 1){
      for(i in 2:nn){
        
        id <- x$IntCoefs[, i]
        
        if(par != TRUE){
          par(ask = TRUE)
        }
        par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, -1))
        plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                          theta = - 45, phi = 20,
                          col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                              range = range(id)),
                                     "stuff")$outputs,
                          ticktype = "detailed", pch = 19,
                          colvar = id,
                          xlab="x",ylab="y",zlab="t",
                          main = names(x$IntCoefs)[i])
        par(mar = c(5, 4, 4, 2) + 0.1)
        if(par != TRUE){
          par(ask = FALSE)
        }
        
      }
    }
    
    
  }
  
  
  if(!inherits(x$CovCoefs, "numeric")){
    
    if(ncol(x$CovCoefs) == 3){
      if(par == T){
        par(mfrow = c(1, 3))
      } else{
        par(mfrow = c(1, 1))
      }
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$sigma)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$sigma,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(sigma))
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$alpha)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$alpha,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(alpha))
      if(par != TRUE){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$beta)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$beta,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(beta))
      par(mar = c(5, 4, 4, 2) + 0.1)
      if(par != TRUE){
        par(ask = FALSE)
      }
      
    } else {
      if(par == TRUE){
        par(mfrow = c(2, 3))
      } else{
        par(mfrow = c(1, 1))
      }
      
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$sigma)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$sigma,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(sigma))
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$alpha)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$alpha,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(alpha))
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$beta)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$beta,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(beta))
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$gamma_s)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$gamma_s,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(gamma))
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$gamma_t)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$gamma_t,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(gamma[t]))
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$delta)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$delta,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(delta))
      par(mar = c(5, 4, 4, 2) + 0.1)
      if(par != TRUE){
        par(ask = FALSE)
      }
    }
  }
  
}

localplot.locstppm <- function(x, par = TRUE){
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  nn <- length(names(x$IntCoefs))
  
  if(par == TRUE){
    par(mfrow = c(1, nn))
  } else{
    par(mfrow = c(1, 1))
    par(ask = TRUE)
  }
  
  par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
  plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                    theta = - 45, phi = 20,
                    col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                        range = range(x$IntCoefs_local[, 1])),
                               "stuff")$outputs,
                    ticktype = "detailed", pch = 19,
                    colvar = x$IntCoefs_local[, 1],
                    xlab="x",ylab="y",zlab="t",
                    main = c("Intercept"))
  par(mar = c(5, 4, 4, 2) + 0.1)
  if(par != TRUE){
    par(ask = FALSE)
  }
  
  if(nn > 1){
    for(i in 2:nn){
      
      id <- x$IntCoefs_local[, i]
      
      if(par != TRUE){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, -1))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(id)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = id,
                        xlab="x",ylab="y",zlab="t",
                        main = names(x$IntCoefs_local)[i])
      par(mar = c(5, 4, 4, 2) + 0.1)
      if(par != TRUE){
        par(ask = FALSE)
      }
      
    }
  }
  
}

localsummary.stlgcppm <- function(x,
                                  scaler = c("silverman", "IQR", "sd", "var"),
                                  do.points = TRUE,
                                  print.bw = FALSE,
                                  zap = 0.00001,
                                  par = TRUE){
  
  if(inherits(x$IntCoefs, "numeric")){
    if(x$formula == "~1"){
      stop("No inhomogeneous intensity to summarise")
    } else {
      stop("No inhomogeneous intensity made by local parameters")
    }
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  mod_Min <- mod_1stQu <- mod_Median <- mod_3rdQu <- mod_Max <- x$mod_global
  
  mod_Min$coefficients <- apply(x$IntCoefs, 2, summary)[1, ]
  mod_1stQu$coefficients <- apply(x$IntCoefs, 2, summary)[2, ]
  mod_Median$coefficients <- apply(x$IntCoefs, 2, summary)[3, ]
  mod_3rdQu$coefficients <- apply(x$IntCoefs, 2, summary)[5, ]
  mod_Max$coefficients <- apply(x$IntCoefs, 2, summary)[6, ]
  
  mark_int <- x$l
  l_Min <- predict(mod_Min, newdata = x$newdata)
  l_1stQu <- predict(mod_1stQu, newdata = x$newdata)
  l_Median <- predict(mod_Median, newdata = x$newdata)
  l_3rdQu <- predict(mod_3rdQu, newdata = x$newdata)
  l_Max <- predict(mod_Max, newdata = x$newdata)
  
  ppx_int <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = mark_int,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_Min <- suppressWarnings(spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Min,
                                                 window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y))))
  ppx_1stQu <- suppressWarnings(spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_1stQu,
                                                   window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y))))
  ppx_Median <- suppressWarnings(spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Median,
                                                    window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y))))
  ppx_3rdQu <- suppressWarnings(spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_3rdQu,
                                                   window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y))))
  ppx_Max <- suppressWarnings(spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Max,
                                                 window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y))))
  
  sig <- sparr::OS(unmark(ppx_int), scaler = scaler)
  
  s_int <- spatstat.explore::Smooth(ppx_int, sigma = sig)
  s_Min <- spatstat.explore::Smooth(ppx_Min, sigma = sig)
  s_1stQu <- spatstat.explore::Smooth(ppx_1stQu, sigma = sig)
  s_Median <- spatstat.explore::Smooth(ppx_Median, sigma = sig)
  s_3rdQu <- spatstat.explore::Smooth(ppx_3rdQu, sigma = sig)
  s_Max <- spatstat.explore::Smooth(ppx_Max, sigma = sig)
  
  g00 <- spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                  range = c(min(s_int, s_Min$v[1], s_1stQu$v[1], s_Median$v[1],
                                                s_3rdQu$v[1], s_Max$v[1]),
                                            max(s_int, s_Min$v[1], s_1stQu$v[1], s_Median$v[1],
                                                s_3rdQu$v[1], s_Max$v[1])))
  
  if(par == T){
    par(mfrow = c(2, 3))
  } else {
    par(mfrow = c(1, 1))
    par(ask = FALSE)
  }
  
  if(x$formula == "~1"){
    if(par != T){
      par(ask = TRUE)
    }
    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
    plot(s_int,
         col = g00,
         main = c("First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Min,
         col = g00,
         main = paste("Min. First-order Intensity in space \n lambda = ", round(s_Min$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Min), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_1stQu,
         col = g00,
         main = paste("1st Qu. First-order Intensity in space \n lambda = ",  round(s_1stQu$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_1stQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Median,
         col = g00,
         main = paste("Median First-order Intensity in space \n lambda = ",  round(s_Median$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Median), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_3rdQu,
         col = g00,
         main = paste("3rd Qu. First-order Intensity in space \n lambda = ",  round(s_3rdQu$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_3rdQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Max,
         col = g00,
         main = paste("Max. First-order Intensity in space \n lambda = ",  round(s_Max$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Max), add = T)}
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != T){
      par(ask = FALSE)
    }
  } else {
    if(par != T){
      par(ask = TRUE)
    }
    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
    plot(s_int, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_int)),
         main = c("First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Min, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Min)),
         main = c("Min. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Min), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_1stQu, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_1stQu)),
         main = c("1st Qu. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_1stQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Median, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Median)),
         main = c("Median First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Median), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_3rdQu, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_3rdQu)),
         main = c("3rd Qu. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_3rdQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Max, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Max)),
         main = c("Max. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Max), add = T)}
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != T){
      par(ask = FALSE)
    }
  }
  
  
  if(print.bw == T){print(sig)}
}

localsummary.locstppm <- function(x,
                                  scaler = c("silverman", "IQR", "sd", "var"),
                                  do.points = TRUE,
                                  print.bw = FALSE,
                                  zap = 0.00001,
                                  par = TRUE){
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  mod_Min <- mod_1stQu <- mod_Median <- mod_3rdQu <- mod_Max <- x$mod_global
  
  mod_Min$coefficients <- apply(x$IntCoefs_local, 2, summary)[1, ]
  mod_1stQu$coefficients <- apply(x$IntCoefs_local, 2, summary)[2, ]
  mod_Median$coefficients <- apply(x$IntCoefs_local, 2, summary)[3, ]
  mod_3rdQu$coefficients <- apply(x$IntCoefs_local, 2, summary)[5, ]
  mod_Max$coefficients <- apply(x$IntCoefs_local, 2, summary)[6, ]
  
  mark_int <- x$l_local
  l_Min <- predict(mod_Min, newdata = x$newdata)
  l_1stQu <- predict(mod_1stQu, newdata = x$newdata)
  l_Median <- predict(mod_Median, newdata = x$newdata)
  l_3rdQu <- predict(mod_3rdQu, newdata = x$newdata)
  l_Max <- predict(mod_Max, newdata = x$newdata)
  
  ppx_int <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = mark_int,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_Min <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Min,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_1stQu <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_1stQu,
                                  window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_Median <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Median,
                                   window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_3rdQu <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_3rdQu,
                                  window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_Max <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Max,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  
  sig <- sparr::OS(unmark(ppx_int), scaler = scaler)
  
  s_int <- spatstat.explore::Smooth(ppx_int, sigma = sig)
  s_Min <- spatstat.explore::Smooth(ppx_Min, sigma = sig)
  s_1stQu <- spatstat.explore::Smooth(ppx_1stQu, sigma = sig)
  s_Median <- spatstat.explore::Smooth(ppx_Median, sigma = sig)
  s_3rdQu <- spatstat.explore::Smooth(ppx_3rdQu, sigma = sig)
  s_Max <- spatstat.explore::Smooth(ppx_Max, sigma = sig)
  
  g00 <- spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                  range = c(min(s_int, s_Min$v[1], s_1stQu$v[1], s_Median$v[1],
                                                s_3rdQu$v[1], s_Max$v[1]),
                                            max(s_int, s_Min$v[1], s_1stQu$v[1], s_Median$v[1],
                                                s_3rdQu$v[1], s_Max$v[1])))
  
  if(par == T){
    par(mfrow = c(2, 3))
  } else {
    par(mfrow = c(1, 1))
    par(ask = FALSE)
  }
  
  if(x$formula == "~1"){
    if(par != T){
      par(ask = TRUE)
    }
    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
    plot(s_int,
         col = g00,
         main = c("First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Min,
         col = g00,
         main = paste("Min. First-order Intensity in space \n lambda = ", round(s_Min$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Min), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_1stQu,
         col = g00,
         main = paste("1st Qu. First-order Intensity in space \n lambda = ",  round(s_1stQu$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_1stQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Median,
         col = g00,
         main = paste("Median First-order Intensity in space \n lambda = ",  round(s_Median$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Median), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_3rdQu,
         col = g00,
         main = paste("3rd Qu. First-order Intensity in space \n lambda = ",  round(s_3rdQu$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_3rdQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Max,
         col = g00,
         main = paste("Max. First-order Intensity in space \n lambda = ",  round(s_Max$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Max), add = T)}
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != T){
      par(ask = FALSE)
    }
  } else {
    if(par != T){
      par(ask = TRUE)
    }
    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
    plot(s_int, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_int)),
         main = c("First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Min, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Min)),
         main = c("Min. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Min), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_1stQu, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_1stQu)),
         main = c("1st Qu. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_1stQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Median, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Median)),
         main = c("Median First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Median), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_3rdQu, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_3rdQu)),
         main = c("3rd Qu. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_3rdQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Max, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Max)),
         main = c("Max. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Max), add = T)}
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != T){
      par(ask = FALSE)
    }
  }
  
  
  if(print.bw == T){print(sig)}
}