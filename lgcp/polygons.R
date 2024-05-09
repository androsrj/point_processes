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

setAs("tess", "SpatialPolygons", function(from) as.SpatialPolygons.tess(from))


as.SpatialPolygons.owin <- function(x) {
  #check_spatstat("spatstat.geom")
  stopifnot(spatstat.geom::is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

setAs("owin", "SpatialPolygons", function(from) as.SpatialPolygons.owin(from))


