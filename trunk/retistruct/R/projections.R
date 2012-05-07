##' @title Sinusoidal projection
##' @param r Lattitude-longitude coordinates in a matrix with columns
##' labelled \code{phi} (lattitude) and \code{lambda} (longitude)
##' @param lambda0 Coordinate of central meridian
##' @param lambdalim 
##' @param lines 
##' @return Two-column matrix with columns labelled \code{x} and
##' \code{y} of locations of projection of coordinates on plane 
##' @author David Sterratt
##' @export
sinusoidalproj <- function(r, lambda0=0,
                           lambdalim=NULL, lines=FALSE) {
  x <- (r[,"lambda"] - lambda0)*cos(r[,"phi"])
  y <- r[,"phi"]
  rc <- cbind(x=x, y=y)
  if (!is.null(lambdalim)) {
    rc[r[,"lambda"] > lambdalim[2],] <- NA
    rc[r[,"lambda"] < lambdalim[1],] <- NA
    if (lines) {
      inds <- which(abs(diff(r[,"lambda"])) > 0.5*(lambdalim[2] - lambdalim[1]))
      rc[inds,] <- NA
    }
  }
  return(rc)
}

##' This is the inverse of \code{\link{polar.cart.to.sphere.spherical}}
##'
##' @title Convert spherical coordinates on sphere to  polar
##' projection in Cartesian coordinates
##' @param r 2-column Matrix of spherical coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @return 2-column Matrix of Cartesian coordinates of points on polar
##' projection. Column names should be \code{x} and \code{y}
##' @author David Sterratt
##' @export
lambertproj <- function(r, ...) {
  rho <- sqrt(2*(1 + sin(r[,"phi"])))
  x <- rho*cos(r[,"lambda"])
  y <- rho*sin(r[,"lambda"])
  return(cbind(x=x, y=y))
}
