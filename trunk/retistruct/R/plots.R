##
## Utilities for plotting functions
## 

##' @title Generate colours for strain plots
##' @param x Vector of values of log strain
##' @return Vector of colours corresponding to strains
##' @author David Sterratt
strain.colours <- function(x) {
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  col <- x/log(0.75)*35 + 35
  col[col<1] <- 1
  col[col>70] <- 70
  return(col)
}

##' Place text at bottom right of \code{\link{plot.polar}}
##'
##' @title Put text on the polar plot
##' @param text Test to place
##' @author David Sterratt
##' @export
text.polar <- function(text) {
  mtext(text, 1, adj=1, line=-1)
}

##' This is a helper function for \code{\link{plot.polar}}
##'
##' @title Convert lattitude to radial variable in polar plot
##' @param phi Lattitude
##' @param phi0 Lattitude of top of curtailed sphere
##' @param area.preserving Whether the conversion should preserve area
##' @return Radial variable
##' @author David Sterratt
##' @export
phi.to.rho <- function(phi, phi0, area.preserving) {
  if (area.preserving) {
    phi0d <- phi0*180/pi
    rho0 <- spherical.to.polar.area(phi0)
    return((phi0d + 90)/rho0*spherical.to.polar.area(phi))
  } else {
    return(phi*180/pi + 90)      
  }
}

## FIXME: it would be nice to tidy up the profusion of conversion
## functions

##' This is a helper function for \code{\link{plot.polar}}
##'
##' @title Convert rho variable into normalised radial coordinate in degrees
##' @param rho The rho variable. This may be in radians or the
##' area-preserving coordinate
##' @param phi0 Lattitude of top of curtailed sphere
##' @param area.preserving Whether the conversion should preserve area
##' @return Radial variable
##' @author David Sterratt
##' @export
rho.to.degrees <- function(rho, phi0, area.preserving) {
  if (area.preserving) {
    phi0d <- phi0*180/pi
    rho0 <- spherical.to.polar.area(phi0)
    return((phi0d + 90)/rho0*rho)
  } else {
    return(rho*180/pi)
  }
}


##' @title Plot a Voronoi mosaic in a region with a circular boundary
##' @param x \code{voronoi} object
##' @param R Radius at which to trim circle
##' @param ... Additional low-level graphics parameters
##' @author David Sterratt
plot.voronoi.circular <- function (x, R=110, ...) {
  chop.boundary <- function(x0, x1, y0, y1) {
    R12 <- x1^2 + y1^2
    if (R12 > R^2) {
      R02 <- x0^2 + y0^2
      ##                lambda <- sqrt((R^2 - R02)/(R02 + R12))
      l2 <- (x1 - x0)^2 + (y1 - y0)^2
      R0R1 <- x1*x0 + y0*y1
      lambda <- 1/l2*(-(R0R1 - R02) + sqrt(R0R1^2 - R12*R02 + R^2*l2))
      
      x1 <- (1 - lambda)*x0 + lambda*x1
      y1 <- (1 - lambda)*y0 + lambda*y1
    }
    return(c(x1, y1))
  }
  n <- length(x$x)
  x0 <-  matrix(NA, n, 3)
  x1 <-  matrix(NA, n, 3)
  y0 <-  matrix(NA, n, 3)
  y1 <-  matrix(NA, n, 3)
  lty <- matrix(NA, n, 3)
  for (i in 1:n) {
    if (x$node[i]) {
      if (x$x[i]^2 + x$y[i]^2 <  R^2) {
        tns <- sort(c(x$n1[i], x$n2[i], x$n3[i]))
        for (j in 1:3) {
          if (tns[j] > 0) {
            if (x$node[tns[j]]) {
              x0[i, j] <- x$x[i]
              x1[i, j] <- x$x[tns[j]]
              y0[i, j] <- x$y[i]
              y1[i, j] <- x$y[tns[j]]
              ## If the node is outside the boudary, chop at the boundary
              r1 <- chop.boundary(x0[i, j], x1[i, j], y0[i, j], y1[i, j])
              x1[i, j] <- r1[1]
              y1[i, j] <- r1[2]
              lty[i, j] <- 1
            }
          } else {
            if (tns[j] < 0) {
              x0[i, j] <- x$x[i]
              x1[i, j] <- x$dummy.x[-tns[j]]
              y0[i, j] <- x$y[i]
              y1[i, j] <- x$dummy.y[-tns[j]]
              r1 <- chop.boundary(x0[i, j], x1[i, j], y0[i, j], y1[i, j])
              x1[i, j] <- r1[1]
              y1[i, j] <- r1[2]
              lty[i, j] <- 1
            }
          }
        }
      }
    }
  }
  x0 <- na.omit(as.vector(x0))
  x1 <- na.omit(as.vector(x1))
  y0 <- na.omit(as.vector(y0))
  y1 <- na.omit(as.vector(y1))
  lty <- na.omit(as.vector(lty))
  segments(x0, y0, x1, y1, lty=lty, ...)
}

