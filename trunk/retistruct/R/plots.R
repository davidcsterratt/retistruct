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
rho.to.degrees <- function(rho, phi0, area.preserving) {
  if (area.preserving) {
    phi0d <- phi0*180/pi
    rho0 <- spherical.to.polar.area(phi0)
    return((phi0d + 90)/rho0*rho)
  } else {
    return(rho*180/pi)
  }
}

