##' The identity transformation
##' @param r Coordinates of points in spherical coordinates
##' represented as  2 column matrix with column names \code{phi}
##' (latitude) and \code{lambda} (longitude).
##' @param ... Other arguments
##' @return Identical matrix
##' @author David Sterratt
##' @export
identity.transform <- function(r, ...) {
  return(r)
}

##' Invert sphere about its centre
##' @param r Coordinates of points in spherical coordinates
##' represented as  2 column matrix with column names \code{phi}
##' (latitude) and \code{lambda} (longitude).
##' @param ... Other arguments
##' @return Matrix in same format, but with \code{pi} added to lambda
##' and \code{phi} negated.
##' @author David Sterratt
##' @export
invert.sphere <- function(r, ...) {
  r[,"phi"] <- -r[,"phi"]
  r[,"lambda"] <- (r[,"lambda"] + pi) %% (2*pi)
  return(r)
}

##' Invert image of a partial sphere and scale the longitude so that
##' points at latitude \code{phi0} is projected onto a longitude of 0
##' degrees (the equator).
##' @title Invert sphere to hemisphere
##' @param r Coordinates of points in spherical coordinates
##' represented as  2 column matrix with column names \code{phi}
##' (latitude) and \code{lambda} (longitude).
##' @param phi0 The latitude to map onto the equator
##' @param ... Other arguments
##' @return Matrix in same format, but with \code{pi} added to lambda
##' and \code{phi} negated and scaled so that the longitude
##' \code{phi0} is projected to 0 degrees (the equator)
##' @author David Sterratt
##' @export
invert.sphere.to.hemisphere <- function(r, phi0, ...) {
  r[,"phi"] <- -((r[,"phi"]+pi/2)*(pi/2)/(phi0+pi/2)-pi/2)
  r[,"lambda"] <- (r[,"lambda"] + pi) %% (2*pi)
  return(r)
}
