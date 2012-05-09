##' @title Invert sphere about its centre
##' @param r Coordinates of points in spherical coordinates
##' represented as  2 column matrix with column names \code{phi}
##' (lattitude) and \code{lambda} (longitude).
##' @return Matrix in same format, but with \code{pi} added to lambda
##' and \code{phi} negated.
##' @author David Sterratt
##' @export
invert.sphere <- function(r) {
  r[,"phi"] <- -r[,"phi"]
  r[,"lambda"] <- (r[,"lambda"] + pi) %% (2*pi)
  return(r)
}
