##' Create an object that is specific to retinal datasets. This
##' contains methods that return datapoint and landmark coordinates
##' that have been transformed according to the values of
##' \code{DVflip} and \code{side}.
##'
##' @title RetinalReconstructedOutline constructor
##' @param r Object that inherits \code{ReconstructedOutline}
##' @return \code{RetinalReconstructedOutline} object. This does not
##' contain any extra fields, but there are extra mthods dthat apply
##' to it.
##' @author David Sterratt
RetinalReconstructedOutline <- function(r) {
  if (!(inherits(r, "reconstructedOutline"))) {
    stop("Argument needs to inherit reconstructedOutline")
  }
  class(r) <- addClass("retinalReconstructedOutline", r)
  return(r)
}

getIms.retinalReconstructedOutline <- function(r) {
  ims <- NextMethod()
  if (r$DVflip) {
    if (!is.null(ims)) {
      ims[,"lambda"] <- -ims[,"lambda"]
    }
  }
  if (r$side=="Left") {
    if (!is.null(ims)) {
      ims[,"lambda"] <- 2*pi - ims[,"lambda"]
    }
  }
  return(ims)
}
