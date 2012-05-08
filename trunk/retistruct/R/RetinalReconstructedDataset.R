##' Create an object that is specific to retinal datasets. This
##' contains methods that return datapoint and landmark coordinates
##' that have been transformed according to the values of
##' \code{DVflip} and \code{side}.
##'
##' @title RetinalReconstructedDataset constructor
##' @param r Object that inherits both \code{reconstructedDataset} and
##' \code{dataset}.
##' @return \code{\link{RetinalReconstructedDataset}} object. This does not
##' contain any extra fields, but there are extra mthods dthat apply
##' to it.
##' @author David Sterratt
##' @export
RetinalReconstructedDataset <- function(r) {
  if (!(inherits(r, "reconstructedDataset") &
        inherits(r, "retinalDataset"))) {
    stop("Argument needs to inherit reconstructedDataset and retinalDataset")
  }
  class(r) <- addClass("retinalReconstructedDataset", r)
  r$KDE <- getKDE(r)
  r$KR <-  getKR(r)
  return(r)
}

##' Get spherical coordinates of datapoints, transformed according to
##' the values of \code{DVflip} and \code{side}.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss}
##' @method getDss retinalReconstructedDataset
##' @author David Sterratt
getDss.retinalReconstructedDataset <- function(r) {
  Dss <- NextMethod()
  if (length(Dss) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Dss)) {
        Dss[[i]][,"lambda"] <- -Dss[[i]][,"lambda"]
      }
    }
  }
  return(Dss)
}


##' @title Get grouped variable with locations in spherical coordinates.
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Gss}
##' @method getGss retinalReconstructedDataset
##' @author David Sterratt
##' @export
getGss.retinalReconstructedDataset <- function(r) {
  Gss <- NextMethod()
  if (length(Gss) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Gss)) {
        Gss[[i]][,"lambda"] <- -Gss[[i]][,"lambda"]
      }
    }
  }
  return(Gss)
}

##' Get Karcher mean of datapoints in spherical coordinates,
##' transformed according to the values of \code{DVflip} and
##' \code{side}.
##'
##' @title Get transformed spherical coordinates of Karcher mean of
##' datapoints
##' @param r \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss.mean}
##' @method getDssMean retinalReconstructedDataset
##' @author David Sterratt
getDssMean.retinalReconstructedDataset <- function(r) {
  Dss.mean <- NextMethod()
  Dss.mean[["OD"]] <- NULL
  if (length(Dss.mean) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Dss.mean)) {
        Dss.mean[[i]][,"lambda"] <- -Dss.mean[[i]][,"lambda"]
      }
    }
  }
  return(Dss.mean)
}

##' Get spherical coordinates of landmarks, transformed according to
##' the values of \code{DVflip} and \code{side}.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss}
##' @method getSss retinalReconstructedDataset
##' @author David Sterratt
getSss.retinalReconstructedDataset <- function(r) {
  Sss <- NextMethod()
  if (length(Sss) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Sss)) {
        Sss[[i]][,"lambda"] <- -Sss[[i]][,"lambda"]
      }
    }
  }
  return(Sss)
}

##' This lablels the poles N, D, T and V
##' 
##' @title Polar plot of reconstructed dataset
##' @param r \code{RetinalReconstructedDataset} object
##' @param show.grid Whether or not to show the grid lines of
##' lattitude and longitude
##' @param grid.col Colour of the minor grid lines
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param ... Other graphics parameters.
##' @method polarplot retinalReconstructedDataset
##' @author David Sterratt
##' @export
polarplot.retinalReconstructedDataset <- function(r, show.grid=TRUE,
                                                   grid.col="gray",
                                                   grid.bg="transparent", 
                                                   grid.int.minor=15,
                                                   grid.int.major=45,  ...) {
  ## This will call polarplot.reconstructedDataset()
  NextMethod(flip.horiz=(r$side=="Left"),
             labels=c("N", "D", "T", "V"))
}


