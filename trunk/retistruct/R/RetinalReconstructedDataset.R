##' Create an object that is specific to retinal datasets. This
##' contains methods that return datapoint and landmark coordinates
##' that have been transformed according to the values of
##' \code{DVflip} and \code{side}.
##'
##' @title RetinalReonstructedDataset constructor
##' @param r Object that inherits both \code{reconstructedDataset} and
##' \code{dataset}.
##' @return \code{retinalReonstructedDataset} object. This does not
##' contain any extra fields, but there are extra mthods dthat apply
##' to it.
##' @author David Sterratt
RetinalReconstructedDataset <- function(r) {
  if (!(inherits(r, "reconstructedDataset") &
        inherits(r, "retinalDataset"))) {
    stop("Argument needs to inherit reconstructedDataset and retinalDataset")
  }
  class(r) <- addClass("retinalReconstructedDataset", r)
  return(r)
}

##' Get spherical coordinates of datapoints, transformed according to
##' the values of \code{DVflip} and \code{side}.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{retinalReonstructedDataset} object.
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
    if (r$side=="Left") {
      for (i in 1:length(Dss)) {
        Dss[[i]][,"lambda"] <- 2*pi - Dss[[i]][,"lambda"]
      }
    }
  }
  return(Dss)
}

##' Get Karcher mean of datapoints in spherical coordinates,
##' transformed according to the values of \code{DVflip} and
##' \code{side}.
##'
##' @title Get transformed spherical coordinates of Karcher mean of
##' datapoints
##' @param r \code{retinalReonstructedDataset} object.
##' @return \code{Dss.mean}
##' @method getDss.mean retinalReconstructedDataset
##' @author David Sterratt
getDss.mean.retinalReconstructedDataset <- function(r) {
  Dss.mean <- NextMethod()
  Dss.mean[["OD"]] <- NULL
  if (length(Dss.mean) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Dss.mean)) {
        Dss.mean[[i]][,"lambda"] <- -Dss.mean[[i]][,"lambda"]
      }
    }
    if (r$side=="Left") {
      for (i in 1:length(Dss.mean)) {
        Dss.mean[[i]][,"lambda"] <- 2*pi - Dss.mean[[i]][,"lambda"]
      }
    }
  }
  return(Dss.mean)
}

##' Get spherical coordinates of landmarks, transformed according to
##' the values of \code{DVflip} and \code{side}.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{retinalReonstructedDataset} object.
##' @return \code{Dss}
##' @method getDss retinalReconstructedDataset
##' @author David Sterratt
getSss.retinalReconstructedDataset <- function(r) {
  Sss <- NextMethod()
  if (length(Sss) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Sss)) {
        Sss[[i]][,"lambda"] <- -Sss[[i]][,"lambda"]
      }
    }
    if (r$side=="Left") {
      for (i in 1:length(Sss)) {
        Sss[[i]][,"lambda"] <- 2*pi - Sss[[i]][,"lambda"]
      }
    }
  }
  return(Sss)
}

##' Get spherical coordinates of datapoints, transformed according to
##' the values of \code{DVflip} and \code{side}.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{retinalReonstructedDataset} object.
##' @return \code{Dss}
##' @method getDss retinalReconstructedDataset
##' @author David Sterratt
getTss.retinalReconstructedDataset <- function(r) {
  Tss <- NextMethod()
  if (r$DVflip) {
    for (i in 1:length(Tss)) {
      Tss[[i]][,"lambda"] <- -Tss[[i]][,"lambda"]
    }
  }
  if (r$side=="Left") {
    for (i in 1:length(Tss)) {
      Tss[[i]][,"lambda"] <- 2*pi - Tss[[i]][,"lambda"]
    }
  }
  return(Tss)
}

plot.polar.retinalReconstructedDataset <- function(r, show.grid=TRUE,
                                                   grid.col="gray",
                                                   grid.bg="transparent", 
                                                   grid.int.minor=15,
                                                   grid.int.major=45,
                                                   flip.horiz=FALSE, ...) {
  ## This will call plot.polar.reconstructedDataset()
  NextMethod(flip.horiz=(r$side=="Left"),
             labels=c("N", "D", "T", "V"))
}


