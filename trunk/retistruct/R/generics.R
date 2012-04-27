##' This adds a new class to the class vector. If the class is dataset
##' type it is prepended at the start of the class vector. If it is an
##' outline type, it is put after all the dataset classes, but before
##' the other outline ones. This is needed for the plotting functions
##' to work properly.
##'
##' @title Add new class to class vector of object
##' @param newclass New class to add
##' @param obj Object to which to add it
##' @return New class vector
##' @author David Sterratt
addClass <- function(newclass, obj) {
  if (newclass %in% class(obj)) {
    return(class(obj))
  }
  cl <- unique(c(newclass, class(obj)))
  cld <- cl[grep("[Dd]ataset", cl)]
  clo <- cl[setdiff(1:length(cl), grep("[Dd]ataset", cl))]
  return(c(cld, clo))
}

##' Plot flat representation of object
##'
##' @title Flat plot of object
##' @param x \code{\link{Outline}}, \code{\link{Dataset}} \&c object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
flatplot <- function(x, axt="n", ylim=NULL, ...) {
  UseMethod("flatplot")
}

##' @export
flatplot.default <- function(x, axt="n", ylim=NULL, ...) {
}

##' Plot polar representation of object
##'
##' @title Polar plot of object
##' @param r \code{\link{ReconstructedOutline}},
##' \code{\link{ReconstructedDataset}}  object
##' @param show.grid Whether or not to show the grid lines of lattitude and longitude
##' @param grid.col Colour of the minor grid lines
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param flip.horiz Wether to flip about a horizontal axis
##' @param labels Vector of 4 labels to plot at 0, 90, 180 and 270 degrees 
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
polarplot <- function(r, show.grid=TRUE,
                       grid.col="gray", grid.bg="transparent", 
                       grid.int.minor=15, grid.int.major=45,
                       flip.horiz=FALSE,
                       labels=c(0, 90, 180, 270),...) {
  UseMethod("polarplot")
}

##' @export
polarplot.default <- function(r, show.grid=TRUE,
                               grid.col="gray", grid.bg="transparent", 
                               grid.int.minor=15, grid.int.major=45,
                               flip.horiz=FALSE, labels=c(0, 90, 180, 270),
                               ...) {
  plot.new()
}

##' Draw a sinusoidal plot of a reconstructed outline. This method just sets
##' up the grid lines and the angular labels. 
##'
##' @title Sinusoidal projection
##' @param r \code{ReconstructedOutline} object
##' @param show.grid Whether or not to show the grid lines of lattitude and longitude
##' @param grid.col Colour of the minor grid lines
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param flip.horiz Wether to flip about a horizontal axis
##' @param labels Vector of 4 labels to plot at 0, 90, 180 and 270 degrees 
##' @param ... Other parameters, including graphics ones. The option
##' \code{image} causes an image to be plotted if \code{TRUE}
##' (default \code{TRUE}).  The option \code{preserve.area} creates an
##' area-preserving plot (default \code{FALSE}).
##' @author David Sterratt
##' @export
sinusoidalplot <- function(r, show.grid=TRUE,
                           grid.col="gray",
                           grid.bg="transparent", 
                           grid.int.minor=15,
                           grid.int.major=45,
                           flip.horiz=FALSE,
                           labels=c(0, 90, 180, 270), ...) {
  UseMethod("sinusoidalplot")
}

##' @export
sinusoidalplot.default <- function(r, show.grid=TRUE,
                           grid.col="gray",
                           grid.bg="transparent", 
                           grid.int.minor=15,
                           grid.int.major=45,
                           flip.horiz=FALSE,
                           labels=c(0, 90, 180, 270), ...) {
  plot.new()
}


##' @export
sphericalplot <- function(r, ...) {
  UseMethod("sphericalplot")
}

##' @export
sphericalplot.default <- function(r, ...) {
  rgl.clear()
  rgl.bg(color="white")
}

##' @export
lvsLplot <- function(r) {
  UseMethod("lvsLplot")
}

##' @export
lvsLplot.default <- function(r) {
}

##' Get spherical coordinates of datapoints.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss}
##' @author David Sterratt
##' @export
getDss <- function(r) {
  UseMethod("getDss")
}

##' @export
getDss.default <- function(r) {
  return(NULL)
}

##' Get Karcher mean of datapoints in spherical coordinates.
##'
##' @title Karcher mean of datapoints in spherical coordinates
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss.mean}
##' @author David Sterratt
##' @export
getDssMean <- function(r) {
  UseMethod("getDssMean")
}

##' @export
getDssMean.default <- function(r) {
  return(NULL)
}

##' @title Get grouped variable with locations in spherical coordinates.
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Gss}
##' @author David Sterratt
##' @export
getGss <- function(r) {
  UseMethod("getGss")
}

##' @export
getGss.default <- function(r) {
  return(NULL)
}


##' Get spherical coordinates of landmarks.
##'
##' @title Get transformed spherical coordinates of landmarks.
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Sss}
##' @author David Sterratt
##' @export
getSss <- function(r) {
  UseMethod("getSss")
}

##' @export
getSss.default <- function(r) {
  return(NULL)
}

##' @export
getTss <- function(r) {
  UseMethod("getTss")
}

##' @export
getTss.default <- function(r) {
  return(NULL)
}

##' @title Get coordinates of corners of pixels of image in spherical coordinates 
##' @param r \code{\link{ReconstructedOutline}} or
##' \code{\link{RetinalReconstructedOutline}} object
##' @return Coordinates of corners of pixels in spherical coordinates 
##' @author David Sterratt
##' @export
getIms <- function(r) {
  UseMethod("getIms")
}

##' @export
getIms.default <- function(r) {
  return(NULL)
}

