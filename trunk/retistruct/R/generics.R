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
  cl <- unique(c(newclass, class(obj)))
  cld <- cl[grep("[Dd]ataset", cl)]
  clo <- cl[setdiff(1:length(cl), grep("[Dd]ataset", cl))]
  return(c(cld, clo))
}

##' Plot flat object
##'
##' @title Flat plot of object
##' @param o \code{\link{Outline}}, \code{\link{Dataset}} \&c object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
plot.flat <- function(x, axt="n", ylim=NULL, ...) {
  UseMethod("plot.flat")
}

plot.flat.default <- function(x, axt="n", ylim=NULL, ...) {
}

plot.polar <- function(r, show.grid=TRUE,
                       grid.col="gray", grid.bg="transparent", 
                       grid.int.minor=15, grid.int.major=45, ...) {
  UseMethod("plot.polar")
}

plot.polar.default <- function(r, show.grid=TRUE,
                               grid.col="gray", grid.bg="transparent", 
                               grid.int.minor=15, grid.int.major=45,
                               flip.horiz=FALSE, ...) {
  plot.new()
}

plot.spherical <- function(r, ...) {
  UseMethod("plot.spherical")
}

plot.spherical.default <- function(r, ...) {
  rgl.clear()
  rgl.bg(color="white")
}

plot.l.vs.L <- function(r) {
  UseMethod("plot.l.vs.L")
}

plot.l.vs.L.default <- function(r) {
}

##' Get spherical coordinates of datapoints.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{\link{reconstructedDataset}} or \code{\link{retinalReconstructedDataset}} object.
##' @return \code{Dss}
##' @author David Sterratt
##' @export
getDss <- function(r) {
  UseMethod("getDss")
}

getDss.default <- function(r) {
  return(NULL)
}

##' Get Karcher mean of datapoints in spherical coordinates.
##'
##' @title Karcher mean of datapoints in spherical coordinates
##' @param r \code{\link{reconstructedDataset}} or \code{\link{retinalReconstructedDataset}} object.
##' @return \code{Dss.mean}
##' @author David Sterratt
##' @export
getDss.mean <- function(r) {
  UseMethod("getDss.mean")
}

getDss.mean.default <- function(r) {
  return(NULL)
}

##' Get spherical coordinates of landmarks.
##'
##' @title Get transformed spherical coordinates of landmarks.
##' @param r \code{\link{reconstructedDataset}} or \code{\link{retinalReconstructedDataset}} object.
##' @return \code{Sss}
##' @author David Sterratt
##' @export
getSss <- function(r) {
  UseMethod("getSss")
}

getSss.default <- function(r) {
  return(NULL)
}

getTss <- function(r) {
  UseMethod("getTss")
}

getTss.default <- function(r) {
  return(NULL)
}

##' @title Get coordinates of corners of pixels of image in spherical coordinates 
##' @param r \code{\link{reconstructedOutline}} or \code{\link{retinalReconstructedOutline}} object
##' @return Coordinates of corners of pixels in spherical coordinates 
##' @author David Sterratt
##' @export
getIms <- function(r) {
  UseMethod("getIms")
}

getIms.default <- function(r) {
  return(NULL)
}

