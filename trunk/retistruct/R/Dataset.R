##' Constructor for a \code{dataset} object.
##'
##' @title Constructor for a \code{dataset} object.
##' @param o An \code{outline} object.
##' @param dataset The name of the dataset
##' @param Ds A list of data point sets, with each set being a 2
##' column matrix of X and Y coordinates of data point locations. Each
##' item in the list should be named. Elements with these names should
##' also be in the \code{cols} argument (see below).
##' @param Ss A list of landmarks. These do not need to be named. If
##' any elements are  named, the names should map onto an element in
##' the \code{cols} argument. Any elements that are named \code{""}
##' will be plotted using the default colour.
##' @param cols A list of colours in which to plot datapoints and landmarks.
##' @param raw A place to put raw data in whatever format is desired.
##' @param Gs  A list of grouped point sets, with each set being a 3
##' column matrix of X and Y coordinates and the value Z of the
##' variable at that point.  Each item in the list should be
##' named. Elements with these names should also be in the \code{cols}
##' argument.
##' @return A \code{dataset} object.
##' @author David Sterratt
Dataset <- function(o, dataset, Ds, Ss, cols, raw, Gs=NULL) {
  d <- o
  class(d) <- addClass("dataset", o)
  d$dataset <- dataset
  d$Ds <- Ds
  d$Ss <- Ss
  if (is.null(names(Ss))) {
    names(d$Ss) <- rep("", length(d$Ss))
  }
  d$cols <- cols
  d$raw <- raw
  d$Gs <- Gs
  return(d)
}

##' Name a landmark in a \code{dataset}. The name of element \code{i}
##' of \code{Ss} is set to \code{name}, the name of any element that
##' bore the name is set to "" and all other elements are unaltered.
##'
##' @title Name a landmark in a Dataset
##' @param d \code{dataset} object
##' @param i index of landmark to name
##' @param name name to give landmark
##' @return New \code{dataset} object in which landmark is named
##' @author David Sterratt
##' @export
nameLandmark <- function(d, i, name) {
  if (!is.na(i)) {
    new.names <- names(d$Ss)
    ## If this name already exists, replace it with ""
    j <- getLandmarkID(d, name)
    if (!is.na(j)) {
      new.names[j] <- ""
    }
    new.names[i] <- name
    names(d$Ss) <- new.names
  }
  return(d)
}

getLandmarkID <- function(d, name) {
  id <- which(names(d$Ss) == name)
  if (length(id) == 1) {
    return(id)
  } else {
    return(NA)
  }
}

##' @title Flat plot of Dataset
##' @param x \code{\link{Dataset}} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param datapoints If \code{TRUE}, display data points.
##' @param grouped If \code{TRUE}, dipslay grouped data.
##' @param landmarks If \code{TRUE}, dipslay landmarks.
##' @param ... Graphical parameters to pass to plotting functions
##' @method flatplot dataset
##' @author David Sterratt
##' @export
flatplot.dataset <- function(x, axt="n", ylim=NULL,
                        datapoints=TRUE,
                        grouped=FALSE,
                        landmarks=TRUE,
                        ...) {
  NextMethod()
  ## flatoutline(d, axt=axt, ylim=ylim, ...)
  if (datapoints) {
    with(x, {
      for (name in names(Ds)) {
        suppressWarnings(points(Ds[[name]][,1], Ds[[name]][,2],
                                col=cols[[name]], pch=20, ...))
      }
    })
  }
  
  if (grouped) {
    with(x, {
      for (name in names(Gs)) {
        suppressWarnings(text(Gs[[name]][,1], Gs[[name]][,2], Gs[[name]][,3],
                              col=cols[[name]],  ...))
      }
    })
  }

  if (landmarks) {
    with(x, {
      if (length(Ss) > 0) {
        for (i in 1:length(Ss)) {
          name <- names(Ss)[i]
          col <- ifelse(is.null(name) || (name==""), "default", name)
          suppressWarnings(lines(Ss[[i]][,1], Ss[[i]][,2],
                                 col=cols[[col]], ...))
        }
      }
    })
  }
}


