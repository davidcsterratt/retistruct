##' Constructor for a \code{dataset} object.
##'
##' @title Constructor for a \code{dataset} object.
##' @param o An \code{outline} object.
##' @param dataset The name of the dataset
##' @param Ds A list of datapoints. Each item in the list should be
##' named. Elements with these names should also be in the \code{cols}
##' argument (see below).
##' @param Ss A list of landmarks. These do not need to be named. If
##' any elements are  named, the names should map onto an element in
##' the \code{cols} argument. Any elements that are named \code{""}
##' will be plotted using the default colour.
##' @param cols A list of colours in which to plot datapoints and landmarks.
##' @param raw A place to put raw data in whatever format is desired.
##' @return A \code{dataset} object.
##' @author David Sterratt
Dataset <- function(o, dataset, Ds, Ss, cols, raw) {
  d <- o
  class(d) <- c("dataset", class(o))
  d$dataset <- dataset
  d$Ds <- Ds
  d$Ss <- Ss
  if (is.null(names(Ss))) {
    names(d$Ss) <- rep("", length(d$Ss))
  }
  d$cols <- cols
  d$raw <- raw
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

plot.flat.dataset <- function(d, axt="n", ylim=NULL, ...) {
  args <- list(...)
  plot.datapoints <- is.null(args$datapoints) || args$datapoints
  plot.landmarks <- is.null(args$landmarks) || args$landmarks

  NextMethod()
  ## plot.flat.outline(d, axt=axt, ylim=ylim, ...)
  if (plot.datapoints) {
    with(d, {
      for (col in names(Ds)) {
        suppressWarnings(points(Ds[[col]][,1], Ds[[col]][,2],
                                col=cols[[col]], pch=20, ...))
      }
    })
  }
  if (plot.landmarks) {
    with(d, {
      if (length(Ss) > 0) {
        for (i in 1:length(Ss)) {
          name <- names(Ss)[i]
          col <- ifelse(is.null(name) || (name==""), "default", name)
          suppressWarnings(lines(Ss[[i]][,1], Ss[[i]][,2], col=cols[[col]], ...))
        }
      }
    })
  }
}


