##' Construct an AnnotatedDataset that contains information specific
##' to the dataset in question. 
##'
##' @title AnnotatedDataset constructor
##' @param d A \code{dataset} object
##' @return An \code{annotatedDataset} object. This contains all the
##' information from \code{d} plus:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' @author David Sterratt
AnnotatedDataset <- function(d) {
  a <- d
  class(a) <- c("annotatedDataset", class(a))
  a$DVflip <- FALSE
  a$side <- "Right"
  return(a)
}

##' Plot an annotated dataset. This basically is equivalent to
##' plotting a \code{dataset}, but may perform some transformations to
##' the date or plotting parameters. At present, if \code{DVflip} is
##' \code{TRUE}, it flips the \emph{y}-axis.
##'
##' @title Flat plot of annotated dataset
##' @param a \code{annotatedDataset} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param ... Other plotting parameters
##' @method plot.flat annotatedDataset
##' @author David Sterratt
plot.flat.annotatedDataset <- function(a, axt="n", ylim=NULL, ...) {
  if (a$DVflip) {
    if (is.null(ylim)) {
      ylim <- c(max(a$P[,2]), min(a$P[,2]))
    } else {
      ylim <- sort(ylim, TRUE)
    }
  }
  NextMethod(ylim=ylim)
}
