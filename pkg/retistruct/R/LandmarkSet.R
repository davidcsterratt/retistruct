##' Subclass of \code{\link{FeatureSet}} to represent points
##'
##' @description A LandmarkSet contains information about points
##'   located on \code{\link{Outline}}s. Each LandmarkSet contains a
##'   list of matrices, each of which has columns labelled \code{X}
##'   and \code{Y} describing the cartesian coordinates (in the
##'   unscaled coordinate frame) of points in landmarks on the
##'   Outline.
##'
##' @author David Sterratt
##' @export
LandmarkSet <- R6Class("LandmarkSet",
  inherit = FeatureSet,
  public = list(
    ##' @description Constructor
    ##' @param data List of matrices describing data. Each matrix
    ##'   should have columns named \code{X} and \code{Y}
    ##' @param cols Named vector of colours for each data set. The name is
    ##'   used as the ID (label) for the data set. The colours should be names
    ##'   present in the output of the \code{\link{colors}} function
    initialize = function(data=NULL, cols=NULL) {
      super$initialize(data, cols, "LandmarkSet")
    },
    ##' @description Map the LandmarkSet to a \code{\link{ReconstructedOutline}}
    ##' @param ro The \code{\link{ReconstructedOutline}}
    reconstruct = function(ro) {
      return(ReconstructedLandmarkSet$new(self, ro))
    }
   
  )
)

flatplot.LandmarkSet <- function(x, ids=x$getIDs(), ...) {
  for (id in ids) {
    if (!is.null(x$getFeature(id))) {
      lines(x$getFeature(id)[,"X"], x$getFeature(id)[,"Y"],
            col=x$getCol(id))
    }
  }
}
