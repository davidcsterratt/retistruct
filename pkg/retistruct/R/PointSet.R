##' Subclass of \code{\link{FeatureSet}} to represent points
##'
##' @description A PointSet contains information about points located
##'   on \code{\link{Outline}}s. Each PointSet contains a list of
##'   matrices, each of which has columns labelled \code{X} and
##'   \code{Y} describing the cartesian coordinates (in the unscaled
##'   coordinate frame) of points on the Outline.
##'
##' @author David Sterratt
##' @export
PointSet <- R6Class("PointSet",
  inherit = FeatureSet,
  public = list(
    ##' @description Constructor
    ##' @param data List of matrices describing data. Each matrix
    ##'   should have columns named \code{X} and \code{Y}
    ##' @param cols Named vector of colours for each data set. The name is
    ##'   used as the ID (label) for the data set. The colours should be names
    ##'   present in the output of the \code{\link{colors}} function
    initialize = function(data=NULL, cols=NULL) {
      if (!is.null(data)) {
        if (!all(sapply(data, function(d) (ncol(d) == 2)))) {
          stop("Data must have 2 columns")
        }
        super$initialize(data, cols, "PointSet")
      }
    },
    ##' @description Map the PointSet to a \code{\link{ReconstructedOutline}}
    ##' @param ro The \code{\link{ReconstructedOutline}}
    reconstruct = function(ro) {
      return(ReconstructedPointSet$new(self, ro))
    }
  )
)

flatplot.PointSet <- function(x, ids=x$getIDs(), ...) {
  for (id in ids) {
    if (!is.null(x$data[[id]])) {
      points(x$data[[id]][,1], x$data[[id]][,2],
             col=x$cols[[id]], pch=20)
    }
  }
}
