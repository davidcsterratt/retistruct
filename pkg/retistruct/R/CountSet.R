##' Subclass of \code{\link{FeatureSet}} to represent counts centred
##' on points
##'
##' @description A CountSet contains information about points located
##'   on \code{\link{Outline}}s. Each CountSet contains a list of
##'   matrices, each of which has columns labelled \code{X} and
##'   \code{Y} describing the cartesian coordinates (in the unscaled
##'   coordinate frame) of the centres of boxes in the Outline, and a
##'   column \code{C} representing the counts in those boxes.
##'
##' @author David Sterratt
##' @export
CountSet <- R6Class("CountSet",
  inherit = FeatureSet,
  public = list(
    ##' @description Constructor
    ##' @param data List of matrices describing data. Each matrix
    ##'   should have columns named \code{X}, \code{Y} and \code{C}
    ##' @param cols Named vector of colours for each data set. The name is
    ##'   used as the ID (label) for the data set. The colours should be names
    ##'   present in the output of the \code{\link{colors}} function
    initialize = function(data=NULL, cols=NULL) {
      if (!is.null(data)) {
        if (!all(sapply(data, function(d) (all(c("C") %in% colnames(d)))))) {
          stop("data argument to CountSet needs column marked C")
        }
        if (!all(sapply(data, function(d) (ncol(d) == 3)))) {
          stop("Data must have 3 columns")
        }
        super$initialize(data, cols, "CountSet")
      }
    },
    ##' @description Map the CountSet to a \code{\link{ReconstructedOutline}}
    ##' @param ro The \code{\link{ReconstructedOutline}}
    reconstruct = function(ro) {
      return(ReconstructedCountSet$new(self, ro))
    }
  )
)

##' @method flatplot CountSet
##' @export
flatplot.CountSet <- function(x, ids=x$getIDs(), ...) {
  for (id in ids) {
    if (!is.null(x$data[[id]])) {
      if (nrow(x$data[[id]]) > 0) {
        text(x$data[[id]][,1], x$data[[id]][,2], x$data[[id]][,3],
             col=x$cols[[id]])
      }
    }
  }
}
