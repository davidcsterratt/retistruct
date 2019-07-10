##' CountSet class
##' @return An \code{CountSet} object. This contains the following fields:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' \item{\code{dataset}}{File system path to dataset}
##' @author David Sterratt
##' @export
CountSet <- R6Class("CountSet",
  inherit = FeatureSet,
  public = list(
    initialize = function(data, cols) {
      if (!all(sapply(data, function(d) (all(c("C") %in% colnames(d)))))) {
        stop("data argument to CountSet needs column marked C")
      }
      super$initialize(data, cols, "CountSet")
    },
    reconstruct = function(ro) {
      return(ReconstructedCountSet$new(self, ro))
    }
  )
)

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
