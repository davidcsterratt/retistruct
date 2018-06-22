##' PointSet class
##' @return An \code{PointSet} object. This contains the following fields:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' \item{\code{dataset}}{File system path to dataset}
##' @author David Sterratt
##' @export
PointSet <- R6Class("PointSet",
  inherit = FeatureSet,
  public = list(
    initialize = function(data, cols) {
      super$initialize(data, cols, "PointSet")
    },
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
