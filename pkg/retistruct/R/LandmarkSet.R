##' LandmarkSet class
##' @return An \code{LandmarkSet} object. This contains the following fields:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' \item{\code{dataset}}{File system path to dataset}
##' @author David Sterratt
##' @export
LandmarkSet <- R6Class("LandmarkSet",
  inherit = FeatureSet,
  public = list(
    initialize = function(data, cols) {
      super$initialize(data, cols, "LandmarkSet")
    },
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
