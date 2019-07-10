##' RetinalOutline class
##' @return An \code{retinalOutline} object. This contains the following fields:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' \item{\code{dataset}}{File system path to dataset}
##' @author David Sterratt
##' @export
RetinalOutline <- R6Class("RetinalOutline",
  inherit = StitchedOutline,
  public = list(
    DVflip = FALSE,
    side = "Right",
    dataset = NULL,
    initialize = function(..., dataset=NULL) {
      super$initialize(...)
      self$dataset <- dataset
    }
  )
)

flatplot.RetinalOutline <- function(x, axt="n", ylim=NULL,
                             datapoints=TRUE,
                             grouped=FALSE,
                             landmarks=TRUE,
                             ...) {
  ## This will call projection.reconstructedOutline(), and hence
  ## Outline(), but without drawing a grid.  The grid will be drawn
  ## later, after all the data has appeared.
  if (x$DVflip) {
    if (is.null(ylim)) {
      ylim <- c(max(x$getPoints()[,"Y"]), min(x$getPoints()[,"Y"]))
    } else {
      ylim <- sort(ylim, TRUE)
    }
  }
  NextMethod(grid=FALSE, ylim=ylim)

  ## Plot feature sets
  if (datapoints) {
    fs <- x$getFeatureSet("PointSet")
    if (!is.null(fs)) {
      flatplot(fs, ...)
    }
  }
  
  if (grouped) {
    fs <- x$getFeatureSet("CountSet")
    if (!is.null(fs)) {
      flatplot(fs, ...)
    }
  }

  if (landmarks) {
    fs <- x$getFeatureSet("LandmarkSet")
    if (!is.null(fs)) {
      flatplot(fs, ...)
    }
  }
}

