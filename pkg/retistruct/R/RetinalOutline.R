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
  NextMethod(grid=FALSE)

  ## Plot feature sets
  if (datapoints) {
    fs <- x$getFeatureSet("PointSet")
    if (!is.null(fs)) {
      flatplot(fs, ...)
    }
  }
  
  ## if (grouped) {
  ##   Gs <- x$Gs
  ##   for (id in ids) {
  ##     if (!is.null(Gs[[id]])) {
  ##       if (nrow(Gs[[id]]) > 0) {
  ##         text(Gs[[id]][,1], Gs[[id]][,2], Gs[[id]][,3],
  ##              col=x$cols[[id]])
  ##       }
  ##     }
  ##   }
  ## }

  if (landmarks) {
    fs <- x$getFeatureSet("LandmarkSet")
    if (!is.null(fs)) {
      flatplot(fs, ...)
    }
  }
}

