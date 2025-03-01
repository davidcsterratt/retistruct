##' Class containing functions and data relating to retinal outlines
##'
##' @description In addition to fields inherited from
##'   \link{StitchedOutline}, a RetinalOutline contains a
##'   \code{dataset} field, describing the system path to dataset
##'   directory and metadata specific to retinae and some formats of
##'   retinae
##'
##' @description An \code{retinalOutline} object. This contains the following fields:
##' @author David Sterratt
##' @export
RetinalOutline <- R6Class("RetinalOutline",
  inherit = StitchedOutline,
  public = list(
    ##' @field DVflip \code{TRUE} if the raw data is flipped in
    ##'   the dorsoventral direction
    DVflip = FALSE,
    ##' @field side The side of the eye (\dQuote{Left} or \dQuote{Right})
    side = "Right",
    ##' @field dataset File system path to dataset directory
    dataset = NULL,
    ##' @description Constructor
    ##' @param ... Parameters to superclass constructors
    ##' @param dataset File system path to dataset directory
    initialize = function(..., dataset=NULL) {
      super$initialize(...)
      self$dataset <- dataset
    }
  )
)

##' @method flatplot RetinalOutline
##' @export
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
