##' Class containing functions and data to map \link{CountSet}s to
##' \link{ReconstructedOutline}s
##'
##' @description A ReconstructedCountSet contains information about
##'   features located on \code{\link{ReconstructedOutline}}s. Each
##'   ReconstructedCountSet contains a list of matrices, each of which
##'   has columns labelled \code{phi} (latitude) and \code{lambda}
##'   (longitude) describing the spherical coordinates of points on
##'   the ReconstructedOutline, and a column \code{C} representing the
##'   counts at these points.
##'
##' @author David Sterratt
##' @importFrom geometry delaunayn
##' @export
ReconstructedCountSet <- R6Class("ReconstructedCountSet",
  inherit = ReconstructedFeatureSet,
  public = list(
    ##' @field KR Kernel regression
    KR = NULL,
    ##' @description Constructor
    ##' @param fs \code{\link{FeatureSet}} to reconstruct
    ##' @param ro \code{\link{ReconstructedOutline}} to which feature
    ##'   set should be mapped
    initialize = function(fs=NULL, ro=NULL) {
      if (!is.null(fs)) {
        super$initialize(fs, ro)
        if (!is.null(fs$data) & (length(fs$data) > 0)) {
          for (name in names(fs$data)) {
            self$data[[name]] <- cbind(self$data[[name]], C=fs$data[[name]][,"C"])
          }
        }
      }
    },
    ##' @description Get kernel regression estimate of grouped data points
    ##' @return Kernel regression computed using
    ##'    \code{\link{compute.kernel.estimate}}
    getKR = function() {
      if (is.null(self$KR)) {
        yhat <- function(r, mu, kappa) {
          kr.yhat(r, mu[,1:2], mu[,3], kappa)
        }
        compute.conc <- function(mu) {
          kr.compute.concentration(mu[,1:2], mu[,3])
        }
        Gss <- list()
        for (id in self$getIDs()) {
          Gss[[id]] <- self$getFeature(id)
        }
        ## Get rid of datasets with very little data
        for (n in names(Gss)) {
          if (sum(Gss[[n]][,"C"]) <= 2) {
            Gss[[n]] <- NULL
          }
        }
        return(compute.kernel.estimate(Gss, self$ro$phi0, yhat, compute.conc))
      }
      return(self$KR)
    }
  )
)

##' @method projection ReconstructedCountSet
##' @export
projection.ReconstructedCountSet <-
  function(r,
           phi0,
           transform=identity.transform,
           ids=r$getIDs(),
           axisdir=cbind(phi=90, lambda=0),
           projection=azimuthal.equalarea,
           proj.centre=cbind(phi=0, lambda=0),
           markup=NULL, # Not used in this function; hides from ...
           max.proj.dim=NULL, # Not used in this function; hides from ...
           ...)
{
  ## This will call projection.reconstructedOutline(), but without
  ## drawing a grid.  The grid will be drawn later, after all the data
  ## as appeared..
  ## Datapoints
  for (id in ids) {
    if (!is.null(r$getFeature(id))) {
      if (nrow(r$getFeature(id)) > 0) {
        rc <- projection(rotate.axis(transform(r$getFeature(id)[,c("phi", "lambda")],
                                               phi0=r$phi0),
                                     axisdir*pi/180),
                         proj.centre=pi/180*proj.centre)

        text(rc[,"x"], rc[,"y"], r$getFeature(id)[,"C"],
             col=r$cols[[id]],
             ...)
      }
    }
  }
}
