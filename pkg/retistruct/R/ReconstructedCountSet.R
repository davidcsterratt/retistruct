##' ReconstructedCountSet class
##' @return An \code{ReconstructedCountSet} object. This contains the following fields:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' \item{\code{dataset}}{File system path to dataset}
##' @author David Sterratt
##' @importFrom geometry delaunayn
##' @export
ReconstructedCountSet <- R6Class("ReconstructedCountSet",
  inherit = ReconstructedFeatureSet,
  public = list(
    KR = NULL,
    initialize = function(fs, ro) {
      super$initialize(fs, ro)
      if (!is.null(fs$data) & (length(fs$data) > 0)) {
        for (name in names(fs$data)) {
          self$Ps[[name]] <- cbind(self$Ps[[name]], C=fs$data[[name]][,"C"])
        }
      }
    },
    ## Get kernel regression estimate of grouped data points
    ## @param r \code{\link{ReconstructedDataset}} object
    ## @return See \code{\link{compute.kernel.estimate}}
    ## @author David Sterratt
    ## @export
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

projection.ReconstructedCountSet <-
  function(r,
           phi0,
           transform=identity.transform,
           ids=r$getIDs(),
           axisdir=cbind(phi=90, lambda=0),
           projection=azimuthal.equalarea,
           proj.centre=cbind(phi=0, lambda=0),
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


