##' ReconstructedLandmarkSet class
##' @return An \code{ReconstructedLandmarkSet} object, which inherits
##'   from ReconstructedFeatureSet
##' @author David Sterratt
##' @export
ReconstructedLandmarkSet <- R6Class("ReconstructedLandmarkSet",
  inherit = ReconstructedFeatureSet
)

projection.ReconstructedLandmarkSet <-
  function(r,
           phi0,
           transform=identity.transform,
           ids=r$getIDs(),
           axisdir=cbind(phi=90, lambda=0),
           projection=azimuthal.equalarea,
           proj.centre=cbind(phi=0, lambda=0),
           lambdalim=c(-180, 180),
           ...)
{
  for (id in ids) {
    if (!is.null(r$getFeature(id))) {
      lines(projection(rotate.axis(transform(r$getFeature(id),
                                             phi0=phi0),
                                   axisdir*pi/180),
                       lines=TRUE,
                       lambdalim=lambdalim*pi/180,
                       proj.centre=pi/180*proj.centre),
            col=r$getCol(id), ...)
    }
  }
}


