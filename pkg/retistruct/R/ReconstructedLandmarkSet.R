##' Class containing functions and data to map \link{LandmarkSet}s to
##' \link{ReconstructedOutline}s
##'
##' @description A ReconstructedLandmarkSet contains information about
##'   features located on \code{\link{ReconstructedOutline}}s. Each
##'   ReconstructedLandmarkSet contains a list of matrices, each of
##'   which has columns labelled \code{phi} (latitude) and
##'   \code{lambda} (longitude) describing the spherical coordinates
##'   of points on the ReconstructedOutline.
##' 
##' @author David Sterratt
##' @importFrom geometry delaunayn
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
           markup=NULL, # Not used in this function; hides from ...
           max.proj.dim=NULL, # Not used in this function; hides from ...
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


