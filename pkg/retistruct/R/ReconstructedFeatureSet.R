##' Class containing functions and data to map \link{FeatureSet}s to
##' \link{ReconstructedOutline}s
##'
##' @description A ReconstructedFeatureSet contains information about
##'   features located on \code{\link{ReconstructedOutline}}s. Each
##'   ReconstructedFeatureSet contains a list of matrices, each of
##'   which has columns labelled \code{phi} (latitude) and
##'   \code{lambda} (longitude) describing the spherical coordinates
##'   of points on the ReconstructedOutline. Derived classes, e.g. a
##'   \code{\link{ReconstructedCountSet}}, may have extra columns.
##'   Each matrix in the list has an associated label and colour,
##'   which is used by plotting functions.
##'
##' @author David Sterratt
##' @export
ReconstructedFeatureSet <- R6Class("ReconstructedFeatureSet",
  inherit = FeatureSetCommon,
  public = list(
    ##' @description Constructor
    ##' @param fs \code{\link{FeatureSet}} to reconstruct
    ##' @param ro \code{\link{ReconstructedOutline}} to which feature
    ##'   set should be mapped
    initialize = function(fs=NULL, ro=NULL) {
      if (!is.null(fs)) {
        self$cols <- fs$cols
        self$type <- paste0("Reconstructed", fs$type)
        report(paste("Inferring coordinates of", fs$type))
        if (!is.null(fs$data) & (length(fs$data) > 0)) {
          for (name in names(fs$data)) {
            self$data[[name]] <- ro$mapFlatToSpherical(fs$data[[name]])
          }
        }
      }
    }
  )
)
