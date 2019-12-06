##' ReconstructedFeatureSet class
##' Initialised with a FeatureSet (\code{fs}) and an ReconstructedOutline (\code{ro}) object
##' @return A \code{ReconstructedFeatureSet} object.
##' @author David Sterratt
##' @export
ReconstructedFeatureSet <- R6Class("ReconstructedFeatureSet",
  public = list(
    Ps = NULL,
    cols = NA,
    type = NA,
    fs = NULL,
    initialize = function(fs=NULL, ro=NULL) {
      if (!is.null(fs)) {
        self$cols <- fs$cols
        self$type <- paste0("Reconstructed", fs$type)
        report(paste("Inferring coordinates of", fs$type))
        if (!is.null(fs$data) & (length(fs$data) > 0)) {
          for (name in names(fs$data)) {
            self$Ps[[name]] <- ro$mapFlatToSpherical(fs$data[[name]])
          }
        }
      }
    },
    getID = function(name) {
      id <- which(names(self$Ps) == name)
      if (length(id) == 1) {
        return(id)
      } else {
        return(NA)
      }
    },
    getIDs = function() {
      return(names(self$Ps))
    },
    getFeature = function(name) {
      if (is.na(self$getID(name))) {
        return(NULL)
      }
      return(self$Ps[[name]])
    },
    getCol = function(id) {
      if (id %in% names(self$cols)) {
        return(self$cols[[id]])
      }
      return(self$cols[["default"]])
    }
  )
)
