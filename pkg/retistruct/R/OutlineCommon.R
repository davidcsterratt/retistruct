##' Class containing functionality common to flat and reconstructed outlines
##'
##' @description An OutlineCommon has functionality for retrieving
##'   sets of features (e.g. points or landmarks associated with an
##'   outline)
OutlineCommon <- R6Class("OutlineCommon",
  public = list(
    ##' @field version Version of reconstruction file data format
    version = 6,
    ##' @field featureSets List of feature sets associated with the outline, which may be of various types, e.g. a \link{PointSet} or \link{LandmarkSet}
    featureSets = list(),
    ##' @description Get all the feature sets
    ##' @return List of \link{FeatureSet}s associated with the outline
    getFeatureSets = function() {
      return(self$featureSets)
    },
    ##' @description Get all feature sets of a particular type, e.g. \link{PointSet} or \link{LandmarkSet}
    ##' @param type The type of the feature set as a string
    ##' @return All \link{FeatureSet}s of that type
    getFeatureSet = function(type) {
      if (!(type %in% self$getFeatureSetTypes())) {
        return(NULL)
      }
      ind <- which(sapply(self$featureSets, function(x) {inherits(x, type)}))
      if (length(ind) > 1) {
        stop(paste("More than 1", type, "attached to Outline"))
      }
      return(self$featureSets[[ind]])
    },
    ##' @description Clear all feature sets from the outline
    clearFeatureSets = function() {
      self$featureSets = list()
    },
    ##' @description Get all the distinct IDs contained in the \link{FeatureSet}s
    ##' @return Vector of IDs
    getIDs = function() {
      return(unique(unlist(sapply(self$getFeatureSets(), function(fs) { fs$getIDs() }))))
    },
    ##' @description Get all the distinct types of \link{FeatureSet}s
    ##' @return Vector of types as strings, e.g. \var{PointSet}, \var{LandmarkSet}
    getFeatureSetTypes = function() {
      return(sapply(self$getFeatureSets(), function(fs) return(fs$type)))
    }
  )
)
  
