OutlineCommon <- R6Class("OutlineCommon",
  public = list(
    ## Version of reconstruction file data format
    version = 6,
    featureSets = list(),
    getFeatureSets = function() {
      return(self$featureSets)
    },
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
    clearFeatureSets = function() {
      self$featureSets = list()
    },
    getIDs = function() {
      return(unique(unlist(sapply(self$getFeatureSets(), function(fs) { fs$getIDs() }))))
    },
    getFeatureSetTypes = function() {
      return(sapply(self$getFeatureSets(), function(fs) return(fs$type)))
    }
  )
)
  
