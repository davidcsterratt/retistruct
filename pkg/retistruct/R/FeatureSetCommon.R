##' Class containing functionality common to \code{\link{FeatureSet}}s and
##'   \code{\link{ReconstructedFeatureSet}}s
##'
##' @description An FeatureSetCommon has functionality for retrieving
##'   sets of features (e.g. points or landmarks associated with an
##'   outline)
##' @author David Sterratt
##' @export
FeatureSetCommon <- R6Class("FeatureSetCommon",
  public = list(
    ##' @field data List of matrices describing data
    data = NULL,
    ##' @field cols Vector of colours for each data set
    cols = NA,
    ##' @field type String giving type of feature set
    type = NA,
    ##' @description Get numeric index of features
    ##' @param fid Feature ID (string)
    getIndex = function(fid) {
      i <- which(names(self$data) == fid)
      if (length(i) == 1) {
        return(i)
      } else {
        return(NA)
      }
    },
    ##' @description Get IDs of features
    ##' @return Vector of IDs of features
    getIDs = function() {
      return(names(self$data))
    },
    ##' @description Set name
    ##' @param i Numeric index of feature
    ##' @param fid Feature ID (string)
    setID = function(i, fid) {
      if (!is.na(i)) {
        new.fids <- names(self$data)
        ## If this name already exists, replace it with ""
        j <- self$getIndex(fid)
        if (!is.na(j)) {
          new.fids[j] <- ""
        }
        new.fids[i] <- fid
        names(self$data) <- new.fids
      }
    },
    ##' @description Get feature by feature ID
    ##' @param fid Feature ID string
    ##' @return Matrix describing feature
    getFeature = function(fid) {
      if (is.na(self$getIndex(fid))) {
        return(NULL)
      }
      return(self$data[[fid]])
    },
    ##' @description Get all features
    getFeatures = function() {
      return(self$data)
    },
    ##' @description Get colour in which to plot feature ID
    ##' @param fid Feature ID string
    getCol = function(fid) {
      if (fid %in% names(self$cols)) {
        return(self$cols[[fid]])
      }
      return(self$cols[["default"]])
    }
  )
)
