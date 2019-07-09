##' FeatureSet class
##' @return An \code{FeatureSet} object. This contains the following fields:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' \item{\code{dataset}}{File system path to dataset}
##' @author David Sterratt
##' @export
FeatureSet <- R6Class("FeatureSet",
  public = list(
    data = NA,
    cols = NA,
    type = NA,
    initialize = function(data, cols, type) {
      if (!is.list(data)) {
        stop("data should be a list")
      }
      if (!all(sapply(data, function(d) (all(c("X", "Y") %in% colnames(d)))))) {
        stop("data argument to FeatureSet needs columns marked X and Y")
      }
      self$data <- data
      self$cols <- cols
      if (!("default" %in% names(cols))) {
        cols[["default"]] <- "orange"
      }
      self$type <- type
      if (is.null(names(self$data)) & length(self$data) > 0) {
        names(self$data) <- paste0(type, 1:length(self$data))
      }
      missing_cols <- setdiff(cols, colors())
      if (length(missing_cols)) {
        warning(paste("Colour(s)", missing_cols, "missing from", type,
                      names(cols)[which(missing_cols == cols)]))
      }

    },
    getID = function(name) {
      id <- which(names(self$data) == name)
      if (length(id) == 1) {
        return(id)
      } else {
        return(NA)
      }
    },
    getIDs = function() {
      return(names(self$data))
    },
    setName = function(i, name) {
      if (!is.na(i)) {
        new.names <- names(self$data)
        ## If this name already exists, replace it with ""
        j <- self$getID(name)
        if (!is.na(j)) {
          new.names[j] <- ""
        }
        new.names[i] <- name
        names(self$data) <- new.names
      }
    },
    getFeature = function(name) {
      if (is.na(self$getID(name))) {
        return(NULL)
      }
      return(self$data[[name]][,c("X", "Y")])
    },
    getFeatures = function() {
      return(self$data)
    },
    getCol = function(id) {
      if (id %in% names(self$cols)) {
        return(self$cols[[id]])
      }
      return(self$cols[["default"]])
    }
  )
)
