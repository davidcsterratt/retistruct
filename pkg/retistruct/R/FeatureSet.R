##' Superclass containing functions and data relating to sets of
##' features in flat \code{\link{Outline}}s
##'
##' @description A FeatureSet contains information about features
##'   located on \code{\link{Outline}}s. Each FeatureSet contains a
##'   list of matrices, each of which has columns labelled \code{X}
##'   and \code{Y} describing the cartesian coordinates of points on
##'   the Outline, in the unscaled coordinate frame. Derived classes,
##'   e.g. a \code{\link{CountSet}}, may have extra columns. Each matrix
##'   in the list has an associated label and colour, which is used by
##'   plotting functions.
##'
##' @author David Sterratt
##' @export
FeatureSet <- R6Class("FeatureSet",
  inherit = FeatureSetCommon,
  public = list(
    ##' @description Constructor
    ##' @param data List of matrices describing data. Each matrix
    ##'   should have columns named \code{X} and \code{Y}
    ##' @param cols Named vector of colours for each data set. The name is
    ##'   used as the ID (label) for the data set. The colours should be names
    ##'   present in the output of the \code{\link{colors}} function
    ##' @param type String 
    initialize = function(data=NULL, cols=NULL, type=NULL) {
      if (!is.null(data)) {
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
      }
    }
  )
)
