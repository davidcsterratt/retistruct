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
    ro = NULL,
    fs = NULL,
    initialize = function(fs, ro) {
      self$cols <- fs$cols
      self$ro <- ro
      self$type <- paste0("Reconstructed", fs$type)
      ro$report(paste("Inferring coordinates of", fs$type))
      if (!is.null(fs$data) & (length(fs$data) > 0)) {

        ## Meshpoints in Cartesian coordinates
        Ptc <- sph2cart(theta=ro$lambda, phi=ro$phi, r=1)

        for (name in names(fs$data)) {
          Pb <- tsearch(ro$ol$getPoints()[,"X"],
                        ro$ol$getPoints()[,"Y"],
                        ro$ol$T,
                        fs$data[[name]][,"X"],
                        fs$data[[name]][,"Y"], bary=TRUE)
          oo <- is.na(Pb$idx)           # Points outwith outline
          if (any(oo)) {
            warning(paste(sum(oo), name, "datapoints outwith the outline will be ignored."))
          }
          Pb$p   <- Pb$p[!oo,,drop=FALSE]
          Pb$idx <- Pb$idx[!oo]
          self$Ps[[name]] <- bary2sph(Pb, ro$Tt, Ptc)
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
      return(self$ro$featureSetTransform(self$Ps[[name]]))
    },
    getCol = function(id) {
      if (id %in% names(self$cols)) {
        return(self$cols[[id]])
      }
      return(self$cols[["default"]])
    }
  )
)
