##' Create an object that is specific to retinal datasets. This
##' contains methods that return data point and landmark coordinates
##' that have been transformed according to the values of
##' \code{DVflip} and \code{side}.
##'
##' @title RetinalReconstructedOutline constructor
##' @return \code{RetinalReconstructedOutline} object. This does not
##' contain any extra fields, but there are extra methods that apply
##' to it.
##' @author David Sterratt
##' @export
RetinalReconstructedOutline <- R6Class("RetinalReconstructedOutline",
  inherit = ReconstructedOutline,
  public = list(
    EOD = NULL,
    fst = NULL,          # Transformed feature set
    ## @method getIms retinalReconstructedOutline
    ## @export
    getIms = function() {
      ims <- super$getIms()
      if (self$ol$DVflip) {
        if (!is.null(ims)) {
          ims[,"lambda"] <- -ims[,"lambda"]
        }
      }
      return(ims)
    },
    getTearCoords = function() {
      Tss <- super$getTearCoords()
      if (self$ol$DVflip) {
        for (i in 1:length(Tss)) {
          Tss[[i]][,"lambda"] <- -Tss[[i]][,"lambda"]
        }
      }
      return(Tss)
    },
    reconstruct = function(...) {
      super$reconstruct(...)
      OD <- self$getFeatureSet("LandmarkSet")$getFeature("OD")
      if (!is.null(OD)) {
        ODmean <- karcher.mean.sphere(OD)
        self$EOD <- 90 + ODmean["phi"]*180/pi
      }
    },
    getFeatureSet = function(type) {
      fs <- super$getFeatureSet(type)
      if (self$ol$DVflip) {
        if (is.null(self$fst)) {
          fst <- fs$clone()
          fst$Ps <-
            lapply(fs$Ps,
                   function(x) {
                     x[,"lambda"] <- -x[,"lambda"]
                     return(x)
                   })
        }
        return(fst)
      }
      return(fs)
    }
  )
)

##' Plot projection of reconstructed dataset
##' @param r \code{\link{RetinalReconstructedOutline}} object
##' @param transform Transform function to apply to spherical coordinates
##' before rotation
##' @param projection Projection in which to display object,
##' e.g. \code{\link{azimuthal.equalarea}} or \code{\link{sinusoidal}}
##' @param axisdir Direction of axis (North pole) of sphere in external space
##' @param proj.centre Location of centre of projection as matrix with
##' column names \code{phi} (elevation) and \code{lambda} (longitude).
##' @param lambdalim Limits of longitude (in degrees) to display
##' @param datapoints If \code{TRUE}, display data points
##' @param datapoint.means If \code{TRUE}, display Karcher mean of data points.
##' @param datapoint.contours If \code{TRUE}, display contours around
##' the data points generated using Kernel Density Estimation.
##' @param grouped If \code{TRUE}, display grouped data.
##' @param grouped.contours If \code{TRUE}, display contours around
##' the grouped data generated using Kernel Regression.
##' @param landmarks If \code{TRUE}, display landmarks.
##' @param mesh If \code{TRUE}, display the triangular mesh used in reconstruction
##' @param grid If \code{TRUE}, show grid lines
##' @param image If \code{TRUE}, show the reconstructed image
##' @param ids IDs of groups of data within a dataset, returned using
##' \code{getIDs}.
##' @param ... Graphical parameters to pass to plotting functions
##' @method projection RetinalReconstructedOutline
##' @export
projection.RetinalReconstructedOutline <-
  function(r,
           transform=identity.transform,
           projection=azimuthal.equalarea,
           axisdir=cbind(phi=90, lambda=0),
           proj.centre=cbind(phi=0, lambda=0),
           lambdalim=c(-180, 180),
           datapoints=TRUE,
           datapoint.means=TRUE,
           datapoint.contours=FALSE,
           grouped=FALSE,
           grouped.contours=FALSE,
           landmarks=TRUE,
           mesh=FALSE,
           grid=TRUE,
           image=TRUE,
           ids=r$getIDs(),
           ...) {
    philim <- c(-90, 90)
    colatitude <- FALSE
    pole <- TRUE
    if (!(identical(projection, sinusoidal) |
          identical(projection, orthographic))) {
      philim <- c(-90, r$ol$phi0*180/pi)
      colatitude <- TRUE
      pole <- FALSE
    }
    if (r$ol$side=="Right") {
      labels=c("N", "D", "T", "V")
    } else {
      labels=c("T", "D", "N", "V")
    }
    NextMethod(projection=projection,
               philim=philim,
               labels=labels,
               colatitude=TRUE,
               grid=FALSE,
               mesh=FALSE,
               image=image)


    ## Plot FeatureSets
    
    ## Datapoints
    if (datapoints) {
      message("Plotting points")
      fs <- r$getFeatureSet("PointSet")
      if (!is.null(fs)) {
        projection.ReconstructedPointSet(fs,
                                         phi0=r$phi0,
                                         ids=ids,
                                         transform=transform,
                                         axisdir=axisdir,
                                         projection=projection,
                                         proj.centre=proj.centre,
                                         ...)
      }
    }

    ## Mean datapoints
    if (datapoint.means) {
      message("Plotting point means")
      fs <- r$getFeatureSet("PointSet")
      if (!is.null(fs)) {
        Dss.mean <- fs$getMean()
        for (id in ids) {
          if (!is.null(Dss.mean[[id]])) {
            points(projection(rotate.axis(transform(Dss.mean[[id]],
                                                    phi0=r$phi0),
                                          axisdir*pi/180),
                              proj.centre=pi/180*proj.centre),
                   bg=fs$cols[[id]], col="black",
                   pch=23, cex=1.5)
          }
        }
      }
    }

    ## Count sets, formerly known as groups 
    if (grouped) {
      message("Plotting counts")
      fs <- r$getFeatureSet("CountSet")
      if (!is.null(fs)) {
        projection.ReconstructedCountSet(fs,
                                         phi0=r$phi0,
                                         ids=ids,
                                         transform=transform,
                                         axisdir=axisdir,
                                         projection=projection,
                                         proj.centre=proj.centre,
                                         ...)
      }
    }
    
    ## KDE
    if (datapoint.contours) {
      message("Plotting point contours")
      fs <- r$getFeatureSet("PointSet")
      if (!is.null(fs)) {
        k <- fs$getKDE()
        for (id in ids) {
          if (!is.null(k[[id]])) {
            css <- k[[id]]$contours
            for(cs in css) {
              suppressWarnings(lines(projection(rotate.axis(transform(cs,
                                                                      phi0=r$phi0),
                                                            axisdir*pi/180),
                                                lambdalim=lambdalim*pi/180,
                                                lines=TRUE,
                                                proj.centre=pi/180*proj.centre),
                                     col=fs$cols[[id]]))
            }
            ## FIXME: contours need to be labelled
          }
        }

        ## Plot locations of highest contours
        for (id in ids) {
          if (!is.null(k[[id]])) {
            suppressWarnings(points(projection(rotate.axis(transform(k[[id]]$maxs,
                                                                     phi0=r$phi0),
                                                           axisdir*pi/180),
                                               proj.centre=pi/180*proj.centre),
                                    pch=22, cex=1, lwd=1,
                                    col="black", bg=fs$cols[[id]]))
          }
        }
      }
    }

    ## KR
    if (grouped.contours) {
      message("Plotting count contours")
      fs <- r$getFeatureSet("CountSet")
      if (!is.null(fs)) {
        k <- fs$getKR()
        for (id in ids) {
          if (!is.null(k[[id]])) {
            css <- k[[id]]$contours
            for(cs in css) {
              lines(projection(rotate.axis(transform(cs,
                                                     phi0=r$phi0),
                                           axisdir*pi/180),
                               lambdalim=lambdalim*pi/180,
                               lines=TRUE,
                               proj.centre=pi/180*proj.centre),
                    col=fs$cols[[id]])
            }
            ## FIXME: contours need to be labelled
          }
        }
        ## Plot locations of highest contours
        for (id in ids) {
          if (!is.null(k[[id]])) {
            points(projection(rotate.axis(transform(k[[id]]$maxs,
                                                    phi0=r$phi0),
                                          axisdir*pi/180),
                              proj.centre=pi/180*proj.centre),
                   pch=23, cex=1, lwd=1,
                   col="black", bg=fs$cols[[id]])
          }
        }
      }
    }

    ## Landmarks
    if (landmarks) {
      message("Plotting landmarks")
      fs <- r$getFeatureSet("LandmarkSet")
      if (!is.null(fs)) {
        projection.ReconstructedLandmarkSet(fs,
                                            phi0=r$phi0,
                                            ids=ids,
                                            transform=transform,
                                            axisdir=axisdir,
                                            projection=projection,
                                            proj.centre=proj.centre,
                                            ...)
      }
    }

    
    NextMethod(projection=projection,
               philim=philim,
               labels=labels,
               colatitude=TRUE,
               grid=grid,
               add=TRUE,
               image=FALSE,
               mesh=mesh)

  }

##' @method projection RetinalReconstructedOutline
sphericalplot.RetinalReconstructedOutline <- function(r,
                                                      datapoints=TRUE,
                                                      ids=r$getIDs(), ...) {
  NextMethod()

  if (datapoints) {
    message("Plotting points")
      sphericalplot.ReconstructedPointSet(r,
                                          projection=projection,
                                          ids=ids, ...)
  }
}
