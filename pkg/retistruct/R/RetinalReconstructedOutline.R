##' Create an object that is specific to retinal datasets. This
##' contains methods that return datapoint and landmark coordinates
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
        self$EOD <- 90 + ODmean["phi"] * 180/pi
      }
    }
  )
)

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


    ## Plot feature sets
    if (datapoints) {
      message("Plotting datapoints")
      fs <- r$getFeatureSet("PointSet")
      if (!is.null(fs)) {
        projection.ReconstructedPointSet(fs,
                                         projection=projection,
                                         phi0=r$phi0, ids=ids, ...)
      }
    }

    ## Mean datapoints
    if (datapoint.means) {
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
    
    if (landmarks) {
      message("Plotting landmarks")
      fs <- r$getFeatureSet("LandmarkSet")
      if (!is.null(fs)) {
        projection.ReconstructedLandmarkSet(fs,
                                            projection=projection,
                                            phi0=r$phi0, ids=ids, ...)
      }
    }

    if (grouped) {
      message("Plotting counts")
      fs <- r$getFeatureSet("CountSet")
      if (!is.null(fs)) {
        projection.ReconstructedCountSet(fs,
                                         projection=projection,
                                         phi0=r$phi0, ids=ids, ...)
      }
    }
    
    ## KDE
    if (datapoint.contours) {
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
    
    NextMethod(projection=projection,
               philim=philim,
               labels=labels,
               colatitude=TRUE,
               grid=grid,
               add=TRUE,
               image=FALSE,
               mesh=mesh)

  }
