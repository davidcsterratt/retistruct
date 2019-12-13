##' Class containing functions and data to map \link{PointSet}s to
##' \link{ReconstructedOutline}s
##'
##' @description A ReconstructedPointSet contains information about
##'   features located on \code{\link{ReconstructedOutline}}s. Each
##'   ReconstructedPointSet contains a list of matrices, each of
##'   which has columns labelled \code{phi} (latitude) and
##'   \code{lambda} (longitude) describing the spherical coordinates
##'   of points on the ReconstructedOutline.
##' 
##' @author David Sterratt
##' @importFrom geometry delaunayn
##' @export
ReconstructedPointSet <- R6Class("ReconstructedPointSet",
  inherit = ReconstructedFeatureSet,
  public = list(
    ##' @field KDE Kernel density estimate, computed using
    ##'   \code{\link{compute.kernel.estimate}} in \code{getKDE}
    KDE = NULL,
    ##' @description Get Karcher mean of datapoints in spherical coordinates
    ##' @return Karcher mean of datapoints in spherical coordinates
    getMean = function() {
      Dss.mean <- list()
      for (id in self$getIDs()) {
        km <- karcher.mean.sphere(self$getFeature(id), na.rm=TRUE)
        Dss.mean[[id]] <- cbind(phi=km["phi"], lambda=km["lambda"])
      }
      return(Dss.mean)
    },
    ##' @description Get area of convex hull around data points on sphere
    ##' @return Area in degrees squared
    getHullarea = function() {
      Dss.hullarea <- list()
      for (id in self$getIDs()) {
        if (nrow(self$getFeature(id)) >= 3) {
          Dsp <- sphere.spherical.to.polar.cart(self$getFeature(id), pa=TRUE)
          Dspt <- suppressMessages(delaunayn(Dsp))
          Dss.hullarea[[id]] <- sum(sphere.tri.area(self$getFeature(id), Dspt))*(180/pi)^2
        } else {
          Dss.hullarea[[id]] <- NA
        }
      }
      return(Dss.hullarea)
    },
    ##' @description Get kernel density estimate of data points
    ##' @return See \code{\link{compute.kernel.estimate}}
    getKDE = function() {
      if (is.null(self$KDE)) {
        Dss <- list()
        for (id in self$getIDs()) {
          Dss[[id]] <- self$getFeature(id)
        }
        self$KDE <- compute.kernel.estimate(Dss, self$ro$phi0, kde.fhat, kde.compute.concentration)
      }
      return(self$KDE)
    }
  )
)

##' Compute a kernel estimate over a grid and do a contour analysis
##' of this estimate. The contour heights the determined by finding
##' heights that exclude a certain fraction of the probability. For
##' example, the 95% contour is excludes 95% of the probability mass,
##' and it should enclose about 5% of the points. The contour levels
##' are specified by  the \code{contour.levels} option; by default
##' they are \code{c(5, 25, 50, 75, 95)}.
##' 
##' @title Kernel estimate over grid
##' @param Dss List of datasets. The first two columns of each datasets
##' are coordinates of points on the sphere in spherical polar
##' (latitude, \code{phi}, and longitude, \code{lambda})
##' coordinates. In the case kernel smoothing, there is a third column
##' of values of dependent variables at those points.
##' @param phi0 Rim angle in radians
##' @param fhat Function such as \code{\link{kde.fhat}} or
##' \code{\link{kr.yhat}} to compute the density given data and a
##' value of the concentration parameter \code{kappa} of the Fisher
##' density.
##' @param compute.conc Function to return the optimal value of the
##' concentration parameter kappa given the data.
##' @return A list containing
##' \item{\code{kappa}}{The concentration parameter}
##' \item{\code{h}}{A pseudo-bandwidth parameter, the inverse of the square root of \code{kappa}. Units of degrees.}
##' \item{\code{flevels}}{Contour levels.}
##' \item{\code{labels}}{Labels of the contours.}
##' \item{\code{g}}{Raw density estimate drawn on non-area-preserving projection. Comprises locations of gridlines in Cartesian coordinates (\code{xs} and \code{ys}), density estimates at these points, \code{f} and location of maximum in Cartesian coordinates (\code{max}).}
##' \item{\code{gpa}}{Raw density estimate drawn on area-preserving projection. Comprises same elements as above.}
##' \item{\code{contour.areas}}{Area of each individual contour. One level may have more than one contour; this shows the areas of all such contours.}
##' \item{\code{tot.contour.areas}}{Data frame containing the total area within the contours at each level.}
##' @author David Sterratt
##' @export
compute.kernel.estimate <- function(Dss, phi0, fhat, compute.conc) {
  vols <- getOption("contour.levels")
  res <- 100

  ## Helper function to get kde as locations gs in spherical coordinates
  get.kde <- function(gs, mu, kappa, res, crop=TRUE) {
    ## Make space for the kernel density estimates
    gk <- fhat(gs, mu, kappa)

    ## If we're cropping, were going to eliminate pixels whos centres
    ## lie outwith the outline altogether. If crop isn't set, there's
    ## a wee margin, that ought to allow all contours to lie outwith
    ## the outline.
    if (crop) {
      gk[gs[,"phi"] > phi0] <- NA
    }
    ## else {
    ##  gk[gs[,"phi"] > phi0 + 2/res*(phi0+pi/2)] <- 0
    ## }

    ## Put the estimates back into a matrix. The matrix is filled up
    ## column-wise, so the matrix elements should match the elements of
    ## gxs and gys
    k <- matrix(gk, res, res)
    k[is.na(k)] <- 0
    return(k)
  }
  
  ## Get data points
  KDE <- list()
  if (length(Dss) > 0) {
    ## First create a grid in Cartesian coordinates with
    ## area-preserving coords. The extra margin on the grid is needed
    ## so that the contour lines are continuous round the edge.
    gpa <- create.polar.cart.grid(TRUE, res, min(phi0 + 0.2*(phi0+pi/2), pi/2))
    ## And one without area-preserving coords
    g   <- create.polar.cart.grid(FALSE, res, min(phi0 + 0.2*(phi0+pi/2), pi/2))

    ## Check conversion
    ## gcb <- sphere.spherical.to.polar.cart(gs, pa)
    ## points(rho.to.degrees(gcb, phi0, pa), pch='.')
    
    for (i in names(Dss)) {
      if (nrow(Dss[[i]]) > 2) {
        ## Find the optimal concentration of the kernel density
        ## estimator
        kappa <- compute.conc(Dss[[i]])
        
        ## Now we've found the concentration, let's try to estimate
        ## and display the density at grid points on an azimuthal
        ## equidistant projection (f) and on an aziumuthal equal-area
        ## projection (fpa).
        fpa <- get.kde(gpa$s, Dss[[i]], kappa, res)
        f  <-  get.kde(g$s,   Dss[[i]], kappa, res)
        maxs   <- gpa$s[which.max(fpa),,drop=FALSE]

        ## The above estimates are set to NA outwith the outline. For
        ## the purposes of computing smooth contours, this causes
        ## jagged edges, so we also get uncropped versions there the
        ## density spreads outwith the outline.
        fpau <- get.kde(gpa$s, Dss[[i]], kappa, res, crop=FALSE)
        fu  <-  get.kde(g$s,   Dss[[i]], kappa, res, crop=FALSE)

        ## Determine the value of gk that encloses 0.95 of the
        ## density.  To compute the density, we need to know the
        ## area of each little square, which is why we have used the
        ## are-preserving projection. FIXME: I think this method of
        ## "integration" could be improved.
        vol.contours <- TRUE
        if (vol.contours) {
          f.sort <- sort(as.vector(fpa))
          js <- findInterval(vols/100, cumsum(f.sort)/sum(f.sort))
          flevels <- f.sort[js]
          if (length(unique(flevels)) < length(flevels)) {
            warning("The volume contours method has found duplicated levels - there is probably something wrong with the data.")
            KDE[[i]] <- NULL
            next
          }
        } else {
          flevels <- vols/100*max(fpa)
        }

        ## Store full kde matrices
        KDE[[i]] <- list(kappa=kappa,
                         h=180/pi/sqrt(kappa),
                         flevels=flevels,
                         maxs=maxs,
                         g=  list(xs=g$xs,   ys=g$ys,   f=f  , fu=fu),
                         gpa=list(xs=gpa$xs, ys=gpa$ys, f=fpa, fu=fpau))

        ## Get contours in Cartesian space
        cc <- grDevices::contourLines(gpa$xs, gpa$ys, fpau, levels=flevels)
        cs <- list()
        ## Must be careful, as there is a core function called labels
        labels <- rep(NA, length(cc))
        contour.areas <- rep(NA, length(cc))
        if (length(cc) > 0) {
          for (j in 1:length(cc)) {
            cs[[j]] <- list()
            contour.areas[j] <- geometry::polyarea(cc[[j]]$x, cc[[j]]$y) * (180/pi)^2
            ccj <- cbind(x=cc[[j]]$x, y=cc[[j]]$y)
            ## Convert back to Spherical coordinates
            cs[[j]] <- polar.cart.to.sphere.spherical(ccj, TRUE)
            ## Push any points outwith the outline back into it
            cs[[j]][cs[[j]][,"phi"] > phi0, "phi"] <- phi0
            labels[j] <- vols[which(flevels==cc[[j]]$level)]
          }
        }
        KDE[[i]]$contours <- cs
        KDE[[i]]$labels <- labels
        names(contour.areas) <- c()
        KDE[[i]]$contour.areas <- contour.areas
        KDE[[i]]$tot.contour.areas <- stats::aggregate(contour.areas ~ labels,
                                                       data.frame(labels, contour.areas), sum)
      }
    }
  }
  return(KDE)
}


projection.ReconstructedPointSet <-
  function(r,
           phi0,
           transform=identity.transform,
           ids=r$getIDs(),
           axisdir=cbind(phi=90, lambda=0),
           projection=azimuthal.equalarea,
           proj.centre=cbind(phi=0, lambda=0),
           markup=NULL, # Not used in this function; hides from ...
           max.proj.dim=NULL, # Not used in this function; hides from ...
           ...)
{
  ## This will call projection.reconstructedOutline(), but without
  ## drawing a grid.  The grid will be drawn later, after all the data
  ## as appeared..
  ## Datapoints
  for (id in ids) {
    if (!is.null(r$getFeature(id))) {
      points(projection(rotate.axis(transform(r$getFeature(id), phi0=phi0),
                                    axisdir*pi/180),
                        proj.centre=pi/180*proj.centre),
             col=r$cols[[id]],
             pch=20, ...)
    }
  }
}

##' @method sphericalplot ReconstructedPointSet
sphericalplot.ReconstructedPointSet <- function(r,
                                                datapoints=TRUE,
                                                ids=r$getIDs(), ...) {

  fs <- r$getFeatureSet("PointSet")
  if (!is.null(fs)) {
    size <- 1/10
    for (id in ids) {
      if (!is.null(fs$getFeature(id))) {
        Dc <- sphere.spherical.to.sphere.cart(fs$getFeature(id))
        if (nrow(Dc) > 0) {
          
          ## Find axis in z=0 plane that is orthogonal to projection of
          ## datapoint onto that plane
          ax1 <- 1/sqrt(apply(Dc[,1:2,drop=FALSE]^2, 1, sum)) * cbind(-Dc[,2], Dc[,1], 0)

          ax1 <- as.matrix(ax1, ncol=3)

          ## Find axis that is orthogonal to the plane of axis 1 and the
          ## datapoint
          ax2 <- extprod3d(Dc, ax1)
          ax2 <- matrix(ax2, ncol=3)

          ax2 <- ax2/sqrt(apply(ax2^2, 1, sum))

          ## Create the vertices of an equillateral triangle to plot
          v1 <- Dc + size *  ax1/2
          v2 <- Dc + size * (-ax1/4 + sqrt(3)/4*ax2)
          v3 <- Dc + size * (-ax1/4 - sqrt(3)/4*ax2)

          ## Plot the triangle inside and outside the sphere
          inmag <- 0.98
          outmag <- 1.02
          
          x <- rbind(v2[,1], v1[,1], v3[,1])
          y <- rbind(v2[,2], v1[,2], v3[,2])
          z <- rbind(v2[,3], v1[,3], v3[,3])
          triangles3d(inmag*x, inmag*y, inmag*z, color=fs$cols[[id]])

          x <- rbind(v1[,1], v2[,1], v3[,1])
          y <- rbind(v1[,2], v2[,2], v3[,2])
          z <- rbind(v1[,3], v2[,3], v3[,3])
          triangles3d(outmag*x, outmag*y, outmag*z, color=fs$cols[[id]],
                      pch=20)
        }
      }
    }
  }
}
  
