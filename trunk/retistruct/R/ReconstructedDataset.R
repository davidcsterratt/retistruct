##' This function infers the coordinates of datapoints \code{Ds }and
##' landmarks \code{Ss} in  spherical coordinates.
##'
##' @title Constructor for RecontructedDataset object
##' @param r Object that of clases \code{reconstructedOutline} and \code{dataset}.
##' @param report Function used to report progress.
##' @return \code{reconstructedDataset} object containing the input
##' information and the following modified and extra information:
##' \item{\code{Dsb}}{Datapoints in barycentric coordinates}
##' \item{\code{Dsc}}{Datapoints on reconstructed sphere in cartesian coordinates}
##' \item{\code{Dss}}{Datapoints on reconstructed sphere in spherical coordinates}
##' \item{\code{Ssb}}{Landmarks in barycentric coordinates}
##' \item{\code{Ssc}}{Landmarks on reconstructed sphere in cartesian coordinates}
##' \item{\code{Sss}}{Landmarks on reconstructed sphere in spherical coordinates}
##' @author David Sterratt
##' @import geometry
##' @export
ReconstructedDataset <- function(r, report=message) {
  report("Inferring coordinates of datapoints")
  Dsb <- list() # Datapoints in barycentric coordinates
  Dsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Dss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ds) & (length(r$Ds) > 0)) {
    for (name in names(r$Ds)) {
      Dsb[[name]] <- tsearchn(r$P, r$T, r$Ds[[name]])
      oo <- is.na(Dsb[[name]]$idx)     # Points outwith outline
      if (any(oo)) {
        warning(paste(sum(oo), name, "datapoints outwith the outline will be ignored."))
      }
      Dsb[[name]]$p   <- Dsb[[name]]$p[!oo,,drop=FALSE]
      Dsb[[name]]$idx <- Dsb[[name]]$idx[!oo]
      Dsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Dsb[[name]])
      Dss[[name]] <- sphere.cart.to.sphere.spherical(Dsc[[name]], r$R)
    }
  }
  
  report("Inferring coordinates of grouped datapoints")
  Gsb <- list() # Datapoints in barycentric coordinates
  Gsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Gss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Gs) & (length(r$Gs) > 0)) {
    for (name in names(r$Gs)) {
      Gsb[[name]] <- tsearchn(r$P, r$T, r$Gs[[name]][,1:2,drop=FALSE])
      oo <- is.na(Gsb[[name]]$idx)     # Points outwith outline
      if (any(oo)) {
        warning(paste(sum(oo), name, "datapoints outwith the outline will be ignored."))
      }
      Gsb[[name]]$p   <- Gsb[[name]]$p[!oo,,drop=FALSE]
      Gsb[[name]]$idx <- Gsb[[name]]$idx[!oo]
      Gsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Gsb[[name]])
      Gss[[name]] <- cbind(sphere.cart.to.sphere.spherical(Gsc[[name]], r$R),
                           r$Gs[[name]][!oo,3])
    }
  }
  
  report("Inferring coordinates of landmarks")
  Ssb <- list() # Landmarks in barycentric coordinates
  Ssc <- list() # Landmarks on reconstructed sphere in cartesian coordinates
  Sss <- list() # Landmarks on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ss) & (length(r$Ss) > 0)) {
    for (i in 1:length(r$Ss)) {
      Ssb[[i]] <- with(r, tsearchn(P, T, r$Ss[[i]]))
      Ssc[[i]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Ssb[[i]])
      Sss[[i]] <- sphere.cart.to.sphere.spherical(Ssc[[i]], r$R)
    }
    names(Ssb) <- names(r$Ss)
    names(Ssc) <- names(r$Ss)
    names(Sss) <- names(r$Ss)
  }

  d <- merge(list(Dsb=Dsb, Dsc=Dsc, Dss=Dss,
                  Gsb=Gsb, Gsc=Gsc, Gss=Gss,
                  Ssb=Ssb, Ssc=Ssc, Sss=Sss), r)
  ## Trigger recomputation of Kernel density estimates
  d$KDE <- NULL
  d$KR <- NULL

  class(d) <- addClass("reconstructedDataset", r)
  return(d)
}

##' Get spherical coordinates of datapoints.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{reonstructedDataset} object.
##' @return \code{Dss}
##' @method getDss reconstructedDataset
##' @author David Sterratt
##' @export
getDss.reconstructedDataset <- function(r) {
  return(r$Dss)
}

##' @title Get grouped variable with locations in spherical coordinates.
##' @param r \code{reonstructedDataset} object.
##' @return \code{Gss}
##' @method getGss reconstructedDataset
##' @author David Sterratt
##' @export
getGss.reconstructedDataset <- function(r) {
  return(r$Gss)
}

##' Get Karcher mean of datapoints in spherical coordinates.
##'
##' @title Karcher mean of datapoints in spherical coordinates
##' @param r \code{\link{reconstructedDataset}} or \code{\link{retinalReconstructedDataset}} object.
##' @return \code{Dss.mean}
##' @method getDss.mean reconstructedDataset
##' @author David Sterratt
##' @export
getDss.mean.reconstructedDataset <- function(r) {
  Dss.mean <- list()
  if (length(r$Dss)) {
    for (i in 1:length(r$Dss)) {
      km <- karcher.mean.sphere(r$Dss[[i]], na.rm=TRUE)
      Dss.mean[[i]] <- cbind(phi=km["phi"], lambda=km["lambda"])
    }
  }
  names(Dss.mean) <- names(r$Dss)
  return(Dss.mean)
}

##' Get spherical coordinates of landmarks.
##'
##' @title Get transformed spherical coordinates of landmarks.
##' @param r \code{reonstructedDataset} object.
##' @return \code{Sss}
##' @method getSss reconstructedDataset
##' @author David Sterratt
##' @export
getSss.reconstructedDataset <- function(r) {
  return(r$Sss)
}

##' @title Get kernel density estimate of data points
##' @param r \code{reconstructedDataset} object
##' @return See \code{\link{compute.kernel.estimate}}
##' @author David Sterratt
##' @export
getKDE <- function(r) {
  if (is.null(r$KDE)) {
    return(compute.kernel.estimate(getDss(r), r$phi0, kde.fhat, kde.compute.concentration))
  }
  return(r$KDE)
}

##' Compute a kernel estimate over a grid and do a contour analsysis
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
##' (lattitude, \code{phi}, and longitude, \code{lambda})
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
##' \item{\code{g}}{Raw density estimate drawn on non-area-preserving projection. Comprises locations of gridlines in Cartesian coordinates (\code{xs} and \code{ys}) and density estimates at these points, \code{f}.}
##' \item{\code{gpa}}{Raw density estimate drawn on area-preserving projection. Comprises same elements as above.}
##' \item{\code{contour.areas}}{Area of each individual contour. One level may ahave more than one contour; this shows the areas of all such contours.}
##' \item{\code{tot.contour.areas}}{Data frame containing the total area within the contours at each level.}
##' @author David Sterratt
##' @export
compute.kernel.estimate <- function(Dss, phi0, fhat, compute.conc) {
  vols <- getOption("contour.levels")
  res <- 100

  ## Helper function to get kde as locations gs in spherical coordinates
  get.kde <- function(gs, mu, kappa, res) {
    ## Make space for the kernel density estimates
    gk <- fhat(gs, mu, kappa)

    gk[gs[,"phi"] > phi0] <- NA
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
    ## area-preserving coords
    gpa <- create.polar.cart.grid(TRUE, res, phi0)
    ## And one without area-preserving coords
    g   <- create.polar.cart.grid(FALSE, res, phi0)

    ## Check conversion
    ## gcb <- sphere.spherical.to.polar.cart(gs, pa)
    ## points(rho.to.degrees(gcb, phi0, pa), pch='.')
    
    for (i in names(Dss)) {
      if (nrow(Dss[[i]]) > 2) {
        ## Find the optimal concentration of the kernel density
        ## estimator
        kappa <- compute.conc(Dss[[i]])
        
        ## Now we've found the concentration, let's try to estimate
        ## and display the density over our polar representation of
        ## the data points
        fpa <- get.kde(gpa$s, Dss[[i]], kappa, res)
        f  <-  get.kde(g$s,   Dss[[i]], kappa, res)
        
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
        } else {
          flevels <- vols/100*max(fpa)
        }

        ## Store full kde matrices
        KDE[[i]] <- list(kappa=kappa,
                         h=180/pi/sqrt(kappa),
                         flevels=flevels,
                         labels=vols,
                         g=  list(xs=g$xs,   ys=g$ys,   f=f),
                         gpa=list(xs=gpa$xs, ys=gpa$ys, f=fpa))

        ## Get contours in Cartesian space
        cc <- contourLines(gpa$xs, gpa$ys, fpa, levels=flevels)
        cs <- list()
        ## Must be careful, as there is a core function called labels
        labels <- rep(NA, length(cc))
        contour.areas <- rep(NA, length(cc))
        if (length(cc) > 0) {
          for (j in 1:length(cc)) {
            cs[[j]] <- list()
            contour.areas[j] <- polyarea(cc[[j]]$x, cc[[j]]$y) * (180/pi)^2
            ccj <- cbind(x=cc[[j]]$x, y=cc[[j]]$y)
            ## Convert back to Spherical coordinates
            cs[[j]] <- polar.cart.to.sphere.spherical(ccj, TRUE)
            labels[j] <- vols[which(flevels==cc[[j]]$level)]
          }
        }
        KDE[[i]]$contours <- cs
        KDE[[i]]$labels <- labels
        names(contour.areas) <- labels
        KDE[[i]]$contour.areas <- contour.areas
        KDE[[i]]$tot.contour.areas <- aggregate(contour.areas ~ labels,
                                                data.frame(labels, contour.areas), sum)
      }
    }
  }
  return(KDE)
}

##' @title Get kernel regression estimate of grouped data points
##' @param r \code{reconstructedDataset} object
##' @return See \code{\link{compute.kernel.estimate}}
##' @author David Sterratt
##' @export
getKR <- function(r) {
  if (is.null(r$KR)) {
    yhat <- function(r, mu, kappa) {
      kr.yhat(r, mu[,1:2], mu[,3], kappa)
    }
    compute.conc <- function(mu) {
      kr.compute.concentration(mu[,1:2], mu[,3])
    }
    Gss <- getGss(r)
    ## Get rid of datasets with very little data
    for (n in names(Gss)) {
      if (sum(Gss[[n]][,3]) <= 2) {
        Gss[[n]] <- NULL
      }
    }
    return(compute.kernel.estimate(Gss, r$phi0, yhat, compute.conc))
  }
  return(r$KR)
}

##' Plot datapoints in polar plot
##'
##' @title Polar plot of reconstructed dataset
##' @param r \code{reconstructedDataset} object
##' @param show.grid Whether or not to show the grid lines of lattitude and longitude
##' @param grid.col Colour of the minor grid lines
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param flip.horiz Wether to flip about a horizontal axis
##' @param labels Vector of 4 labels to plot at 0, 90, 180 and 270 degrees 
##' @param ... Other graphics parameters.  The option
##' \code{preserve.area} creates an area-preserving plot (default
##' \code{FALSE}). The option \code{datapoints} causes datapoints to
##' be plotted (default \code{TRUE}).  The option
##' \code{datapoint.means} causes datapoint means to be plotted
##' (default \code{TRUE}).  The option \code{landmakrs} causes
##' landmarks to be plotted (default \code{TRUE}). 
##' @method plot.polar reconstructedDataset
##' @author David Sterratt
##' @export
plot.polar.reconstructedDataset <- function(r, show.grid=TRUE,
                                            grid.col="gray",
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45,
                                            flip.horiz=FALSE,
                                            labels=c(0, 90, 180, 270), ...) {
  NextMethod()

  args <- list(...)
  plot.datapoints <- is.null(args$datapoints) || args$datapoints
  plot.datapoint.means <- is.null(args$datapoint.means) || args$datapoint.means
  plot.datapoint.contours <- is.null(args$datapoint.contours) || args$datapoint.contours
  plot.grouped.contours <- is.null(args$grouped.contours) || args$grouped.contours
  plot.landmarks <- is.null(args$landmarks) || args$landmarks
  plot.preserve.area <- !is.null(args$preserve.area) && args$preserve.area
  plot.mosaic <- !is.null(args$mosaic) && args$mosaic
  plot.kde <- !is.null(args$kde) && args$kde
  pa <- plot.preserve.area

  if (plot.kde) {
    KDE <- getKDE(r)
    image(rho.to.degrees(KDE$red$g$xs, r$phi0, pa),
          rho.to.degrees(KDE$red$g$ys, r$phi0, pa),
          KDE$red$g$f, col=gray((0:100)/100))
    Dss <- getDss(r)
    pos <- sphere.spherical.to.polar.cart(Dss[["red"]], pa)
        suppressWarnings(points(rho.to.degrees(pos, r$phi0, pa),
                                col="red",
                                pch=20, cex=0.2, ...))

  }
  
  ## Datapoints
  if (plot.datapoints) {
    Dss <- getDss(r)
    if (length(Dss)) {
      for (i in 1:length(Dss)) {
        pos <- sphere.spherical.to.polar.cart(Dss[[i]], pa)
        suppressWarnings(points(rho.to.degrees(pos, r$phi0, pa),
                                col=r$cols[[names(Dss)[i]]],
                                pch=20, ...))
      }
    }
  }

  ## Mean datapoints
  if (plot.datapoint.means) {
    Dss.mean <- getDss.mean(r)
    if (length(Dss.mean)) {
      for (i in 1:length(Dss.mean)) {
        pos <- sphere.spherical.to.polar.cart(Dss.mean[[i]], pa)
        suppressWarnings(points(rho.to.degrees(pos, r$phi0, pa),
                                bg=r$cols[[names(Dss.mean)[i]]], col="black",
                                pch=23, cex=1.5, ...))
      }
    }
  }
  
  ## KDE
  if (plot.datapoint.contours) {
    k <- getKDE(r)
    Dss <- getDss(r)
    if (length(k)) {
      for (i in 1:length(k)) {
        if (pa) {
          g <- k[[i]]$gpa
        } else {
          g <- k[[i]]$g
        }
        ## Plot contours
        contour(rho.to.degrees(g$xs, r$phi0, pa),
                rho.to.degrees(g$ys, r$phi0, pa),
                g$f, add=TRUE, levels=k[[i]]$flevels,
                col=r$cols[[names(Dss)[i]]],
                ## drawlabels=FALSE,
                labels=k[[i]]$labels)
      }
    }
  }

  ## KDE
  if (plot.grouped.contours) {
    k <- getKR(r)
    Gss <- getGss(r)
    if (length(k)) {
      for (i in 1:length(k)) {
        if (pa) {
          g <- k[[i]]$gpa
        } else {
          g <- k[[i]]$g
        }
        ## Plot contours
        contour(rho.to.degrees(g$xs, r$phi0, pa),
                rho.to.degrees(g$ys, r$phi0, pa),
                g$f, add=TRUE, levels=k[[i]]$flevels,
                col=r$cols[[names(Gss)[i]]],
                ## drawlabels=FALSE,
                labels=k[[i]]$labels)
      }
    }
  }
  
  ## Voroni
  if (plot.mosaic) {
    Dss <- getDss(r)
    if (length(Dss)) {
      for (i in 1:length(Dss)) {
        if (nrow(Dss[[i]]) >= 2) {
          ## Convert to angle-preserving coordinates
          pos <- sphere.spherical.to.polar.cart(Dss[[i]], preserve="angle")
          ## Create Voronoi mosaic in area-preserving coordinates
          vm <- voronoi.mosaic(pos[,1], pos[,2])

          ## Convert back to polar coords
          vmxy.polar <- polar.cart.to.sphere.spherical(cbind(x=vm$x, y=vm$y),
                                                       preserve="angle")
          ## Convert into whichever representation we're using
          ## and put back on to voroini object for plotting
          pos <- rho.to.degrees(sphere.spherical.to.polar.cart(vmxy.polar, pa), r$phi, pa)
          vm$x <- pos[,1]
          vm$y <- pos[,2]
          ## Do the same for the dummy coordinates
          vmxy.polar <- polar.cart.to.sphere.spherical(cbind(x=vm$dummy.x,
                                                             y=vm$dummy.y),
                                                       preserve="angle")
          ## Convert into whichever representation we're using
          ## and put back on to voroini object for plotting
          pos <- rho.to.degrees(sphere.spherical.to.polar.cart(vmxy.polar, pa), r$phi, pa)
          vm$dummy.x <- pos[,1]
          vm$dummy.y <- pos[,2]

          plot.voronoi.circular(vm, R=r$phi0*180/pi + 90, col=r$cols[[names(Dss)[i]]])
        }
      }
    }
  }
  
  ## Landmarks
  if (plot.landmarks) {
    Sss <- getSss(r)
    if (length(Sss)) {
      for (i in 1:length(Sss)) {
        name <- names(Sss)[i]
        col <- ifelse(is.null(name) || (name==""), "default", name)
        pos <- sphere.spherical.to.polar.cart(Sss[[i]], pa)
        suppressWarnings(lines(rho.to.degrees(pos, r$phi, pa),
                               col=r$cols[[col]], ...))
      }
    }
  }
}

##' Draw a spherical plot of datapoints.
##'
##' @title Spherical plot of reconstructed outline
##' @param r \code{reconstructedOutline} object
##' @param ... Other graphics parameters -- not used at present
##' @method plot.spherical reconstructedDataset
##' @author David Sterratt
##' @export
plot.spherical.reconstructedDataset <- function(r, ...) {
  NextMethod()

  args <- list(...)
  plot.datapoints <- is.null(args$datapoints) || args$datapoints
  size <- r$R/10
  
  if (plot.datapoints) {
    Dsc <- r$Dsc

    for (i in 1:length(Dsc)) {
      Dc <- Dsc[[i]]
      
      ## Find axis in z=0 plane that is orthogonal to projection of
      ## datapoint onto that plane
      ax1 <- 1/sqrt(apply(Dc[,1:2,drop=FALSE]^2, 1, sum)) * cbind(-Dc[,2], Dc[,1], 0)

      ax1 <- as.matrix(ax1, ncol=3)
      print("ax1")
      print(ax1)
      ## Find axis that is orthogonal to the plane of axis 1 and the
      ## datapoint
      ax2 <- extprod3d(Dc, ax1)
      ax2 <- matrix(ax2, ncol=3)
      print("ax2")
      print(ax2)

      ax2 <- ax2/sqrt(apply(ax2^2, 1, sum))

      ## Create the verticies of an equillateral triangle to plot
      v1 <- Dc + size *  ax1/2
      v2 <- Dc + size * (-ax1/4 + sqrt(3)/4*ax2)
      v3 <- Dc + size * (-ax1/4 - sqrt(3)/4*ax2)

      ## Plot the triangle inside and outside the sphere
      inmag <- 0.99
      outmag <- 1.02
      
      x <- rbind(v2[,1], v1[,1], v3[,1])
      y <- rbind(v2[,2], v1[,2], v3[,2])
      z <- rbind(v2[,3], v1[,3], v3[,3])
      triangles3d(inmag*x, inmag*y, inmag*z, color=r$cols[[names(Dsc)[i]]])

      x <- rbind(v1[,1], v2[,1], v3[,1])
      y <- rbind(v1[,2], v2[,2], v3[,2])
      z <- rbind(v1[,3], v2[,3], v3[,3])
      triangles3d(outmag*x, outmag*y, outmag*z, color=r$cols[[names(Dsc)[i]]],
                  pch=20, ...)
    }
  }
}  

