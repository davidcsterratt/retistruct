##' This function infers the coordinates of datapoints \code{Ds }and
##' landmarks \code{Ss} in  spherical coordinates.
##'
##' @title Constructor for RecontructedDataset object
##' @param r Object that of clases \code{reconstructedOutline} and \code{dataset}.
##' @param report Function used to report progress.
##' @return \code{reconstructedDataset} object containing the input
##' information and the following modified and extra information:
##' \item{\code{Dsb}} Datapoints in barycentric coordinates
##' \item{\code{Dsc}} Datapoints on reconstructed sphere in cartesian coordinates
##' \item{\code{Dss}} Datapoints on reconstructed sphere in spherical coordinates
##' \item{\code{Ssb}} Landmarks in barycentric coordinates
##' \item{\code{Ssc}} Landmarks on reconstructed sphere in cartesian coordinates
##' \item{\code{Sss}} Landmarks on reconstructed sphere in spherical coordinates
##' @author David Sterratt
ReconstructedDataset <- function(r, report=message) {
  report("Inferring coordinates of datapoints")
  Dsb <- list() # Datapoints in barycentric coordinates
  Dsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Dss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ds) & (length(r$Ds) > 0)) {
    for (name in names(r$Ds)) {
      Dsb[[name]] <- tsearchn(r$P, r$T, r$Ds[[name]])
      Dsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Dsb[[name]])
      Dss[[name]] <- sphere.cart.to.sphere.spherical(Dsc[[name]], r$R)
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
                  Ssb=Ssb, Ssc=Ssc, Sss=Sss), r)
  class(d) <- unique(c("reconstructedDataset", class(r)))
  return(d)
}

##' Get spherical coordinates of datapoints.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{reonstructedDataset} object.
##' @return \code{Dss}
##' @method getDss reconstructedDataset
##' @author David Sterratt
getDss.reconstructedDataset <- function(r) {
  return(r$Dss)
}

##' Get Karcher mean of datapoints in spherical coordinates.
##'
##' @title Karcher mean of datapoints in spherical coordinates
##' @param r \code{reonstructedDataset} object.
##' @return \code{Dss.mean}
##' @method getDss.mean reconstructedDataset
##' @author David Sterratt
getDss.mean.reconstructedDataset <- function(r) {
  Dss.mean <- list()
  if (length(r$Dss)) {
    for (i in 1:length(r$Dss)) {
      Dss.mean[[i]] <- karcher.mean.sphere(r$Dss[[i]], na.rm=TRUE)
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
getSss.reconstructedDataset <- function(r) {
  return(r$Sss)
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
##' @method plot.spherical reconstructedOutline
##' @author David Sterratt
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
  plot.landmarks <- is.null(args$landmarks) || args$landmarks
  plot.preserve.area <- !is.null(args$preserve.area) && args$preserve.area
  pa <- plot.preserve.area
  
  ## Datapoints
  if (plot.datapoints) {
    Dss <- getDss(r)
    if (length(Dss)) {
      for (i in 1:length(Dss)) {
        phis    <- Dss[[i]][,"phi"]
        lambdas <- Dss[[i]][,"lambda"]
        xpos <- cos(lambdas)*phi.to.rho(phis, r$phi0, pa)
        ypos <- sin(lambdas)*phi.to.rho(phis, r$phi0, pa)
        suppressWarnings(points(xpos, ypos, col=r$cols[[names(Dss)[i]]],
                                pch=20, ...))
      }
    }
  }

  ## Mean datapoints
  if (plot.datapoint.means) {
    Dss.mean <- getDss.mean(r)
    if (length(Dss.mean)) {
      for (i in 1:length(Dss.mean)) {
        phis    <- Dss.mean[[i]]["phi"]
        lambdas <- Dss.mean[[i]]["lambda"]
        xpos <- cos(lambdas)*phi.to.rho(phis, r$phi0, pa)
        ypos <- sin(lambdas)*phi.to.rho(phis, r$phi0, pa)
        suppressWarnings(points(xpos, ypos,
                                bg=r$cols[[names(Dss.mean)[i]]], col="black",
                                pch=23, cex=1.5, ...))
      }
    }
  }

  ## Contours
  vols <- 0.5
  res <- 100
  if (plot.datapoint.contours) {
    Dss <- getDss(r)
    if (length(Dss)) {
      ## First create a grid in Cartesian coordinates
      lim <- sphere.spherical.to.polar.cart(cbind(phi=r$phi0, lambda=0), pa)[1,"x"]
      xs <- seq(-lim, lim, len=res)
      ys <- seq(-lim, lim, len=res)

      ## Create grid
      gxs <- outer(xs, ys*0, "+")
      gys <- outer(xs*0, ys, "+")

      ## gxs and gys are both 101 by 100 matrixes We now combine both
      ## matrices as a 101*100 by 2 matrix. The conversion as.vector() goes
      ## down the columns of the matrices gxs and gys
      gc <- cbind(x=as.vector(gxs), y=as.vector(gys))

      ## Now convert the cartesian coordinates to polar coordinates
      gs <- polar.cart.to.sphere.spherical(gc, pa)

      ## Check conversion
      ## gcb <- sphere.spherical.to.polar.cart(gs, pa)
      ## points(180/pi*gcb[,"x"], 180/pi*gcb[,"y"], pch='.')
      
      for (i in 1:length(Dss)) {
        mu <- cbind(phi=Dss[[i]][,"phi"], lambda=Dss[[i]][,"lambda"])
        if (nrow(mu) > 2) {
          ## Find the optimal bandwidth of the kernel density estimator
          sigma <- compute.bandwidth(mu, K)
          message(paste("sigma=", sigma))
          
          ## points(sphere.spherical.to.polar.cart(mu, pa)*180/pi)
          
          ## Now we've found sigma, let's try to estimate and display the
          ## density over our polar representation of the data points

          ## Make space for the kernel density estimates
          gk <- rep(0, nrow(gs))
          for (j in 1:nrow(gs)) {
            gk[j] <- K(gs[j,], mu, sigma)
          }

          gk[gs[,"phi"] > r$phi0] <- NA
          ## Put the estimates back into a matrix. The matrix is filled up
          ## column-wise, so the matrix elements should match the elements of
          ## gxs and gys
          k <- matrix(gk, res, res)
          k[is.na(k)] <- 0
          
          ## Determine the value of gk that encloses 0.95 of the
          ## density. FIXME: But of course to compute the density, we need to
          ## know the area of each little square...
          k.vec <- as.vector(k)
          js <- findInterval(1 - vols, cumsum(sort(k.vec))/sum(k.vec))
          klevels <- k.vec[js]
          klevels <- (1-vols)*max(k)
          message(paste("klevels=", klevels))
          

          ## Plot contours
          contour(180/pi*xs, 180/pi*ys, k, add=TRUE, levels=klevels,
                  col=r$cols[[names(Dss)[i]]])
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
        phi    <- Sss[[i]][,"phi"]
        lambda <- Sss[[i]][,"lambda"]
        x <- cos(lambda)*phi.to.rho(phis, r$phi0, pa)
        y <- sin(lambda)*phi.to.rho(phis, r$phi0, pa)
        suppressWarnings(lines(x, y, col=r$cols[[col]], ...))
      }
    }
  }
}

##' Draw a spherical plot of datapoints.
##'
##' @title Spherical plot of reconstructed outline
##' @param r \code{reconstructedOutline} object
##' @param ... Other graphics parameters -- not used at present
##' @method plot.spherical reconstructedOutline
##' @author David Sterratt
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

