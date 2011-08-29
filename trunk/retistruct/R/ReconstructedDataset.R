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
  plot.landmarks <- is.null(args$landmarks) || args$landmarks
  
  ## Datapoints
  if (plot.datapoints) {
    Dss <- getDss(r)
    for (i in 1:length(Dss)) {
      phis    <- Dss[[i]][,"phi"]
      lambdas <- Dss[[i]][,"lambda"]
      xpos <- cos(lambdas) * ((phis * 180/pi) + 90)
      ypos <- sin(lambdas) * ((phis * 180/pi) + 90)
      suppressWarnings(points(xpos, ypos, col=r$cols[[names(Dss)[i]]],
                              pch=20, ...))
    }
  }

  ## Mean datapoints
  if (plot.datapoint.means) {
    Dss.mean <- getDss.mean(r)
    for (i in 1:length(Dss.mean)) {
      phis    <- Dss.mean[[i]]["phi"]
      lambdas <- Dss.mean[[i]]["lambda"]
      xpos <- cos(lambdas) * ((phis * 180/pi) + 90)
      ypos <- sin(lambdas) * ((phis * 180/pi) + 90)
      suppressWarnings(points(xpos, ypos,
                              bg=r$cols[[names(Dss.mean)[i]]], col="black",
                              pch=23, cex=1.5, ...))
    }
  }
  
  ## Landmarks
  if (plot.landmarks) {
    Sss <- getSss(r)
    if (length(Sss) > 0) {
      for (i in 1:length(Sss)) {
        name <- names(Sss)[i]
        col <- ifelse(is.null(name) || (name==""), "default", name)
        phi    <- Sss[[i]][,"phi"]
        lambda <- Sss[[i]][,"lambda"]
        x <- cos(lambda) * ((phi * 180/pi) + 90)
        y <- sin(lambda) * ((phi * 180/pi) + 90)
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
