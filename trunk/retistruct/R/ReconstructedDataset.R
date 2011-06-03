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
  if (!is.null(r$Ds)) {
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
  if (!is.null(r$Ss)) {
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
  NextMethod(r)

  args <- list(...)
  plot.datapoints <- is.null(args$datapoints) || args$datapoints
  plot.landmarks <- is.null(args$landmarks) || args$landmarks

  phi0d <- r$phi0*180/pi
  grid.pos <- c(seq(-90, phi0d, by=grid.int.minor), phi0d)
  maxlength <- diff(range(grid.pos))

  ## Radial Labels
  if (!is.null(labels)) {
    angles <- seq(0, by=2*pi/length(labels), len=length(labels))
    xpos <- cos(angles) * maxlength * 1.05
    ypos <- sin(angles) * maxlength * 1.05
    text(xpos, ypos, labels)
  }

  ## FIXME: need to think about how we do this plotting, depending on
  ## what we decide about whether the spherical coordinates ought to
  ## be flipped
  ## ylim <- range(ypos)
  ## if (r$DVflip) {
  ##   if (is.null(ylim)) {
  ##     ylim <- range(ypos)
  ##   }
  ##   ylim <- sort(ylim, TRUE)
  ## }
  ## NextMethod(ylim=ylim)
  
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
