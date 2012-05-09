## This is the old polarplot code, including the links to the ray
## tracing. There may be some useful stuff here, including the
## projection of aziumuthal and elevation lines onto the polar plots.

## FIXME: it would be nice to tidy up the profusion of conversion
## functions

##' This is a helper function for \code{\link{polarplot}}
##'
##' @title Convert rho variable into normalised radial coordinate in degrees
##' @param rho The rho variable. This may be in radians or the
##' area-preserving coordinate
##' @param phi0 Lattitude of top of curtailed sphere
##' @param area.preserving Whether the conversion should preserve area
##' @return Radial variable
##' @author David Sterratt
##' @export
rho.to.degrees <- function(rho, phi0, area.preserving) {
  if (area.preserving) {
    phi0d <- phi0*180/pi
    rho0 <- spherical.to.polar.area(phi0)
    return((phi0d + 90)/rho0*rho)
  } else {
    return(rho*180/pi)
  }
}


##' Plot polar representation of object
##'
##' @title Polar plot of object
##' @param r \code{\link{ReconstructedOutline}},
##' \code{\link{ReconstructedDataset}}  object
##' @param show.grid Whether or not to show the grid lines of lattitude and longitude
##' @param grid.col Colour of the minor grid lines
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param flip.horiz Wether to flip about a horizontal axis
##' @param labels Vector of 4 labels to plot at 0, 90, 180 and 270 degrees 
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
polarplot <- function(r, show.grid=TRUE,
                       grid.col="gray", grid.bg="transparent", 
                       grid.int.minor=15, grid.int.major=45,
                       flip.horiz=FALSE,
                       labels=c(0, 90, 180, 270),...) {
  UseMethod("polarplot")
}

##' @export
polarplot.default <- function(r, show.grid=TRUE,
                               grid.col="gray", grid.bg="transparent", 
                               grid.int.minor=15, grid.int.major=45,
                               flip.horiz=FALSE, labels=c(0, 90, 180, 270),
                               ...) {
  plot.new()
}


##' Draw a polar plot of reconstructed outline. This method just sets
##' up the grid lines and the angular labels. 
##'
##' @title Polar plot of reconstructed outline
##' @param r \code{reconstructedOutline} object
##' @param show.grid Whether or not to show the grid lines of lattitude and longitude
##' @param grid.col Colour of the minor grid lines
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param flip.horiz Wether to flip about a horizontal axis
##' @param labels Vector of 4 labels to plot at 0, 90, 180 and 270 degrees 
##' @param ... Other parameters, including graphics ones. The option
##' \code{image} causes an image to be plotted if \code{TRUE}
##' (default \code{TRUE}).  The option \code{preserve.area} creates an
##' area-preserving plot (default \code{FALSE}).
##' @method polarplot reconstructedOutline
##' @author David Sterratt
##' @export
polarplot.reconstructedOutline <- function(r, show.grid=TRUE,
                                            grid.col="gray",
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45,
                                            flip.horiz=FALSE,
                                            labels=c(0, 90, 180, 270), ...) {
  args <- list(...)
  plot.image <- is.null(args$image) || args$image
  plot.contours <- FALSE # Contour plotting happens only with image plotting at present
  plot.space <- !is.null(args$space) && args$space
  plot.preserve.area <- !is.null(args$preserve.area) && args$preserve.area
  pa <- plot.preserve.area
  
  phi0d <- r$phi0*180/pi
  
  grid.pos <- c(seq(-90, phi0d, by=grid.int.minor), phi0d)
  maxlength <- diff(range(grid.pos))
  if (flip.horiz) {
    xlim <- c(maxlength, -maxlength)
  } else {
    xlim <- c(-maxlength, maxlength)
  }
  plot(NA, NA, xlim=xlim, ylim=c(-maxlength, maxlength), 
       type = "n", axes = FALSE, xlab = "", ylab = "", asp=1)
  ims <- getIms(r)
  if (plot.image && !is.null(ims)) {
    ## Reconstitute image from stored values of phi and lambda
    ## coordinates of corners of pixels
    N <- ncol(r$im)
    M <- nrow(r$im)

    ## Compute x and y positions of corners of pixels
    xpos <- matrix(cos(ims[,"lambda"])*phi.to.rho(ims[,"phi"], r$phi0, pa), M+1, N+1)
    ypos <- matrix(sin(ims[,"lambda"])*phi.to.rho(ims[,"phi"], r$phi0, pa), M+1, N+1)
    
    ## Convert these to format suitable for polygon()
    impx <- rbind(as.vector(xpos[1:M    , 1:N    ]),
                  as.vector(xpos[1:M    , 2:(N+1)]),
                  as.vector(xpos[2:(M+1), 2:(N+1)]),
                  as.vector(xpos[2:(M+1), 1:N    ]),
                  NA)
    impy <- rbind(as.vector(ypos[1:M    , 1:N    ]),
                  as.vector(ypos[1:M    , 2:(N+1)]),
                  as.vector(ypos[2:(M+1), 2:(N+1)]),
                  as.vector(ypos[2:(M+1), 1:N    ]),
                  NA)

    ## Plot the polygon, masking as we go
    with(r, polygon(impx[,immask], impy[,immask],
                    col=im[immask], border=im[immask]))

    if (plot.contours) {
      ## Find centre locations of polygons
      xposc <- 0.25*(xpos[1:M    , 1:N] +
                     xpos[1:M    , 2:(N+1)] +
                     xpos[2:(M+1), 1:N] +
                     xpos[2:(M+1), 2:(N+1)])
      yposc <- 0.25*(ypos[1:M,     1:N] +
                     ypos[1:M,     2:(N+1)] +
                     ypos[2:(M+1), 1:N] +
                     ypos[2:(M+1), 2:(N+1)])

      ## Compute intensity of image
      imrgb <- col2rgb(r$im)  
      imin <- matrix(0.3*imrgb[1,] + 0.59*imrgb[2,] + 0.11*imrgb[3,],
                     nrow(r$im), ncol(r$im), byrow=TRUE)

      ## Mask all data
      xposm <- xposc[r$immask]
      yposm <- yposc[r$immask]
      imm <- imin[r$immask]

      ## Interporlate. Length 20 (rather than the default 40) gets
      ## rid of some noise, though this is quite a crude way of
      ## doing things
      im.smooth <- interp(xposm, yposm, imm,
                          xo=seq(min(xposm), max(xposm), len=20),
                          yo=seq(min(yposm), max(yposm), len=20))  
      contour(im.smooth, add=TRUE, nlevels=5)
    }
  }
  
  ## Plot the grid
  ## Tangnential lines
  angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
  if (show.grid) {
    for (i in seq(length(grid.pos), 1, by = -1)) {
      xpos <- cos(angles)*phi.to.rho(grid.pos[i]*pi/180, r$phi0, pa)
      ypos <- sin(angles)*phi.to.rho(grid.pos[i]*pi/180, r$phi0, pa)
      if (((grid.pos[i] %% grid.int.major) == 0) || (i == length(grid.pos))) {
        col <- "black"
      } else {
        col <- grid.col
      }
      polygon(xpos, ypos, border = col, col = grid.bg)
    }
  }

  ## Radial lines
  angles <- seq(0, 180-grid.int.minor, by = grid.int.minor)
  col <- rep(grid.col, length(angles))
  col[(angles %% grid.int.major) == 0] <- "black"
  angles <- angles * pi/180
  xpos <- cos(angles) * maxlength
  ypos <- sin(angles) * maxlength
  segments(xpos, ypos, -xpos, -ypos, col=col)

  ## Radial Labels
  rlabels <- c(seq(-90, phi0d, by=grid.int.major), phi0d)
  label.pos <- phi.to.rho(rlabels*pi/180, r$phi0, pa)
  text(label.pos, -maxlength/15, rlabels, xpd=TRUE, pos=2, offset=0)

  ## Tangential Labels
  if (!is.null(labels)) {
    angles <- seq(0, by=2*pi/length(labels), len=length(labels))
    xpos <- cos(angles) * maxlength * 1.1
    ypos <- sin(angles) * maxlength * 1.1
    text(xpos, ypos, labels, xpd=TRUE)
  }
  
  ## Plot outline
  Tss <- getTss(r)
  for (Ts in Tss) {
    ## Plot
    x <- with(r, cos(Ts[,"lambda"])*phi.to.rho(Ts[,"phi"], phi0, pa))
    y <- with(r, sin(Ts[,"lambda"])*phi.to.rho(Ts[,"phi"], phi0, pa))
    suppressWarnings(lines(x, y, col=getOption("TF.col"), ...))
  }

  ## FIXME: This should probably go in RetinalReconstructedOutline.R
  ## Plot the projection of visual space onto the retina
  if (plot.space) {
    alphads <- seq(-90, 90, by = 10)    # elevations
    alphas <- alphads*pi/180
    thetads <- seq(-180, 180, by = 10) # azimuths
    thetas <- thetads*pi/180 # azimuths
    s0azel <- r$s0azel
    if (is.null(s0azel)) {
      ## s0azel <- cbind(alpha=0, theta=0)
      s0azel = cbind(alpha=42, theta=62)*pi/180
    }

    ray.trace <- FALSE
    ## Lines of constant Elevation 
    for (alphad in alphads) {
      ## Coordinates to represent
      sazel <- cbind(alpha=alphad*pi/180, theta=thetas)
      ## Convert to spherical colattitude
      ssc <- azel.to.sphere.colattitude(sazel, s0azel)
      
      ## Weed out points likely not to fall on retina
      ##sazel <- sazel[ssc[,"psi"] < r$phi0 + pi/2,, drop=FALSE]
      ## ssc   <- ssc[ssc[,"psi"] < r$phi0 + pi/2,, drop=FALSE]
      ssc[ssc[,"psi"] > r$phi0 + pi/2,] <- c(NA, NA)
      ## Trace through lens
      if (ray.trace) {
        ssc[,"psi"] <- incident.to.lens.angle(ssc[,"psi"], 0.2, S.mouse.ss.23)
      }
      ## Invert and convert to spherical coordinates
      sss <- cbind(ssc[,"psi"] - pi/2, ssc[,"lambda"] - pi)
      colnames(sss) <- c("phi", "lambda")
      ## Project onto polar plot
      pos <- sphere.spherical.to.polar.cart(sss, pa=pa)
      ##      suppressWarnings(lines(rho.to.degrees(pos[sss[,"phi"] < r$phi0,], r$phi, pa),
      ## col="black", ...))
      suppressWarnings(lines(rho.to.degrees(pos, r$phi, pa),
                             col="red", ...))

      
      tpos <- rho.to.degrees(pos[which(sazel[,"theta"]==0),], r$phi, pa)
      text(tpos[1], tpos[2], alphad, cex=0.8, col="red")
    }

    ## Lines of constant Aziumth
    for (thetad in thetads) {
      ## Coordinates to represent
      sazel <- cbind(alpha=alphas, theta=thetad*pi/180)
      ## Convert to spherical colattitude
      ssc <- azel.to.sphere.colattitude(sazel, s0azel)
      ## Weed out points likely not to fall on retina
      ssc[ssc[,"psi"] > r$phi0 + pi/2,] <- c(NA, NA)
      ## Trace through lens
      if (ray.trace) {
        ssc[,"psi"] <- incident.to.lens.angle(ssc[,"psi"], 0.2, S.mouse.ss.23)
      }
      ## Convert to spherical coordinates
      sss <- cbind(ssc[,"psi"] - pi/2, ssc[,"lambda"] -pi)
      colnames(sss) <- c("phi", "lambda")
      ## Project onto poloar plot
      pos <- sphere.spherical.to.polar.cart(sss, pa=pa)
      suppressWarnings(lines(rho.to.degrees(pos, r$phi, pa),
                             col="red", lty=2, ...))

      tpos <- rho.to.degrees(pos[which(sazel[,"alpha"]==0),], r$phi, pa)
      text(tpos[1], tpos[2], thetad, cex=0.8, col="red")

    }
  }
}


##' Plot datapoints in polar plot
##'
##' @title Polar plot of reconstructed dataset
##' @param r \code{\link{ReconstructedDataset}} object
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
##' @method polarplot reconstructedDataset
##' @author David Sterratt
##' @export
polarplot.reconstructedDataset <- function(r, show.grid=TRUE,
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
    Dss.mean <- getDssMean(r)
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
                col=r$cols[[names(k)[i]]],
                ## drawlabels=FALSE,
                labels=k[[i]]$labels)
        pos <- sphere.spherical.to.polar.cart(k[[i]]$maxs, pa)
        suppressWarnings(points(rho.to.degrees(pos, r$phi0, pa),
                                bg=r$cols[[names(k)[i]]], col="black", 
                                pch=22, cex=1.5, ...))
      }
    }
  }

  ## KDE
  if (plot.grouped.contours) {
    k <- getKR(r)
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
                col=r$cols[[names(k)[i]]],
                ## drawlabels=FALSE,
                labels=k[[i]]$labels)
        points(rho.to.degrees(g$max[1], r$phi0, pa),
               rho.to.degrees(g$max[2], r$phi0, pa),
               pch=23, cex=1, lwd=1, col="black", bg=r$cols[[names(k)[i]]])
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


##' This lablels the poles N, D, T and V
##' 
##' @title Polar plot of reconstructed dataset
##' @param r \code{RetinalReconstructedDataset} object
##' @param show.grid Whether or not to show the grid lines of
##' lattitude and longitude
##' @param grid.col Colour of the minor grid lines
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param ... Other graphics parameters.
##' @method polarplot retinalReconstructedDataset
##' @author David Sterratt
##' @export
polarplot.retinalReconstructedDataset <- function(r, show.grid=TRUE,
                                                   grid.col="gray",
                                                   grid.bg="transparent", 
                                                   grid.int.minor=15,
                                                   grid.int.major=45,  ...) {
  ## This will call polarplot.reconstructedDataset()
  NextMethod(flip.horiz=(r$side=="Left"),
             labels=c("N", "D", "T", "V"))
}



