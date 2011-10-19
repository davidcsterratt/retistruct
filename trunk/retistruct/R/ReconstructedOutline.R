##' Transform an image into the reconstructed space. The four corner
##' coordinates of each pixel are transformed into spherical
##' coordinates and a mask matrix with the same dimensions as
##' \code{im} is created. This has \code{TRUE} for pixels that should
##' be displayed and \code{FALSE} for ones that should not.
##'
##' @title Transform an image into the reconstructed space
##' @param r \code{reconstructedOutline} object
##' @return \code{reconstructedOutline} object with extra elements
##' \item{\code{ims}}{Coordinates of corners of pixes in spherical coordinates}
##' \item{\code{immask}}{Mask matrix with same dimensions as image \code{im}}
##' @author David Sterratt
transform.image.reconstructedOutline <- function(r) {
  if (!is.null(r$im)) {
    ## Need to find the *boundaries* of pixels
    N <- ncol(r$im)
    M <- nrow(r$im)

    ## Create grid coords of corners of pixels.  These run from the
    ## top left of the image down each column of the image.
    xs <- 0:N
    ys <- M:0
    ## x-coords of pixel corners, arranged in (N+1) by (M+1) grid 
    Ix <- outer(ys*0, xs, FUN="+")
    ## Ditto for y-coords
    Iy <- outer(ys, xs*0, FUN="+")
    ## Join to give (x, y) coordinates of all corners
    I  <- cbind(as.vector(Ix), as.vector(Iy))
    
    ## Find Barycentric coordinates of corners of pixels
    Ib <- tsearchn(r$P, r$T, I)
    
    ## Create mask depending on whether corners are in outline
    idx <- matrix(Ib$idx, M+1, N+1)
    r$immask <- (!is.na(idx[1:M    , 1:N    ]) & 
                 !is.na(idx[1:M    , 2:(N+1)]) &
                 !is.na(idx[2:(M+1), 1:N    ]) &
                 !is.na(idx[2:(M+1), 2:(N+1)]))
    
    ## Convert to spherical coordinates
    Ic <- with(r, bary.to.sphere.cart(phi, lambda, R, Tt, Ib))
    r$ims <- with(r, sphere.cart.to.sphere.spherical(Ic, R))
  }
  return(r)
}


##' Get spherical coordinates of tears.
##'
##' @title Get spherical coordinates of tears.
##' @param r \code{reconstructedOutline} object.
##' @return \code{Tss}
##' @method getTss reconstructedOutline
##' @author David Sterratt
getTss.reconstructedOutline <- function(r) {
  Tss <- list()
  for (TF in r$TFset) {
    ## Convert indicies to the spherical frame of reference
    j <- r$ht[TF]
    Tss <- with(r, c(Tss, list(cbind(phi=phi[j], lambda=lambda[j]))))
  }
  return(Tss)
}

##' Plot a mesh of gridlines from the spherical retina (described by
##' points \code{phi}, \code{lambda} and triangulation \code{Tt} and
##' cutoff point \code{phi0}) onto a flattened retina (described by
##' points \code{P} and triangulation \code{T}).
##'
##' @title Flat plot of reconstructed outline
##' @param r \code{reconstructedOutline} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param ... Other plotting parameters
##' @method plot.flat reconstructedOutline
##' @author David Sterratt
plot.flat.reconstructedOutline <- function(r, axt="n", ylim=NULL, ...) {
  NextMethod()

  args <- list(...)
  plot.grid <-   is.null(args$grid) || args$grid
  plot.strain <- !is.null(args$strain) && args$strain

  if (plot.strain) {
    o <- getStrains(r)
    palette(rainbow(100))
    scols <- strain.colours(o$flat$logstrain)
    with(r, 
         segments(P[Cu[,1],1], P[Cu[,1],2],
                  P[Cu[,2],1], P[Cu[,2],2], col=round(scols)))
  }
  
  ## Plot a gridline from the spherical retina (described by points phi,
  ## lambda and triangulation Tt) onto a flattened retina (described by
  ## points P and triangulation T). The gridline is described by a
  ## normal n to a plane and a distance to the plane. The intersection of
  ## the plane and the spehere is the gridline.
  plot.gridline.flat <- function(P, T, phi, lambda, Tt, n, d, ...) {
    mu <- compute.intersections.sphere(phi, lambda, Tt, n, d)

    ## Take out rows that are not intersections. If a plane intersects
    ## one side of a triangle and the opposing vertex, in the row
    ## corresponding to the triangle, there will be a 0, a 1 and a
    ## value between 0 and 1. We get rid of the 1 in the
    ## following. Triangles in which one line is in the plane have mu
    ## values 0, 1 and NaN; we want to include these.
    tri.int <- ((rowSums((mu >= 0) & (mu < 1)) == 2) |
                apply(mu, 1, function(x) setequal(x, c(0, 1, NaN))))

    if (any(tri.int)) {
      T  <- T[tri.int,,drop=FALSE]
      mu <- mu[tri.int,,drop=FALSE]

      ## Create a logical matrix of which points are involved in lines
      ## that interscect the plane.
      line.int <- (mu >= 0) & (mu < 1)
      ## If any element of mu contained a NaN, due to a line being in
      ## the plane, this should be set to false as the point opposite
      ## the NaN is not in the plane
      line.int[is.na(line.int)] <- FALSE

      ## Order rows so that the false indicator is in the third column
      T[!line.int[,2] ,] <- T[!line.int[,2], c(3,1,2)]
      mu[!line.int[,2],] <- mu[!line.int[,2],c(3,1,2)]
      T[!line.int[,1] ,] <- T[!line.int[,1], c(2,3,1)]
      mu[!line.int[,1],] <- mu[!line.int[,1],c(2,3,1)]

      P1 <- mu[,1] * P[T[,3],] + (1-mu[,1]) * P[T[,2],]
      P2 <- mu[,2] * P[T[,1],] + (1-mu[,2]) * P[T[,3],]
      suppressWarnings(segments(P1[,1], P1[,2], P2[,1], P2[,2], ...))
    }
  }
  
  if (plot.grid) {
    grid.int.minor <- 15
    grid.int.major <- 45
    grid.maj.col <- getOption("grid.maj.col")
    grid.min.col <- getOption("grid.min.col")

    phi0d <- r$phi0 * 180/pi
    
    Phis <- seq(-90, phi0d, by=grid.int.minor)
    Lambdas <- seq(0, 180-grid.int.minor, by=grid.int.minor)
    for (Phi in Phis) {
      if ((!(Phi %% grid.int.major) || Phi == phi0d)) {
        col <- grid.maj.col
      } else {
        col <- grid.min.col
      }
      with(r, plot.gridline.flat(P, T, phi, lambda, Tt, c(0,0,1), sin(Phi*pi/180), col=col, ...))
    }
    for (Lambda in Lambdas) {
      if (!(Lambda %% grid.int.major)) {
        col <- grid.maj.col
      } else {
        col <- grid.min.col
      }
      Lambda <- Lambda * pi/180
      with(r, plot.gridline.flat(P, T, phi, lambda, Tt, c(sin(Lambda),cos(Lambda),0), 0, col=col, ...))
    }
  }
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
##' @param ... Other graphics parameters -- not used at present
##' @method plot.polar reconstructedOutline
##' @author David Sterratt
plot.polar.reconstructedOutline <- function(r, show.grid=TRUE,
                                            grid.col="gray",
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45,
                                            flip.horiz=FALSE,
                                            labels=c(0, 90, 180, 270), ...) {
  args <- list(...)
  plot.image <- is.null(args$image) || args$image
  plot.contours <- TRUE # Contour plotting happens only with image plotting at present
      
  phi0d <- r$phi0*180/pi
  par(mar=c(0.5, 0.5, 0.5, 0.5))
  grid.pos <- c(seq(-90, phi0d, by=grid.int.minor), phi0d)
  maxlength <- diff(range(grid.pos))
  if (flip.horiz) {
    xlim <- c(maxlength, -maxlength)
  } else {
    xlim <- c(-maxlength, maxlength)
  }
  plot(NA, NA, xlim=xlim, ylim=c(-maxlength, maxlength), 
       type = "n", axes = FALSE, xlab = "", ylab = "", asp=1)
  if (plot.image && !is.null(r$im)) {
    ## Reconstitute image from stored values of phi and lambda
    ## coordinates of corners of pixels
    with(r, {
      N <- ncol(r$im)
      M <- nrow(r$im)

      ## Compute x and y positions of corners of pixels
      xpos <- matrix(cos(ims[,"lambda"]) * ((ims[,"phi"] * 180/pi) + 90), M+1, N+1)
      ypos <- matrix(sin(ims[,"lambda"]) * ((ims[,"phi"] * 180/pi) + 90), M+1, N+1)

      ## Convert these to format suitable to polygon
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
      polygon(impx[,immask], impy[,immask],  col=im[immask], border=NA)

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

        ## Alternative ways of doing this: use Kriging (from spatial
        ## package). May be better, but slower.
        ##
        ## im.kr <- surf.ls(6, data.frame(x=xposm, y=yposm, z=imm))
        ## im.kr <- surf.gls(2, expcov, data.frame(x=xposm, y=yposm, z=imm),
        ## d=1)
        ## trsurf=trmat(im.kr, -124, 124, -124, 124, 100)
        ## contour(trsurf, add=TRUE, nlevels=20)
      }
    })
  }
  
  ## Plot the grid
  ## Tangnential lines
  angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
  if (show.grid) {
    for (i in seq(length(grid.pos), 1, by = -1)) {
      xpos <- cos(angles) * (grid.pos[i] - grid.pos[1])
      ypos <- sin(angles) * (grid.pos[i] - grid.pos[1])
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
  label.pos <- rlabels - min(grid.pos)
  text(label.pos, -maxlength/15, rlabels)

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
    x <- with(r, cos(Ts[,"lambda"]) * ((Ts[,"phi"] * 180/pi) + 90))
    y <- with(r, sin(Ts[,"lambda"]) * ((Ts[,"phi"] * 180/pi) + 90))
    suppressWarnings(lines(x, y, col=getOption("TF.col"), ...))
  }
}

##' Draw a spherical plot of reconstructed outline. This method just
##' draws the mesh.
##'
##' @title Spherical plot of reconstructed outline
##' @param r \code{reconstructedOutline} object
##' @param ... Other graphics parameters -- not used at present
##' @method plot.spherical reconstructedOutline
##' @author David Sterratt
plot.spherical.reconstructedOutline <- function(r, ...) {
  NextMethod()
  
  args <- list(...)
  plot.strain <- !is.null(args$strain) && args$strain
  plot.surf   <- is.null(args$surf)    ||  args$surf
  
  ## FIXME: This needs to be looked at with a view to replacing
  ## functions in plots.R
  with(r, {
    ## Obtain Cartesian coordinates of points
    P <- sphere.spherical.to.sphere.cart(phi, lambda, R)

    if (plot.surf) {
      ## Outer triangles
      fac <- 1.005
      triangles3d(matrix(fac*P[t(Tt[,c(2,1,3)]),1], nrow=3),
                  matrix(fac*P[t(Tt[,c(2,1,3)]),2], nrow=3),
                  matrix(fac*P[t(Tt[,c(2,1,3)]),3], nrow=3),
                  color="darkgrey", alpha=1)
      
      ## Inner triangles
      triangles3d(matrix(P[t(Tt),1], nrow=3),
                  matrix(P[t(Tt),2], nrow=3),
                  matrix(P[t(Tt),3], nrow=3),
                  color="white", alpha=1)
    }
    
    ## Plot any flipped triangles
    ft <- flipped.triangles(phi, lambda, Tt, R)
    with(ft, points3d(cents[flipped,1], cents[flipped,2], cents[flipped,3],
                      col="blue", size=5))

    ## Shrink so that they appear inside the hemisphere
    fac <- 0.997
    rgl.lines(fac*rbind(P[ht[gb[gb]],1], P[ht[gb],1]),
              fac*rbind(P[ht[gb[gb]],2], P[ht[gb],2]),
              fac*rbind(P[ht[gb[gb]],3], P[ht[gb],3]),
              lwd=3, color=getOption("TF.col"))
    
    fac <- 1.006
    rgl.lines(fac*rbind(P[ht[gb[gb]],1], P[ht[gb],1]),
              fac*rbind(P[ht[gb[gb]],2], P[ht[gb],2]),
              fac*rbind(P[ht[gb[gb]],3], P[ht[gb],3]),
              lwd=3, color=getOption("TF.col"))

    if (plot.strain) {
      o <- getStrains(r)
      palette(rainbow(100))
      scols <- strain.colours(o$spherical$logstrain)

      fac <- 0.999
      P1 <- fac*P[Cut[,1],]
      P2 <- fac*P[Cut[,2],]
      
      width <- 40
      ## Compute displacement vector to make sure that strips are
      ## parallel to surface of sphere
      d <- extprod3d(P1, P2-P1)
      d <- width/2*d/vecnorm(d)
      PA <- P1 - d
      PB <- P1 + d
      PC <- P2 + d
      PD <- P2 - d

      ## This is a ridiculously inefficient way of drawing the strain,
      ## but if you try presenting a color vector, it makes each line
      ## multi-coloured. It has taking HOURS of fiddling round to
      ## discover this! GRRRRRRRRRRRRRRRRR!
      for (i in 1:nrow(PA)) {
        quads3d(rbind(PA[i,1], PB[i,1], PC[i,1], PD[i,1]),
                rbind(PA[i,2], PB[i,2], PC[i,2], PD[i,2]),
                rbind(PA[i,3], PB[i,3], PC[i,3], PD[i,3]),
                color=round(scols[i]), alpha=1)
      }
    }
  })
}

##' Plot the fractional change in length of mesh edges. The length of
##' each edge in the mesh in the reconstructed object is plotted
##' against each edge in the spherical object. The points are
##' colour-coded according to the amount of log strain each edge is
##' under.
##'
##' @title Plot the fractional change in length of mesh edges
##' @param r \code{reconstructedOutline} object
##' @method plot.l.vs.L.reconstructedOutline
##' @author David Sterratt
plot.l.vs.L.reconstructedOutline <- function(r) {
  o <- getStrains(r)$spherical
  op <- par()["mar"]
  par(mar=c(4.5, 4.5, 0.5,0.5))
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  cols <- strain.colours(o$logstrain)
  with(o, plot(L, l, col=cols, pch=20,
               xlim=c(0, max(L, l)), ylim=c(0, max(L, l)),
               xlab="Length on flattened object",
               ylab="Length on reconstructed object", asp=1))
  par(xpd=FALSE)
  abline(0, 1)
  abline(0, 0.75, col="blue")
  abline(0, 1.25, col="red")
  with(o, text(0.7*max(L), 0.7*max(L)*0.75, "25% compressed", col="blue",
               pos=4))
  with(o, text(0.75*max(L), 0.75*max(L)*1.25, "25% expanded", col="red",
               pos=2))
  par(op)
}
