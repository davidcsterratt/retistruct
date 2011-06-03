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
    grid.col <- "gray"

    phi0d <- r$phi0 * 180/pi
    
    Phis <- seq(-90, phi0d, by=grid.int.minor)
    Lambdas <- seq(0, 180-grid.int.minor, by=grid.int.minor)
    for (Phi in Phis) {
      if ((!(Phi %% grid.int.major) || Phi == phi0d)) {
        col <- "black"
      } else {
        col <- grid.col
      }
      with(r, plot.gridline.flat(P, T, phi, lambda, Tt, c(0,0,1), sin(Phi*pi/180), col=col, ...))
    }
    for (Lambda in Lambdas) {
      if (!(Lambda %% grid.int.major)) {
        col <- "black"
      } else {
        col <- grid.col
      }
      Lambda <- Lambda * pi/180
      with(r, plot.gridline.flat(P, T, phi, lambda, Tt, c(sin(Lambda),cos(Lambda),0), 0, col=col, ...))
    }
  }
  if (plot.strain) {
    o <- compute.strain(r)
    palette(rainbow(100))
    scols <- strain.colours(log(o$strain))
    with(r, 
         segments(P[Cu[,1],1], P[Cu[,1],2],
                  P[Cu[,2],1], P[Cu[,2],2], col=round(scols)))
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
##' @param ... Other graphics parameters -- not used at present
##' @method plot.polar reconstructedOutline
##' @author David Sterratt
plot.polar.reconstructedOutline <- function(r, show.grid=TRUE,
                                            grid.col="gray",
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45,
                                            flip.horiz=FALSE, ...) {

  phi0d <- r$phi0*180/pi
  par(mar=c(1, 1, 1, 1))
  grid.pos <- c(seq(-90, phi0d, by=grid.int.minor), phi0d)
  maxlength <- diff(range(grid.pos))
  if (flip.horiz) {
    xlim <- c(maxlength, -maxlength)
  } else {
    xlim <- c(-maxlength, maxlength)
  }
  plot(NA, NA, xlim=xlim, ylim=c(-maxlength, maxlength), 
       type = "n", axes = FALSE, xlab = "", ylab = "", asp=1)

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
  labels <- c(seq(-90, phi0d, by=grid.int.major), phi0d)
  label.pos <- labels - min(grid.pos)
  text(label.pos, -maxlength/15, labels)

  ## Plot outline
  Tss <- getTss(r)
  for (Ts in Tss) {
    ## Plot
    x <- with(r, cos(Ts[,"lambda"]) * ((Ts[,"phi"] * 180/pi) + 90))
    y <- with(r, sin(Ts[,"lambda"]) * ((Ts[,"phi"] * 180/pi) + 90))
    suppressWarnings(lines(x, y, ...))
  }
}
