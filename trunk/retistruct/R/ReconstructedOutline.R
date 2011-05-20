## plot.gridline.flat(P, T, phi, lambda, Tt, n, d)
##
## Plot a gridline from the spherical retina (described by points phi,
## lambda and triangulation Tt) onto a flattened retina (described by
## points P and triangulation T). The gridline is described by a
## normal n to a plane and a distance to the plane. The intersection of
## the plane and the spehere is the gridline.
plot.gridline.flat <- function(P, T, phi, lambda, Tt, n, d, ...) {
  mu <- compute.intersections.sphere(phi, lambda, Tt, n, d)

  ## Take out rows that are not intersections
  tri.int <- (rowSums((mu >=0) & (mu <=1)) == 2)

  if (any(tri.int)) {
    T  <- T[tri.int,,drop=FALSE]
    mu <- mu[tri.int,,drop=FALSE]

    line.int <- (mu >=0) & (mu <=1)
    
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

## plot.gridlines.flat(P, T, phi, lambda, Tt, phi0)
##
## Plot a mesh of gridlines from the spherical retina (described by
## points phi, lambda and triangulation Tt and cutoff point phi0) onto
## a flattened retina (described by points P and triangulation T).
plot.flat.reconstructedOutline <- function(r, axt="n", ylim=NULL, ...) {
  NextMethod()

  args <- list(...)
  plot.grid <-   is.null(args$grid) || args$grid
  plot.strain <- !is.null(args$strain) && args$strain

  if (plot.grid) {
    Phis=(-8:9)*pi/18
    Lambdas=(0:17)*pi/18
    grid.int.minor=15
    grid.int.major=45
    grid.col="gray"

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


plot.polar.reconstructedOutline <- function(r, show.grid=TRUE,
                                            grid.col="gray",
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45, ...) {

  phi0d <- r$phi0*180/pi
  par(mar=c(1, 1, 1, 1))
  grid.pos <- c(seq(-90, phi0d, by=grid.int.minor), phi0d)
  maxlength <- diff(range(grid.pos))
  plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
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
}
