## This code was used to produce a Voronoi tesselation in projection.ReconstructedDataset()

plot.mosaic <- !is.null(args$mosaic) && args$mosaic

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

##' @title Plot a Voronoi mosaic in a region with a circular boundary
##' @param x \code{voronoi} object
##' @param R Radius at which to trim circle
##' @param ... Additional low-level graphics parameters
##' @author David Sterratt
plot.voronoi.circular <- function (x, R=110, ...) {
  chop.boundary <- function(x0, x1, y0, y1) {
    R12 <- x1^2 + y1^2
    if (R12 > R^2) {
      R02 <- x0^2 + y0^2
      ##                lambda <- sqrt((R^2 - R02)/(R02 + R12))
      l2 <- (x1 - x0)^2 + (y1 - y0)^2
      R0R1 <- x1*x0 + y0*y1
      lambda <- 1/l2*(-(R0R1 - R02) + sqrt(R0R1^2 - R12*R02 + R^2*l2))
      
      x1 <- (1 - lambda)*x0 + lambda*x1
      y1 <- (1 - lambda)*y0 + lambda*y1
    }
    return(c(x1, y1))
  }
  n <- length(x$x)
  x0 <-  matrix(NA, n, 3)
  x1 <-  matrix(NA, n, 3)
  y0 <-  matrix(NA, n, 3)
  y1 <-  matrix(NA, n, 3)
  lty <- matrix(NA, n, 3)
  for (i in 1:n) {
    if (x$node[i]) {
      if (x$x[i]^2 + x$y[i]^2 <  R^2) {
        tns <- sort(c(x$n1[i], x$n2[i], x$n3[i]))
        for (j in 1:3) {
          if (tns[j] > 0) {
            if (x$node[tns[j]]) {
              x0[i, j] <- x$x[i]
              x1[i, j] <- x$x[tns[j]]
              y0[i, j] <- x$y[i]
              y1[i, j] <- x$y[tns[j]]
              ## If the node is outside the boudary, chop at the boundary
              r1 <- chop.boundary(x0[i, j], x1[i, j], y0[i, j], y1[i, j])
              x1[i, j] <- r1[1]
              y1[i, j] <- r1[2]
              lty[i, j] <- 1
            }
          } else {
            if (tns[j] < 0) {
              x0[i, j] <- x$x[i]
              x1[i, j] <- x$dummy.x[-tns[j]]
              y0[i, j] <- x$y[i]
              y1[i, j] <- x$dummy.y[-tns[j]]
              r1 <- chop.boundary(x0[i, j], x1[i, j], y0[i, j], y1[i, j])
              x1[i, j] <- r1[1]
              y1[i, j] <- r1[2]
              lty[i, j] <- 1
            }
          }
        }
      }
    }
  }
  x0 <- na.omit(as.vector(x0))
  x1 <- na.omit(as.vector(x1))
  y0 <- na.omit(as.vector(y0))
  y1 <- na.omit(as.vector(y1))
  lty <- na.omit(as.vector(lty))
  segments(x0, y0, x1, y1, lty=lty, ...)
}

