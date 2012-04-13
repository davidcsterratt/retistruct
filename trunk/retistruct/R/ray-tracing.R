##' @title Compute intersection of ray and arc
##' @param Ptheta Vector comprising cartesian coordinates of incident
##' ray origin and angle with x-axis of incident ray
##' @param A Left-edge of circle on x-axis
##' @param R Radius of circle
##' @param C Angle at which arc cut off
##' @return Vector of x, y and lambda
##' @author David Sterratt
ray.arc.intersection <- function(Ptheta, A, R, C) {
  P <- Ptheta[1:2]
  theta <- Ptheta[3]
  X <- A + R                           # X-co-ord of centre of circle
  x <- P[1]
  y <- P[2]

  ## Solve quadratic
  bh <- (x - X)*cos(theta) + y*sin(theta)
  ac <- ((x - X)^2 + y^2 - R^2)

  ## The ray doesn't intersect the circle at all
  if ((bh^2 - ac) < 0) {
    return(NA)
  }

  ## The ray does intersect the circle; compute the distances along
  ## the ray and the intersection points
  lambda <- c(-bh - sign(R)*sqrt(bh^2 - ac),
              -bh + sign(R)*sqrt(bh^2 - ac))
  x <- x + lambda*cos(theta)
  y <- y + lambda*sin(theta)

  ## Find out if the intersections are on the arc
  cosgamma <- (X - x)/R
  if (any(cosgamma > 1)) {
    warning("cosgamma above 1")
    cosgamma[cosgamma > 1] <- 1
  }
    if (any(cosgamma < -1)) {
    warning("cosgamma below -1")
    cosgamma[cosgamma < -1] <- -1
  }

  gamma <- acos(cosgamma)
  ## print(gamma)
  ##  print(C)
  ## If not, return NA
  if (all(gamma > C)) {
    return(NA)
  }
  out <- cbind(x, y, lambda)
  ## print(out)
  out <- out[gamma <= C,,drop=FALSE]
  ## print(out)
  out <- out[out[,"lambda"]>1e-10,,drop=FALSE]
  ## print(out)
  if (nrow(out) == 0) {
    return(NA)
  }
  ## print(out)
  out <- out[which.min(out[,"lambda"]),]
  ## print(out)
  
  ## Otherwise return the first intersection
  return(out)
}

##' @title Compute origin and deflection of new ray
##' @param r0 Vector comprising Cartesian coordinates of incident ray
##' origin, angle of incident ray with x-axis and index of system
##' surface that origin lies on, \code{NA} if the origin is not
##' associated with any surface.
##' @param S S System defined by data frame containing left margin of
##' surfaces (A), radii of surfaces (R), cut-offs (C) and refracive
##' indices (n).
##' @return Vector comprinsing new origin and angle
##' @author David Sterratt
new.ray <- function(r0, S) {
  P0 <- r0[1:2]
  theta0 <- r0[3]
  j <- r0[4]
  ## Compute the next intersection with a surface
  lambda.min <- Inf                     # Minmum distance along ray
  P1 <- NA

  ## Indicies of surfaces to search trhough. Ignore the one associated
  ## with the origin
  iarc <- 1:nrow(S)

  ## Do the seraching
  ## print(iarc)
  for (i in iarc) {
    P <- with(S, ray.arc.intersection(r0, A[i], R[i], C[i]))
    ## print(lambda.min)
    ## print(P)
    if (!any(is.na(P))) {
      if (P[3] < lambda.min) {
        lambda.min <- P[3]
        P1 <- P
        j <- i
      }
    }
  }
  if (all(is.na(P1))) {
    return(NA)
  }

  alpha <- NA
  ## Compute incident angle to lens surface
  ## if (j != nrow(S)) {
    ## Angle of normal for spherical arc
    alpha <- -asin(P1[2]/S$R[j])
    ## Incident angle of beam, mapped onto the range -pi/2, pi/2
    I0 <- theta0 - alpha
    I0 <- ((I0 + pi/2) %% pi) - pi/2
    ## Snell's law
    I1 <- with(S, asin(n[j]*sin(I0)/n[j+1]))
    ## theta1 <- I1 + alpha
    theta1 <- theta0 + I1 - I0
    message(sprintf("theta0 = % 03.1f; alpha = % 03.1f; I0 = % 03.1f; I1 = % 03.1f; theta1 = % 03.1f", theta0*180/pi, alpha*180/pi, I0*180/pi,  I1*180/pi, theta1*180/pi))
 ##  } else {
    ## theta1 <- NA
##   }
  out <- c(P1[1:2],  theta1, j, alpha)
  names(out) <- c("x", "y", "theta", "j", "alpha")
  return(out)
}

##' @title Draw a arc
##' @param A Left edge of arc
##' @param R Radius of arc
##' @param C Cutoff angle
##' @param ... Graphics parameters
##' @author David Sterratt
draw.arc <- function(A, R, C=pi/2, ...) {
  X <- A + R
  angles <- seq(-C, C, len=100)
  P <- R*circle(100)
  x <- X - R*cos(angles)
  y <- R*sin(angles)
  lines(x, y, ...)
}

##' @title Trace ray through system
##' @param r0 Initial ray vector
##' @param S System defined by data frame containing left margin of
##' surfaces (A), radii of surfaces (R), cut-offs (C) and refracive
##' indices (n).
##' @return 5-column matrix with Cartesian coordinates of ray origin
##' in first two columns (\code{x} and \code{y}), angles of rays to
##' x-axis in third column (\code{theta}), index of surface at which
##' intersection occurs (\code{j}) and angle to x-axis at which
##' intersection occurs on surface (\code{alpha}). 
##' @author David Sterratt
trace.ray <- function(r0, S) {
  r0 <- c(r0, NA, NA)
  r <- matrix(r0, ncol=5)
  colnames(r) <- c("x", "y", "theta", "j", "alpha")

  iter <- 10
  while(!is.na(r0[3]) & iter) {
    ## print(r0)
    r1 <- new.ray(r0, S)
    if (!is.na(r1[1])) {
      r <- rbind(r, r1)
      r0 <- r1
    }
    iter <- iter - 1
  }
  lines(r[,1], r[,2])
  return(r)
}

##' @title Draw a system of lenses
##' @param S System defined  by data frame containing left margin of
##' surfaces (A), radii of surfaces (R)  and refracive indices (n).
##' @param xoff x offset of left edge
##' @param ypad y padding
##' @param add If \code{TRUE} add to an existing plot
##' @param ... Graphics parameters
##' @author David Sterratt
draw.system <- function(S, xoff=-10, ypad=1, add=FALSE, ...) {
  with(S, {
    xlim <- range(A)
    xlim[1] <- xlim[1] + xoff
    ylim <- c(min(-R*sin(C) - ypad), max(R*sin(C) + ypad))
    print(ylim)
    if (!add)
      plot(NA, NA, xlim=xlim, ylim=ylim, asp=1)
    abline(0, 0)
    for (i in 1:length(A)) {
      draw.arc(A[i], R[i], C[i], ...)
    }})
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Trace ray at an angle
##' @param psi Angle from axis
##' @param anp Nodal point (hopefully anterior nodal point)
##' @param S System
##' @return Ray - see \code{\link{trace.ray}}
##' @author David Sterratt
trace.ray.angle <- function(psi, anp, S) {
  return(trace.ray(c(anp-2*cos(psi), -2*sin(psi), psi), S))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Convert incident angle to angle on lens
##' @param psi Incident angles
##' @param anp Anterior nodal point on x-axis
##' @param S Lens system
##' @return Angles on lens
##' @examples
##' S <- S.mouse.ss.23
##' draw.system(S, xoff=-2)
##' imfile <- "eye23.png"
##' im <- as.raster(readPNG(imfile))
##' scale <- 0.0086
##' xoff <- -0.52
##' yoff <- 0.31
##' rasterImage(im, xoff, yoff-scale*nrow(im)/2, xoff+scale*ncol(im), yoff+scale*nrow(im)/2)

##' draw.system(S, add=TRUE, col="yellow")

##' anp <- 1
##' psi <- 0.1
##' psi <- 0.2

##' ## Trace ray angles using nodal point at 0.2 - not too bad an approximation
##' ## trace.ray.angle(0.9, 0.2, S)

##' ## All that's needed is to get the angle at which the ray hits the final lens.
##' ## Then we have a mapping from psi to phi.

##' ## r <- trace.ray(c(-500, 101.5, -0.2), S)

##' angles <- seq(0, 90, by=10)
##' lens.angles <- 180/pi*incident.to.lens.angle(angles*pi/180, 0.2, S)
##' @export
##' @author David Sterratt
incident.to.lens.angle <- function(psi, anp, S) {
  return(sapply(psi, function(psi) {
    r <- trace.ray.angle(psi, anp, S)
    r0 <- r[nrow(r),]
    S0 <- S[r0["j"],]
    return(with(S0, atan2(r0["y"], r0["x"] - A - R)))
  }))
}



## Define system by data frame containing left margin of surfaces (A),
## radii of surfaces (R) and refracive indices (n).


## Actual one
## S.mouse.ss.23 <- data.frame(A=c(0,      0.06,  0.263,  2.039,  2.863, 3.037),
##                             R=c(1.43,  1.38,  0.95, -1.0, -1.47, -1.61),
##                             C=c(0.8*pi/2,   0.8*pi/2,   0.78*pi/2 ,1.05*pi/2,   1.2*pi/2,  1.1*pi/2),
##                             n=c(1,      1.4060, 1.3376, 1.6778, 1.3365, 1.3365))

##' @export
S.mouse.ss.23 <- data.frame(A=c(0,        0.06,    0.263,  2.039,  2.863, 3.037),
                            R=c(1.43,     1.38,     0.95, -1.0, -1.35, -1.5),
                            C=c(0.8*pi/2, 0.8*pi/2, 1.09*pi/2, 0.8*pi/2,   1.1*pi/2,  1.1*pi/2),
                            n=c(1,      1.4015, 1.3336, 0.0005*23+1.557, 1.3329, 1.351))


