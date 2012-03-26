source("../retistruct/R/geometry.R")

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
  gamma <- acos((X - x)/R)
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
##' @param r0 Vector comprising cartesian coordinates of incident ray
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

  ## Compute incident angle to lens surface
  if (j != nrow(S)) {
    alpha1 <- atan2(P1[2], P1[1] - (S$A[j] + S$R[j]))
    I0 <- pi + theta0 - alpha1
    if (I0 >= pi) I0 <- I0 - pi
    ## Snell's law
    I1 <- with(S, asin(n[j]*sin(I0)/n[j+1]))
    theta1 <- I1 + alpha1 - pi
    if (theta1 <=  -pi) theta1 <- theta1 + pi
    message(paste("alpha =", alpha1*180/pi, "; I0 =", I0*180/pi, "; I1 =", I1*180/pi, "; theta0 =", theta0*180/pi, ";theta1 =", theta1*180/pi))
  } else {
    theta1 <- NA
  }
  out <- c(P1[1:2],  theta1, j)
  names(out) <- c("x", "y", "theta", "j")
  return(out)
}

##' @title Draw a arc
##' @param A Left edge of arc
##' @param R Radius of arc
##' @param C Cutoff angle
##' @author David Sterratt
draw.arc <- function(A, R, C=pi/2) {
  X <- A + R
  angles <- seq(-C, C, len=100)
  P <- R*circle(100)
  x <- X - R*cos(angles)
  y <- R*sin(angles)
  lines(x, y)
}

##' @title Trace ray through system
##' @param r0 Initial ray vector
##' @param S System defined by data frame containing left margin of
##' surfaces (A), radii of surfaces (R), cut-offs (C) and refracive
##' indices (n).
##' @return Matrix with x & y coordinates in first column and angles
##' in final colum.
##' @author David Sterratt
trace.ray <- function(r0, S) {
  r0 <- c(r0, NA)
  r <- matrix(r0, ncol=4)
  colnames(r) <- c("x", "y", "theta", "j")

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
##' @author David Sterratt
draw.system <- function(S) {
  plot(NA, NA, xlim=c(-4, 6), ylim=c(-5, 5), asp=1)
  abline(0, 0)
  with(S, {
    for (i in 1:length(A)) {
      draw.arc(A[i], R[i])
    }})
}

## Define system by data frame containing left margin of surfaces (A),
## radii of surfaces (R) and refracive indices (n).
S <- data.frame(A=c(3,  4,  8),
                R=c(1, -1, -4),
                C=c(pi/2, pi/2, pi/2),
                n=c(1, 1.4, 1))

S.hughes <- data.frame(A=c(0,     0.26,  0.881, 1.778,  3.695,  4.591,  5.981),
                       R=c(2.965, 2.705, 2.340, 0.958, -0.958, -2.340, -2.813),
                       C=c( pi/2,  pi/2,  pi/2,  pi/2,  pi/2,   pi/2,   pi/2),
                       n=c(1,     1.38,  1.337, 1.390,  1.50,   1.390,  1.337))

S <- S.hughes
draw.system(S)
#trace.ray(c(0, 0.5, 0), S)
#trace.ray(c(0, 0.5, -0.05), S)
#trace.ray(c(0, 0.5, 0.05), S)
#trace.ray(c(0, 0.5, 0.15), S)
## r <- trace.ray(c(-500, 0.5, -0.002), S)
## r <- trace.ray(c(-500, 0.5, -0.001), S)
## r <- trace.ray(c(-500, 0.5, -0), S)

## r <- trace.ray(c(-500, 10, -0.02), S)

## r <- trace.ray(c(-500, 50, -0.1), S)
## r <- trace.ray(c(-500, 50.5, -0.1), S)

## r <- trace.ray(c(-500, 101, -0.2), S)

## r <- trace.ray(c(-500, 102.5, -0.2), S)
## r <- trace.ray(c(-500, 102, -0.2), S)
## r <- trace.ray(c(-500, 100, -0.2), S)
## r <- trace.ray(c(-500, 100.8, -0.2), S)
## r <- trace.ray(c(-500, 99, -0.2), S)

## r <- trace.ray(c(1, 4, -pi/2), S)
r <- trace.ray(c(1.5, 4, -pi/2), S)
r <- trace.ray(c(-500, 101.5, -0.2), S)

