source("../retistruct/R/geometry.R")

##' @title Compute intersection of ray and cirle
##' @param Ptheta0 Vector comprising cartesian coordinates of incident
##' ray origin and angle with x-axis of incident ray
##' @param A Left-edge of circle on x-axis
##' @param R Radius of circle
##' @return 
##' @author David Sterratt
ray.lens.intersection <- function(Ptheta, A, R) {
  P <- Ptheta[1:2]
  theta <- Ptheta[3]
  X <- A + R                           # X-co-ord of centre of circle
  x <- P[1]
  y <- P[2]

  ## Solve quadratic

  bh <- (x - X)*cos(theta) + y*sin(theta)
  ac <- ((x - X)^2 + y^2 - R^2)
  if ((bh^2 - ac) < 0) {
    return(NA)
  }
  lambda <- -bh -sign(R)*sqrt(bh^2 - ac)
  return(c(x + lambda*cos(theta), y + lambda*sin(theta)))
}

##' @title Compute origin and deflection of new ray
##' @param Ptheta0 Vector comprising cartesian coordinates of incident
##' ray origin and angle with x-axis of incident ray
##' @param A1 Left-edge of circle on x-axis
##' @param R1 Radius of circle
##' @param n0 Refractive index of medium to left of circle
##' @param n1 Refractive index of medium to right of circle
##' @return Vector comprinsing new origin and angle
##' @author David Sterratt
new.ray <- function(Ptheta0, A1, R1, n0, n1) {
  P0 <- Ptheta0[1:2]
  theta0 <- Ptheta0[3]
  ## Compute intersection
  P1 <- ray.lens.intersection(Ptheta0, A1, R1)
  if (any(is.na(P1))) {
    return(NA)
  }
  ## Compute incident angle to lens surface
  alpha1 <- sin(P1[2]/R1)
  I0 <- theta0 + alpha1
  ## Snell's law
  I1 <- asin(n0*sin(I0)/n1)
  theta1 <- I1 - alpha1
  return(c(P1, theta1))
}

##' @title Draw a lens
##' @param A Left edge of lens
##' @param R Radius of lens
##' @param C Cutoff angle
##' @author David Sterratt
draw.lens <- function(A, R, C=pi/2) {
  X <- A + R
  angles <- seq(-C, C, len=100)
  P <- R*circle(100)
  x <- X - R*cos(angles)
  y <- R*sin(angles)
  lines(x, y)
}

##' @title Trace ray through system
##' @param r0 Initial ray vector
##' @param S System defined  by data frame containing left margin of
##' surfaces (A), radii of surfaces (R)  and refracive indices (n).
##' @return Matrix with x & y coordinates in first column and angles
##' in final colum.
##' @author David Sterratt
trace.ray <- function(r0, S) {
  r <- matrix(r0, ncol=3)
  with(S, {
    for (i in 1:(nrow(S))) {
      r1 <- new.ray(r0, A[i], R[i], n[i], n[i+1])
      if (!is.na(r1[1])) {
        r <- rbind(r, r1)
        r0 <- r1
      }
    }
    ## lines(c(r[-nrow(r),1], r[-1,1]), c(r[-nrow(r),2], r[-1,2]))
    lines(r[,1], r[,2])
    return(r)
  })
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
      draw.lens(A[i], R[i])
    }})
}

## Define system by data frame containing left margin of surfaces (A),
## radii of surfaces (R) and refracive indices (n).
S <- data.frame(A=c(3,  4,  8),
                R=c(1, -1, -4),
                n=c(1, 1.4, 1))

S.hughes <- data.frame(A=c(0,     0.26,  0.881, 1.778,  3.695,  4.591,  5.981),
                       R=c(2.965, 2.705, 2.340, 0.958, -0.958, -2.340, -2.813),
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

r <- trace.ray(c(-500, 10, -0.02), S)

r <- trace.ray(c(-500, 50, -0.1), S)
r <- trace.ray(c(-500, 50.5, -0.1), S)

r <- trace.ray(c(-500, 101, -0.2), S)
r <- trace.ray(c(-500, 101.5, -0.2), S)
r <- trace.ray(c(-500, 102, -0.2), S)
r <- trace.ray(c(-500, 103, -0.2), S)

r <- trace.ray(c(1, 4, -pi/2), S)







6.291 - 1.7 - 1.778
