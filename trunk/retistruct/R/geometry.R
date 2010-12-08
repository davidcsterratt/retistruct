##
## Geometry functions
## 

geometry.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev$"))))
}

## scalar product of two column matricies
dot <- function(x, y) {
  return(rowSums(x * y))
}

## Distance norm
vecnorm <- function(X) {
  if (is.vector(X)) {
    return(sqrt(sum(X^2)))
  } else {
    return(sqrt(rowSums(X^2)))
  }
}

## Function to return "signed area" of triangles on a plane
## given points P and a triangulation Pt. Positive sign
## indicates points are anticlockwise direction; negative indicates
## clockwise
tri.area.signed <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]
  AB <- cbind(B-A, 0)
  BC <- cbind(C-B, 0)
  return(0.5 * extprod3d(AB, BC)[,3])
}

## Function to return area of triangles on a plane
## given points P and a triangulation Pt. Positive sign
## indicates points are anticlockwise direction; negative indicates
## clockwise
tri.area <- function(P, Pt) {
  return(abs(tri.area.signed(P, Pt)))
}

##' Determine the intersection of two lines L1 and L2 in two dimensions,
##' using the formula described by Weisstein.
##' 
##' @title Determine intersection between two lines 
##' @usage P <- line.line.intersection(P1, P2, P3, P4)
##' @param P1 vector containing x,y coordinates of one end of L1
##' @param P2 vector containing x,y coordinates of other end of L1
##' @param P3 vector containing x,y coordinates of one end of L2
##' @param P4 vector containing x,y coordinates of other end of L2
##' @return Vector containing x,y coordinates of intersection of L1 and L2.
##' If L1 and L2 are parallel, this is infinite-valued.
##' @source Weisstein, Eric W. "Line-Line Intersection."
##' From MathWorld--A Wolfram Web Resource.
##' \url{http://mathworld.wolfram.com/Line-LineIntersection.html}
##' @author David Sterratt
##' @examples
##' ## Intersection of two intersecting lines
##' line.line.intersection(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
##'
##' ## Two lines that don't intersect
##' line.line.intersection(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
line.line.intersection <- function(P1, P2, P3, P4) {
  P1 <- as.vector(P1)
  P2 <- as.vector(P2)
  P3 <- as.vector(P3)
  P4 <- as.vector(P4)

  dx1 <- P1[1] - P2[1]
  dx2 <- P3[1] - P4[1]
  dy1 <- P1[2] - P2[2]
  dy2 <- P3[2] - P4[2]

  D <- det(rbind(c(dx1, dy1),
                 c(dx2, dy2)))
  if (D==0) {
    return(c(Inf, Inf))
  }
  D1 <- det(rbind(P1, P2))
  D2 <- det(rbind(P3, P4))
  
  X <- det(rbind(c(D1, dx1),
                 c(D2, dx2)))/D
  Y <- det(rbind(c(D1, dy1),
                 c(D2, dy2)))/D
  return(c(X, Y))
}

## Remove identical consecutive rows from a matrix
##
## This is simlar to unique(), but spares rows which are duplicated, but 
## at different points in the matrix
##
remove.identical.consecutive.rows <- function(P) {
  for (i in 2:nrow(P)) {
    if (identical(P[i-1,], P[i,])) {
      return(remove.identical.consecutive.rows(P[-i,]))
    }
  }
  return(P)
}

## Remove intersections between adjacent segements in a circular path
##
## Suppose segments AB and CD intersect.  Point B is replaced by the
## intersection point, defined B'.  Point C is replaced by a point C'
## on the line B'D. The maxium distance of B'C' is given by the
## parameter d. If the distance l B'D is less than 2d, the distance
## B'C' is l/2.
## 
## Arguments:
## P   - The points, as a 2-column matrix
## d   - Criterion for maximum distance when points are inser
##
remove.intersections <- function(P, d=50) {
  N <- nrow(P)
  for (i in 1:N) {
    R <- line.line.intersection(P[i,],            P[mod1(i+1, N),],
                            P[mod1(i+2, N),], P[mod1(i+3, N),])
    if (identical(P[mod1(i+1, N),], P[mod1(i+2, N),])) {
      R <- P[mod1(i+1, N),]
    }
    if (is.finite(R[1])) {
      print("Intersection found. Old points:")
      print(P[i,])
      print(P[mod1(i+1, N),])
      print(P[mod1(i+2, N),])
      print(P[mod1(i+3, N),])

      P[mod1(i+1, N),] <- R
      print("Point i+1 has been changed:")
      print(P[i,])
      print(P[mod1(i+1, N),])
      print(P[mod1(i+2, N),])
      print(P[mod1(i+3, N),])

      l <- vecnorm(P[mod1(i+1, N),] - P[mod1(i+3, N),])
      if (l > 2*d) {
        a <- d/l
      } else {
        a <- 0.5
      }
      print(paste("a=", a))
      print(paste("l=", l))
      P[mod1(i+2, N),] <- a*P[mod1(i+1, N),] + (1-a)*P[mod1(i+3, N),]
      print("New points:")
      print(P[i,])
      print(P[mod1(i+1, N),])
      print(P[mod1(i+2, N),])
      print(P[mod1(i+3, N),])
    }
  }
  return(P)
}

## Create n uniformly spaced points on the unit circle in an anti-clockwise
## direction
circle <- function(n) {
  return(cbind(cos((1:n)/n*2*pi), sin((1:n)/n*2*pi)))
}

## compute.intersections.sphere(phi, lambda, T, n, d)
##
## Find the interections of the plane defined by the normal n and the
## distance d expressed as a fractional distance along the side of
## each triangle.
compute.intersections.sphere <- function(phi, lambda, T, n, d) {
  P <- cbind(cos(phi)*cos(lambda),
             cos(phi)*sin(lambda),
             sin(phi))
  return(cbind((d - P[T[,2],] %*% n)/((P[T[,3],] - P[T[,2],]) %*% n),
               (d - P[T[,3],] %*% n)/((P[T[,1],] - P[T[,3],]) %*% n),
               (d - P[T[,1],] %*% n)/((P[T[,2],] - P[T[,1],]) %*% n)))
}

## Function to convert locations of points on sphere in spherical
## coordinates to points in 3D cartesian space
##
## Arguemnts:
## phi    - lattitude of mesh points
## lambda - longitude of mesh points
## R      - radius of sphere
sphere.spherical.to.sphere.cart <- function(phi, lambda, R) {
  P <- cbind(R*cos(phi)*cos(lambda),
             R*cos(phi)*sin(lambda),
             R*sin(phi))
  colnames(P) <- c("X", "Y", "Z")
  return(P)
}

## Function to determine the locations of points on the reconstructed surface
## in Cartesian (X, Y, Z) coordinates
## phi    - lattitude of mesh points
## lambda - longitude of mesh points
## R      - radius of sphere
## Tt     - triagulation
## cb     - object returned by tsearch containing information on the
##          triangle in which a point occurs and the barycentric coordinates
##          within that triangle
##
bary.to.sphere.cart <- function(phi, lambda, R, Tt, cb) {
  ## Initialise output
  cc <- matrix(0, 0, 3)
  colnames(cc) <- c("X", "Y", "Z")

  ## If there are no points, exit
  if (nrow(cb$p) == 0) {
    return(cc)
  }

  ## Find 3D coordinates of mesh points
  P <- sphere.spherical.to.sphere.cart(phi, lambda, R)

  ## Now find locations cc of datapoints in Cartesian coordinates  
  for(i in 1:nrow(cb$p)) {
    cc <- rbind(cc, bary2cart(P[Tt[cb$idx[i],],], cb$p[i,]))
  }
  return(cc)
}

## Function to convert locations on the surface of a sphere in cartesian
## (X, Y, Z) coordinates to spherical (phi, lambda) coordinates
##
## Arguments:
## Dsc    - locations of points on sphere
## R      - radius of sphere
##
## Returns:
## Matrix wtih columns ("phi" and "lambda") of new locations
##
sphere.cart.to.sphere.spherical <- function(Dsc, R) {
  return(cbind(phi   =asin(Dsc[,"Z"]/R),
               lambda=atan2(Dsc[,"Y"], Dsc[,"X"])))
}

## Convert elevation in spherical coordinates into radius in polar
## coordinates in an area-preserving projection
spherical.to.polar.area <- function(phi) { return(sqrt(2*(1 +
  sin(phi)))) }

## Convert polar coordinates to cartesian coordinates
polar.to.cart <- function(r, theta) {
  return(cbind(x=r*cos(theta), y=r*sin(theta)))   
}

## Compute mean on sphere
sphere.mean.sphere <- function(phi, lambda) {
  ## First estimate of mean
  P <- rbind(cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi))
  P.mean <- apply(P, 1, mean)
  phi.mean <-    asin(P.mean[3])
  lambda.mean <- atan2(P.mean[2], P.mean[1])

  print(c(phi.mean, lambda.mean))
  opt <- optim(c(phi.mean, lambda.mean),
        function(p) { sum((central.angle(phi, lambda, p[1], p[2]))^2) })
  return(opt$par)
  ## return(c(phi.mean, lambda.mean))
}
