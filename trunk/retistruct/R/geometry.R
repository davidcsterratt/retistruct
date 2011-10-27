##
## Geometry functions
## 

## scalar product of two column matricies
dot <- function(x, y) {
  return(rowSums(x * y))
}

##' @title Vector norm
##' @param X Vector or matrix. 
##' @return If a vector, returns the 2-norm  of the
##' vector. If a matrix, returns the 2-norm of each row of the matrix
##' @author David Sterratt
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
##' @usage P <- line.line.intersection(P1, P2, P3, P4, interior.only=TRUE)
##' @param P1 vector containing x,y coordinates of one end of L1
##' @param P2 vector containing x,y coordinates of other end of L1
##' @param P3 vector containing x,y coordinates of one end of L2
##' @param P4 vector containing x,y coordinates of other end of L2
##' @param interior.only boolean flag indicating whether only
##' intersections inside L1 and L2 should be returned.
##' @return Vector containing x,y coordinates of intersection of L1
##' and L2.  If L1 and L2 are parallel, this is infinite-valued.  If
##' \code{interior.only} is \code{TRUE}, then when the intersection
##' does not occur between P1 and P2 and P3 and P4, a vector
##' containing \code{NA}s is returned.
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
line.line.intersection <- function(P1, P2, P3, P4, interior.only=FALSE) {
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
  
  if (interior.only) {
    ## Compute the fractions of L1 and L2 at which the intersection
    ## occurs
    lambda1 <- -((X-P1[1])*dx1 + (Y-P1[2])*dy1)/(dx1^2 + dy1^2)
    lambda2 <- -((X-P3[1])*dx2 + (Y-P3[2])*dy2)/(dx2^2 + dy2^2)
    if (!((lambda1>0) & (lambda1<1) &
          (lambda2>0) & (lambda2<1))) {
      return(c(NA, NA))
    }
  }
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
                                P[mod1(i+2, N),], P[mod1(i+3, N),],
                                interior.only=TRUE)
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


##' Return points on the unit circle in an anti-clockwise
##' direction. If \code{L} is not specified \code{n} points are
##' returned. If \code{L} is specified, the same number of points are
##' returned as there are elements in \code{L}, the interval between
##' successive points being proportional to \code{L}.
##'
##' @title Return points on the unit circle
##' @param n Number of points
##' @param L Intervals between points
##' @return The cartesian coordinates of the points
##' @author David Sterratt
circle <- function(n=12, L=NULL) {
  if (is.null(L)) {
    angles <- (0:(n-1))/n*2*pi
  } else {
    angles <- (cumsum(L)-L[1])/sum(L)*2*pi
  }
  return(cbind(cos(angles), sin(angles)))
}

##' Find the interections of the plane defined by the normal \code{n} and the
##' distance \code{d} expressed as a fractional distance along the side of
##' each triangle.
##'
##' @title Find the intersection of a plane with triangles on a sphere
##' @param phi Lattitude of grid points on sphere centred on origin.
##' @param lambda Longitude of grid points on sphere centred on origin.
##' @param T Triangulation
##' @param n Normal of plane
##' @param d Distance of plane along normal from origin.
##' @return Matrix with same dimensions as \code{T}. Each row gives
##' the intersection of the plane  with the corresponding triangle in
##' \code{T}. Column 1 gives the fractional distance from vertex 2 to
##' vertex 3. Column 2 gives the fractional distance from vertex 3 to
##' vertex 1. Column 2 gives the fractional distance from vertex 1 to
##' vertex 2. A value of \code{NaN} indicates that the corresponding
##' side lies in the plane. A value of \code{Inf} indicates that the
##' side lies parallel to the plane but outside it.
##' @author David Sterratt
compute.intersections.sphere <- function(phi, lambda, T, n, d) {
  P <- cbind(cos(phi)*cos(lambda),
             cos(phi)*sin(lambda),
             sin(phi))
  return(cbind((d - P[T[,2],] %*% n)/((P[T[,3],] - P[T[,2],]) %*% n),
               (d - P[T[,3],] %*% n)/((P[T[,1],] - P[T[,3],]) %*% n),
               (d - P[T[,1],] %*% n)/((P[T[,2],] - P[T[,1],]) %*% n)))
}


##' Convert locations of points on sphere in spherical coordinates to
##' points in 3D cartesian space
##'
##' @title Convert from spherical to Cartesian coordinates
##' @param phi vector of lattitudes of N points
##' @param lambda vector of longitudes of N points
##' @param R radius of sphere 
##' @return An N-by-3 matrix in which each row is the cartesian (X, Y,
##' Z) coordinates of each point
##' @author David Sterratt
sphere.spherical.to.sphere.cart <- function(phi, lambda, R=1) {
  P <- cbind(R*cos(phi)*cos(lambda),
             R*cos(phi)*sin(lambda),
             R*sin(phi))
  colnames(P) <- c("X", "Y", "Z")
  return(P)
}

##' Given a triangular mesh on a sphere described by mesh locations
##' (\code{phi}, \code{lambda}), a radius \code{R} and a triangulation
##' \code{Tt}, determine the Cartesian coordinates of points \code{cb}
##' given in barycentric coordinates with respect to the mesh.
##'
##' @title Convert barycentric coordinates of points in mesh on sphere
##' to cartesian coordinates 
##' @param phi Lattitudes of mesh points
##' @param lambda Longitudes of mesh points
##' @param R Radius of sphere 
##' @param Tt Triagulation
##' @param cb Object returned by tsearch containing information on the
##' triangle in which a point occurs and the barycentric coordinates
##' within that triangle
##' @return An N-by-3 matrix of the Cartesian coordinates of the points
##' @author David Sterratt
bary.to.sphere.cart <- function(phi, lambda, R, Tt, cb) {
  ## Initialise output
  cc <- matrix(NA, nrow(cb$p), 3)
  colnames(cc) <- c("X", "Y", "Z")

  ## If there are no points, exit
  if (nrow(cb$p) == 0) {
    return(cc)
  }

  ## Find 3D coordinates of mesh points
  P <- sphere.spherical.to.sphere.cart(phi, lambda, R)

  ## Now find locations cc of datapoints in Cartesian coordinates  
  for(i in 1:nrow(cb$p)) {
    cc[i,] <- bary2cart(P[Tt[cb$idx[i],],], cb$p[i,])
  }
  return(cc)
}

##' Convert locations on the surface of a sphere in cartesian
##' (X, Y, Z) coordinates to spherical (phi, lambda) coordinates. 
##'
##' It is assumed that all points are lying on the surface of a sphere
##' of radius R.
##' @title Convert from Cartesian to spherical coordinates
##' @param P locations of points on sphere as N-by-3 matrix with labelled columns "X", "Y" and "Z"
##' @param R radius of sphere 
##' @return N-by-2 Matrix wtih columns ("phi" and "lambda") of locations of points in spherical coordinates 
##' @author David Sterratt
sphere.cart.to.sphere.spherical <- function(P, R=1) {
  return(cbind(phi   =asin(P[,"Z"]/R),
               lambda=atan2(P[,"Y"], P[,"X"])))
}

##' Project spherical coordinate system \eqn{(\phi, \lambda)} to a polar
##' coordinate system \eqn{(\rho, \lambda)} such that the area of each
##' small region is preserved.
##'
##' This requires \deqn{R^2\delta\phi\cos\phi\delta\lambda =
##' \rho\delta\rho\delta\lambda}.  Hence \deqn{R^2\int^{\phi}_{-\pi/2}
##' \cos\phi' d\phi' = \int_0^{\rho} \rho' d\rho'}.  Solving gives
##' \eqn{\rho^2/2=R^2(\sin\phi+1)} and hence
##' \deqn{\rho=R\sqrt{2(\sin\phi+1)}}.
##' 
##' As a check, consider that total area needs to be preserved.  If
##' \eqn{\rho_0} is maximum value of new variable then
##' \eqn{A=2\pi R^2(\sin(\phi_0)+1)=\pi\rho_0^2}. So
##' \eqn{\rho_0=R\sqrt{2(\sin\phi_0+1)}}, which agrees with the formula
##' above.
##' @title Convert lattitude on sphere to radial variable in
##' area-preserving projection
##' @param phi Lattitude
##' @return Coordinate \code{rho} that has the dimensions of length
##' @author David Sterratt
spherical.to.polar.area <- function(phi, R=1) {
  return(R*sqrt(2*(1 + sin(phi))))
}

## Convert polar coordinates to cartesian coordinates
polar.to.cart <- function(r, theta) {
  return(cbind(x=r*cos(theta), y=r*sin(theta)))   
}

##' On a sphere the central angle between two points is defined as the
##' angle whose vertex is the centre of the sphere and that subtends
##' the arc formed by the great circle between the points. This
##' function computes the central angle for two points \eqn{(\phi_1,
##' \lambda_1)}{(phi1, lambda1)} and \eqn{(\phi_2,\lambda_2)}{(phi2,
##' lambda2)}.
##'
##' @title Central angle between two points on a sphere
##' @param phi1 Lattitude of first point
##' @param lambda1 Longitude of first point
##' @param phi2 Lattitude of second point
##' @param lambda2 Longitude of secone point
##' @return Central angle
##' @source Wikipedia \url{http://en.wikipedia.org/wiki/Central_angle}
##' @author David Sterratt
central.angle <- function(phi1, lambda1, phi2, lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)))
}

##' The Karcher mean of a set of points on a manifold is defined as
##' the point whose sum of squared Riemmann distances to the points is
##' minimal. On a sphere using sphereical coordinates this distance
##' can be computed using the formula for central angle.
##'
##' @title Karcher mean on the sphere
##' @param x Matrix of points on sphere as N-by-2 matrix with labelled
##' columns \\code{phi} (lattitude) and \code{lambda} (longitude)
##' @param a logical value indicating whether \code{NA} values should
##' be stripped before the computation proceeds.
##' @param var logical value indicating whether variance should be
##' returned too.
##' @return Vector of means with components named \code{phi} and
##' \code{lambda}. If \code{var} is \code{TRUE}, a list containing
##' mean and varience in elements \code{mean} and \code{var}.
##' @references Heo, G. and Small, C. G. (2006). Form representations
##' and means for landmarks: A survey and comparative
##' study. \emph{Computer Vision and Image Understanding},
##' 102:188-203.
##' @seealso \code{\link{central.angle}}
##' @author David Sterratt
karcher.mean.sphere <- function(x, na.rm=FALSE, var=FALSE) {
  if (na.rm) {
    x <- na.omit(x)
  }
  if (nrow(x) == 0) {
    return(x)
  }
  ## Compute first estimate of mean by computing centroid in 3D and
  ## then finding angle to this
  P <- cbind(cos(x[,"phi"])*cos(x[,"lambda"]),
             cos(x[,"phi"])*sin(x[,"lambda"]),
             sin(x[,"phi"]))
  N <- nrow(P)
  P.mean <- apply(P, 2, mean)
  phi.mean <-    asin(P.mean[3])
  lambda.mean <- atan2(P.mean[2], P.mean[1])

  ## Now minimise sum of squared distances
  if (all(!is.nan(c(phi.mean, lambda.mean)))) {
    opt <- optim(c(phi.mean, lambda.mean),
                 function(p) { sum((central.angle(x[,"phi"], x[,"lambda"], p[1], p[2]))^2) })
    mu <- opt$par
    names(mu) <- c("phi", "lambda")
    if (var) {
      X <- list(mean=mu, var=opt$value/N)
    } else {
      X <- mu
    }
  } else {
    X <- cbind(phi=NaN, lambda=NaN)
  }
  return(X)
}

## karcher.variance.sphere <- function()
