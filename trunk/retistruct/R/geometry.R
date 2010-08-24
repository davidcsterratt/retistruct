##
## Geometry functions
## 

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

## Check whether line from P1 to Q1 and from P2 to Q2 intersect
find.intersection <- function(P1, Q1, P2, Q2) {
  if ((max(P1[1], Q1[1]) < min(P2[1], Q2[1])) ||
      (min(P1[1], Q1[1]) > max(P2[1], Q2[1])) ||
      (max(P1[2], Q1[2]) < min(P2[2], Q2[2])) ||
      (min(P1[2], Q1[2]) > max(P2[2], Q2[2]))) {
    return(FALSE)
  }
  M <- cbind(Q1-P1, -Q2+P2)
  if (!(det(M) == 0)) {
    lambda <- solve(M) %*% (P2-P1)
    if (all((lambda<1) & (lambda>0))) {
      return(list(lambda=lambda, R=(1-lambda[1])*P1 + lambda[1]*Q1))
    }
  }
  return(FALSE)
}

## Remove intersections between adjacent segements in a circular path
remove.intersections <- function(P) {
  N <- nrow(P)
  for(i in 1:N) {
    fi <- find.intersection(P[i,],            P[mod1(i+1, N),],
                            P[mod1(i+2, N),], P[mod1(i+3, N),])
    if (is.list(fi)) {
      print(fi)
      P[mod1(i+1, N),] <- fi$R
      P <- P[-mod1(i+2, N),] 
      return(remove.intersections(P))
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
bary.to.sphere.cart <- function(phi, lambda, R, Tt, cb) {
  ## Find 3D coordinates of mesh points
  P <- sphere.spherical.to.sphere.cart(phi, lambda, R)
  
  ## Now find locations cc of datapoints in Cartesian coordinates
  cc <- matrix(0, 0, 3)
  colnames(cc) <- c("X", "Y", "Z")
  for(i in 1:(dim(cb$p)[1])) {
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
