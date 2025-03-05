##' Stretch the mesh in the flat retina to a circular outline
##'
##' @title Stretch mesh
##' @param Cu Edge matrix
##' @param L Lengths in flat outline
##' @param i.fix Indices of fixed points
##' @param P.fix Coordinates of fixed points
##' @return New matrix of 2D point locations
##' @author David Sterratt
stretchMesh <- function(Cu, L, i.fix, P.fix) {
  N <- max(Cu)
  M <- length(L)
  C <- matrix(0, 2*N, 2*N)
  for (i in 1:M) {
    C[2*Cu[i,1]-1:0,2*Cu[i,2]-1:0] <- diag(2) / L[i]
    C[2*Cu[i,2]-1:0,2*Cu[i,1]-1:0] <- diag(2) / L[i]
  }

  dupC <- duplicated(C)
  if (any(dupC)) {
    i <- which(dupC)
    Ci <- C[i,]
    dups <- which(apply(C, 1, function(x) {identical(x, Ci)}))
    message(paste("dups", dups))
    message(paste("Ci", Ci))
    for (d in dups) {
      message(paste("d", d, ":", which(C[d,]==1)))
    }
  }

  ind <- as.vector(rbind(2*i.fix-1, 2*i.fix))
  A <- C[-ind, -ind]
  B <- C[-ind,  ind]
  P <- matrix(t(P.fix), ncol(B), 1)
  D <- diag(apply(cbind(A, 2*B), 1, sum))
  # if (is.infinite(det(D))) stop ("det(D) is infinite")
  Q <- 2 * solve(D - A) %*% B %*% P
  Q <- matrix(Q, nrow(Q)/2, 2, byrow=TRUE)
  R <- matrix(0, nrow(Q) + length(i.fix), 2)
  R[i.fix,] <- P.fix
  R[-i.fix,] <- Q
  return(R)
}


##
## Energy/error functions
##

## Calculate lengths of connections on sphere
compute.lengths <- function(phi, lambda, Cu, R) {
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)

  return(l)
}

## Calculate lengths of connections on sphere
##' @importFrom geometry dot
compute.areas <- function(phi, lambda, Tr, R) {
  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Find areas of all triangles
  areas <- -0.5/R * dot(P[Tr[,1],], extprod3d(P[Tr[,2],], P[Tr[,3],]), 2)

  return(areas)
}

##' Piecewise, smooth function that increases linearly with negative arguments.
##' \deqn{    f(x) = \left\{
##'        \begin{array}{ll}
##'          -(x - x_0/2) & x < 0 \\
##'          \frac{1}{2x_0}(x - x_0)^2 & 0 < x <x_0 \\
##'          0 & x \ge x_0
##'          \end{array} \right.
##' }
##'
##' @title Piecewise smooth function used in area penalty
##' @param x Main argument
##' @param x0 The cut-off parameter. Above this value the function is zero.
##' @return The value of the function.
##' @author David Sterratt
f <- function(x, x0) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -(x[c1] - x0/2)
  y[c2] <- 1/2/x0*(x0 - x[c2])^2
  y[c3] <- 0

  return(y)
}

##' Derivative of \code{\link{f}}
##'
##' @title Piecewise smooth function used in area penalty
##' @param x Main argument
##' @param x0 The cut-off parameter. Above this value the function is zero.
##' @return The value of the function.
##' @author David Sterratt
fp <- function(x, x0) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -1
  y[c2] <- -1/x0*(x0 - x[c2])
  y[c3] <- 0

  return(y)
}

##' The function that computes the energy (or error) of the
##' deformation of the mesh from the flat outline to the sphere. This
##' depends on the locations of the points given in spherical
##' coordinates. The function is designed to take these as a vector
##' that is received from the \code{optim} function.
##'
##' @title The deformation energy function
##' @param p Parameter vector of \code{phi} and \code{lambda}
##' @param Cu The upper part of the connectivity matrix
##' @param C The connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param B Connectivity matrix
##' @param Tr Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of the sphere
##' @param Rset Indices of points on the rim
##' @param i0 Index of fixed point on rim
##' @param phi0 Latitude at which sphere curtailed
##' @param lambda0 Longitude of fixed points
##' @param Nphi Number of free values of \code{phi}
##' @param N Number of points in sphere
##' @param alpha Area scaling coefficient
##' @param x0 Area cut-off coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A single value, representing the energy of this particular
##' configuration
##' @author David Sterratt
E <- function(p, Cu, C, L, B, Tr, A, R, Rset, i0, phi0, lambda0, Nphi, N,
              alpha=1, x0,  nu=1, verbose=FALSE) {
  ## Extract phis and lambdas from parameter vector
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ## Find cartesian coordinates of points
  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Compute elastic energy
  return(Ecart(P, Cu, L, Tr, A, R,
               alpha, x0, nu, verbose))

}

##' The function that computes the gradient of the  energy (or error)
##' of the deformation of the mesh from the flat outline to the
##' sphere. This depends on the locations of the points given in
##' spherical coordinates. The function is designed to take these as a
##' vector that is received from the \code{optim} function.
##'
##' @title The deformation energy gradient function
##' @param p Parameter vector of \code{phi} and \code{lambda}
##' @param Cu The upper part of the connectivity matrix
##' @param C The connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param B Connectivity matrix
##' @param Tr Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of the sphere
##' @param Rset Indices of points on the rim
##' @param i0 Index of fixed point on rim
##' @param phi0 Latitude at which sphere curtailed
##' @param lambda0 Longitude of fixed points
##' @param Nphi Number of free values of \code{phi}
##' @param N Number of points in sphere
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cut-off coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A vector representing the derivative of the energy of this
##' particular configuration with respect to the parameter vector
##' @author David Sterratt
dE <- function(p, Cu, C, L, B, Tr, A, R, Rset, i0, phi0, lambda0, Nphi, N,
               alpha=1, x0, nu=1, verbose=FALSE) {
  ## Extract phis and lambdas from parameter vector
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  cosp <- cos(phi)
  cosl <- cos(lambda)
  sinl <- sin(lambda)
  sinp <- sin(phi)
  ## Find cartesian coordinates of points
  P <- R * cbind(cosp*cosl,
                 cosp*sinl,
                 sinp)

  ## Compute force in Cartesian coordinates
  dE.dp <- -Fcart(P, C, L, Tr, A, R,
                  alpha, x0, nu, verbose)

  ## Convert to Spherical coordinates
  dp.dphi <- R * cbind(-sinp * cosl,
                       -sinp * sinl,
                       cosp)
  dp.dlambda <- R * cbind(-cosp * sinl,
                          cosp * cosl,
                          0)

  dE.dphi    <- rowSums(dE.dp * dp.dphi)
  dE.dlambda <- rowSums(dE.dp * dp.dlambda)

  ## Return, omitting uncessary indices
  return(c(dE.dphi[-Rset], dE.dlambda[-i0]))
}

##' The function that computes the energy (or error) of the
##' deformation of the mesh from the flat outline to the sphere. This
##' depends on the locations of the points given in spherical
##' coordinates. The function is designed to take these as a vector
##' that is received from the \code{optim} function.
##'
##' @title The deformation energy function
##' @param P N-by-3 matrix of point coordinates
##' @param Cu The upper part of the connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param Tr Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of sphere
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cut-off coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A single value, representing the energy of this particular
##' configuration
##' @author David Sterratt
Ecart <- function(P, Cu, L, Tr, A, R,
                  alpha=1, x0, nu=1, verbose=FALSE) {
  ## Compute elastic energy

  ## using the upper triagular part of the
  ## connectivity matrix Cu to extract coordinates of end-points of
  ## each edge
  ## Compute lengths of edges
  ## l <- vecnorm(P2 - P1)
  l <- 2*R*asin(vecnorm(P[Cu[,2],] - P[Cu[,1],])/2/R)
  if (verbose==2) { report(l) }

  ## Compute spring energy
  E.E <- 0.5/sum(L)*sum((l - L)^2/L)
  if (verbose>=1) { report(E.E) }

  ## Compute areal penalty term if alpha is nonzero
  E.A <- 0
  if (alpha) {
    ## Find signed areas of all triangles
    a <- -0.5/R * dot(P[Tr[,1],], extprod3d(P[Tr[,2],], P[Tr[,3],]), 2)

    ## Now compute area energy
    E.A <- sum((A/mean(A))^nu*f(a/A, x0=x0))
    ## E.A <- sum(f(a/A, x0=x0))
  }
  return(E.E + alpha*E.A)
}

##' The function that computes the gradient of the  energy (or error)
##' of the deformation of the mesh from the flat outline to the
##' sphere. This depends on the locations of the points given in
##' spherical coordinates. The function is designed to take these as a
##' vector that is received from the \code{optim} function.
##'
##' @title The deformation energy gradient function
##' @param P N-by-3 matrix of point coordinates
##' @param C The connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param Tr Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of sphere
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cut-off coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A vector representing the derivative of the energy of this
##' particular configuration with respect to the parameter vector
##' @author David Sterratt
##' @useDynLib retistruct
Fcart <- function(P, C, L, Tr, A, R,
                  alpha=1, x0, nu=1, verbose=FALSE) {
  ## Compute derivative of elastic energy

  ## Lengths of springs
  dP <- P[C[,2],] - P[C[,1],]
  l <- 2*R*asin(vecnorm(dP)/2/R)
  if (verbose==2) { report(l) }

  ## Compute general scaling factor
  fac <- 1/sum(L)*(l - c(L, L))/c(L, L)/c(L, L) #sqrt(1-(d/2/R)^2)/d

  ## Now compute the derivatives
  ## This method is slower than using dense matrix multiplication
  ## dF <- fac * dP

  ## F.E <- matrix(0, nrow(P), 3)
  ## for (i in 1:nrow(C)) {
  ##   F.E[C[i,1],] <- F.E[C[i,1],] + dF[i,]
  ## }

  ## Now compute the derivatives
  ## F.E <- as.matrix(B %*% (fac * dP))

  F.E <- matrix(0, nrow(P), 3)
  F.E <- .Call("sum_force_components", fac*dP, C, F.E, PACKAGE="retistruct")

  ## Compute the derivative of the area component if alpha is nonzero
  if (alpha) {
    dEdpi <- matrix(0, nrow(P), 3)
    ## Here follows computation of the derivative - it's a bit
    ## complicated!

    ## Expand triangulation so that every point is the first point
    ## once. The number of points is effectively tripled.
    Tr <- rbind(Tr, Tr[,c(2,3,1)], Tr[,c(3,1,2)])
    A <- c(A, A, A)

    ## Compute the derivative of area with respect to the first points
    dAdPt1 <- -0.5/R * extprod3d(P[Tr[,2],], P[Tr[,3],])

    ## Find areas of all triangles
    a <- dot(P[Tr[,1],], dAdPt1, 2)

    ## Now convert area derivative to energy derivative
    dEdPt1 <- -(A/mean(A))^nu*fp(a/A, x0=x0)/A*dAdPt1
    ## dEdPt1 <- -fp(a/A, x0=x0)/A*dAdPt1

    ## Map contribution of first vertex of each triangle back onto the
    ## points
    ## for(m in 1:nrow(T)) {
    ##   dEdpi[Tr[m,1],] <- dEdpi[Tr[m,1],] - dEdPt1[m,]
    ## }
    dEdpi <- .Call("sum_force_components", -dEdPt1, Tr, dEdpi , PACKAGE="retistruct")

    F.E <- F.E - alpha*dEdpi
  }
  return(F.E)
}

##' Restore points to spherical manifold after an update of the
##' Lagrange integration rule
##'
##' @title Restore points to spherical manifold
##' @param P Point positions as N-by-3 matrix
##' @param R Radius of sphere
##' @param Rset Indices of points on rim
##' @param i0 Index of fixed point
##' @param phi0 FullCut-off of curtailed sphere in radians
##' @param lambda0 Longitude of fixed point on rim
##' @return Points projected back onto sphere
##' @author David Sterratt
Rcart <- function(P, R, Rset, i0, phi0, lambda0) {

  ## Now ensure that Lagrange constraint is obeyed

  ## Points on rim
  P[Rset,1:2] <- R*cos(phi0)*P[Rset,1:2]/vecnorm(P[Rset,1:2])
  P[Rset,3]   <- R*sin(phi0)

  ## All points lie on sphere
  P[-Rset,] <- R*P[-Rset,]/vecnorm(P[-Rset,])

  ## Fixed point is set
  P[i0,] <- R*c(cos(phi0)*cos(lambda0),
                cos(phi0)*sin(lambda0),
                sin(phi0))

  return(P)
}

##' In the projection of points onto the sphere, some triangles maybe
##' flipped, i.e. in the wrong orientation.  This function determines
##' which triangles are flipped by computing the vector pointing to
##' the centre of each triangle and comparing this direction to vector
##' product of two sides of the triangle.
##'
##' @title Determine indices of triangles that are flipped
##' @param P Points in Cartesian coordinates
##' @param Trt Triangulation of points
##' @param R Radius of sphere
##' @return List containing:
##' \item{\code{flipped}}{Indices of in rows of \code{Trt} of flipped triangles.}
##' \item{\code{cents}}{Vectors of centres.}
##' \item{\code{areas}}{Areas of triangles.}
##' @author David Sterratt
flipped.triangles.cart <- function(P, Trt, R) {
  ## Plot any flipped triangles
  ## First find vertices and find centres and normals of the triangles
  P1 <- P[Trt[,1],]
  P2 <- P[Trt[,2],]
  P3 <- P[Trt[,3],]
  cents <- (P1 + P2 + P3)/3
  normals <- 0.5 * extprod3d(P2 - P1, P3 - P2)

  areas <- -0.5/R * dot(P1, extprod3d(P2, P3), 2)

  flipped <- (-dot(cents, normals, 2) < 0)
  return(list(flipped=flipped, cents=cents, areas=areas))
}

##' In the projection of points onto the sphere, some triangles maybe
##' flipped, i.e. in the wrong orientation.  This functions determines
##' which triangles are flipped by computing the vector pointing to
##' the centre of each triangle and comparing this direction to vector
##' product of two sides of the triangle.
##'
##' @title Determine indices of triangles that are flipped
##' @param Ps N-by-2 matrix with columns containing latitudes
##'   (\code{phi}) and longitudes (\code{lambda}) of N points
##' @param Trt Triangulation of points
##' @param R Radius of sphere
##' @return List containing:
##' \item{\code{flipped}}{Indices of in rows of \code{Trt} of flipped triangles.}
##' \item{\code{cents}}{Vectors of centres.}
##' \item{\code{areas}}{Areas of triangles.}
##' @author David Sterratt
flipped.triangles <- function(Ps, Trt, R=1) {
  return(flipped.triangles.cart(sphere.spherical.to.sphere.cart(Ps, R), Trt, R))
}
