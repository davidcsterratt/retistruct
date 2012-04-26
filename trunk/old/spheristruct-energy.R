## Now for the dreaded elastic error function....
E2 <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
              E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  cosp <- cos(phi)
  sinp <- sin(phi)
  cosl <- cos(lambda)
  sinl <- sin(lambda)

  x    <-  R*cosp*cosl
  y    <-  R*cosp*sinl
  z    <-  R*sinp

  E.E <- 0
  for (i in 1:nrow(Cu)) {
    j <- Cu[i,1]
    k <- Cu[i,2]
    l <- sqrt((x[j]-x[k])^2 + 
              (y[j]-y[k])^2 +
              (z[j]-z[k])^2)
    E.E <- E.E + 0.5 * (l - L[i])^2/L[i]
  }

  E.A <- 0
  if (E0.A) {
  }
  return(E.E + E0.A*E.A)
}


## Alternative, sequential version of the gradient
dE2 <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
                E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]
  
  cosp <- cos(phi)
  sinp <- sin(phi)
  cosl <- cos(lambda)
  sinl <- sin(lambda)

  x    <-  R*cosp*cosl
  y    <-  R*cosp*sinl
  z    <-  R*sinp

  dxdp <- -R*sinp*cosl
  dydp <- -R*sinp*sinl
  dzdp <-  R*cosp

  dxdl <- -R*cosp*sinl
  dydl <-  R*cosp*cosl
  dzdl <-  rep(0, N)

  dEdp <- rep(0, N)
  dEdl <- rep(0, N)

  if (verbose) print(paste(x[10], y[10], z[10]))
  
  for (i in 1:nrow(Cu)) {
    j <- Cu[i,1]
    k <- Cu[i,2]
    # print(paste("i", i))
    # print(j)
    # print(k)

    dx <- x[j]-x[k]
    dy <- y[j]-y[k]
    dz <- z[j]-z[k]
    
    l <- sqrt(dx^2 + dy^2 + dz^2)
    
    dlidpj <-  1/l*(dx*dxdp[j] + dy*dydp[j] + dz*dzdp[j])
    dlidpk <- -1/l*(dx*dxdp[k] + dy*dydp[k] + dz*dzdp[k])
    dlidlj <-  1/l*(dx*dxdl[j] + dy*dydl[j] + dz*dzdl[j])
    dlidlk <- -1/l*(dx*dxdl[k] + dy*dydl[k] + dz*dzdl[k])

    #print(dlidpj)
#    print(dlidpk)
#    print(dlidlj)
#    print(dxdl[j])
#    print(dydl[j])
#    print(dzdl[j])
  #  print(dlidlk)
    dll <- (l - L[i])/L[i]
      
    dEdp[j] <- dEdp[j] + dlidpj * dll
    dEdp[k] <- dEdp[k] + dlidpk * dll
    dEdl[j] <- dEdl[j] + dlidlj * dll
    dEdl[k] <- dEdl[k] + dlidlk * dll
  }
#  print(dEdl)
  dEAdp <- rep(0, N)
  dEAdl <- rep(0, N)


  if (E0.A) {
  }
  
  return(c(dEdp[-Rset]  + E0.A * dEAdp[-Rset],
           dEdl[-i0]    + E0.A * dEAdl[-i0]))
  ## return(c(dEdp, dEdl))
}

## Alternative, sequential version of the gradient
dE2c <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
                E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  dE <- .Call("spheristruct_E2",
              phi, lambda, Cu, L, R)

  dEdp <- dE[1:N]
  dEdl <- dE[N+(1:N)]
  return(c(dEdp[-Rset], dEdl[-i0]))
  #return(dE[-c(Rset, i0+Nphi)])
  ## return(dE)
}

## FIXME: E.ca and dE.ca should go if next round of testing works.

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
##' @param T Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of the sphere
##' @param Rset Indicies of points on the rim
##' @param i0 Index of fixed point on rim
##' @param phi0 Lattitude at which sphere curtailed
##' @param lambda0 Longitude of fixed points
##' @param Nphi Number of free values of \code{phi}
##' @param N Number of points in sphere
##' @param alpha Area scaling coefficient
##' @param x0 Area cutoff coefficient
##' @param verbose How much information to report
##' @return A single value, representing the energy of this particular
##' configuration
##' @author David Sterratt
E.ca <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi, N,
              alpha=1, x0, verbose=FALSE) {
  ## Extract phis and lambdas from parameter vector
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ## Compute elastic energy

  ## using the upper triagular part of the
  ## connectivity matrix Cu to extract coordinates of end-points of
  ## each edge
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]

  ## Compute lengths of edges
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)
  if (verbose==2) { print(l) }

  ## Compute spring energy
  E.E <- 0.5/sum(L)*sum((l - L)^2/L)
  if (verbose>=1) { print(E.E) }

  ## Compute areal penalty term if alpha is nonzero
  E.A <- 0
  if (alpha) {
    ## Compute areas by first finding cartesian coordinates of points
    P <- R * cbind(cos(phi)*cos(lambda),
                   cos(phi)*sin(lambda),
                   sin(phi))

    ## Find signed areas of all triangles
    a <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]), 2)

    ## Now compute area energy
    E.A <- sum(f(a/A, x0=x0))
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
##' @param p Parameter vector of \code{phi} and \code{lambda}
##' @param Cu The upper part of the connectivity matrix
##' @param C The connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param B Connectivity matrix
##' @param T Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of the sphere
##' @param Rset Indicies of points on the rim
##' @param i0 Index of fixed point on rim
##' @param phi0 Lattitude at which sphere curtailed
##' @param lambda0 Longitude of fixed points
##' @param Nphi Number of free values of \code{phi}
##' @param N Number of points in sphere
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cutoff coefficient
##' @param verbose How much information to report
##' @return A vector representing the derivative of the energy of this
##' particular configuration with respect to the parameter vector
##' @author David Sterratt
dE.ca <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi, N,
               alpha=1, x0, verbose=FALSE) {
  ## Extract phis and lambdas from parameter vector
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ## Compute derivative of elastic energy
  phii   <- phi[C[,1]]
  lambdai <- lambda[C[,1]]
  phij    <- phi[C[,2]]
  lambdaj <- lambda[C[,2]]

  ## x is the argument of the acos in the central angle
  u <- sin(phii)*sin(phij) + cos(phii)*cos(phij)*cos(lambdai-lambdaj)
  ## the central angle
  l <- R * acos(u)
  if (verbose==2) { print(l) }

  ## Compute general scaling factor
  fac <- 1/sum(L)*R*(l - c(L, L))/(c(L, L)*sqrt(1 - u^2))

  ## Now compute the derivatives using the B matrix
  dE.E.phii     <- B %*% (fac * (sin(phii)*cos(phij)*cos(lambdai-lambdaj)
                               - cos(phii)*sin(phij)))
  dE.E.dlambdai <- B %*% (fac * cos(phii)*cos(phij)*sin(lambdai-lambdaj))

  ## Compute the derivative of the area component if alpha is nonzero
  dE.A.dphi <- rep(0, N)
  dE.A.dlambda <- 0

  if (alpha) {
    ## Find cartesian coordinates of points
    P <- R * cbind(cos(phi)*cos(lambda),
                   cos(phi)*sin(lambda),
                   sin(phi))

    ## Here follows computation of the derivative - it's a bit
    ## complicated!
    
    ## Expand triangulation so that every point is the first point
    ## once. The number of points is effectively tripled.
    T <- rbind(T, T[,c(2,3,1)], T[,c(3,1,2)])
    A  <- c(A, A, A)

    ## Compute the derivative of area with respect to the first points
    dAdPt1 <- -0.5/R * extprod3d(P[T[,2],], P[T[,3],])

    ## Find areas of all triangles
    a <- dot(P[T[,1],], dAdPt1, 2)
    
    ## Now convert area derivative to energy derivative
    dEdPt1 <- fp(a/A, x0=x0)/A * dAdPt1

    ## Now map back onto phis and lambdas
    Pt1topi <- matrix(0, length(phi), nrow(T))
    for(m in 1:nrow(T)) {
      Pt1topi[T[m,1],m] <- 1
    }
    dEdpi <- Pt1topi %*% dEdPt1
    dpidphi <- R * cbind(-sin(phi) * cos(lambda),
                         -sin(phi) * sin(lambda),
                         cos(phi))
    dpidlambda <- R * cbind(-cos(phi) * sin(lambda),
                            cos(phi) * cos(lambda),
                            0)

    dE.A.dphi    <- rowSums(dEdpi * dpidphi)
    dE.A.dlambda <- rowSums(dEdpi * dpidlambda)
  }
  
  return(c(dE.E.phii[-Rset]      + alpha * dE.A.dphi[-Rset],
           dE.E.dlambdai[-i0]    + alpha * dE.A.dlambda[-i0]))
}

