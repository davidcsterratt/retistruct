## Return sequence of indicies in path between i and j, governed by
## pointer vector p
path <- function(i, j, g, h) {
  if (i == j) {
    return(i)
  } else {
    if (h[i] == i) {
      return(c(i, path(g[i], j, g, h)))
    } else {
      return(c(i, path(h[i], j, g, h)))
    }
  }
}

## Return sequence of indicies in path between i and j, governed by
## pointer vector p
path.length <- function(i, j, g, h, P) {
  if (any(is.na(c(i, j)))) {
    stop("i or j contains NA")
  }
  if (i == j) {
    return(0)
  } else {
    if (h[i] == i) {
      return(sqrt(sum((P[i,] - P[g[i],])^2)) + path.length(g[i], j, g, h, P))
    } else {
      return(path.length(h[i], j, g, h, P))
    }
  }
}

 ## order.Rset(Rset, gf, hf)
 ##
 ## It is nice to create Rset as an ordered set
 order.Rset <- function(Rset, gf, hf) {
   ## To to this, join the path from the first two members of the set.
   R12 <- path(Rset[1], Rset[2], gf, hf)
   R21 <- path(Rset[2], Rset[1], gf, hf)
   Rset <- c(R12[-1], R21[-1])
   return(Rset)
 }

## Convert a matrix containing on each line the indicies of the points
## forming a segment, and convert this to two sets of ordered pointers
segments2pointers <- function(S) {
  g <- c()
  j <- 1                                # Row of S
  k <- 1                                # Column of S
  while(nrow(S) > 0) {
    i <- S[j,3-k]                       # i is index of the next point
    g[S[j,k]] <- i                      # Set the pointer to i
    S <- S[-j,,drop=FALSE]              # We have used this row of S
    if (nrow(S) == 0) {
      return(g)
    }
    j <- which(S[,1] == i)            # Is i in the first column of S?
    if (length(j) > 1) {
      stop("The segment list is not valid as it contains an element more than twice.")
    }
    if (length(j)) {              # If so, set the current column to 1
      k <- 1
    } else {
      j <- which(S[,2] == i) # Otherwise, look for i in the second column
      k <- 2
      if (!length(j)) {
        stop(paste("No matching index for point", i, "in S."))
        return(NULL)
      }
    }
  }
  return(g)
}

## Convert a set of ordered pointers to a matrix containing on each
## line the indicies of the points forming a segment
pointers2segments <- function(g) {
  S1 <- which(!is.na(g))
  S2 <- g[S1]
  return(cbind(S1, S2))
}

## Merge points and edges before optimising the elastic energy
## The information to be merged is contained in the
## list t, and the index of the landmark i0
##
## The information includes 
## h    - the point correspondence mapping
## Rset - the set of points on the rim
## i0   - the index of the landmark
## T    - the triangulation
## Cu   - the edge matrix
## L    - the edge lengths
## P    - the point locations
##
## The function returns merged and transformed versions of all these
## objects (all suffixed with t), as well as a matrix Bt, which maps a
## binary vector representation of edge indicies onto a binary vector
## representation of the indicies of the points linked by the edge
merge.points.edges <- function(t) {
  h <- t$h
  T <- t$T
  Cu <- t$Cu
  L <- t$L
  P <- t$P
  gf <- t$gf
  
  ## Form the mapping from a new set of consecutive indicies
  ## the existing indicies onto the existing indicies
  u <- unique(h)

  ## Transform the point set into the new indicies
  Pt  <- P[u,]

  ## Transform the point correspondance mapping to the new index space  
  ht <- c()
  for (i in 1:length(h)) {
    ht[i] <- which(u == h[i])
  }

  ## DOESN'T WORK
  ## Form the inverse mapping from the existing indicies to the new
  ## set of consecutive indicies
  ## uinv <- c()
  ## uinv[u] <- 1:length(u)
  ## ht <- uinv[h[u]]

  ## Transform the triangulation to the new index space
  Tt  <- matrix(ht[T], ncol=3)

  ## Tansform the forward pointer into the new indicies
  gft <- ht[gf]

  ## Determine H, the mapping from edges onto corresponding edges
  Cut <- matrix(ht[Cu], ncol=2)
  Cut <- t(apply(Cut, 1, sort))
  M <- nrow(Cut)
  H <- rep(0, M)
  for (i in 1:M) {
    if (!H[i]) {
      H[i] <- i
      for (j in i:M) {
        if (identical(Cut[i,], Cut[j,])) {
          H[j] <- i
        }
      }
    }
  }

  ## Form the mapping from a new set of consecutive edge indicies
  ## onto the existing edge indicies
  U <- unique(H)

  ## Transform the edge set into the new indicies
  Cut <- Cut[U,]

  ## Transform the edge correspondance mapping to the new index space  
  Ht <- c()
  for (i in 1:length(H)) {
    Ht[i] <- which(U == H[i])
  }

  ## Create the lengths of the merge edges by averaging
  Lt <- c()
  for (k in 1:length(U)) {
    is <- which(Ht == k)
    ## if (length(is)>1) {
    ##   print(L[is])
    ## }
    Lt[k] <- mean(L[is])
  }

  ## Transform the rim set
  Rset <- order.Rset(t$Rset, t$gf, t$hf)
  Rsett <- unique(ht[Rset])
  i0t <- ht[t$i0]

  ## Create the symmetric connection set
  Ct <- rbind(Cut, Cut[,2:1])

  ## Matrix to map line segments onto the points they link
  Bt <- matrix(0, nrow(Pt), nrow(Ct))
  for (i in 1:nrow(Ct)) {
    Bt[Ct[i,1],i] <- 1
  }

  m <- merge(list(Pt=Pt, Tt=Tt, Ct=Ct, Cut=Cut, Bt=Bt, Lt=Lt, ht=ht,
                  Rset=Rset, Rsett=Rsett, i0t=i0t, P=P, H=H, Ht=Ht), t)
  class(m) <- class(t)
  return(m)
}

## Stretch the mesh in the flat retina to a circular outline
##
## Arguments:
## Cu - connection matrix
## L  - lengths
## i.fix - indicies of fixed points
## P.fix - coordinates of fixed points
##
## Returns:
## New matrix of point locations
stretch.mesh <- function(Cu, L, i.fix, P.fix) {
  N <- max(Cu)
  M <- length(L)
  C <- matrix(0, 2*N, 2*N)
  for (i in 1:M) {
    C[2*Cu[i,1]-1:0,2*Cu[i,2]-1:0] <- diag(2) / L[i]
    C[2*Cu[i,2]-1:0,2*Cu[i,1]-1:0] <- diag(2) / L[i]
  }

  ## FIXME: This debugging output should probably go once the "Lapack
  ## routine dgesv: system is exactly singular" bug has been
  ## vanquished
  message(paste("det(C) =", det(C)))
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
  message(paste("det(A) =", det(A)))
  B <- C[-ind,  ind]
  P <- matrix(t(P.fix), ncol(B), 1)
  D <- diag(apply(cbind(A, 2*B), 1, sum))
  message(paste("det(D) =", det(D)))
  if (is.infinite(det(D))) stop ("det(D) is infinite")
  Q <- 2 * solve(D - A) %*% B %*% P
  Q <- matrix(Q, nrow(Q)/2, 2, byrow=TRUE)
  R <- matrix(0, nrow(Q) + length(i.fix), 2)
  R[i.fix,] <- P.fix
  R[-i.fix,] <- Q
  return(R)
}

## Project mesh points in the flat outline onto a sphere
##
## The information to be merged is contained in the list r and the
## lattitude at which the sphere is cut off is given in phi0
##
## The information includes:
## Pt     - the mesh point coordinates
## Rsett  - the set of points on the rim
## A.tot  - the area of the flat outline
##
## The function returns a list with the following members:
## phi    - lattitude of mesh points
## lambda - longitude of mesh points
## R      - radius of sphere
## phi0   - lattitude at which sphere is cut off (from input)
project.to.sphere <- function(r) {
  Pt <- r$Pt
  Rsett <- r$Rsett
  i0t <- r$i0t
  A.tot <- r$A.tot
  Cut <- r$Cut
  Lt <- r$Lt
  phi0 <- r$phi0
  lambda0 <- r$lambda0
    
  Nt <- nrow(Pt)
  Nphi <- Nt - length(Rsett)

  ## From this we can infer what the radius should be from the formula
  ## for the area of a sphere which is cut off at a lattitude of phi0
  ## area = 2 * PI * R^2 * (sin(phi0)+1)
  R <- sqrt(A.tot/(2*pi*(sin(phi0)+1)))
  
  ## Stretch mesh points to circle
  Ps <- stretch.mesh(Cut, Lt, Rsett, circle(length(Rsett)))
  x <- Ps[,1]
  y <- Ps[,2]
  phi <- -pi/2 + sqrt(x^2 + y^2)*(phi0+pi/2)
  phi[Rsett] <- phi0
  lambda <- atan2(y, x)
  lambda <- lambda - lambda[i0t] + lambda0

  p <- merge(list(phi=phi, lambda=lambda, R=R,
                  phi0=phi0, lambda0=lambda0, Ps=Ps),
             r)
  class(p) <- unique(c("reconstructedOutline", class(r)))
  return(p)
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
compute.areas <- function(phi, lambda, T, R) {
  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Find areas of all triangles
  areas <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))
  
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
##' @param x0 The cutoff parameter. Above this value the function is zero.
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
##' @param x0 The cutoff parameter. Above this value the function is zero.
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
E <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi, N,
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
    a <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))

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
dE <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi, N,
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
    a <- dot(P[T[,1],], dAdPt1)
    
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

##' In the projection of points onto the sphere, some triangles maybe
##' flipped, i.e. in the wrong orientation.  This functions determines
##' which triangles are flipped by computing the vector pointing to
##' the centre of each triangle and comparing this direction to vector
##' product of two sides of the triangle.
##'
##' @title Determine indicies of triangles that are flipped
##' @param phi Vector of lattitudes of points
##' @param lambda Vector of longitudes of points
##' @param Tt Triangulation of points
##' @param R Radius of sphere
##' @return List containing:
##' \item{\code{flipped}}{Indicies of in rows of \code{Tt} of flipped triangles.}
##' \item{\code{cents}}{Vectors of centres.}
##' \item{\code{areas}}{Areas of triangles.}
##' @author David Sterratt
flipped.triangles <- function(phi, lambda, Tt, R) {
  P <- sphere.spherical.to.sphere.cart(phi, lambda, R)
  ## Plot any flipped triangles
  ## First find verticies and find centres and normals of the triangles
  P1 <- P[Tt[,1],]
  P2 <- P[Tt[,2],]
  P3 <- P[Tt[,3],]
  cents <- (P1 + P2 + P3)/3
  normals <- 0.5 * extprod3d(P2 - P1, P3 - P2)

  areas <- -0.5/R * dot(P1, extprod3d(P2, P3))
  
  flipped <- (-dot(cents, normals) < 0)
  return(list(flipped=flipped, cents=cents, areas=areas))
}
  
##' Optimise the mapping from the flat outline to the sphere
##'
##' @title Optimise mapping
##' @param r Reconstruction object
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cutoff coefficient
##' @param method Method to pass to \code{optim}
##' @param plot.3d If \code{TRUE} make a 3D plot in an RGL window
##' @param dev.grid Device handle for plotting grid to
##' @param dev.polar Device handle for plotting ploar plot to
##' @return Reconstruction object
##' @author David Sterratt
optimise.mapping <- function(r, alpha=4, x0=0.5, method="BFGS",
                             plot.3d=FALSE, dev.grid=NA, dev.polar=NA) {
  phi <- r$phi
  lambda <- r$lambda
  R <- r$R
  phi0 <- r$phi0
  lambda0 <- r$lambda0
  Tt <- r$Tt
  A <- r$A
  Cut <- r$Cut
  Ct <- r$Ct
  Pt <- r$Pt
  Lt <- r$Lt
  Bt <- r$Bt
  Rsett <- r$Rsett
  i0t <- r$i0t
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)
  
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi[-Rsett], lambda[-i0t])
  opt$conv <- 1
  
  while (opt$conv) {
    ## Optimise
    opt <- optim(opt$p, E, gr=dE,
                 method=method,
                 T=Tt, A=A, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R,
                 alpha=alpha,  N=Nt, x0=x0,
                 Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi,
                 verbose=FALSE)

    ## Report
    E.tot <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               alpha=alpha,  N=Nt, x0=x0,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
    E.l <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               alpha=0,  N=Nt, x0=x0,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)

    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)
    message(sprintf("E = %8.5f | E_L = %8.5f | E_A = %8.5f | %3d flippped triangles", E.tot, E.l, E.tot - E.l,  nflip))
    if (nflip) {
      print(data.frame(rbind(id=which(ft$flipped),
                             A=A[ft$flipped],
                             a=ft$areas[ft$flipped])))
    }
    
    ## Decode p vector
    phi          <- rep(phi0, Nt)
    phi[-Rsett]  <- opt$p[1:Nphi]
    lambda       <- rep(lambda0, Nt)
    lambda[-i0t] <- opt$p[Nphi+1:(Nt-1)]

    ## Plot
    if (plot.3d) {
      plot.sphere.spherical(phi, lambda, R, Tt, Rsett)
      plot.outline.spherical(phi, lambda, R, r$gb, r$ht)
    }

    if (!is.na(dev.grid)) {
      dev.set(dev.grid)
      plot.flat(r, grid=TRUE, 
                datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      r$phi <- phi
      r$lambda <- lambda
      plot.polar(r)
    }
  }
  o <- merge(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
                  E.tot=E.tot, E.l=E.l),
             r)
  class(o) <- class(r)
  return(o)
}

## Function to plot the fractional change in length of connections 
compute.strain <- function(r) {
  ## Original lengths in flattened outline is a vector with
  ## M elements, the number of rows of Cu
  L <- r$L
  ## New lenghts in reconstructed object is a vector wtih Mt < M
  ## elements, the number of rows of Cut
  ls <- compute.lengths(r$phi, r$lambda, r$Cut, r$R)
  ## For each connection in the flattened object, we want the length of
  ## the corresponding connection in the reconstructed object
  ## The mapping Ht achieves this@
  l <- ls[r$Ht]
  strain <- l/L
  logstrainnorm <- -log(strain)/max(abs(log(strain)))
  return(data.frame(L=L, l=l, strain=strain, logstrainnorm=logstrainnorm))
}

##' Reconstruct outline into spherical surface. Reconstruction
##' proceeds in a number of stages:
##'
##' \enumerate{
##' 
##' \item The flat object is triangulated with at least \code{n}
##' triangles. This can introduce new vertices in the rim. 
##'
##' \item The triangulated object is stitched.
##'
##' \item The stitched object is triangulated again, but this time it
##' is not permitted to add extra vertices to the rim.
##'
##' \item The corresponding points determined by the stitching process
##' are merged to form a new set of merged points and a new
##' triangulation.
##'
##' \item The merged points are projected roughly to a sphere.
##'
##' \item The locations of the points on the sphere are moved so as to
##' minimise the energy function.
##' }
##'
##' @title Reconstruct outline into spherical surface
##' @param o Outline list, containing the following information:\describe{
##' \item{\code{P}}{outline points as N-by-2 matrix}
##' \item{\code{V0}}{indicies of the apex of each tear}
##' \item{\code{VF}}{indicies of the forward vertex of each tear}
##' \item{\code{VB}}{indicies of the backward vertex of each tear}
##' \item{\code{i0}}{index of the landmark on the rim}
##' \item{\code{phi0}}{lattitude of rim of partial sphere}
##' \item{\code{lambda0}}{longitude of landmark on rim}
##' }
##' @param n Number of points in triangulation.
##' @param alpha Area scaling coefficient
##' @param x0 Area cutoff coefficient
##' @param report Function used to report progress.
##' @param plot.3d Whether to show 3D picture during optimisation.
##' @param dev.grid Device to plot grid onto. Value of \code{NA} (default)
##' means no plotting.
##' @param dev.polar Device to plot polar plot onto. Value of NA
##' (default) means no plotting.
##' @return \code{reconstructedOutline} object containing the input
##' information and the following modified and extra information:
##' \item{\code{P}}{New set of points in flattened object}
##' \item{\code{gf}}{New set of forward pointers in flattened object}
##' \item{\code{gb}}{New set of backward pointers in flattened object}
##' \item{\code{phi}}{lattitude of new points on sphere}
##' \item{\code{lambda}}{longitude of new points on sphere}
##' \item{\code{Tt}}{New triangulation}
##' @author David Sterratt
ReconstructedOutline <- function(o, 
                                 n=500, alpha=4, x0=0.5,
                                 report=print,
                                 plot.3d=FALSE, dev.grid=NA, dev.polar=NA) {
  ## Clear polar plot, if it's required
  if (!is.na(dev.polar)) {
    dev.set(dev.polar)
    plot.polar(o$phi0)
  }
  
  report("Triangulating...")
  t <- triangulate.outline(o, n=n)
  if (!is.na(dev.grid)) {
    dev.set(dev.grid)
    plot.flat(t)
  }
    
  report("Stitching...")
  s <- stitch.outline(t)
  if (!is.na(dev.grid)) {
    dev.set(dev.grid)
    plot.flat(s, datapoints=FALSE)
  }

  report("Triangulating...")  
  r <- triangulate.outline(s, n=n,
                           suppress.external.steiner=TRUE)
  
  if (!is.na(dev.grid)) {
    dev.set(dev.grid)
    plot.flat(r, datapoints=FALSE)
  }

  report("Merging points...")
  r <- merge.points.edges(r)
  
  report("Projecting to sphere...")
  r <- project.to.sphere(r)
  
  if (!is.na(dev.grid)) {
    ## Plot of initial gridlines
    dev.set(dev.grid)
      plot.flat(r, grid=TRUE, 
                datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE)
    
    ## Initial plot in 3D space
    plot.sphere.spherical(r$phi, r$lambda, r$R, r$Tt, r$Rsett)
    plot.outline.spherical(r$phi, r$lambda, r$R, r$gb, r$ht)
  }

  report("Optimising mapping...")
  r <- optimise.mapping(r, alpha=alpha, x0=x0,
                        plot.3d=plot.3d,
                        dev.grid=dev.grid, dev.polar=dev.polar)
  
  report(paste("Mapping optimised. Error:", format(r$opt$value, 5),
               ";", r$nflip, "flipped triangles."))
  class(r) <- unique(c("reconstructedOutline", class(r)))
  return(r)
}
