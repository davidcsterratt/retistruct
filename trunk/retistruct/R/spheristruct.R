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
  Rsett <- unique(ht[t$Rset])
  i0t <- ht[t$i0]

  ## Create the symmetric connection set
  Ct <- rbind(Cut, Cut[,2:1])

  ## Matrix to map line segments onto the points they link
  Bt <- matrix(0, nrow(Pt), nrow(Ct))
  for (i in 1:nrow(Ct)) {
    Bt[Ct[i,1],i] <- 1
  }
  
  return(list(Pt=Pt, Tt=Tt, Ct=Ct, Cut=Cut, Bt=Bt, Lt=Lt, ht=ht,
              Rsett=Rsett, i0t=i0t, P=P, H=H, Ht=Ht))
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
  print(det(C))
  dupC <- duplicated(C)
  if (any(dupC)) {
    i <- which(dupC)
    Ci <- C[i,]
    dups <- which(apply(C, 1, function(x) {identical(x, Ci)}))
    print(dups)
    print(Ci)
    for (d in dups) {
      print(d)
      print(which(C[d,]==1))
    }
  }
  
  ind <- as.vector(rbind(2*i.fix-1, 2*i.fix))
  A <- C[-ind, -ind]
  print(det(A))
  B <- C[-ind,  ind]
  P <- matrix(t(P.fix), ncol(B), 1)
  D <- diag(apply(cbind(A, 2*B), 1, sum))
  print(det(D))

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

  return(list(phi=phi, lambda=lambda, R=R, phi0=phi0, lambda0=lambda0, Ps=Ps))
}

##
## Energy/error functions
## 

## Formula for central angle
central.angle <- function(phi1, lambda1, phi2, lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)))
}

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

f <- function(x, x0=0.1, k=1) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -k*(x[c1] - x0/2)
  y[c2] <- k/2/x0*(x0 - x[c2])^2
  y[c3] <- 0

  return(y)
}

fp <- function(x, x0=0.1, k=1) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -k
  y[c2] <- -k/x0*(x0 - x[c2])
  y[c3] <- 0

  return(y)
}

## Now for the dreaded elastic error function....
E <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
              E0.A=1, k.A=1, x0=0.1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ##
  ## Compute derivative of elastic energy
  ##
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)
  if (verbose==2) {
    print(l)
  }
  E.E <- 0.5 * sum((l - L)^2/L)
  ## E.E <- 0.5 * sum((l - L)^2)
  ## E.E <- 0.5*sum((l)^2/L)
  if (verbose>=1) {
    print(E.E)
  }

  E.A <- 0
  if (E0.A) {
    ##
    ## Compute areas
    ##
    P <- R * cbind(cos(phi)*cos(lambda),
                   cos(phi)*sin(lambda),
                   sin(phi))

    ## Find areas of all triangles
    a <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))
    ##E.A <- 0.5 * sum((areas - A)^2/A)
    ## E.A <- sum(exp(-k.A*areas/A))
    E.A <- sum(sqrt(A)*f(a/A, x0=x0))
  }
  
  return(E.E + E0.A*E.A)
}

## ... and the even more dreaded gradient of the elastic error
dE <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
               E0.A=1, k.A=1, x0=0.1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ##
  ## Compute derivative of elastic energy
  ##
  phii   <- phi[C[,1]]
  lambdai <- lambda[C[,1]]
  phij    <- phi[C[,2]]
  lambdaj <- lambda[C[,2]]
  ## x is the argument of the acos in the central angle
  u <- sin(phii)*sin(phij) + cos(phii)*cos(phij)*cos(lambdai-lambdaj)
  ## the central angle
  l <- R * acos(u)
  if (verbose==2) {
    print(l)
  }
  fac <- R * (l - c(L, L))/(c(L, L) * sqrt(1-u^2))
  ## fac <- R * (l - c(L, L))/(sqrt(1-u^2))
  ## fac <- R * (l)/(c(L, L) * sqrt(1-u^2))
  dE.E.phii     <- B %*% (fac * (sin(phii)*cos(phij)*cos(lambdai-lambdaj)
                               - cos(phii)*sin(phij)))
  dE.E.dlambdai <- B %*% (fac * cos(phii)*cos(phij)*sin(lambdai-lambdaj))

  dE.A.dphi <- rep(0, N)
  dE.A.dlambda <- 0
  if (E0.A) {
    ##
    ## Compute derivative of areas
    ##
    P <- R * cbind(cos(phi)*cos(lambda),
                   cos(phi)*sin(lambda),
                   sin(phi))

    ## expand triangulation
    T <- rbind(T, T[,c(2,3,1)], T[,c(3,1,2)])
    A  <- c(A, A, A)

    ## Slow way of computing gradient
    dAdPt1 <- -0.5/R * extprod3d(P[T[,2],], P[T[,3],])

    ## Find areas of all triangles
    a <- dot(P[T[,1],], dAdPt1)

##     dEdPt1 <- (areas - A)/A * dAdPt1
##    dEdPt1 <- -k.A/A*exp(-k.A*areas/A) * dAdPt1
    dEdPt1 <- fp(a/A, x0=x0)/sqrt(A) * dAdPt1
    
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
  
  return(c(dE.E.phii[-Rset]      + E0.A * dE.A.dphi[-Rset],
           dE.E.dlambdai[-i0]    + E0.A * dE.A.dlambda[-i0]))
}

## flipped.triangles - Determine indicies of triangles that are flipped
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
  print(length(areas))
  
  flipped <- (-dot(cents, normals) < 0)
  return(list(flipped=flipped, cents=cents, areas=areas))
}
  
## Grand optimisation function
optimise.mapping <- function(r, E0.A=10, k.A=1, x0=0.5, method="BFGS",
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

  ## iE <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)

  ## idE <- dE(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
  ## idE2 <- dE2(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi, verbose=TRUE)

  ## idE2c <- dE2c(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
  ## print(iE)
  
  ## print(idE[1:10])
  ## print("R dE2")
  ## print(idE2[1:10])
  ## print("C dE2")
  ## print(idE2c[1:10])
  ## print("Ratio")
  ## print(idE2c/idE2)
  
  while (opt$conv) {
    ## Optimise
    opt <- optim(opt$p, E, gr=dE,
                 method=method,
                 T=Tt, A=A, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R,
                 E0.A=E0.A, k.A=k.A, N=Nt, x0=x0,
                 Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi,
                 verbose=FALSE)

    ## Report
    E.tot <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               E0.A=E0.A, k.A=k.A, N=Nt, x0=x0,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
    E.l <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               E0.A=0, k.A=k.A, N=Nt, x0=x0,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
    print(paste("Total error:", E.tot, opt$value, "; Length error:", E.l))
    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)
    print(paste(nflip, "flipped triangles:"))
    print(which(ft$flipped))
    print("Areas")
    print(ft$areas[ft$flipped])
    print(A[ft$flipped])
    print(length(A))
    
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
      with(r, plot.outline.flat(P, gb))
      plot.gridlines.flat(r$P, r$T, phi, lambda, Tt, phi0*180/pi)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      plot.polar(phi0 * 180/pi)
      r$phi <- phi
      r$lambda <- lambda
      plot.outline.polar(r)
    }
  }
  return(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
              E.tot=E.tot, E.l=E.l))
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
                         n=500,
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
  r$Rset <- order.Rset(r$Rset, r$gf, r$hf)

  report("Merging points...")
  m <- merge.points.edges(r)
  r <- merge(m, r)

  report("Projecting to sphere...")
  p <- project.to.sphere(r)
  r <- merge(p, r)
  
  if (!is.na(dev.grid)) {
    ## Plot of initial gridlines
    dev.set(dev.grid)
    with(r, plot.outline.flat(P, gb))
    with(r, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0*180/pi))
    
    ## Initial plot in 3D space
    plot.sphere.spherical(r$phi, r$lambda, r$R, r$Tt, r$Rsett)
    plot.outline.spherical(r$phi, r$lambda, r$R, r$gb, r$ht)
  }

  report("Optimising mapping...")
  ## This pass is experimental - ideally it wouldn't be here, but until
  ## we can fix problems with fixed triangles, it must stay.
  ##o <- optimise.mapping(r, E0.A=exp(1), k.A=2,
  ##                        plot.3d=plot.3d,
  ## dev.grid=dev.grid, dev.polar=dev.polar)
  ## r <- merge.lists(r, o)

  ## This pass is original
  ## o <- list(nflip=1)
  E0.A <- 32
  ##while(o$nflip>0) {
  o <- optimise.mapping(r, E0.A=E0.A,
                        plot.3d=plot.3d,
                        dev.grid=dev.grid, dev.polar=dev.polar)
  r <- merge(o, r)
  ## o <- fem.optimise.mapping(r, nu=0.45,
  ##                           plot.3d=plot.3d,
  ##                           dev.grid=dev.grid, dev.polar=dev.polar)
  ## r <- merge.lists(r, o)
  ##  E0.A <- E0.A * 2;
  ##  }
  
  report(paste("Mapping optimised. Error:", format(r$opt$value,5),
               ";", r$nflip, "flipped triangles."))
  class(r) <- c("reconstructedOutline", class(r))
  return(r)
}

## Infer coordinates of datapoints
## Arguments:
## r    - object returned by fold.outline
## Ds   - list of sets of datapoints to plot
## Ss   - list of sets of landmarks to plot
##
## Returns:
## New object r, with new objects attached:
## 
## Dsb  - Datapoints in barycentric coordinates
## Dsc  - Datapoints on reconstructed sphere in cartesian coordinates
## Dss  - Datapoints on reconstructed sphere in spherical coordinates
## Ssb  - Landmarks in barycentric coordinates
## Ssc  - Landmarks on reconstructed sphere in cartesian coordinates
## Sss  - Landmarks on reconstructed sphere in spherical coordinates
##
infer.datapoint.landmark.coordinates <- function(r, report=print) {
  report("Inferring coordinates of datapoints")
  Dsb <- list() # Datapoints in barycentric coordinates
  Dsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Dss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ds)) {
    for (name in names(r$Ds)) {
      Dsb[[name]] <- tsearchn(r$P, r$T, r$Ds[[name]])
      Dsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Dsb[[name]])
      Dss[[name]] <- sphere.cart.to.sphere.spherical(Dsc[[name]], r$R)
    }
  }

  report("Inferring coordinates of landmarks")
  Ssb <- list() # Landmarks in barycentric coordinates
  Ssc <- list() # Landmarks on reconstructed sphere in cartesian coordinates
  Sss <- list() # Landmarks on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ss)) {
    for (name in names(r$Ss)) {
      Ssb[[name]] <- with(r, tsearchn(P, T, r$Ss[[name]]))
      Ssc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Ssb[[name]])
      Sss[[name]] <- sphere.cart.to.sphere.spherical(Ssc[[name]], r$R)
    }
  }

  r <- merge(list(Dsb=Dsb, Dsc=Dsc, Dss=Dss,
                  Ssb=Ssb, Ssc=Ssc, Sss=Sss), r)
}

## Infer coordinates of tears
## Arguments:
## r    - object returned by fold.outline
##
## Returns:
## New object r, with new object attached:
## 
## Tss  - Tear coordinates on reconstructed sphere in spherical coordinates
##
infer.tear.coordinates <- function(r,
                                   report=print) {
  report("Inferring coordinates of tears")
  Tss <- list() # Landmarks on reconstructed sphere in spherical coordinates
  if (!is.null(r$TFset)) {
    for (TF in r$TFset) {
      ## Convert indicies to the spherical frame of reference
      j <- r$ht[TF]
      Tss <- c(Tss, list(cbind(phi=r$phi[j], lambda=r$lambda[j])))
    }
  }
  r <- merge(list(Tss=Tss), r)
}
