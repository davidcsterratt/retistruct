## Now for the dreaded elastic error function....
E.alabama <- function(p, Cu, C, L, B, T, A, 
                      E0.A=0, x0=0.1, verbose=FALSE, ...) {
  ## Extract points from parameters
  P <- matrix(p, length(p)/3, 3)

  ##
  ## Compute the elastic energy
  ##
  ## Use the upper triagular part of the connectivity matrix Cu
  P1    <- P[Cu[,1],]
  P2    <- P[Cu[,2],]
  l <- sqrt(rowSums(P1 - P2)^2)
  E.E <- 0.5 * sum((l - L)^2/L)
  if (verbose>=1) {
    print(E.E)
  }

  if (E0.A) {
    ## Find areas of all triangles
    a <- -0.5 * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))
    E.A <- sum(sqrt(A)*f(a/A, x0=x0))
  } else {
    E.A <- 0
  }
  
  return(E.E + E0.A*E.A)
}

## ... and the even more dreaded gradient of the elastic error
dE.alabama <- function(p, Cu, C, L, B, T, A, R, P, E0.A=1, k.A=1, x0=0.1, N, verbose=FALSE, ...) {
  ## Extract points from parameters
  P <- matrix(p, length(p)/3, 3)

  ##
  ## Compute derivative of elastic energy
  ##
  Pi   <- P[C[,1],]
  Pj   <- P[C[,2],]

  l <- sqrt(rowSums(Pi - Pj)^2)
  if (verbose==2) {
    print(l)
  }
  fac <- (l - c(L, L))/(c(L, L) * l)
  dE.E     <- B %*% (fac * (Pi - Pj))

  dE.A.dphi <- rep(0, N)
  dE.A.dlambda <- 0
  if (0) {
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
  
  return(as.vector(dE.E))
}


## heq.alabama <- function(p,  P0, Rset, i0, R, ...) {
##   ## Extract points from parameters
##   P <- matrix(p, length(p)/3, 3)
##   heq <- c(P[i0,1] - P0[1],
##            P[i0,2] - P0[2],
##            P[Rset, 3] - P0[3],
##            rowSums(P^2) - R^2)
##   return(heq)
## }

heq.alabama <- function(p,  P0, Rset, i0, R, ...) {
  ## Extract points from parameters
  P <- matrix(p, length(p)/3, 3)
  Rsetp <- setdiff(Rset, i0)
  heq <- c(P[i0, 1] - P0[1],
           P[i0, 2] - P0[2],
           P[i0, 3] - P0[3],
           rowSums(P[Rsetp, 1:2]^2) - P0[1]^2 -P0[2]^2,
           P[Rsetp, 3] - P0[3],
           rowSums(P[-Rset,]^2) - R^2)
  return(heq)
}

heq.jac.alabama <- function(p,  P0, Rset, i0, R, ...) {
  ## Extract points from parameters
  N <- length(p)
  NR <- length(Rset)
  M <- 3 + 2*(NR - 1) + (N - NR)

  P <- matrix(p, N, 3)

  jac <- matrix(0, M, N)
  jac[1:3, 1:3] <- diag(1 , 3)
  jac[4:NR
  return(jac)
}


## Grand optimisation function
optimise.mapping.alabama <- function(r, E0.A=10, k.A=1, x0=0.5, method="BFGS",
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
  P <- sphere.spherical.to.sphere.cart(phi, lambda, R)
  opt <- list()
  opt$p <- as.vector(P)
  opt$conv <- 1
  ## Location of fixed point
  P0 <- sphere.spherical.to.sphere.cart(phi0, lambda0, R)
  
  while (opt$conv) {
    ## Optimise
    opt <- auglag(opt$p, E.alabama, gr=dE.alabama,
                  heq=heq.alabama, #heq.jac=heq.jac.alabama,
                  T=Tt, A=A, Cu=Cut, C=Ct, L=Lt, B=Bt,
                  E0.A=E0.A, k.A=k.A, x0=x0,
                  Rset=Rsett, i0=i0t, P0=P0, R=R, 
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
    Ps <- sphere.cart.to.sphere.spherical(matrix(opt$p, 3, Nt), R)
    phi <- Ps[,"phi"]
    lambda <- Ps[,"lambda"]

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
      polarplot(phi0 * 180/pi)
      r$phi <- phi
      r$lambda <- lambda
      plot.outline.polar(r)
    }
  }
  return(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
              E.tot=E.tot, E.l=E.l))
}
