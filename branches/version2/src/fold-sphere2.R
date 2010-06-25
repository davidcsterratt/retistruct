source("tsearch.R")
source("triangulate.R")                 # for tri.area and tri.area.signed
source("geometry.R")                    # for dot() and norm()
source("misc.R")                        # For Mod
require("rgl")
require("plotrix")                      # For polar plots

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

## Stitch together tears in an outline
##
## Input arguments:
## P     - the coordinates of points in a mesh, including those the outline
## gf    - the forward pointer list
## gb    - the backward pointer list
## T     - the tear matrix
## i0    - the index of the landmark; this needs to be in the rim
##
## The function returns a list contatining:
## Rset  - the set of points on the rim
## i0    - the index of the landmark
## P     - a new set of meshpoints
## A     - indicies of the apex of each tear
## VF    - indicies of the forward vertex of each tear
## VB    - indicies of the backward vertex of each tear
## TFset - list containing indicies of points in each foward tear
## TBset - list containing indicies of points in each backward tear
## gf    - new forward pointer list
## gb    - new backward pointer list
## h     - correspondence mapping
## hf    - correspondence mapping in forward direction for
##         points on boundary
## hb    - correspondence mapping in backward direction for
##         points on boundary
##
stitch.outline <- function(P, gf, gb, T, i0=NA) {
  ## Extract information from tear matrix
  A  <- T[,1]                            # apicies of tears
  VB <- T[,2]                           # forward verticies
  VF <- T[,3]                           # backward verticies

  ## Create initial sets of correspondances
  N <- nrow(P)                          # Number of points
  hf <- 1:N
  hb <- 1:N
  hf[VB] <- VF
  hb[VF] <- VB
  h <- hf

  ## Initialise the set of points in the rim
  ## We don't assume that P is the entire set of points; instead
  ## get this information from the pointer list.
  Rset <- na.omit(gf)
  
  ## Create lists of forward and backward tears
  TFset <- list()
  TBset <- list()
  
  ## Iterate through the tears to create tear sets and rim set
  for (j in 1:nrow(T)) {
    ## Create sets of points for each tear and remove these points from
    ## the rim set
    ##print(paste("Forward tear", j))
    TFset[[j]] <- Mod(path(A[j], VF[j], gf, hf), N)
    ##print(paste("Backward tear", j))
    TBset[[j]] <- Mod(path(A[j], VB[j], gb, hb), N)
    Rset <- setdiff(Rset, setdiff(TFset[[j]], VF[j]))
    Rset <- setdiff(Rset, setdiff(TBset[[j]], VB[j]))
  }

  ## If not set, set the landmark marker index. Otherwise
  ## check it
  if (is.na(i0)) {
    i0 <- Rset[1]
  } else {
    if (!(i0 %in% Rset)) {
      return(NULL)
    }
  }
  
  ## Iterate through tears to insert new points
  for (j in 1:nrow(T)) {
    ## Compute the total path length along each side of the tear
    Sf <- path.length(A[j], VF[j], gf, hf, P)
    Sb <- path.length(A[j], VB[j], gb, hb, P)
    ## print(paste("Sf", Sf))
    ## print(paste("Sb", Sb))

    ## For each point in the forward path, create one in the backwards
    ## path at the same fractional location
    for (i in setdiff(TFset[[j]], c(A[j], VF[j]))) {
      sf <- path.length(A[j], i, gf, hf, P)
      ## print(paste("sf", sf/Sf))
      ## print(TBset[[j]])
      for (k in TBset[[j]]) {
        sb <- path.length(A[j], k, gb, hb, P)
        ## print(paste("k", k, "; sb/Sb", sb/Sb))
        if (sb/Sb > sf/Sf) {
          break;
        }
        k0 <- k
        sb0 <- sb
      }

      ## If a new point hasn't already been created for a
      ## corresponding point, Create new point
      if (hb[i] == i) {
        f <- (sf/Sf*Sb-sb0)/(sb-sb0)
        p <- (1-f) * P[k0,] + f * P[k,]
        P <- rbind(P, p)

        ## Update forward and backward pointers
        n <- nrow(P)
        gb[n]     <- k
        gf[n]     <- gf[k]
        gb[gf[k]] <- n
        gf[k]     <- n

        ## Update correspondences
        hf[n] <- n
        hb[n] <- n
        h[i] <- n
      } else {
        h[i] <- h[hb[i]]
      }
      
      ## print(paste("n =", n, "; k =", k, "; k0 =", k0,
      ##            "; gf[", n, "] =", gf[n], "; gb[", n,  "] =", gb[n],
      ##             "; gf[", k, "] =", gf[k], "; gb[", k0, "] =", gb[k0]))
    }

    ## plot.outline(P, gb)
    ## print(paste("Forwards", j))
    ## readline("Press <Enter> to continue")
      
    ## Go along backward path
    for (i in setdiff(TBset[[j]], c(A[j], VB[j]))) {
      sb <- path.length(A[j], i, gb, hb, P)
      ## print(paste("i", i, "sb", sb/Sb))
      ## print(TFset[[j]])
      for (k in TFset[[j]]) {
        sf <- path.length(A[j], k, gf, hf, P)
        ## print(paste("k", k, "; sf/Sf", sf/Sf))
        if (sf/Sf > sb/Sb) {
          break;
        }
        k0 <- k
        sf0 <- sf
      }

      ## If a new point hasn't already been created for a
      ## corresponding point, Create new point
      if (hf[i] == i) {
        f <- (sb/Sb*Sf-sf0)/(sf-sf0)
        p <- (1-f) * P[k0,] + f * P[k,]
        P <- rbind(P, p)
        
        ## Update forward and backward pointers
        n <- nrow(P)
        gf[n]  <- k
        gb[n]  <- gb[k]
        gf[gb[k]] <- n
        gb[k]  <- n
        
        ## Update correspondences
        hf[n] <- n
        hb[n] <- n
        h[i] <- n
      } else {
        h[i] <- h[hf[i]]
      }
      
      ## print(paste("n =", n, "; k =", k, "; k0 =", k0,
      ##            "; gf[", n,  "] =", gf[n], "; gb[", n,  "] =", gb[n],
      ##            "; gf[", k0, "] =", gf[k0], "; gb[", k, "] =", gb[k]))
    }
    
    ## plot.outline(P, gb)
    ## print(paste("Backwards", j))
    ## readline("Press <Enter> to continue")
  }

  ## Make sure that there are no chains of correspondences
  h <- c(h, (length(h)+1):nrow(P))
  while (!all(h==h[h])) {
    h <- h[h]
  }
  
  return(list(Rset=Rset, i0=i0,
              VF=VF, VB=VB, A=A,
              TFset=TFset, TBset=TBset,
              P=P, h=h, hf=hf, hb=hb,
              gf=gf, gb=gb))
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

## Now for the dreaded elastic error function....
E <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, Nphi, E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(0, N)
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
    areas <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))
    ##E.A <- 0.5 * sum((areas - A)^2/A)
    E.A <- sum(exp(-k.A*areas/A))
  }
  
  return(E.E + E0.A*E.A)
}

## ... and the even more dreaded gradient of the elastic error
dE <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, Nphi, E0.A=0.1, k.A=1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(0, N)
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
    areas <- dot(P[T[,1],], dAdPt1)

##     dEdPt1 <- (areas - A)/A * dAdPt1
    dEdPt1 <- -k.A/A*exp(-k.A*areas/A) * dAdPt1
    
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
  
  return(c(dE.E.phii[-Rset]           + E0.A * dE.A.dphi[-Rset],
           dE.E.dlambdai[-i0]    + E0.A * dE.A.dlambda[-i0]))
}

E.dE <- function(p, Cu, C, L, B, T, A, R, Rset, phi0, Nphi, E0.A=0.1, N) {
  E <- E(p, Cu, C, L, B, T, A, R, Rset, phi0, Nphi, E0.A=E0.A, N)
  attr(E, "gradient") <- dE(p, Cu, C, L, B, T, A, R, Rset, phi0, Nphi, E0.A=E0.A, N)
  return(E)
}

## Combined energy function
##E <- function(p, Cu, C, L, B, Pt, A, R,
##              E0.E=1, E0.A=1,
##              Rset, phi0, verbose=FALSE) {
##  E <- E0.E * E.E(p, Cu, C, L, B, R, Rset, phi0, verbose=verbose) 
##  if (E0.A) {
##    E <- E + E0.A * E.A(p, Pt, A, R, Rset, phi0, verbose=verbose)
##  }
##  return(E) 
##}

## Combined gradient
##dE <- function(p, Cu, C, L, B, Pt, A, R,
##               E0.E=1, E0.A=1,
##               Rset, phi0, verbose=FALSE) {
##  dE <- E0.E * dE.E(p, Cu, C, L, B, R, Rset, phi0, verbose=verbose)
##  if (E0.A) {
##    dE <- dE + E0.A * dE.A(p, Pt, A, R, Rset, phi0, verbose=verbose)
##  }
##  return(dE)
##}

## Grand optimisation function
optimise.mapping <- function(p, m, t, s, E0.A=1, k.A=1, method="BFGS") {
  phi <- p$phi
  lambda <- p$lambda
  R <- p$R
  phi0 <- p$phi0
  Tt <- m$Tt
  a <- t$a
  Cut <- m$Cut
  Ct <- m$Ct
  Pt <- m$Pt
  Lt <- m$Lt
  Bt <- m$Bt
  Rsett <- m$Rsett
  i0t <- m$i0t
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)

  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi[-Rsett], lambda[-i0t])
  opt$conv <- 1
  while (opt$conv) {
    opt <- optim(opt$p, E, gr=dE,
                 method=method,
                 T=Tt, A=a, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R,
                 E0.A=E0.A, k.A=k.A, N=Nt, 
                 Rset=Rsett, i0=i0t, phi0=phi0, Nphi=Nphi, verbose=FALSE)
    ## print(opt)
    ##               control=list(maxit=200))
    print(E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=a,
            E0.A=E0.A, N=Nt,
            Rset=Rsett, i0=i0t, phi0=phi0, Nphi=Nphi))
    phi            <- rep(phi0, Nt)
    phi[-Rsett]    <- opt$p[1:Nphi]
    lambda         <- rep(0, Nt)
    lambda[-i0t] <- opt$p[Nphi+1:(Nt-1)]

    lt <- compute.lengths(phi, lambda, Cut, R)
    ## lt <- R*central.angle(phi1, lambda1, phi2, lambda2)
    plot.retina(phi, lambda, R, Tt, Rsett) ## , ts.red, ts.green, edge.inds)
    with(s, plot.outline(P, gb))
    with(t, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0))
  }
  return(list(phi=phi, lambda=lambda))
}

## Try to simulate the mapping using Euler integation
solve.mapping <- function(p, m, t, s, E0.A=0, dt=1E-12, nstep=100, Rexp=1, verbose=FALSE) {
  phi <- p$phi
  lambda <- p$lambda
  R <- p$R
  phi0 <- p$phi0
  Tt <- m$Tt
  a <- t$a
  Cut <- m$Cut
  Ct <- m$Ct
  Pt <- m$Pt
  Lt <- m$Lt
  Bt <- m$Bt
  Rsett <- m$Rsett
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)

  p <- c(phi[-Rsett], lambda[-Rsett[1]])
  ## Optimisation and plotting
  for (i in 0:nstep) {
    dEbydp <- dE(p, T=Tt, A=a, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R*Rexp, E0.A=E0.A, N=Nt, Rset=Rsett, phi0=phi0, Nphi=Nphi)
    p <- p - dEbydp * dt

    p1 <- p - dEbydp * dt/2
    dEbydp <- dE(p1, T=Tt, A=a, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R*Rexp, E0.A=E0.A, N=Nt, Rset=Rsett, phi0=phi0, Nphi=Nphi)
    p1 <- p1 - dEbydp * dt/2

    Delta <- max(abs(p-p1))

  ##  print(dt)
          lt <- compute.lengths(phi, lambda, Cut, R)
      print(c(E(p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=a,
            E0.A=E0.A, N=Nt,
            Rset=Rsett, phi0=phi0, Nphi=Nphi, verbose=verbose), cor(lt, Lt)))
    
    dt <- dt*(0.001/Delta)^(1/2)
    print(dt)
    
    phi            <- rep(phi0, Nt)
    phi[-Rsett]    <- p[1:Nphi]
    lambda         <- rep(0, Nt)
    lambda[-Rsett[1]] <- p[Nphi+1:(Nt-1)]

    ## Output
    if (!(i %% 100)) {
      lt <- compute.lengths(phi, lambda, Cut, R)
      print(c(E(p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=a,
            E0.A=E0.A, N=Nt,
            Rset=Rsett, phi0=phi0, Nphi=Nphi, verbose=verbose), cor(lt, Lt)))
      plot.retina(phi, lambda, R, Tt, Rsett) ## , ts.red, ts.green, edge.inds)
      with(s, plot.outline(P, gb))
      with(t, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0))
    }
  }
  return(list(phi=phi, lambda=lambda))
}

## Try to simulate the mapping using Euler integation
solve.mapping.momentum <- function(p, m, t, s, E0.A=0, dt=1E-6, nstep=100, Rexp=1, verbose=FALSE) {
  phi <- p$phi
  lambda <- p$lambda
  R <- p$R
  phi0 <- p$phi0
  Tt <- m$Tt
  a <- t$a
  Cut <- m$Cut
  Ct <- m$Ct
  Pt <- m$Pt
  Lt <- m$Lt
  Bt <- m$Bt
  Rsett <- m$Rsett
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)

  p <- c(phi[-Rsett], lambda[-Rsett[1]])
  p1 <- p
  p2 <- p
  ## Optimisation and plotting
  gamma <- 0.0002
  ## gamma <- 0
  mu <- 0.005
  for (i in 0:nstep) {
    dEbydp <- dE(p, T=Tt, A=a, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R*Rexp, E0.A=E0.A, N=Nt, Rset=Rsett, phi0=phi0, Nphi=Nphi)
    p <- 1/(mu/dt^2 + gamma/2/dt)*(-dEbydp + 2*mu/dt^2 * p - (mu/dt^2 - gamma/2/dt)*p1)
    p2 <- p1
    p1 <- p
    
    phi            <- rep(phi0, Nt)
    phi[-Rsett]    <- p[1:Nphi]
    lambda         <- rep(0, Nt)
    lambda[-Rsett[1]] <- p[Nphi+1:(Nt-1)]


    ## Output
    if (!(((i*dt) / 1E-6) %% 1000)) {
       lt <- compute.lengths(phi, lambda, Cut, R)
       print(c(E(p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=a,
             E0.A=E0.A, N=Nt,
                 Rset=Rsett, phi0=phi0, Nphi=Nphi, verbose=verbose), cor(lt, Lt)))
       plot.retina(phi, lambda, R, Tt, Rsett) ## , ts.red, ts.green, edge.inds)
       with(s, plot.outline(P, gb))
       with(t, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0))
    }
  }
  return(list(phi=phi, lambda=lambda))
}

## Method from Wang et al (2002)
solve.mapping.momentum2 <- function(p, m, t, s, E0.A=0, dt=1E-6, nstep=100, Rexp=1, verbose=FALSE) {
  phi <- p$phi
  lambda <- p$lambda
  R <- p$R
  phi0 <- p$phi0
  Tt <- m$Tt
  a <- t$a
  Cut <- m$Cut
  Ct <- m$Ct
  Pt <- m$Pt
  Lt <- m$Lt
  Bt <- m$Bt
  Rsett <- m$Rsett
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)

  p <- c(phi[-Rsett], lambda[-Rsett[1]])
  v <- rep(0, length(p))

  for (i in 0:nstep) {
    f <- -dE(p, T=Tt, A=a, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R*Rexp, E0.A=E0.A, N=Nt, Rset=Rsett, phi0=phi0, Nphi=Nphi)
    p <- p + dt*v + dt^2/2*f
    v <- 0.9*v + dt*f
    
    phi            <- rep(phi0, Nt)
    phi[-Rsett]    <- p[1:Nphi]
    lambda         <- rep(0, Nt)
    lambda[-Rsett[1]] <- p[Nphi+1:(Nt-1)]

    ## Output
    if (!(((i*dt) / 1E-6) %% 100)) {
       lt <- compute.lengths(phi, lambda, Cut, R)
       print(c(E(p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=a,
             E0.A=E0.A, N=Nt,
                 Rset=Rsett, phi0=phi0, Nphi=Nphi, verbose=verbose), cor(lt, Lt)))
       plot.retina(phi, lambda, R, Tt, Rsett) ## , ts.red, ts.green, edge.inds)
       with(s, plot.outline(P, gb))
       with(t, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0))
    }
  }
  return(list(phi=phi, lambda=lambda))
}


optimise.mapping.nlm <- function(p, m, t, s, E0.A=0, Rexp=1, verbose=FALSE,
                                 iterlim=1000) {
  phi <- p$phi
  lambda <- p$lambda
  R <- p$R
  phi0 <- p$phi0
  Tt <- m$Tt
  a <- t$a
  Cut <- m$Cut
  Ct <- m$Ct
  Pt <- m$Pt
  Lt <- m$Lt
  Bt <- m$Bt
  Rsett <- m$Rsett
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)

  p <- c(phi[-Rsett], lambda[-Rsett[1]])
  ## Optimisation and plotting
    opt <- nlm(E.dE, p, T=Tt, A=a, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R*Rexp, E0.A=E0.A, N=Nt, Rset=Rsett, phi0=phi0, Nphi=Nphi,
               print.level=1, check.analyticals=TRUE,
               fscale=1000, gradtol=1e-4, iterlim=iterlim, stepmax=10)

  
    phi        <- rep(phi0, Nt)
    phi[-Rsett] <- opt$est[1:Nphi]
    lambda         <- rep(0, Nt)
    lambda[-Rsett[1]] <- p[Nphi+1:(Nt-1)]

    ## Output

      lt <- compute.lengths(phi, lambda, Cut, R)
      print(c(E(opt$est, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=a,
            E0.A=E0.A, N=Nt,
            Rset=Rsett, phi0=phi0, Nphi=Nphi, verbose=verbose), cor(lt, Lt)))
      plot.retina(phi, lambda, R, Tt, Rsett) ## , ts.red, ts.green, edge.inds)
      with(s, plot.outline(P, gb))
      with(t, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0))


  return(list(phi=phi, lambda=lambda))
}


## Stuff for doing nstiff - not yet tried
##   ## Convert phis and lambdas to Carteisan coordinates
##   q0 <-  c(cos(phi)*cos(lambda),
##            cos(phi)*sin(lambda),
##            sin(phi))

##   ## Now define matricies needed for nstiff
##   gamma <- 1                            # friction coef
##   M <- diag(rep(1, Nt))
##   Q <- function(q, v) { - gamma * v }   # forces - friction and tension
  

##
## Geometry functions
## 

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

## plot.gridline.flat(P, T, phi, lambda, Tt, n, d)
##
## Plot a gridline from the spherical retina (described by points phi,
## lambda and triangulation Tt) onto a flattened retina (described by
## points P and triangulation T). The gridline is described by a
## normal n to a plane and a distance to the plane. The intersection of
## the plane and the spehere is the gridline.
plot.gridline.flat <- function(P, T, phi, lambda, Tt, n, d, ...) {
  mu <- compute.intersections.sphere(phi, lambda, Tt, n, d)

  ## Take out rows that are not intersections
  tri.int <- (rowSums((mu >=0) & (mu <=1)) == 2)

  if (any(tri.int)) {
    T  <- T[tri.int,,drop=FALSE]
    mu <- mu[tri.int,,drop=FALSE]

    line.int <- (mu >=0) & (mu <=1)
    
    ## Order rows so that the false indicator is in the third column
    T[!line.int[,2] ,] <- T[!line.int[,2], c(3,1,2)]
    mu[!line.int[,2],] <- mu[!line.int[,2],c(3,1,2)]
    T[!line.int[,1] ,] <- T[!line.int[,1], c(2,3,1)]
    mu[!line.int[,1],] <- mu[!line.int[,1],c(2,3,1)]

    P1 <- mu[,1] * P[T[,3],] + (1-mu[,1]) * P[T[,2],]
    P2 <- mu[,2] * P[T[,1],] + (1-mu[,2]) * P[T[,3],]
    segments(P1[,1], P1[,2], P2[,1], P2[,2], ...)
  }
}

## plot.gridlines.flat(P, T, phi, lambda, Tt, phi0)
##
## Plot a mesh of gridlines from the spherical retina (described by
## points phi, lambda and triangulation Tt and cutoff point phi0) onto
## a flattened retina (described by points P and triangulation T).
plot.gridlines.flat <- function(P, T, phi, lambda, Tt, phi0,
                                Phis=(-5:6)*pi/12, Lambdas=(0:23)*pi/24, ...) {
  Phis <- Phis[Phis<phi0]
  for (Phi in Phis) {
    plot.gridline.flat(P, T, phi, lambda, Tt, c(0,0,1), sin(Phi), ...)
  }
  for (Lambda in Lambdas) {
    plot.gridline.flat(P, T, phi, lambda, Tt, c(sin(Lambda),cos(Lambda),0), 0, ...)
  }
}

##
## Plotting functions
## 

## plot.outline(P, gb)
##
## Plot outline of retina given set of outline points P and backwards
## pointer gb
plot.outline <- function(P, gb, add=FALSE, ...) {
  if (!add) {
    par(mar=c(1, 1, 1, 1))
    plot(P, pch=".", xaxt="n", yaxt="n", xlab="", ylab="",
         bty="n")
  }
  segments(P[,1], P[,2], P[gb,1], P[gb, 2], ...)
}

## plot.stitch(P, s)
##
## Plot stitch given set of outline points stitch information s
plot.stitch <- function(s, ...) {
  with(s, {
    plot.outline(P, gb, ...)
    points(P[VF,], col="red", pch="+")
    points(P[VB,], col="orange", pch="+")
    points(P[A, ], col="cyan", pch="+")
    for (TF in TFset) {
      lines(P[TF,], col="red", ...)
    }
    for (TB in TBset) {
      lines(P[TB,], col="orange", ...)
    }
    for (j in 1:length(h)) {
      if (h[j] != j) {
        lines(P[c(j,h[j]),], col="blue", ...)
      }
    }
    
    for (j in 1:length(hf)) {
      if (hf[j] != j) {
        lines(P[c(j,s$hf[j]),], col="green", ...)
      }
    }
  })
  ## points(P[s$Rset,], col="red")
}

## Function to plot the retina in spherical coordinates
## phi - lattitude of points
## lambda - longitude of points
## R - radius of sphere
## Tt - triagulation
## Rsett - members of rim set
plot.retina <- function(phi, lambda, R, Tt, Rsett) {
  ## Now plot this in 3D space....
  x <- R*cos(phi)*cos(lambda) 
  y <- R*cos(phi)*sin(lambda)
  z <- R*sin(phi)
  P <- cbind(x, y, z)
  rgl.clear()
  rgl.bg(color="white")

  ## Outer triangles
  triangles3d(matrix(1.01*x[t(Tt[,c(2,1,3)])], nrow=3),
              matrix(1.01*y[t(Tt[,c(2,1,3)])], nrow=3),
              matrix(1.01*z[t(Tt[,c(2,1,3)])], nrow=3),
              color="darkgrey", alpha=1)
  
  ## Inner triangles
  triangles3d(matrix(x[t(Tt)], nrow=3),
              matrix(y[t(Tt)], nrow=3),
              matrix(z[t(Tt)], nrow=3),
              color="white", alpha=1)

  ## Plot any flipped triangles
  ## First find verticies and find centres and normals of the triangles
  P1 <- P[Tt[,1],]
  P2 <- P[Tt[,2],]
  P3 <- P[Tt[,3],]
  cents <- (P1 + P2 + P3)/3
  normals <- 0.5 * extprod3d(P2 - P1, P3 - P2)
  areas <- apply(normals^2, 1, sum)
  ##  print(cents)
  ##  print(areas)
  flipped <- (-dot(cents, normals) < 0)
  print(paste(sum(flipped), "flipped triangles:"))
  print(which(flipped))
  points3d(cents[flipped,1], cents[flipped,2], cents[flipped,3], col="blue", size=5)
}

## Function to determine the locations of cell bodies on the folded
## retina in Cartesian (X, Y, Z) coordinates
## phi    - lattitude of mesh points
## lambda - longitude of mesh points
## R      - radius of sphere
## Tt     - triagulation
## cb     - object returned by tsearch containing information on the
cell.bodies.folded.cart <- function(phi, lambda, R, Tt, cb) {
  ## Obtain Cartesian coordinates of points
  P <- cbind(R*cos(phi)*cos(lambda),
             R*cos(phi)*sin(lambda),
             R*sin(phi))

  ## Now find locations cc of cell bodies in Cartesian coordinates
  cc <- matrix(0, 0, 3)
  colnames(cc) <- c("X", "Y", "Z")
  for(i in 1:(dim(cb$p)[1])) {
    cc <- rbind(cc, bary2cart(P[Tt[cb$idx[i],],], cb$p[i,]))
  }
  return(cc)
}

## Function to determine the locations of cell bodies on the folded
## retina in spherical (lambda, phi) coordinates
## phi    - lattitude of mesh points
## lambda - longitude of mesh points
## R      - radius of sphere
## Tt     - triagulation
## cb     - object returned by tsearch containing information on the
cell.bodies.folded.sphere <- function(phi, lambda, R, Tt, cb) {
  ## Get locations in Cartesian coordinates
  cc <- cell.bodies.folded.cart(phi, lambda, R, Tt, cb)
  ## Convert to spherical coordinates
  return(list(phi=asin(cc[,"Z"]/R),
              lambda=atan2(cc[,"Y"], cc[,"X"])))
}

## Function to plot cell bodies on a retina in spherical coordinates
## It assumes that plot.retina has been called already
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## Tt     - triagulation
## cb     - object returned by tsearch containing information on the
##          triangle in which a cell body is found and its location
##          within that triangle in barycentric coordinates
## radius - radius of the spheres to plot
## color  - colour of the spheres to plot
plot.cell.bodies <- function(phi, lambda, R, Tt, cb, size=R/10, color="red") {
  ## Obtain Cartesian coordinates of points
  cc <- cell.bodies.folded.cart(phi, lambda, R, Tt, cb)
  
  ## Plot
  ## shade3d( translate3d( cube3d(col=color), cc[,1], cc[,2], cc[,3]))
  ## rgl.spheres(cc[,1], cc[,2], cc[,3], radius, color=color)
  ## points3d(cc[,1], cc[,2], cc[,3], size=size, color=color)
  ## cc <- cc * 1.01
  ## points3d(cc[,1], cc[,2], cc[,3], size=size, color=color)

  ## Custom code required to plot triangles
  ax1 <- 1/sqrt(apply(cc[,1:2]^2, 1, sum)) * cbind(-cc[,2], cc[,1], 0)
  ## print(ax1)
  
  ax2 <- extprod3d(cc, ax1)
  ax2 <- ax2/sqrt(apply(ax2^2, 1, sum))
  ##print(ax2)

  print(dot(ax1, ax2))
  
  v1 <- cc + size *  ax1/2
  v2 <- cc + size * (-ax1/4 + sqrt(3)/4*ax2)
  v3 <- cc + size * (-ax1/4 - sqrt(3)/4*ax2)

  inmag <- 0.99
  outmag <- 1.02
  
  x <- rbind(v2[,1], v1[,1], v3[,1])
  y <- rbind(v2[,2], v1[,2], v3[,2])
  z <- rbind(v2[,3], v1[,3], v3[,3])
  triangles3d(inmag*x, inmag*y, inmag*z, color=color)

  x <- rbind(v1[,1], v2[,1], v3[,1])
  y <- rbind(v1[,2], v2[,2], v3[,2])
  z <- rbind(v1[,3], v2[,3], v3[,3])
  triangles3d(outmag*x, outmag*y, outmag*z, color=color)
}

## Function to plot cell bodies in spherical coordinates on a polar plot
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## Tt     - triagulation
## cbs    - list of objects returned by tsearch containing information on the
##          triangle in which a cell body is found and its location
##          within that triangle in barycentric coordinates
## phi0   - lattitude of the rim in radians
## cols   - colour of points to plot for each object in cbs
plot.cell.bodies.polar <- function(phi, lambda, R, Tt, cbs, phi0, cols="red",
                                   pch=".", ...) {
  ## Need to organise phis and lambdas into matricies, with
  ## one column per set of data
  phis <- matrix(NA, length(cbs), 0)
  lambdas <- matrix(NA, length(cbs), 0)
  for (i in 1:length(cbs)) {
    cs <- cell.bodies.folded.sphere(phi, lambda, R, Tt, cbs[[i]])
    d <- length(cs$phi) - ncol(phis)
    print(d)
    if (d>0) {
      phis <- cbind(phis, matrix(NA, length(cbs), d))
      lambdas <- cbind(lambdas, matrix(NA, length(cbs), d))
    }
    phis[i,1:length(cs$phi)] <- cs$phi
    lambdas[i,1:length(cs$lambda)] <- cs$lambda
  }
  print(phis)
  radial.lim <- c(seq(-90, phi0*180/pi, by=10), phi0*180/pi)
  radial.labels <- radial.lim
  radial.labels[(radial.lim %% 90) != 0] <- ""
  radial.labels[length(radial.labels)] <- phi0*180/pi
  polar.plot(phis*180/pi, polar.pos=lambdas*180/pi+90,
             rp.type="s", point.col=cols,
             radial.lim=radial.lim,
             radial.labels=radial.labels,
             label.pos=c(0, 90, 180, 270),
             labels=c("N", "D", "T", "V"),
             point.symbols=pch, ...)
}

## Compute mean on sphere
folded.mean.sphere <- function(phi, lambda) {
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

## Convert elevation in spherical coordinates into radius in polar
## coordinates in an area-preserving projection
spherical.to.polar.area <- function(phi) { return(sqrt(2*(1 +
  sin(phi)))) }

## Convert polar coordinates to cartesian coordinates
polar.to.cart <- function(r, theta) {
  return(cbind(x=r*cos(theta), y=r*sin(theta)))   
}

## Function to plot cell bodies in spherical coordinates on a polar plot
## phi    - lattitude of points
## lambda - longitude of points
## R      - radius of sphere
## Tt     - triagulation
## cbs    - list of objects returned by tsearch containing information on the
##          triangle in which a cell body is found and its location
##          within that triangle in barycentric coordinates
## phi0   - lattitude of the rim in radians
## cols   - colour of points to plot for each object in cbs
plot.cell.bodies.polar.area <- function(phi, lambda, R, Tt, cbs, phi0, cols="red",
                                   pch=".", ...) {
  plot(NA, NA, xlim=c(-2,2), ylim=c(-2, 2))
  for (i in 1:length(cbs)) {
    cs <- cell.bodies.folded.sphere(phi, lambda, R, Tt, cbs[[i]])
    ## Turn into polar coordinates, shifting round by 90 degress for plotting
    lambdas <- cs$lambda+pi/2
    p <- polar.to.cart(spherical.to.polar.area(cs$phi), lambdas)
    points(p[,"x"], p[,"y"], pch=pch, col=cols[i], ...)

    ## Compute mean and plot
    m <- folded.mean.sphere(cs$phi, lambdas)
    p <- polar.to.cart(spherical.to.polar.area(m[1]), m[2])
    points(p[,"x"], p[,"y"], col=cols[i], pch="+", ...)
  }
  ## Draw circular grid
  dl <- 2*pi/90
  lambdas <- seq(dl, 2*pi, by=dl)
  phi.degs <- seq(-80, phi0*180/pi, by=10)
  rs <- spherical.to.polar.area(phi.degs*pi/180)
  polygon(rbind(outer(cos(lambdas), rs), NA),
          rbind(outer(sin(lambdas), rs), NA),  col=NA, border="grey")

  ## Draw axes and label
  axis(side=1, pos=0, at=c(-max(rs), 0, 1, max(rs)), labels=c(NA, -90, 0, phi0*180/pi))
  axis(side=2, pos=0, at=c(-max(rs), max(rs)), labels=c(NA, NA))
  text(2, 0, "N")
  text(-2, 0, "T")
  text(0, 2, "D")
  text(0, -2, "V")
  
}

plot.outline.retina <- function(phi, lambda, R, gb, h, ...) {
  ## Obtain Cartesian coordinates of points
  Pc <- cbind(R*cos(phi)*cos(lambda),
             R*cos(phi)*sin(lambda),
             R*sin(phi))

  P <- Pc*0.99
##   segments3d(rbind(P[h[gb[gb]],1], P[h[gb],1]),
##              rbind(P[h[gb[gb]],2], P[h[gb],2]),
##              rbind(P[h[gb[gb]],3], P[h[gb],3]),
##              ...)
  rgl.lines(rbind(P[h[gb[gb]],1], P[h[gb],1]),
            rbind(P[h[gb[gb]],2], P[h[gb],2]),
            rbind(P[h[gb[gb]],3], P[h[gb],3]),
             ...)
  
   P <- Pc*1.001
##   segments3d(rbind(P[h[gb[gb]],1], P[h[gb],1]),
##              rbind(P[h[gb[gb]],2], P[h[gb],2]),
##              rbind(P[h[gb[gb]],3], P[h[gb],3]),
##              ...)

  rgl.lines(rbind(P[h[gb[gb]],1], P[h[gb],1]),
            rbind(P[h[gb[gb]],2], P[h[gb],2]),
            rbind(P[h[gb[gb]],3], P[h[gb],3]),
             ...)
  
}

## make.triagulation(P, n)
##
## Create a triangulation of the outline defined by the points P,
## which are represented as N*2 matrix. There should be at least n
## triangles in the triangulation
## Returns a list comprising:
## P   - The set of new points, with the existing points at the start
## T   - The triangulation
## a   - Array containing area of each triangle
## A   - Total area of outline
## Cu  - Unique set of M connections, as M*2 matrix
## L   - Length of each connection
## h   - Correspondances vector?????
make.triangulation <- function(P, h=1:nrow(P), g=NULL, n=200,
                               suppress.external.steiner=FALSE) {
  S <- NULL
  if (!is.null(g)) {
    S <- pointers2segments(g)
  }
  ## Make initial triangulation to determine area
  out <- triangulate(P, S, Y=TRUE)
  A <- sum(with(out, tri.area(P, T)))
  print(A)
  if (!is.na(n)) {
    out <- triangulate(P, S, a=A/n, q=20, Y=suppress.external.steiner)
  }
  P <- out$P
  T <- out$T

  ## Create pointers from segments

  ## To ensure the correct orientaion, we use the fact that the
  ## triangles are all anticlockwise in orinentation, and that the
  ## orientation of the first row of the segment matrix determines the
  ## orientation of all the other rows.

  ## We therefore find triangle which contains the first segment
  S <- out$S
  T1 <- which(apply(out$T, 1, function(x) {all(S[1,] %in% x)}))

  ## Then find out which of the vertices in the triangle is not the
  ## one we need
  i <- which((out$T[T1,] %in% S[1,]) == FALSE)
  if (i == 3) S[1,] <- out$T[T1,c(1,2)]
  if (i == 2) S[1,] <- out$T[T1,c(3,1)]
  if (i == 1) S[1,] <- out$T[T1,c(2,3)]

  print(out$T[T1,])
  print(S[1,])
  
  gf <- segments2pointers(S)
  gb <- gf
  gb[na.omit(gf)] <- which(!is.na(gf))
  Rset <- na.omit(gf)
  
  ## trimesh(T, P, col="grey", add=TRUE)

  ## Derive edge matrix from triangulation
  Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
  Cu <- Unique(Cu, TRUE)

  ## If we are in the business of refining triangles (i.e. specifying
  ## n), remove lines which join non-ajancent parts of the outline
  if (!is.na(n)) {
    for (i in 1:nrow(Cu)) {
      C1 <- Cu[i,1]
      C2 <- Cu[i,2]
      if (all(Cu[i,] %in% Rset)) {
        if (!((C1 == gf[C2]) ||
              (C2 == gf[C1]))) {
          ## Find triangles containing the line
          ## segments(P[C1,1], P[C1,2], P[C2,1], P[C2,2], col="yellow")
          Tind <- which(apply(T, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
          print(paste("Non-adjacent points in rim connected by line:", C1, C2))
          print(paste("In triangle:", Tind))
          ## Find points T1 & T2 in the two triangles which are not common
          ## with the edge
          T1 <- setdiff(T[Tind[1],], Cu[i,])
          T2 <- setdiff(T[Tind[2],], Cu[i,])
          print(paste("Other points in triangles:", T1, T2))
          ## Create a new point at the centroid of the four verticies
          ## C1, C2, T1, T2
          p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
          P <- rbind(P, p)
          n <- nrow(P)
          ## Remove the two old triangles, and create the four new ones
          T[Tind[1],] <- c(n, C1, T1)
          T[Tind[2],] <- c(n, C1, T2)
          T <- rbind(T,
                     c(n, C2, T1),
                     c(n, C2, T2))
        }
      }
    }

    ## Add the new points to the correspondances vector
    h <- c(h, (length(h)+1):nrow(P))

    ## Create the edge matrix from the triangulation
    Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
    Cu <- Unique(Cu, TRUE)
  }

  ## Swap orientation of triangles which have clockwise orientation
  a.signed <- tri.area.signed(P, T)
  T[a.signed<0,c(2,3)] <- T[a.signed<0,c(3,2)]
  a <- abs(a.signed)
    
  ## Find lengths of connections
  L <- norm(P[Cu[,1],] - P[Cu[,2],])

  ## Check there are no zero-length lines
  if (any(L==0)) {
    print("WARNING: zero-length lines")
  }
  
  return(list(P=P, T=T, Cu=Cu, h=h, a=a, A=A, L=L,
              gf=gf, gb=gb, S=out$S, E=out$E, EB=out$EB))
}

## Convert a matrix containing on each line the indicies of the points
## forming a segment, and convert this to two sets of ordered pointers
segments2pointers <- function(S) {
  g <- c()
  j <- 1                                # Row of S
  k <- 1                                # Column of S
  while(nrow(S) > 0) {
    i <- S[j,3-k]
    g[S[j,k]] <- i
    S <- S[-j,,drop=FALSE]
    if (nrow(S) == 0) {
      return(g)
    }
    j <- which(S[,1] == i)
    if (length(j)) {
      k <- 1
    } else {
      j <- which(S[,2] == i)
      k <- 2
      if (!length(j)) {
        print("No matching index")
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
## triangulation list t and the stitch list s, and the index of the landmark i0
##
## The information includes 
## h    - (in s) the point correspondence mapping
## Rset - (in s) the set of points on the rim
## i0   - (is s) the index of the landmark
## T    - (in t) the triangulation
## Cu   - (in t) the edge matrix
## L    - (in t) the edge lengths
## P    - (in t) the point locations
##
## The function returns merged and transformed versions of all these
## objects (all suffixed with t), as well as a matrix Bt, which maps a
## binary vector representation of edge indicies onto a binary vector
## representation of the indicies of the points linked by the edge
merge.points.edges <- function(t, s) {
  h <- t$h
  T <- t$T
  Cu <- t$Cu
  L <- t$L
  P <- t$P
  
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
  gft <- ht[s$gf]

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
  ## the existing indicies onto the existing edge indicies
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
    if (length(is)>1) {
      print(L[is])
    }
    Lt[k] <- mean(L[is])
  }

  ## Transform the rim set
  Rsett <- unique(ht[s$Rset])
  i0t <- ht[s$i0]

  ## Create the symmetric connection set
  Ct <- rbind(Cut, Cut[,2:1])

  ## Matrix to map line segments onto the points they link
  Bt <- matrix(0, nrow(Pt), nrow(Ct))
  for (i in 1:nrow(Ct)) {
    Bt[Ct[i,1],i] <- 1
  }
  
  return(list(Pt=Pt, Tt=Tt, Ct=Ct, Cut=Cut, Bt=Bt, Lt=Lt, ht=ht,
              Rsett=Rsett, i0t=i0t, P=P))
}

## Project mesh points in the flat outline onto a sphere
##
## The information to be merged is contained in the
## merge structure m and the triangulation t and the lattitude
## at which the sphere is cut off is given in phi0
##
## The information includes:
## Pt     - (in m) the mesh point coordinates
## Rsett  - (in m) the set of points on the rim
## A      - (in t) the area of the flat outline
##
## The function returns a list with the following members:
## phi    - lattitude of mesh points
## lambda - longitude of mesh points
## R      - radius of sphere
## phi0   - lattitude at which sphere is cut off (from input)
project.to.sphere <- function(m, t, phi0=50*pi/180, lambda0=lambda0) {
  Pt <- m$Pt
  Rsett <- m$Rsett
  i0t <- m$i0t
  A <- t$A
  
  Nt <- nrow(Pt)
  Nphi <- Nt - length(Rsett)

  ## From this we can infer what the radius should be from the formula
  ## for the area of a sphere which is cut off at a lattitude of phi0
  ## area = 2 * PI * R^2 * (sin(phi0)+1)
  R <- sqrt(A/(2*pi*(sin(phi0)+1)))

  ## Now assign each point to a location in the phi, lambda coordinates
  ## Shift coordinates to rough centre of grid
  x <- Pt[,1] - mean(Pt[,1]) 
  y <- Pt[,2] - mean(Pt[,2]) 
  phi <- -pi/2 + sqrt(x^2 + y^2)/(R)
  phi[Rsett] <- phi0
  lambda <- atan2(y, x)
  lambda <- lambda-lambda[i0t] + lambda0

  return(list(phi=phi, lambda=lambda, R=R, phi0=phi0, lambda0=lambda0))
}

## Folding routine
##
## Input arguments:
## P         - outline points as N-by-2 matrix
## tearmat   - tear matrix
## phi0      - lattitude of rim of folded object
## i0        - index of the landmark on the rim
## lambda0   - longitude of landmark on rim
## graphical - whether to plot graphs during computation
## 
## Returns list containing:
## t         - triangulation information
## s         - stitching information
## m         - information about merged points and edges
## p         - information about projection onto sphere
## r         - information about locations of gridlines
##
fold.outline <- function(P, tearmat, phi0=50, i0=NA, lambda0=0,
                         graphical=TRUE,
                         report=print) {
  report("Triangulating...")
  t <- make.triangulation(P, h=1:nrow(P), n=200)
  if (graphical) {
    with(t, trimesh(T, P, col="black"))
  }

  report("Stitching...")
  s <- stitch.outline(t$P, t$gf, t$gb, tearmat, i0)
  if (is.null(s)) {
    report("ERROR: Fixed point is not on the rim")
    return(NULL)
  }
  if (graphical) {
    plot.stitch(s)
  }

  report("Triangulating...")  
  t1 <- make.triangulation(s$P, h=s$h, g=s$gf, n=200,
                           suppress.external.steiner=TRUE)
  if (graphical) {
    plot.stitch(s)
    with(t1, trimesh(T, P, col="grey", add=TRUE))

  }

  report("Merging points...")
  m <- merge.points.edges(t1, s)

  ## if (graphical) {
  ##   plot(P)
  ##   with(s, plot.outline(P, gb))
  ## }

  report("Projecting to sphere...")
  p <- project.to.sphere(m, t1, phi0=phi0*pi/180, lambda0=lambda0*pi/180)

  if (graphical) {
    ## Initial plot in 3D space
    plot.retina(p$phi, p$lambda, p$R, m$Tt, m$Rsett)
  }

  report("Optimising mapping...")
  r <- optimise.mapping(p, m, t1, s, E0.A=exp(3), k.A=1)
  p1 <- p
  p1$phi <- r$phi
  p1$lambda <- r$lambda
  r <- optimise.mapping(p1, m, t1, s, E0.A=exp(10), k.A=20)

  report("Mapping optimised.")
  return(list(t=t1, s=s, m=m, p=p, r=r))
}

