source("triangulate.R")                 # for tri.area and tri.area.signed
require("rgl")                 

## Modulus function, that returns
## i , if i <= N
## i - N, if i > N
Mod <- function(i, N) {
  return((i - 1) %% N + 1)
}

## Distance norm
norm <- function(X) {
  if (is.vector(X)) {
    return(sqrt(sum(X^2)))
  } else {
    return(sqrt(apply(X^2, 1, sum)))
  }
}

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

## P is the set of points describing the edge.
## T is the tear matrix
stitch.retina <- function(P, T) {
  ## Extract information from tear matrix
  A  <- T[,1]                            # apicies of tears
  VB <- T[,2]                           # forward verticies
  VF <- T[,3]                           # backward verticies

  ## Create forward and backward pointers
  N <- nrow(P)                          # Number of points
  gf <- Mod((1:N) + 1, N) # c(2:N, 1)
  gb <- Mod((1:N) - 1, N) # c(N, 1:(N-1))
  
  ## Create initial sets of correspondances
  hf <- 1:N
  hb <- 1:N
  hf[VB] <- VF
  hb[VF] <- VB
  h <- hf

  ## Initialise the set of points in the rim
  Rset <- 1:N
  
  ## Create lists of forward and backward tears
  TFset <- list()
  TBset <- list()
  
  ## Iterate through the tears to create tear sets and rim set
  for (j in 1:nrow(T)) {
    ## Create sets of points for each tear and remove these points from
    ## the rim set
    TFset[[j]] <- Mod(path(A[j], VF[j], gf, hf), N)
    TBset[[j]] <- Mod(path(A[j], VB[j], gb, hb), N)
    Rset <- setdiff(Rset, setdiff(TFset[[j]], VF[j]))
    Rset <- setdiff(Rset, setdiff(TBset[[j]], VB[j]))
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

    ## Go along backward path
    for (i in setdiff(TBset[[j]], c(A[j], VB[j]))) {
      sb <- path.length(A[j], i, gb, hb, P)
      ## print(paste("i", i, "sb", sb/Sb))
      ## print(TFset[[j]])
      for (k in TFset[[j]]) {
        sf <- path.length(A[j], k, gf, hb, P)
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
  }

  ## Make sure that there are no chains of correspondences
  h <- c(h, (length(h)+1):nrow(P))
  while (!all(h==h[h])) {
    h <- h[h]
  }
  
  return(list(Rset=Rset,
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
compute.lengths <- function(phi, lamba, Cu, R) {
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)

  return(l)
}

## Now for the dreaded elastic error function....
E.E <- function(p, Cu, C, L, B, R, Rset, phi0, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- p[Nphi+1:N]
  
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)
  if (verbose==2) {
    print(l)
  }
  E <- sum((l - L)^2/L)
  # E <- sum((l - L)^2)
  if (verbose>=1) {
    print(E)
  }
  return(E)
}

## ... and the even more dreaded gradient of the elastic error
dE.E <- function(p, Cu, C, L, B, R, Rset, phi0, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- p[Nphi+1:N]

  phii   <- phi[C[,1]]
  lambdai <- lambda[C[,1]]
  phij    <- phi[C[,2]]
  lambdaj <- lambda[C[,2]]
  ## x is the argument of the acos in the central angle
  x <- sin(phii)*sin(phij) + cos(phii)*cos(phij)*cos(lambdai-lambdaj)
  ## the central angle
  l <- R * acos(x)
  if (verbose==2) {
    print(l)
  }
  fac <- R * (l - c(L, L))/(c(L, L) * sqrt(1-x^2))
  # fac <- R * (l - c(L, L))/sqrt(1-x^2)
  dE.phii     <- B %*% (fac * (sin(phii)*cos(phij)*cos(lambdai-lambdaj)
                               - cos(phii)*sin(phij)))
  dE.dlambdai <- B %*% (fac * cos(phii)*cos(phij)*sin(lambdai-lambdaj))
  return(c(dE.phii[-Rset], dE.dlambdai))
}

## Combined energy function
E <- function(p, Cu, C, L, B, Pt, A, R,
              E0.E=1, E0.A=1,
              Rset, phi0, verbose=FALSE) {
  E <- E0.E * E.E(p, Cu, C, L, B, R, Rset, phi0, verbose=verbose) 
##  if (E0.A) {
##    E <- E + E0.A * E.A(p, Pt, A, R, Rset, phi0, verbose=verbose)
##  }
  return(E) 
}

## Combined gradient
dE <- function(p, Cu, C, L, B, Pt, A, R,
               E0.E=1, E0.A=1,
               Rset, phi0, verbose=FALSE) {
  dE <- E0.E * dE.E(p, Cu, C, L, B, R, Rset, phi0, verbose=verbose)
##  if (E0.A) {
##    dE <- dE + E0.A * dE.A(p, Pt, A, R, Rset, phi0, verbose=verbose)
##  }
  return(dE)
}

## Grand optimisation function
optimise.mapping <- function(E0.E=1, E0.A=1, Rset=Rsett, L=L, method="BFGS") {
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi[-Rset], lambda)
  opt$conv <- 1
  while (opt$conv) {
    opt <- optim(opt$p, E, gr=dE,
                 method=method,
                 Pt=Tt, A=areas, Cu=Cut, C=Ct, L=L, B=Bt, R=R,
                 E0.E=E0.E, E0.A=E0.A,
                 Rset=Rset, phi0=phi0, verbose=FALSE)
    ## print(opt)
    ##               control=list(maxit=200))
    print(E(opt$p, Cut, Ct, L, Bt, Tt, areas, R, 
            E0.E=E0.E, E0.A=E0.A, Rset=Rset, phi0=phi0))
    phi        <- rep(phi0, N)
    phi[-Rset] <- opt$p[1:Nphi]
    lambda     <- opt$p[Nphi+1:N]
    plot.retina(phi, lambda, R, Tt) ## , ts.red, ts.green, edge.inds)
  }
  return(list(phi=phi, lambda=lambda))
}

##
## Plotting functions
## 

plot.stitch <- function(P, s) {
  plot(s$P, pch=".", cex=4, xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  with(s, segments(P[,1], P[,2], P[gb,1], P[gb, 2]))
  points(P[s$VF,], col="purple", pch="+")
  points(P[s$VB,], col="blue", pch="+")
  points(P[s$A, ], col="cyan", pch="+")
  for (TFset in s$TFset) {
    lines(P[TFset,], col="purple")
  }
  for (TBset in s$TBset) {
    lines(P[TBset,], col="blue")
  }
  for (j in 1:length(s$h)) {
    if (s$h[j] != j) {
      lines(s$P[c(j,s$h[j]),], col="orange")
    }
  }

  for (j in 1:length(s$hf)) {
    if (s$hf[j] != j) {
      lines(P[c(j,s$hf[j]),], col="green")
    }
  }
  ## points(P[s$Rset,], col="red")
}

plot.retina <- function(phi, lambda, R, Dt) {
  ## Now plot this in 3D space....
  x <- R*cos(phi)*cos(lambda) 
  y <- R*cos(phi)*sin(lambda)
  z <- R*sin(phi)
  P <- cbind(x, y, z)
  rgl.clear()
  rgl.bg(color="white")

  triangles3d(matrix(x[t(Dt)], nrow=3),
              matrix(y[t(Dt)], nrow=3),
              matrix(z[t(Dt)], nrow=3))
}


## Parameters
d <- 200           # Minimum spacing at which new points are inserted 

## ds <- with(s, apply(cbind(norm(P - P[gf,]), norm(P - P[gb,]))/2, 1, max))

## Read in and curate data to give edge.path and map
source("M634-4.R")
P <- edge.path[-nrow(edge.path),1:2]

s <- stitch.retina(P, tearmat)

plot.stitch(P, s)

## Create ordered version of P for determining outline
Po <- with(s, P[path(1, gb[1], gf, 1:nrow(P)),])

## Attempt to create Nrand points in the retina
P <- s$P
Nrand <- 1000
xmin <- min(P[,1])
xmax <- max(P[,1])
ymin <- min(P[,2])
ymax <- max(P[,2])
for (i in 1:Nrand) {
  C <- c(runif(1, xmin, xmax),
         runif(1, ymin, ymax))
  if (point.in.polygon(C[1], C[2], Po[,1], Po[,2])) {
    if (all(norm(t(t(P) - C)) > d)) {
      P <- rbind(P, C)
      ## ds <- c(ds, d)
    }
  }
}

## Delaunay triangulation and remove triangles outwith the retinal
## outline
T <- delaunayn(P)
Pc <- (P[T[,1],] + P[T[,2],] + P[T[,3],])/3 # Centres
T <- T[point.in.polygon(Pc[,1], Pc[,2], Po[,1], Po[,2])==1,]

trimesh(T, P, col="gray", add=TRUE)

## Find lines which join non-adjacent parts of the outline
Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
Cu <- Unique(Cu, TRUE)

for (i in 1:nrow(Cu)) {
  C1 <- Cu[i,1]
  C2 <- Cu[i,2]
  if (all(Cu[i,] %in% s$Rset)) {
    if (!((C1 == s$gf[C2]) ||
          (C2 == s$gf[C1]))) {
      ## Find triangles containing the line
      segments(P[C1,1], P[C1,2], P[C2,1], P[C2,2], col="red")
      Tind <- which(apply(T, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
      ## print(Cu[i,])
      ## print(Tind)
      T1 <- setdiff(T[Tind[1],], Cu[i,])
      T2 <- setdiff(T[Tind[2],], Cu[i,])
      ## print(T1)
      ## print(paste(C1, C2, T1, T2))
      p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
      points(p[1], p[2], col="red")
      P <- rbind(P, p)
      n <- nrow(P)
      T[Tind[1],] <- c(n, C1, T1)
      T[Tind[2],] <- c(n, C1, T2)
      T <- rbind(T,
                  c(n, C2, T1),
                  c(n, C2, T2))
    }
  }
}


Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
Cu <- Unique(Cu, TRUE)
print(length(which(Cu[,1] == s$gf[Cu[,2]])))
print(length(which(Cu[,2] == s$gf[Cu[,1]])))

l <- norm(P[Cu[,1],] - P[Cu[,2],])
while (max(l) > 2*d) {
  i <- which.max(l)
  ##  print(l[i])

  C1 <- Cu[i,1]
  C2 <- Cu[i,2]
  
  ## Find triangles containing the line
  ## segments(P[C1,1], P[C1,2],
  ## P[C2,1], P[C2,2], col="red")
  Tind <- which(apply(T, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
  ##  print(Cu[i,])
  ##  print(Tind)
  T1 <- setdiff(T[Tind[1],], Cu[i,])
  T2 <- setdiff(T[Tind[2],], Cu[i,])
  ##  print(T1)
  ##  print(paste(C1, C2, T1, T2))

  p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
  ## points(p[1], p[2], col="red")
  P <- rbind(P, p)
  n <- nrow(P)
  T[Tind[1],] <- c(n, C1, T1)
  T[Tind[2],] <- c(n, C1, T2)
  T <- rbind(T,
              c(n, C2, T1),
              c(n, C2, T2))
  
  Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
  Cu <- Unique(Cu, TRUE)

  ## print(length(which(Cu[,1] == s$gf[Cu[,2]])))
  ## print(length(which(Cu[,2] == s$gf[Cu[,1]])))

  ## Exlcude line segments which are on the edge for splitting
  Cu <- Cu[-which(Cu[,1] == s$gf[Cu[,2]]),]
  Cu <- Cu[-which(Cu[,2] == s$gf[Cu[,1]]),]
  ## Cu <- Cu[!((Cu[,1] %in% s$gf) & (Cu[,2] %in% s$gf)),] 
  l <- norm(P[Cu[,1],] - P[Cu[,2],])
}

## Check there are no zero-length lines
if (any(l==0)) {
  print("WARNING: zero-length lines")
}

## Add the new points to the correspondances vector
h <- c(s$h, (length(s$h)+1):nrow(P))

## Swap orientation of triangles which have clockwise orientation
a.signed <- tri.area.signed(P, T)
T[a.signed<0,c(2,3)] <- T[a.signed<0,c(3,2)]

## Estimate the area. It's roughly equal to the number of remaining points
## times the area of the rhomboid.
## area <- nrow(T) * L^2 * sqrt(3)/2
a <- abs(a.signed)
A <- sum(a)

## Find lengths of connections
Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
Cu <- Unique(Cu, TRUE)
L <- sqrt(apply((P[Cu[,1],] - P[Cu[,2],])^2, 1, sum))

## Plotting
plot(P)
trimesh(T, P, col="black")
lines(Po)

## Translation into unique points
u <- unique(h)
ht <- c()
for (i in 1:length(h)) {
  ht[i] <- which(u == h[i])
}
Tt  <- matrix(ht[T], ncol=3)

## Determine correspondance of line edges
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
U <- unique(H)
Cut <- Cut[U,]
Ht <- c()
for (i in 1:length(H)) {
  Ht[i] <- which(U == H[i])
}
Lt <- c()
for (k in 1:length(U)) {
  is <- which(Ht == k)
  Lt[k] <- mean(L[is])
}

## Cut <- rbind(Tt[,1:2], Tt[,2:3], Tt[,c(3,1)])
## Cut <- Unique(Cut, TRUE)
Pt  <- P[u,]
Rsett <- unique(ht[s$Rset])
Ct <- rbind(Cut, Cut[,2:1])
## Matrix to map line segments onto the points they link
Bt <- matrix(0, nrow(Pt), nrow(Ct))
for (i in 1:nrow(Ct)) {
  Bt[Ct[i,1],i] <- 1
}

## Plot stiched retina in 2D (messy)
trimesh(Tt, Pt, col="black")

N <- nrow(Pt)
Nphi <- N - length(Rsett)

## From this we can infer what the radius should be from the formula
## for the area of a sphere which is cut off at a lattitude of phi0
## area = 2 * PI * R^2 * (sin(phi0)+1)
phi0 <- 50*pi/180
R <- sqrt(A/(2*pi*(sin(phi0)+1)))

## Now assign each point to a location in the phi, lambda coordinates
## Shift coordinates to rough centre of grid
x <- Pt[,1] - mean(Pt[,1]) 
y <- Pt[,2] - mean(Pt[,2]) 
phi <- -pi/2 + sqrt(x^2 + y^2)/(R)
phi[Rsett] <- phi0
lambda <- atan2(y, x)

## Initial plot in 3D space
plot.retina(phi, lambda, R, Tt)

