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
      }

      n <- nrow(P)
      h[i] <- n

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
      }
      n <- nrow(P)
      h[i] <- n
      
      ## print(paste("n =", n, "; k =", k, "; k0 =", k0,
      ##            "; gf[", n,  "] =", gf[n], "; gb[", n,  "] =", gb[n],
      ##            "; gf[", k0, "] =", gf[k0], "; gb[", k, "] =", gb[k]))
    }
  }
    
  return(list(Rset=Rset,
              VF=VF, VB=VB, A=A,
              TFset=TFset, TBset=TBset,
              P=P, h=h, hf=hf, hb=hb,
              gf=gf, gb=gb))
}

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
Pt <- delaunayn(P)
Pc <- (P[Pt[,1],] + P[Pt[,2],] + P[Pt[,3],])/3 # Centres
Pt <- Pt[point.in.polygon(Pc[,1], Pc[,2], Po[,1], Po[,2])==1,]

trimesh(Pt, P, col="gray", add=TRUE)

## Find lines which join non-adjacent parts of the outline
Cu <- rbind(Pt[,1:2], Pt[,2:3], Pt[,c(3,1)])
Cu <- Unique(Cu, TRUE)

for (i in 1:nrow(Cu)) {
  C1 <- Cu[i,1]
  C2 <- Cu[i,2]
  if (all(Cu[i,] %in% s$Rset)) {
    if (!((C1 == s$gf[C2]) ||
          (C2 == s$gf[C1]))) {
      ## Find triangles containing the line
      segments(P[C1,1], P[C1,2], P[C2,1], P[C2,2], col="red")
      Tind <- which(apply(Pt, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
      ## print(Cu[i,])
      ## print(Tind)
      T1 <- setdiff(Pt[Tind[1],], Cu[i,])
      T2 <- setdiff(Pt[Tind[2],], Cu[i,])
      ## print(T1)
      ## print(paste(C1, C2, T1, T2))
      p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
      points(p[1], p[2], col="red")
      P <- rbind(P, p)
      n <- nrow(P)
      Pt[Tind[1],] <- c(n, C1, T1)
      Pt[Tind[2],] <- c(n, C1, T2)
      Pt <- rbind(Pt,
                  c(n, C2, T1),
                  c(n, C2, T2))
    }
  }
}


Cu <- rbind(Pt[,1:2], Pt[,2:3], Pt[,c(3,1)])
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
  Tind <- which(apply(Pt, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
  ##  print(Cu[i,])
  ##  print(Tind)
  T1 <- setdiff(Pt[Tind[1],], Cu[i,])
  T2 <- setdiff(Pt[Tind[2],], Cu[i,])
  ##  print(T1)
  ##  print(paste(C1, C2, T1, T2))

  p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
  ## points(p[1], p[2], col="red")
  P <- rbind(P, p)
  n <- nrow(P)
  Pt[Tind[1],] <- c(n, C1, T1)
  Pt[Tind[2],] <- c(n, C1, T2)
  Pt <- rbind(Pt,
              c(n, C2, T1),
              c(n, C2, T2))
  
  Cu <- rbind(Pt[,1:2], Pt[,2:3], Pt[,c(3,1)])
  Cu <- Unique(Cu, TRUE)

  print(length(which(Cu[,1] == s$gf[Cu[,2]])))
  print(length(which(Cu[,2] == s$gf[Cu[,1]])))
  Cu <- Cu[-which(Cu[,1] == s$gf[Cu[,2]]),]
  Cu <- Cu[-which(Cu[,2] == s$gf[Cu[,1]]),]
  ## Cu <- Cu[!((Cu[,1] %in% s$gf) & (Cu[,2] %in% s$gf)),] 
  l <- norm(P[Cu[,1],] - P[Cu[,2],])
}

plot(P)
trimesh(Pt, P, col="black")
lines(Po)

