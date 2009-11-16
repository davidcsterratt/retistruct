## Modulus function, that returns
## i , if i <= N
## i - N, if i > N
Mod <- function(i, N) {
  return((i - 1) %% N + 1)
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
  A <- T[,1]                            # apicies of tears
  VB <- T[,2]                           # forward verticies
  VF <- T[,3]                           # backward verticies

  ## Create forward and backward pointers
  N <- nrow(P)                          # Number of points
  gf <- Mod((1:N) + 1, N) # c(2:N, 1)
  gb <- Mod((1:N) - 1, N) # c(N, 1:(N-1))
  
  ## Adjust forward and backward pointers appropriately
  ## f[VB] <- VF
  ## b[VF] <- VB
  
  ## Create initial set of correspondances - store as
  hf <- 1:N
  hb <- 1:N
  h <- 1:N
  hf[VB] <- VF
  hb[VF] <- VB

  ## Initialise the set of points in the rim
  Rset <- 1:N
  
  ## Do not include the verticies in the tears
  TFset <- list()
  TBset <- list()
  ## Create sets of points for each tear and remove these points from
  ## the rim set
  for (j in 1:nrow(T)) {
    TFset[[j]] <- Mod(path(A[j], VF[j], gf, hf), N)
    TBset[[j]] <- Mod(path(A[j], VB[j], gb, hb), N)
    Rset <- setdiff(Rset, setdiff(TFset[[j]], VF[j]))
    Rset <- setdiff(Rset, setdiff(TBset[[j]], VB[j]))

    Sf <- path.length(A[j], VF[j], gf, hf, P)
    Sb <- path.length(A[j], VB[j], gb, hb, P)
    print(paste("Sf", Sf))
    print(paste("Sb", Sb))
    ## Go along forward path
    print("TF")
    for (i in setdiff(TFset[[j]], c(A[j], VF[j]))) {
      sf <- path.length(A[j], i, gf, hf, P)
      print(paste("sf", sf/Sf))
      print(TBset[[j]])
      for (k in TBset[[j]]) {
        sb <- path.length(A[j], k, gb, hb, P)
        print(paste("k", k, "; sb/Sb", sb/Sb))
        if (sb/Sb > sf/Sf) {
          break;
        }
        sb0 <- sb
      }
      print(paste("sb", sb/Sb))
      print(paste("sb0", sb0/Sb))
      f <- (sf/Sf*Sb-sb0)/(sb-sb0)
      print(paste("f", f))
      p <- (1-f) * P[gf[k],] + f * P[k,]
      P <- rbind(P, p)
      h <- c(h, i)
    }

    ## Go along backward path
    print("TB")
    for (i in setdiff(TBset[[j]], c(A[j], VB[j]))) {
      sb <- path.length(A[j], i, gb, hb, P)
      print(paste("i", i, "sb", sb/Sb))
      print(TFset[[j]])
      for (k in TFset[[j]]) {
        sf <- path.length(A[j], k, gf, hb, P)
        print(paste("k", k, "; sf/Sf", sf/Sf))
        if (sf/Sf > sb/Sb) {
          break;
        }
        sf0 <- sf
      }
      print(paste("sf", sf/Sf))
      print(paste("sf0", sf0/Sf))
      f <- (sb/Sb*Sf-sf0)/(sf-sf0)
      print(paste("f", f))
      p <- (1-f) * P[gb[k],] + f * P[k,]
      P <- rbind(P, p)
      h <- c(h, i)
    }

    ## print(Sb)
  }
    
  return(list(Rset=Rset,
              VF=VF, VB=VB, A=A,
              TFset=TFset, TBset=TBset,
              P=P, h=h, hf=hf, hb=hb))
}

## Read in and curate data to give edge.path and map
source("M634-4.R")

P <- edge.path[-nrow(edge.path),1:2]

s <- stitch.retina(P, tearmat)

plot(s$P, pch=".", cex=4)

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
  if (s$h != j) {
    lines(s$P[c(j,s$h[j]),], col="orange")
  }
}

for (j in 1:length(s$hf)) {
  if (s$hf != j) {
    lines(P[c(j,s$hf[j]),], col="green")
  }
}

points(P[s$Rset,], col="red")
