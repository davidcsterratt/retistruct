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
  A  <- T[,1]                            # apicies of tears
  VB <- T[,2]                           # forward verticies
  VF <- T[,3]                           # backward verticies

  ## Create forward and backward pointers
  N <- nrow(P)                          # Number of points
  gf <- Mod((1:N) + 1, N) # c(2:N, 1)
  gb <- Mod((1:N) - 1, N) # c(N, 1:(N-1))
  
  ## Create initial set of correspondances
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

  ## Create copies of the pointers
  gfp <- gf
  gbp <- gb
  
  ## Iterate through the tears
  for (j in 1:nrow(T)) {
    ## Create sets of points for each tear and remove these points from
    ## the rim set
    TFset[[j]] <- Mod(path(A[j], VF[j], gf, hf), N)
    TBset[[j]] <- Mod(path(A[j], VB[j], gb, hb), N)
    Rset <- setdiff(Rset, setdiff(TFset[[j]], VF[j]))
    Rset <- setdiff(Rset, setdiff(TBset[[j]], VB[j]))

    ## Compute the total path length along each side of the tear
    Sf <- path.length(A[j], VF[j], gf, hf, P)
    Sb <- path.length(A[j], VB[j], gb, hb, P)
    ## print(paste("Sf", Sf))
    ## print(paste("Sb", Sb))

    ## For each point in the forward path, create one in the backwards
    ## path at the same fractional location
    ## print("TF")
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
      ## print(paste("sb", sb/Sb))
      ## print(paste("sb0", sb0/Sb))
      f <- (sf/Sf*Sb-sb0)/(sb-sb0)
      ## print(paste("f", f))
      p <- (1-f) * P[k0,] + f * P[k,]
      print(nrow(P))
      print(p)
      P <- rbind(P, p)
      print(nrow(P))

      ## Update forward and backward pointers
      n <- nrow(P)
      ## k0 <- gfp[k]
      gbp[n]      <- k
      gfp[n]      <- gfp[k] # k0
      gbp[gfp[k]] <- n
      gfp[k]      <- n
      print(paste("n =", n, "; k =", k, "; k0 =", k0,
                  "; gfp[", n, "] =", gfp[n], "; gbp[", n,  "] =", gbp[n],
                  "; gfp[", k, "] =", gfp[k], "; gbp[", k0, "] =", gbp[k0]))
      ## Update correspondences
      h[i] <- n
      hf[n] <- n
      hb[n] <- n
    }

    ## Go along backward path
    print("TB")
    for (i in setdiff(TBset[[j]], c(A[j], VB[j]))) {
      sb <- path.length(A[j], i, gb, hb, P)
      ## print(paste("i", i, "sb", sb/Sb))
      ## print(TFset[[j]])
      for (k in TFset[[j]]) {
        sf <- path.length(A[j], k, gfp, hb, P)
        ## print(paste("k", k, "; sf/Sf", sf/Sf))
        if (sf/Sf > sb/Sb) {
          break;
        }
        sf0 <- sf
        k0 <- k
      }
      ## print(paste("sf", sf/Sf))
      ## print(paste("sf0", sf0/Sf))
      f <- (sb/Sb*Sf-sf0)/(sf-sf0)
      ## print(paste("f", f))
      p <- (1-f) * P[k0,] + f * P[k,]

      print(nrow(P))
      print(p)
      P <- rbind(P, p)
      print(nrow(P))      


      ## Update forward and backward pointers
      n <- nrow(P)
      gfp[n]  <- k
      gbp[n]  <- gbp[k]
      gfp[gbp[k]] <- n
      gbp[k]  <- n

      print(paste("n =", n, "; k =", k, "; k0 =", k0,
                  "; gfp[", n,  "] =", gfp[n], "; gbp[", n,  "] =", gbp[n],
                  "; gfp[", k0, "] =", gfp[k0], "; gbp[", k, "] =", gbp[k]))
      
      ## Update correspondences
      h[i] <- n
      hf[n] <- n
      hb[n] <- n
    }
  }
    
  return(list(Rset=Rset,
              VF=VF, VB=VB, A=A,
              TFset=TFset, TBset=TBset,
              P=P, h=h, hf=hf, hb=hb,
              gf=gfp, gb=gbp))
}

## Read in and curate data to give edge.path and map
source("M634-4.R")

P <- edge.path[-nrow(edge.path),1:2]

s <- stitch.retina(P, tearmat)


plot(s$P, pch=".", cex=4)
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
#    lines(P[c(j,s$hf[j]),], col="green")
  }
}

points(P[s$Rset,], col="red")


