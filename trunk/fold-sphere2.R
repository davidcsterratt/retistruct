## Modulus function, that returns
## i , if i <= N
## i - N, if i > N
Mod <- function(i, N) {
  return((i - 1) %% N + 1)
}

## Return sequence of indicies in path between i and j, governed by
## pointer vector p
path <- function(i, j, p) {
  if (i == j) {
    return(j)
  } else {
    return(c(i, path(p[i], j, p)))
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
  f <- Mod((1:N) + 1, N) # c(2:N, 1)
  b <- Mod((1:N) - 1, N) #c(N, 1:(N-1))
  
  ## Adjust forward and backward pointers appropriately
  f[VB] <- VF
  b[VF] <- VB
  
  ## Create initial set of correspondances - store as matrix
  C <- cbind(VF, VB)

  ## Initialise the set of points in the rim
  Rset <- 1:N
  
  ## Do not include the verticies in the tears
  TFset <- list()
  TBset <- list()
  ## Create sets of points for each tear and remove these points from
  ## the rim set
  for (j in 1:nrow(T)) {
    TFset[[j]] <- Mod(path(A[j], VF[j], f), N)
    TBset[[j]] <- Mod(path(A[j], VB[j], b), N)
    Rset <- setdiff(Rset, setdiff(TFset[[j]], VF[j]))
    Rset <- setdiff(Rset, setdiff(TBset[[j]], VB[j]))
  }

  return(list(Rset=Rset,
              VF=VF, VB=VB, A=A,
              TFset=TFset, TBset=TBset,
              C=C))
}

## Read in and curate data to give edge.path and map
source("M634-4.R")

P <- edge.path[-nrow(edge.path),1:2]

s <- stitch.retina(P, tearmat)

plot(P, pch=".", cex=4)
points(P[s$Rset,], col="red", pch=".", cex=4)
points(P[s$VF,], col="purple", pch="+")
points(P[s$VB,], col="blue", pch="+")
points(P[s$A, ], col="cyan", pch="+")
for (TFset in s$TFset) {
  lines(P[TFset,], col="purple")
}
for (TBset in s$TBset) {
  lines(P[TBset,], col="blue")
}
for (j in 1:nrow(s$C)) {
  lines(P[s$C[j,],], col="orange")
}
