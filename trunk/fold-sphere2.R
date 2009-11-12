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
  N <- nrow(P)                          # Number of points

  ## Create forward and backward pointers
  f <- c(2:N, 1)
  b <- c(N, 1:(N-1))

  ## Extract information from tear matrix
  A <- T[,1]                            # apicies of tears
  VB <- T[,2]                           # forward verticies
  VF <- T[,3]                           # backward verticies

  ## Create initial set of correspondances - store as matrix
  C <- cbind(VF, VB)

  ## Initialise the set of points in the rim
  Rset <- 1:N
  
  ## Create sets of points for each tear and remove these points from
  ## the rim set
  for (j in 1:nrow(T)) {
    print(paste(f[VB[j]], b[VF[j]]))
    Tset[[j]] <- Mod(path(f[VB[j]], b[VF[j]], f), N)
    Rset <- setdiff(Rset, Tset[[j]])
  }

  ## Adjust forward and backward pointers appropriately
  f[VB] <- VF
  b[VF] <- VB

  ## Do not include the verticies in the tears
  Tset <- list()
  
  return(list(Rset=Rset, VF=VF, VB=VB, A=A))
}

## Read in and curate data to give edge.path and map
source("M634-4.R")

P <- edge.path[-nrow(edge.path),1:2]



s <- stitch.retina(P, tearmat)

plot(P, pch=".", cex=4)
points(P[s$Rset,], col="red", pch=".", cex=4)
points(P[s$VF,], col="purple", pch="+")
points(P[s$VB,], col="purple", pch="+")
points(P[s$A, ], col="blue", pch="+")

