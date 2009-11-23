require("sp")
require("geometry")

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

## Read in and curate data to give edge.path and map
source("M634-4.R")

P <- edge.path[-nrow(edge.path),1:2]

d <- 200                                 # minspacing
N <- 200                                 # number of points

xmin <- min(P[,1])
xmax <- max(P[,1])
ymin <- min(P[,2])
ymax <- max(P[,2])

Q <- matrix(0, 0, 2)

for (i in 1:N) {
  ntry <- 100
  while (ntry) {
    ## Generate candidate point
    C <- c(runif(1, xmin, xmax),
           runif(1, ymin, ymax))
    if (point.in.polygon(C[1], C[2], P[,1], P[,2])) {
      if ((nrow(Q) == 0) || all(norm(t(t(rbind(Q, P)) - C)) > d)) {
        Q <- rbind(Q, C)
        print(ntry)
        break;
      }
    }
    ntry <- ntry - 1
  }
}

plot(P, type="l")
points(Q, col="red")

## Triangulate
S <- rbind(P, Q)
St <- delaunayn(S)

## Remove triangles outwith the retinal outline
Sc <- (S[St[,1],] + S[St[,2],] + S[St[,3],])/3 # Centres
St <- St[point.in.polygon(Sc[,1], Sc[,2], P[,1], P[,2])==1,]

trimesh(St, S, col="gray", add=TRUE)
