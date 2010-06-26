source("cluster-analysis.R")
source("common.R")

join.segments <- function(segs) {
  map <- matrix(0, 0, 2)
  for (i in 1:length(segs)) { 
    map <- rbind(map, segs[[i]])
  }
  return(map) 
}

norm2 <- function(v) {
  return(sqrt(as.vector(t(v) %*% v)))
}  

## Find angle between two vectors
find.angle <- function(v1, v2) {
  sp <- (t(v1) %*% v2) / (norm2(v1) * norm2(v2))
  ##  print(sp)
  return(acos(sp))
}

find.corners <- function(map) {
  angle <- rep(NA, nrow(map))
  for (i in 2:(nrow(map)-1)) {
    v1 <- map[i,] - map[i-1,]
    ## print(v1)
    v2 <- map[i+1,] - map[i,]
    ## print(v2)
    angle[i] <-  find.angle(v1, v2) 
  }
  return(angle)
}

plot.triangle <- function(map, inds) {
  i <- inds[c(1:3,1)]
  lines(map[i,1], map[i,2])
}

## Compute area of triangle defined by corners x1, x2 & x3
tri.area <- function(x1, x2, x3) {
  a <- x2-x1
  b <- x3-x1
  area <- 0.5 * sqrt(t(a) %*% a * t(b) %*% b - (t(a) %*% b)^2)
  return(as.vector(area))
}

## Plot a mesh of radial and tangential lines 
plot.mesh <- function(gmap) {
  ## Plot the tangential lines
  for (j in js) {
    inds <- (gmap[,"j"] == j)
    polygon(gmap[inds,"X"], gmap[inds,"Y"])
  }
  ## Plot the radial lines
  for (i in is) {
    inds <- (gmap[,"i"] == i)
    lines(gmap[inds,"X"], gmap[inds,"Y"])
  }
}

## Plot a mesh of radial and tangential lines 
plot.mesh2 <- function(P, L) {
  for (i in 1:nrow(P)) {
    js <- which(!is.na(L[i,]))
    for (j in js) {
      lines(P[c(i,j),1], P[c(i,j),2])
    }
  }
}

## Plot a mesh of radial and tangential lines 
plot.mesh3 <- function(P, L) {
  for (i in 1:ncol(L)) {
    js <- which(L[,i] != 0)
    for (j in js) {
      lines(P[c(i,j),1], P[c(i,j),2])
    }
  }
}


## matrix that has two stripes on either side of the diagonal
stripe.matrix <- function(n) {
  M <- matrix(0, n, n)
  diag(M[2:n,1:(n-1)]) <- 1
  M[1,n] <- 1
  M <- M + t(M)
  return(M)
}

## which.max.matrix - find row and column of matrix
which.max.matrix <- function(M) {
  ## Find the maxium within each row
  row.maxs <-  apply(M, 1, max)
  ## Find which row has the biggest max
  i <- which.max(row.maxs)
  ## Find the location of the max in that row
  j <- which.max(M[i,])
  return(c(i,j))
}

highlight.possible.corners <- function(edge) {
  angles <- find.corners(edge) * 180/pi

  ## Find the distance of each point from the centroid of the data
  cent <- apply(edge, 2, mean)
  dist <- sqrt(apply((t(edge) - cent)^2, 2, sum))

  ## Heuristics for finding corners
  P.corner <- exp(-0.5*((angles-90)/30)^2) * (dist/max(dist))^1

  palette(rainbow(100))

  ## points(edge[,1], edge[,2], col=dist/max(dist)*100, pch=20)
  ## s <- sort(dist, index.return=TRUE, decreasing=TRUE)
  ## inds <- s$ix[1:6]

  points(edge[,1], edge[,2], col=P.corner*50, pch=20)
  s <- sort(P.corner, index.return=TRUE, decreasing=TRUE)
  inds <- s$ix[1:14] + 1

  points(edge[inds,1], edge[inds,2], pch="x")
  text(edge[inds,1], edge[inds,2], labels=format(inds,digits=2))
}

## Find length of all links between points in P, where links
## are specified by not NAs in L
link.lengths <- function(P, L) {
  l <- matrix(0, nrow(L), ncol(L))
  for (i in 1:nrow(L)) {
    js <- which(!is.na(L[i,]))
    l[i, js] <- sqrt(apply((t(P[js,]) - P[i,])^2, 2, sum))
  }
  return(l)
}

link.lengths2 <- function(P, L) {
  l <- matrix(0, nrow(L), ncol(L))
  for (j in 1:ncol(L)) {
    is <- which(L[,j] != 0)
##    print(length(is))
    if (length(is) == 1) {
      l[is, j] <- sqrt(sum((t(P[j,]) - P[is,])^2))
    } else {
      l[is, j] <- sqrt(apply((t(P[j,]) - P[is,])^2, 2, sum))
    }
  }
  return(l)
}

###
### Start of code
### 

###
### Start of the code proper
### 

## Read in and curate data to give edge.path and map
source("M634-4.R")

## Create a hemisphere
if (!exists("hs")) {
  source("hemisphere.R")
}


## Distribute points equally along edges
P <- find.points.in.path(seq(0, by=1/hs$nfix, len=hs$nfix), rim.path)

## In hs$P the first hs$nfix points are fixed, and the rest are
## free
M <- hs$nfix                            # number of fixed points
m <- nrow(hs$P) - M                     # number of free points
Cu <- hs$Cu
Ls <- hs$Ls
PQtri <- hs$Ptri

## Create A and B matricies
A <- matrix(0, m, m)
B <- matrix(0, m, M)

for (k in 1:nrow(Cu)) {
  if ((Cu[k,1] > M) && (Cu[k,2] > M)) {
    ## Connection between two free points
    i <- Cu[k,1] - M
    j <- Cu[k,2] - M
    A[i,j] <- 1/Ls[k]
    A[j,i] <- 1/Ls[k]
  }
  if ((Cu[k,1] >  M) && (Cu[k,2] <= M)) {
    ## Connection between fixed and free points
    i <- Cu[k,1] - M
    j <- Cu[k,2]
    B[i,j] <- 1/Ls[k]
  }
  if ((Cu[k,1] <= M) && (Cu[k,2] > M)) {
    ## Connection between fixed and free points
    i <- Cu[k,2] - M 
    j <- Cu[k,1]
    B[i,j] <- 1/Ls[k]
  }
}
for (iter in 1:1) {
  ## Convert to proximity matricies A and B
  ## Points 1:((N-1)*M+1) are variable points
  ## Points ((N-1)*M+2):(N*M+1) are fixed
  ## Diagonal matrix of row sums of [ A B ]
  D <- diag(apply(cbind(A, 2*B), 1, sum))

  ## The matrix Q is an n by 2 matrix comprising the row vectors of the
  ## solution points
  Q <- 2 * solve(D - A) %*% B %*% P
  
  ## Plot the outline, highlighting the rim in red
  plot.map(map, seginfo=FALSE)
  
  ## lines(rim.path[,"X"], rim.path[,"Y"], col="red", lwd=2)
  
  points(P, pch=16)
  points(Q, pch=16, col="gray")

  ## Plot the triangles
  PQ <- rbind(P, Q)
  ## trimesh(PQtri, PQ, col="gray", add=TRUE)
  
  ls <- sqrt(apply((R[Cu[,1],] - R[Cu[,2],])^2, 1, sum))
  strain <- ls/mean(ls)/(Ls/mean(Ls))
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  segments(PQ[Cu[,1],1], PQ[Cu[,1],2],
           PQ[Cu[,2],1], PQ[Cu[,2],2], col=(-log(strain) * 30 + 35))

  
  ## plot.mesh2(R, L)

  ## plot.mesh3(rbind(Q, P, S), cbind(A, B, C))
  ## Next step: what happens if the longest links are removed?
  ## Need to find the length of all the links
  ## l <- link.lengths(R, L[1:((N-1)*M)+1,])
  ## m.i <- which.max.matrix(l)
  ## L[m.i[1], m.i[2]] <- NA
  ## L[m.i[2], m.i[1]] <- NA
  ## C[which.max(l)] <- 0
}
