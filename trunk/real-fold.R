source("../../data/Anatomy/cluster-analysis.R")

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

## Check whether line from x1 to y1 and from x2 to y2 intersect
check.intersection <- function(x1, y1, x2, y2) {
  m <- cbind(y1-x1, -y2+x2)
  if (!(det(m) == 0)) {
    lambda <- solve(m) %*% (x2-x1)
    if (all((lambda<1) & (lambda>0))) {
      return(lambda)
    }
  }
  return(FALSE)
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

## Read in data
map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))

## Corner analysis
segs <- map.to.segments(map)
## Some curation is required here
segs4 <- segs[[4]][23:1,]
edge <- join.segments(list(segs[[1]], segs[[3]], segs4))
##map <- segs[[1]]
angles <- find.corners(edge) * 180/pi

## Find the distance of each point from the centroid of the data
cent <- apply(edge, 2, mean)
dist <- sqrt(apply((t(edge) - cent)^2, 2, sum))

## Make a grid of indicies
## Specify maxium elevation and grid spacing in degrees
theta.max.deg <- 120
dtheta.deg <- 10 
dphi.deg   <- 10 
r <- 2000                                # radius

## Make some indicies
is <- c(1:(360/dphi.deg))
js <- c(1:(theta.max.deg/dtheta.deg))

M <- length(is)
N <- length(js)

## Convert to radians
dphi <-   dphi.deg * pi/180
dtheta <- dtheta.deg * pi/180
theta.max <- theta.max.deg * pi/180

## Put these into gmap structure
gmap <- rbind(c(i=0, j=0), expand.grid(i=is, j=js))
gmap <- cbind(gmap,
              phi=  gmap[,"i"]*dphi,
              theta=gmap[,"j"]*dtheta)
gmap <- cbind(gmap, X=r*gmap[,"theta"]*cos(gmap[,"phi"]))
gmap <- cbind(gmap, Y=r*gmap[,"theta"]*sin(gmap[,"phi"]))

## Find the matrix L of lengths to neigbours
## No connection is indicated by NA
L <- matrix(NA, 1+N*M, 1+N*M)
## First block: connections from (0,0) to neigbours
L[1, 1+(1:M)] <- dtheta

## Subsequent blocks, grouped by increment in theta
for (j in js) {
  rows <- 1 + (j-1)*M + 1:M
  # connections from (dphi * i, dtheta * j) to neigbours
  L[rows,rows] <- stripe.matrix(M) * dphi * sin(gmap[rows,"theta"])
  if (j == 1) {
    ## connections from the pole
    L[rows,1] <- dtheta
  } else {
    L[rows,rows-M] <-  diag(M) * dtheta
  }
  if (j < N) {
    L[rows,rows+M] <- diag(M) * dtheta
  }
}
L[L==0] <- NA
L <- L * r

## Heuristics for finding corners
## P.corner <- exp(-0.5*((angles-90)/30)^2) * (dist/max(dist))^1

## points(edge[,1], edge[,2], col=dist/max(dist)*100, pch=20)
## s <- sort(dist, index.return=TRUE, decreasing=TRUE)
## inds <- s$ix[1:6]

## points(edge[,1], edge[,2], col=P.corner*50, pch=20)
## s <- sort(P.corner, index.return=TRUE, decreasing=TRUE)
## inds <- s$ix[1:14] + 1

## points(edge[inds,1], edge[inds,2], pch="x")
## text(edge[inds,1], edge[inds,2], labels=format(inds,digits=2))

## Plot the grid
## plot(NA,NA,xlim=c(-pi, pi), ylim=c(-pi, pi))
## plot.mesh(gmap)

## Hand pick corners

c(16, 31)
c(53, 65)
c(83, 102)

## Define rim as list of line segments
rim <- list(edge[16:31,],
            edge[53:65,],
            edge[83:94,])

## Find distance along edges
## Need a matrix to store segment ind and ind with in segment
rim.dist <- matrix(0, 0, 4)
s.cum <- 0
for (i in 1:length(rim)) {
  seg <- rim[[i]]
  v <- diff(seg)
  l <- sqrt(apply(v^2, 1, sum))
  s <- s.cum + cumsum(l)
  print(s)
  s.cum <- s[length(s)]
  rim.dist <- rbind(rim.dist,
                    cbind(s,
                          l,
                          seg=rep(i, length(s)),
                          i=1:length(s)))
}
rim.dist[,"s"] <- c(0,rim.dist[-nrow(rim.dist),"s"])

## Distribute points equally along edges
## Find the distance along the rim at which the points should be
## The last "s" element of rim.dist is actually the distance at
## the *start* of the last line segment, so we need to add on "l".
## This is not ideal, but is necessetated by the behaviour of findInterval()
s.P <- seq(0, by=sum(rim.dist[nrow(rim.dist),c("s","l")])/36, len=36)

## Find the ind of the line segment in which they occur
inds <- findInterval(s.P, rim.dist[,"s"]) 

## Find the distance along the segment
ds.P <- s.P - rim.dist[inds,"s"] 

## Now find the points themselves
P <- matrix(NA, length(s.P), 2)
f <- ds.P/rim.dist[inds,"l"]
for (k in 1:length(s.P)) {
  i <-  rim.dist[inds[k],"i"]
  si <- rim.dist[inds[k],"seg"]
  ends <- rim[[si]][c(i,i+1),]
  P[k,] <- (1 - f[k]) * ends[1,] + f[k] * ends[2,]
##  print(norm2(as.vector(v)))
##  print(ds.P[k]/rim.dist[k,"l"])
}

## Fitting of mesh to data


for (iter in 1:1) {
  ## Convert to proximity matricies A and B
## Points 1:((N-1)*M+1) are variable points
## Points ((N-1)*M+2):(N*M+1) are fixed
C <- cbind(1/L[1:((N-1)*M+1),1:((N-1)*M+1)],
           1/L[1:((N-1)*M+1),((N-1)*M+2):(N*M+1)])

C[is.na(C)] <- 0

A <- C[,1:((N-1)*M+1)]
B <- C[,((N-1)*M+2):(N*M+1)]
## Diagonal matrix of row sums of [ A B ]
D <- diag(apply(C, 1, sum))

## The matrix Q is an n by 2 matrix comprising the row vectors of the
## solution points
Q <- solve(D - A) %*% B %*% P

## Plotting

palette(rainbow(100))

plot.map(map, seginfo=FALSE)

## Highlight rim
for (i in 1:length(rim)) {
  seg <- rim[[i]]
  lines(seg[,1], seg[,2], col="red")
}

points(P[,1], P[,2], pch=16)
#text(P[,1], P[,2], labels=1:M)

R <- rbind(Q, P)
plot.mesh2(R, L)

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

## Next step: what happens if the longest links are removed?
## Need to find the length of all the links
l <- link.lengths(R, L[1:((N-1)*M)+1,])
m.i <- which.max.matrix(l)
L[m.i[1], m.i[2]] <- NA
L[m.i[2], m.i[1]] <- NA
C[which.max(l)] <- 0
}

## Now try improving on this initial guess by using the energy function
## (l - L)^2
E <- function(p, P, L) {
  ## Construct Q from p
  Q <- matrix(p, length(p)/2, 2)
  l <- link.lengths(rbind(Q, P), L[1:nrow(Q),])
  E <- sum(sum((L[1:nrow(Q),] - l)^2, na.rm=TRUE))
  print(E)
  return(E)
}

dE <- function(p, P, L) {
  ## Construct Q from p
  Q <- matrix(p, length(p)/2, 2)

  l <- link.lengths(rbind(Q, P), L[1:nrow(Q),])

  C <- 1 - L[1:nrow(Q),]/l
  C[is.na(C)] <- 0

  ## Diagonal matrix of row sums of [ A B ]
  D <- diag(apply(C, 1, sum))
  dEdQ <- D %*% Q - C %*% rbind(Q, P)
  return(as.vector(dEdQ))
}

opt <- optim(as.vector(Q), E, gr=dE, method="BFGS", L=L, P=P, control=list(maxit=100))
plot.map(map)
plot.mesh2(rbind(Q, P), L)
