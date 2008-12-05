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
    polygon(gmap[inds,"X"], gmap[inds,"Y"])
  }
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
theta.max.deg <- 140 
dtheta.deg <- 10 
dphi.deg   <- 10 
r <- 1000                                # radius

## Make some indicies
is <- 1:(360/dphi.deg)
js <- 1:(theta.max.deg/dtheta.deg)

## Convert to radians
dphi <-   dphi.deg * pi/180
dtheta <- dtheta.deg * pi/180
theta.max <- theta.max.deg * pi/180

## The matricies X and Y specify the mapping from point i,j in polar space
## to cartesian coordinates
X <- matrix(NA, length(is), length(js))
Y <- matrix(NA, length(is), length(js))

## Put these into gmap structure
gmap <- expand.grid(i=is, j=js)
gmap <- cbind(gmap,
              phi=  gmap[,"i"]*dphi,
              theta=gmap[,"j"]*dtheta)
gmap <- cbind(gmap, X=r*gmap[,"theta"]*cos(gmap[,"phi"]))
gmap <- cbind(gmap, Y=r*gmap[,"theta"]*sin(gmap[,"phi"]))

## Find the neigbours and put them in columns
## - nl - tangential left (in direction of increasing phi)
## - nr - tangential right (in direction of decreasing phi)
## - nu - radial up (in direction of increasing theta)
## - nd - radial down (in direction of decreasing theta)
M <- length(is)
N <- length(js)
gmap[,"nl"] <- (1:nrow(gmap))+ 1 - M * (((1:nrow(gmap)) %% M) == 0)
gmap[,"nr"] <- (1:nrow(gmap))- 1 + M * (((1:nrow(gmap)) %% M) == 1)
gmap[,"nu"] <- (1:nrow(gmap))+ M
gmap[gmap[,"nu"]>nrow(gmap),"nu"] <- NA
gmap[,"nd"] <- (1:nrow(gmap))- M
gmap[gmap[,"nd"]<1,"nd"] <- NA

## Computes the distance between grid points if they were
## spread over the surface of a sphere
## Find distance to left neigbours
gmap[,"Ll"] <- r * dphi * sin(gmap[,"theta"])
gmap[,"Lr"] <- r * dphi * sin(gmap[,"theta"])
gmap[,"Lu"] <- r * dtheta
gmap[is.na(gmap[,"nu"]),"Lu"] <- NA
gmap[,"Ld"] <- r * dtheta
gmap[is.na(gmap[,"nd"]),"Ld"] <- NA


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
## Find the distance along the rim at which the points shoudl be
s.P <- seq(0, by=max(rim.dist[,"s"])/36, len=36)

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

## Convert to proximity matricies A and B
## Points 1:((N-1)*M) are variable points
## Points ((N-1)*M+1):(N*M) are fixed
C <- matrix(0, (N-1)*M, N*M)
mgmap <- as.matrix(gmap)
for (k in 1:(((N-1)*M))) {
  C[k,na.omit(mgmap[k,c("nl","nr","nu","nd")])] <- 1/na.omit(mgmap[k,c("Ll","Lr","Lu","Ld")])
}
A <- C[,1:((N-1)*M)]
B <- C[,((N-1)*M+1):(N*M)]
## Diagonal matrix of row sums of [ A B ]
D <- diag(apply(cbind(C), 1, sum))

## Initial test with P lying on circle
## P <-  mgmap[((N-1)*M+1):(N*M),c("X","Y")]

## The matrix Q is an n by 2 matrix comprising the row vectors of the
## solution points
Q <- solve(D - A) %*% B %*% P

## Put new distances back in map for plotting
mgmap[,c("X","Y")] <- rbind(Q, P)

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

plot.mesh(mgmap)

## To fix:
## 1. get evenly spaced points on rim
## 2. get centre point
## 3. get plotting function working
