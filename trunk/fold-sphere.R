source("cluster-analysis.R")
source("triangulate.R")
source("common.R")
source("dipole.R")
library("geometry")
source("tsearch.R")
library("gtools")
library("rgl")

## Function to return "signed area" of triangles on a plane
## given points P and a triangulation Pt. Positive sign
## indicates points are anticlockwise direction; negative indicates
## clockwise
tri.area.signed <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]
  AB <- cbind(B-A, 0)
  BC <- cbind(C-B, 0)
  return(0.5 * extprod3d(AB, BC)[,3])
}

## Given the link matrix C and edge.inds, the indicies
## of points on the edge, classify the links as either being
## edge links or centre links. 
classify.links <- function(C, edge.inds=NA) {
  edge.links <- c()
  cent.links <- 1:nrow(C)
  if(!any(is.na(edge.inds))) {
    edge.links <- which(apply(matrix(C %in% edge.inds, nrow(C), 2), 1, all))
    ## Remove edge.links between nodes that have more than 1 from and to links
    n.to   <- hist(C[edge.links,1], breaks=c(-1, edge.inds+0.5), plot=FALSE)$counts
    n.from <- hist(C[edge.links,2], breaks=c(-1, edge.inds+0.5), plot=FALSE)$counts
    omit.to   <- which(n.to   == 0)
    omit.from <- which(n.from == 0)
    print(omit.to)
    print(omit.from)
    dup.to   <- which(n.to   > 1)
    dup.from <- which(n.from > 1)
    dup.links <- which((C[,1] %in% dup.from) & (C[,2] %in% dup.to))
    print(dup.to)
    print(dup.from)
    print(dup.links)
    edge.links <- setdiff(edge.links, dup.links)
    ##    C[edge.links,1] 
    
    cent.links <- setdiff(cent.links, edge.links)
  }
  return(list(edge=edge.links, cent=cent.links, dup.to=dup.to, dup.from=dup.from))
}

## Function to plot the retina in spherical coordinates
## phi - lattitude of points
## lambda - longitude of points
## R - radius of sphere
## C - connection table
## Ct - triagulated connection table
## sys - data frame containing positions of points
plot.retina <- function(phi, lambda, R, C,
                        Ct=matrix(NA, 0, 3), ts.red, ts.green,
                        edge.inds=NA) {
  ## Now plot this in 3D space....
  x <- R*cos(phi)*cos(lambda) 
  y <- R*cos(phi)*sin(lambda)
  z <- R*sin(phi)
  P <- cbind(x, y, z)
  rgl.clear()
  rgl.bg(color="white")

  ## Plot links, highlighting the edges in black
  links <- classify.links(C, edge.inds)
  ## If edge connections are specified, find the connections
  ## which lie on the edge and those which lie on the centre
  if (!any(is.na(edge.inds))) {

    ## Plot the edge links
    segments3d(rbind(x[C[links$edge,1]],x[C[links$edge,2]]),
               rbind(y[C[links$edge,1]],y[C[links$edge,2]]),
               rbind(z[C[links$edge,1]],z[C[links$edge,2]]),
               xlab="x", color="black", size=2)
  }
  ## Plot the centre links
  segments3d(rbind(x[C[links$cent,1]],x[C[links$cent,2]]),
             rbind(y[C[links$cent,1]],y[C[links$cent,2]]),
             rbind(z[C[links$cent,1]],z[C[links$cent,2]]),
             xlab="x", color="grey")

  ## In order to plot the points, we need to know in which triangle they
  ## lie and where in that triangle they are. The first job is to create
  ## a triangulation, and then we can use tsearch to find the identity of the
  ## triangles and the location in the triangles.

  ## Now plot points
  Pi.cart <- matrix(0, 0, 3)
  for(i in 1:(dim(ts.red$p)[1])) {
    Pi.cart <- rbind(Pi.cart, bary2cart(P[Ct[ts.red$idx[i],],], ts.red$p[i,]))
  }
  points3d(Pi.cart[,1], Pi.cart[,2], Pi.cart[,3], color="red", size=5)
  Pi.cart <- matrix(0, 0, 3)
  for(i in 1:(dim(ts.green$p)[1])) {
    Pi.cart <- rbind(Pi.cart, bary2cart(P[Ct[ts.green$idx[i],],], ts.green$p[i,]))
  }
  points3d(Pi.cart[,1], Pi.cart[,2], Pi.cart[,3], color="green", size=5)

  ## Plot any flipped triangles
  ## First find verticies
  P1 <- P[Ct[,1],]
  P2 <- P[Ct[,2],]
  P3 <- P[Ct[,3],]
  cents <- (P1 + P2 + P3)/3
  normals <- 0.5 * extprod3d(P2 - P1, P3 - P2)
  areas <- apply(normals^2, 1, sum)
##  print(cents)
##  print(areas)
  flipped <- (-dot(cents, normals) < 0)
  points3d(cents[flipped,1], cents[flipped,2], cents[flipped,3], col="blue", size=5)
  h <- 2 * cbind(areas/sqrt(sum((P2 - P3)^2)),
                 areas/sqrt(sum((P3 - P1)^2)),
                 areas/sqrt(sum((P2 - P1)^2)))
  closest <- apply(h, 1, which.min)
  cp <- Ct[(nrow(Ct)*(closest-1))+1:nrow(Ct)] ## closest point
  ## print(Ct[,closest])
##  print(h1[flipped])
##  points3d(cents[flipped,1], cents[flipped,2], cents[flipped,3], col="blue", size=5)
points3d(P[cp[flipped],1], P[cp[flipped],2], P[cp[flipped],3], col="blue", size=5)
##  rgl.triangles(x[Ct[flipped,]], y[Ct[flipped,]], z[Ct[flipped,]], col="gray")
}

## Formula for central angle
central.angle <- function(phi1, lambda1, phi2, lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)))
}

## Given a matrix of points Q in trignometric normalised coordinates,
## return the cartesian coordinates P, given the grid offset P0 and
## the grid spacing L
trig.to.cart <- function(Q, P0, L) {
  M <- L * matrix(c(1, -1/2, 0, sqrt(3)/2), 2, 2)
  P <- Q %*% M + matrix(P0, nrow(Q), 2, byrow=TRUE)
  colnames(P) <- c("X", "Y")
  return(P)
}

## Given a matrix of points P in cartesian coordinates,
## return the trignometric coordinates Q, given the grid offset P0 and
## the grid spacing L
cart.to.trig <- function(P, P0, L) {
  if (is.vector(P)) {
    P <- t(as.matrix(P))
  }
  M <- 1/L * matrix(c(1, 1/sqrt(3), 0, 2/sqrt(3)), 2, 2)  
  return(t(t(P) - P0) %*% M)
}

## Given a matrix of points P in cartesian coordinates,
## return the trignometric coordinates Q, given the grid offset P0 and
## the grid spacing L
xy.to.uw <- function(P, P0, L) {
  if (is.vector(P)) {
    P <- t(as.matrix(P))
  }
  M <- 1/L * matrix(c(1, -1/sqrt(3), 0, -2/sqrt(3)), 2, 2)  
  return(t(t(P) - P0) %*% M)
}

## Given a matrix of points P in cartesian coordinates,
## return the trignometric coordinates Q, given the grid offset P0 and
## the grid spacing L
uw.to.xy <- function(P, P0, L) {
  if (is.vector(P)) {
    P <- t(as.matrix(P))
  }
  M <- 1/L * matrix(c(1, -1/2, 0, -sqrt(3)/2), 2, 2)  
  return(t(t(P) - P0) %*% M)
}

tri.area <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]
  AB <- cbind(B-A, 0)
  AC <- cbind(C-A, 0)
  return(0.5 * abs(extprod3d(AB, AC)[,3]))
}


###
### Start of the code proper
### 

L <- 200                            # spacing between grid points

## Read in and curate data to give edge.path and map
source("M634-4.R")

## All the points on the edge path apart from the first, which is repeated by the last
P <- cbind(X=edge.path[-1,"X"], Y=edge.path[-1,"Y"])

## Find the bottom left corner of the box around the edge
P0   <- c(min(P[,"X"])-1, min(P[,"Y"])-1)
## Find the top right corner of the box around the edge
P1   <- c(max(P[,"X"])+1, max(P[,"Y"])+1)
## Find the top left corner of the box around the edge
P.TL <- c(min(P[,"X"])-1, max(P[,"Y"])+1)
## Find the bottom right corner of the box around the edge
P.BR <- c(max(P[,"X"])+1, min(P[,"Y"])-1)

## Convert these points to the trignometric coordinates, and
## find the maximum coordinates along the u v axes
Q0   <- cart.to.trig(P0, P0, L)
Q1   <- cart.to.trig(P1, P0, L)
R.TL <- xy.to.uw(P.TL, P.TL, L)
R.BR <- xy.to.uw(P.BR, P.TL, L)

umax <- ceiling(Q1[1])
vmax <- ceiling(Q1[2])
wmax <- ceiling(R.BR[2])

## Plot the outline of the flattened retina
plot(edge.path[,"X"], edge.path[,"Y"], xlab="x", ylab="y")
plot.sys.map(sys, map)
lines(edge.path[,"X"], edge.path[,"Y"],lwd=4)

## Plot the entire grid
## segments(G[Cu[,1],"x"], G[Cu[,1],"y"],
##         G[Cu[,2],"x"], G[Cu[,2],"y"])


## Now remove points that are outside the retina
## To do this, the intersections of each of the horizontal
## grid lines with the retina are determined
## Points outwith the range of these connections are discareded
## Points of grid
Q.grid <- matrix(0, 0, 2)
P.edge <- P
colnames(P.edge) <- c("X", "Y")
for (v in 0:vmax) {
  vfix <- create.path(list(trig.to.cart(matrix(c(0, v, umax, v), 2, 2, byrow=TRUE), P0, L)))
  ## This returns a list of intersection points
  ## There should be an even number of them
  ci <- check.intersection.paths(vfix, edge.path)

  if (!is.null(ci)) {
    ## Make sure points are in order of X so that there is no double-couting
    ## of points
    ci <- ci[order(ci[,"X"]),]
    for(i in 1:(nrow(ci)/2)) {
      cii <- ci[(i-1)*2+1:2,]
      ## points(cii[,"X"], cii[,"Y"], col="blue")
      ## Add all points between the intersection points
      Qlim <- cart.to.trig(cii[,c("X","Y")], P0, L)
      print(nrow(Q.grid))
      print(Qlim)
      Q.grid <- rbind(Q.grid, cbind(ceiling(min(Qlim[,1])):floor(max(Qlim[,1])), v))
      ## Add the intersection points of the line to the grid
      P.edge <- rbind(P.edge,cii[,c("X","Y")])
    }
  }
}
P.grid <- trig.to.cart(Q.grid, P0, L)
colnames(P.grid) <- c("X", "Y")

##Do the same thing along the w-axis
for (w in 0:wmax) {
  wfix <- create.path(list(uw.to.xy(matrix(c(w, 0, w, vmax), 2, 2, byrow=TRUE), P.TL, L)))
  ## This returns a list of intersection points
  ## There should be an even number of them
  ci <- check.intersection.paths(wfix, edge.path)

  if (!is.null(ci)) {
    for(i in 1:(nrow(ci)/2)) {
      cii <- ci[(i-1)*2+1:2,]
      ## Add the intersection points of the line to the grid
      P.edge <- rbind(P.edge, cii[,c("X","Y")])
    }
  }
}
## Do the same thing along the u-axis
for (u in 0:umax) {
  ufix <- create.path(list(trig.to.cart(matrix(c(u, 0, u, vmax), 2, 2, byrow=TRUE), P0, L)))
  ## This returns a list of intersection points
  ## There should be an even number of them
  ci <- check.intersection.paths(ufix, edge.path)

  if (!is.null(ci)) {
    for(i in 1:(nrow(ci)/2)) {
      cii <- ci[(i-1)*2+1:2,]
      ## Add the intersection points of the line to the grid
      P.edge <- rbind(P.edge, cii[,c("X","Y")])
    }
  }
}
##points(P.edge[,"X"], P.edge[,"Y"], col="blue")
##points(P.grid[,"X"], P.grid[,"Y"], col="orange")
## Full set of points
N.edge <- nrow(P.edge)
P <- rbind(P.edge, P.grid)
## P <- P.grid

## Find Delaunay triangulation of the points
Pt <- delaunayn(P)
areas.signed <- tri.area.signed(P, Pt)
areas <- abs(areas.signed)
## Swap orientation of triangles which have clockwise orientation
Pt[areas.signed<0,c(2,3)] <- Pt[areas.signed<0,c(3,2)]
## This gets rid of any very small triangles - it is a bit of a hack...
Pt    <- Pt[areas>0.1,]
areas <- areas[areas>0.1]

## Get rid of triangles whose centres are outside the retina
##for (i in 1:nrow(Pt)) {
## for (i in  which(apply(matrix(Pt %in% 1:nrow(P.edge), nrow(Pt), 3), 1, any))) {
for (i in  1:nrow(Pt)) {
  print(i)
  ## Find the centre of each triangle
  C <- apply(P[Pt[i,],], 2, mean)
  ## Does a line between C and P0 intersect the edge
  C.P0 <- create.path(list(rbind(C, P0)))
  ## This returns a list of intersection points
  ## There should be an even number of them
  ci <- check.intersection.paths(C.P0, edge.path)
  ## If there is an even number of intersections, then discard the triangle
  if (is.null(ci) || even(nrow(ci))) {
    Pt[i,] <- c(NA, NA, NA)
    areas[i] <- NA
  }
}
Pt <- na.omit(Pt)
areas <- na.omit(areas)

## Plot the triangles that remain after the pruning
trimesh(Pt, P, col="gray", add=TRUE)

## In order to plot the points, we need to know in which triangle they
## lie and where in that triangle they are. The first job is to create
## a triangulation, and then we can use tsearch to find the identity of the
## triangles and the location in the triangles.
P.red   <- cbind(na.omit(sys[,"XRED"]), na.omit(sys[, "YRED"]))
P.green <- cbind(na.omit(sys[,"XGREEN"]), na.omit(sys[, "YGREEN"]))
ts.red   <- tsearchn(P, Pt, P.red)
ts.green <- tsearchn(P, Pt, P.green)

## Now test this by replotting points using new coordinates
##for(i in 1:(dim(Pi)[1])) {
##  Pi.cart <- bary2cart(P[Pt[ts$idx[i],],], ts$p[i,])
##  points(Pi.cart[1], Pi.cart[2], pch="x")
##}

## Create the asymmetric connectivity matrix
Cu <- rbind(Pt[,1:2], Pt[,2:3], Pt[,c(3,1)])
##Cu <- t(apply(Cu, 1, sort))
##Cu <- Cu[!duplicated(Cu),]
Cu <- Unique(Cu, TRUE)

## C is the symmetric connectivity matrix
C <- rbind(Cu, Cu[,2:1])

## Matrix to map line segments onto the points they link
B <- matrix(0, nrow(P), nrow(C))
for (i in 1:nrow(C)) {
  B[C[i,1],i] <- 1
}

## Find lengths of connections
Ls <- sqrt(apply((P[Cu[,1],] - P[Cu[,2],])^2, 1, sum))

## Estimate the area. It's roughly equal to the number of remaining points
## times the area of the rhomboid.
## area <- nrow(Pt) * L^2 * sqrt(3)/2
area <- sum(areas)

## From this we can infer what the radius should be from the formula
## for the area of a sphere which is cut off at a lattitude of phi0
## area = 2 * PI * R^2 * (sin(phi0)+1)
phi0 <- 50*pi/180
R <- sqrt(area/(2*pi*(sin(phi0)+1)))

## Now assign each point to a location in the phi, lambda coordinates
## Shift coordinates to rough centre of grid
x <- P[,"X"] - mean(P[,"X"]) 
y <- P[,"Y"] - mean(P[,"Y"]) 
phi <- -pi/2 + sqrt(x^2 + y^2)/(R)
lambda <- atan2(y, x)

## Initial plot in 3D space
plot.retina(phi, lambda, R, Cu, Pt, ts.red, ts.green, 1:N.edge)

##
## Energy/error functions
## 

## Now for the dreaded elastic error function....
E.E <- function(p, Cu, C, L, B, R, verbose=FALSE) {
  phi    <- p[1:(length(p)/2)]
  lambda <- p[((length(p)/2)+1):length(p)]
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
  if (verbose>=1) {
    print(E)
  }
  return(E)
}

## ... and the even more dreaded gradient of the elastic error
dE.E <- function(p, Cu, C, L, B, R, verbose=FALSE) {
  phi    <- p[1:(length(p)/2)]
  lambda <- p[((length(p)/2)+1):length(p)]
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
  dE.phii     <- B %*% (fac * (sin(phii)*cos(phij)*cos(lambdai-lambdaj)
                               - cos(phii)*sin(phij)))
  dE.dlambdai <- B %*% (fac * cos(phii)*cos(phij)*sin(lambdai-lambdaj))
  return(c(dE.phii, dE.dlambdai))
}

## Optimisation with just the elastic energy
optimise.mapping <- function() {
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi, lambda)
  opt$conv <- 1
  while (opt$conv) {
    opt <- optim(opt$p, E.E, gr=dE.E,
                 method="BFGS", Cu=Cu, C=C, L=Ls, B=B, R=R, verbose=FALSE)
    ##               control=list(maxit=200))
    phi    <- opt$p[1:(length(opt$p)/2)]
    lambda <- opt$p[((length(opt$p)/2)+1):length(opt$p)]
    plot.retina(phi, lambda, R, Cu, Pt, ts.red, ts.green, 1:N.edge)
  }
  ## CG min: 288037
  ##
  ## BFGS maxit=100 is 261000
  ## with maxit=10 it is 277226
  ## with maxit=200 it is 273882
  ## with maxit=100 and phi0=50 it is 325785
}

## Energy function with areas
E.area <- function(p, Pt, A, R, verbose=FALSE) {
  phi    <- p[1:(length(p)/2)]
  lambda <- p[((length(p)/2)+1):length(p)]

  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Find areas of all triangles
  areas <- -0.5/R * dot(P[Pt[,1],], extprod3d(P[Pt[,2],], P[Pt[,3],]))
  E <- sum(0.5 * (areas - A)^2/A)
  return(E)
}

dE.area <- function(p, Pt, A, R, verbose=FALSE) {
  phi    <- p[1:(length(p)/2)]
  lambda <- p[((length(p)/2)+1):length(p)]

  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## expand triangulation
  Pt <- rbind(Pt, Pt[,c(2,3,1)], Pt[,c(3,1,2)])
  A  <- c(A, A, A)
  

  ## Slow way of computing gradient
  dAdPt1 <- -0.5/R * extprod3d(P[Pt[,2],], P[Pt[,3],])

  ## Find areas of all triangles
  areas <- dot(P[Pt[,1],], dAdPt1)

  dEdPt1 <- (areas - A)/A * dAdPt1

  Pt1topi <- matrix(0, length(phi), nrow(Pt))
  for(m in 1:nrow(Pt)) {
    Pt1topi[Pt[m,1],m] <- 1
  }
  dEdpi <- Pt1topi %*% dEdPt1
  dpidphi <- R * cbind(-sin(phi) * cos(lambda),
                   -sin(phi) * sin(lambda),
                   cos(phi))
  dpidlambda <- R * cbind(-cos(phi) * sin(lambda),
                           cos(phi) * cos(lambda),
                           0)

  dEdphi    <- apply(dEdpi * dpidphi,    1, sum)
  dEdlambda <- apply(dEdpi * dpidlambda, 1, sum)
  return(c(dEdphi, dEdlambda))

}

optimise.mapping.area <- function() {
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi, lambda)
  opt$conv <- 1
  while (opt$conv) {
    opt <- optim(opt$p, E.area, gr=dE.area,
                 method="BFGS", Pt=Pt, A=areas, R=R, verbose=FALSE)
    ##               control=list(maxit=200))
    phi    <- opt$p[1:(length(opt$p)/2)]
    lambda <- opt$p[((length(opt$p)/2)+1):length(opt$p)]
    print(E.area(opt$p, Pt, areas, R))
    plot.retina(phi, lambda, R, Cu, Pt, ts.red, ts.green, 1:N.edge)
  }
}


E.E.area <- function(p, Cu, C, L, B, Pt, A, R, verbose=FALSE) {
  return(E.E(p, Cu, C, L, B, R, verbose=verbose)
         + E.area(p, Pt, A, R, verbose=verbose))
}
    

dE.E.area <- function(p, Cu, C, L, B, Pt, A, R, verbose=FALSE) {
  return(dE.E(p, Cu, C, L, B, R, verbose=verbose)
         + dE.area(p, Pt, A, R, verbose=verbose))
}

optimise.mapping.E.area <- function() {
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi, lambda)
  opt$conv <- 1
  while (opt$conv) {
    opt <- optim(opt$p, E.E.area, gr=dE.E.area,
                 method="BFGS",
                 Cu=Cu, C=C, L=Ls, B=B,
                 Pt=Pt, A=areas, R=R, verbose=FALSE)
    ##               control=list(maxit=200))
    phi    <- opt$p[1:(length(opt$p)/2)]
    lambda <- opt$p[((length(opt$p)/2)+1):length(opt$p)]
    print(E.area(opt$p, Pt, areas, R))
    plot.retina(phi, lambda, R, Cu, Pt, ts.red, ts.green, 1:N.edge)
  }
}

## Energy function with elastic and dipole forces
E.D <- function(p, Cu, R, L, edge.inds, verbose=FALSE) {
  phi    <- p[1:(length(p)/2)]
  lambda <- p[((length(p)/2)+1):length(p)]
  
  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Consider points on the edge only
  P <- P[edge.inds,]

  ## Need to go round rim finding neigbours
  ## First find distances between all pairs of points on the rim
  ## (Maybe this isn't necessary?)
  E <- (outer(P[,1], P[,1], "-")^2 +
        outer(P[,2], P[,2], "-")^2 +
        outer(P[,3], P[,3], "-")^2) < (2*L)^2
  diag(E) <- FALSE

  print(Cu)
  ## Find the vector product of neigbouring pairs of points
  links <- classify.links(Cu, edge.inds)
  
  return(links$edge)
}

dE.D <- function(p, Cu, C, L, B, R, verbose=FALSE) {
  return(dE)
}

optimise.mapping.dipole <- function(dt) {
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi, lambda)
  opt$conv <- 1
  for (epoch in 1:200) {
    for (time in 1:10) {
      opt$p <- opt$p - dt * dE.E(opt$p, Cu, C, L, B, R)
    }
    phi    <- opt$p[1:(length(opt$p)/2)]
    lambda <- opt$p[((length(opt$p)/2)+1):length(opt$p)]
    print(E(opt$p, Cu, C, L, B, R))
    plot.retina(phi, lambda, R, Cu, Pt, ts.red, ts.green)
  }
}





