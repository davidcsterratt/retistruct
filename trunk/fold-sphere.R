source("cluster-analysis.R")
source("triangulate.R")
library("geometry")
source("tsearch.R")
library("gtools")
library("rgl")

## We need a structure in which to store paths, such as the rim, which
## may have gaps in them. The structure should have one row for each
## point. Each row should contain:
## * the X & Y coords of each point
## * the cumulative distance along the edge
## * the length of the edge which starts at that point (the
##   last point will have no information
## Each continuous section of a path is called a segment. At the
## boundary between two segments, there is a row containing NA 

## create.path(segments)
## create a path from a list of line segments
create.path <- function(segs, close=FALSE) {
  p <- matrix(0, 0, 4)                  # the path matrix
  colnames(p) <- c("X", "Y", "s", "l")
  s.cum <- 0
  for (i in 1:length(segs)) {
    seg <- segs[[i]]
    v <- diff(seg)
    l <- sqrt(apply(v^2, 1, sum))
    s <- s.cum + c(0, cumsum(l))
    s.cum <- s[length(s)]
    if ((i > 1) && !close) {
      p <- rbind(p, rep(NA, 4))
    }
    p <- rbind(p,
               cbind(seg,
                     s,
                     c(l, NA)))
  }
  if (close) {
    p <- rbind(p, p[1,]) 
  }
  return(p)
}

## Check whether line from P1 to Q1 and from P2 to Q2 intersect
check.intersection <- function(P1, Q1, P2, Q2) {
  if ((max(P1[1], Q1[1]) < min(P2[1], Q2[1])) ||
      (min(P1[1], Q1[1]) > max(P2[1], Q2[1])) ||
      (max(P1[2], Q1[2]) < min(P2[2], Q2[2])) ||
      (min(P1[2], Q1[2]) > max(P2[2], Q2[2]))) {
    return(FALSE)
  }
  M <- cbind(Q1-P1, -Q2+P2)
  if (!(det(M) == 0)) {
    lambda <- solve(M) %*% (P2-P1)
    if (all((lambda<1) & (lambda>0))) {
      return(list(lambda=lambda, R=(1-lambda[1])*P1 + lambda[1]*Q1))
    }
  }
  return(FALSE)
}

## Test to see if there is an intersection between two paths p1 and p2
## Return a matrix whose columns comprise:
## * P  : the point of intersection
## * s1 : the distance along path 1 of the intersection
## * s2 : the distance along path 2 of the intersection
## * f1 : the fraction along path 1 of the intersection
## * f2 : the fraction along path 2 of the intersection
check.intersection.paths <- function(p1, p2) {
  out <- matrix(0, 0, 14)
  colnames(out) <- c("X", "Y", "s1", "s2", "f1", "f2", "V1X", "V1Y", "R1X", "R1Y", "V2X", "V2Y", "R2X", "R2Y")
  for(i in 1:(nrow(p1)-1)) {
    for(j in 1:(nrow(p2)-1)) {
      if (all(!is.na(p1[c(i,i+1),c("X", "Y")]) &
              !is.na(p2[c(j,j+1),c("X", "Y")]))) {
        ci <- check.intersection(p1[i,c("X", "Y")],p1[i+1,c("X", "Y")],
                                 p2[j,c("X", "Y")],p2[j+1,c("X", "Y")])
        if (is.list(ci)) {
          s1 <- p1[i,"s"]+p1[i,"l"]*ci$lambda[1]
          s2 <- p2[j,"s"]+p2[j,"l"]*ci$lambda[2]
          stot1 <- p1[nrow(p1), "s"]
          stot2 <- p2[nrow(p2), "s"]
          f1 <- s1/stot1
          f2 <- s2/stot2
          v1 <- (p1[i+1,c("X", "Y")] - p1[i,c("X", "Y")]) * stot1/p1[i,"l"]
          v2 <- (p2[i+1,c("X", "Y")] - p2[i,c("X", "Y")]) * stot2/p2[i,"l"]
          r1 <- p1[i,c("X", "Y")] - p1[i,"s"]/stot1 * v1
          r2 <- p2[i,c("X", "Y")] - p2[i,"s"]/stot2 * v2
          out <- rbind(out,
                       c(ci$R,
                         s1, s2,
                         f1, f2,
                         v1, r1,
                         v2, r2))
        }
      }
    }
  }
  if(length(out) > 0) {
    return(out)
  } else {
    return(NULL)
  }
}

## Function to plot the retina in spherical coordinates
## phi - lattitude of points
## lambda - longitude of points
## R - radius of sphere
## C - connection table
## Ct - triagulated connection table
## sys - data frame containing positions of points
plot.retina <- function(phi, lambda, R, C,
                        Ct=matrix(NA, 0, 3), ts.red, ts.green) {
  ## Now plot this in 3D space....
  x <- R*cos(phi)*cos(lambda) 
  y <- R*cos(phi)*sin(lambda)
  z <- R*sin(phi)
  P <- cbind(x, y, z)
  rgl.clear()
  rgl.bg(color="white")

  ## Pt.surf = t(surf.tri(P, Pt))
  ## Not run: 
  ## rgl.triangles(P[Pt.surf,1], P[Pt.surf,2], P[Pt.surf,3],
  ## col="blue", alpha=.2)
  segments3d(rbind(x[C[,1]],x[C[,2]]),
             rbind(y[C[,1]],y[C[,2]]),
             rbind(z[C[,1]],z[C[,2]]),xlab="x", color="black")

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
  
}

## Formula for central angle
central.angle <- function(phi1, lambda1, phi2, lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)))
}

## Now for the dreaded error function....
E <- function(p, Cu, C, L, B, R, verbose=FALSE) {
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

## ... and the even more dreaded gradient of the error
dE <- function(p, Cu, C, L, B, R, verbose=FALSE) {
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

## Read in data
sys <- read.sys("../../data/Anatomy/ALU/M643-4/CONTRA")
map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))

## Corner analysis
segs <- map.to.segments(map)
## Some curation is required here
segs4 <- segs[[4]][23:1,]
edge.path <- create.path(list(segs[[1]], segs[[3]], segs4), close=TRUE)

## All the points on the edge path apart from the first, which is repeated by the last
P <- cbind(X=edge.path[-1,"X"], Y=edge.path[-1,"Y"])

## Find the bottom left corner of the box around the edge
P0 <- c(min(P[,"X"])-1, min(P[,"Y"])-1)
## Find the top right corner of the box around the edge
P1 <- c(max(P[,"X"])+1, max(P[,"Y"])+1)
## Find the top left corner of the box around the edge
P.TL <- c(min(P[,"X"])-1, max(P[,"Y"])+1)
## Find the bottom right corner of the box around the edge
P.BR <- c(max(P[,"X"])+1, min(P[,"Y"])-1)

## Convert these points to the trignometric coordinates, and
## find the maximum coordinates along the u v axes
Q0 <- cart.to.trig(P0, P0, L)
Q1 <- cart.to.trig(P1, P0, L)
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
P <- rbind(P.edge, P.grid)
## P <- P.grid

## Find Delaunay triangulation of the points
Pt <- delaunayn(P)
areas <- tri.area(P, Pt)
Pt <- Pt[areas>0,]

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
  }
}
Pt <- na.omit(Pt)

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
Cu <- rbind(Pt[,1:2], Pt[,2:3], Pt[,c(1,3)])
##Cu <- t(apply(Cu, 1, sort))
##Cu <- Cu[!duplicated(Cu),]
Cu <- Unique(Cu)

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

## Now assign each point to a location in the phi, lambda
## Shift coordinates to rough centre of grid
x <- P[,"X"] - mean(P[,"X"]) 
y <- P[,"Y"] - mean(P[,"Y"]) 
phi <- -pi/2 + sqrt(x^2 + y^2)/(R)
lambda <- atan2(y, x)

## Initial plot in 3D space
plot.retina(phi, lambda, R, Cu, Pt, ts.red, ts.green)

optimise.mapping <- function() {
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi, lambda)
  opt$conv <- 1
  while (opt$conv) {
    opt <- optim(opt$p, E, gr=dE,
                 method="BFGS", Cu=Cu, C=C, L=Ls, B=B, R=R, verbose=1)
    ##               control=list(maxit=200))
    phi    <- opt$p[1:(length(opt$p)/2)]
    lambda <- opt$p[((length(opt$p)/2)+1):length(opt$p)]
    plot.retina(phi, lambda, R, Cu, Pt, ts.red, ts.green)
  }
  ## CG min: 288037
  ##
  ## BFGS maxit=100 is 261000
  ## with maxit=10 it is 277226
  ## with maxit=200 it is 273882
  ## with maxit=100 and phi0=50 it is 325785
}






