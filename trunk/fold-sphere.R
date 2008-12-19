source("cluster-analysis.R")
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
plot.retina <- function(phi, lambda, R, C) {
  ## Now plot this in 3D space....
  x <- R*cos(phi)*cos(lambda) 
  y <- R*cos(phi)*sin(lambda)
  z <- R*sin(phi)
  rgl.clear()
  rgl.bg(color="white")
  segments3d(rbind(x[C[,1]],x[C[,2]]),
             rbind(y[C[,1]],y[C[,2]]),
             rbind(z[C[,1]],z[C[,2]]),xlab="x", color="black")
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
  E <- sum((l - L)^2)
  if (verbose>=1) {
    print(E)
  }
  return(E)
}

## ... and the even more dreaded gradient of the error
dE <- function(p, C, Cu, L, B, R, verbose=FALSE) {
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
  fac <- R * (l - L)/sqrt(1-x^2)
  dE.phii     <- B %*% (fac * (sin(phii)*cos(phij)*cos(lambdai-lambdaj)
                               - cos(phii)*sin(phij)))
  dE.dlambdai <- B %*% (fac * cos(phii)*cos(phij)*sin(lambdai-lambdaj))
  return(c(dE.phii, dE.dlambdai))
}

###
### Start of the code proper
### 

M <- 60                               # Number of grid points in x-dir
N <- 40                               # Number of grid points in y-dir
x0 <-  1000                             # minium x
y0 <- -1000                             # minium y
L <- 200                            # spacing between grid points

## G is contains the grid; ix and iy are x & y indicies and x and y
## are the actual x and y positions.  In order to make the grid
## triangular, for each successive y position, the x position of the
## row is shifted by L/2.
G <- expand.grid(ix=1:M, iy=1:N)
G[,"x"] <- x0 + L*(-M/2+G[,"ix"]-G[,"iy"]/2)
G[,"y"] <- y0 + L*(-N/2+G[,"iy"])*sqrt(3)/2

## Cu is the part of the connectivity matrix above the leading
## diagonal: One row for each connection, which contains indicies of
## pairs of connected points
## Connect (x,y) to (x,y+1)
## Connect (x,y) to (x+1,y+1)
## Connect (x,y) to (x,y+1)
m <- M*N                              # total number of points in grid
## Start points & endpoints
sp <- (1:m)[(1:m %% M) != 0]
ep <- sp + 1
Cu <- cbind(sp, ep)
sp <- 1:(m-M)
ep <- sp + M
Cu <- rbind(Cu, cbind(sp, ep))
sp <- (1:(m-M))[(1:(m-M) %% M) != 0]
ep <- sp + M + 1
Cu <- rbind(Cu, cbind(sp, ep))

## Read in data
map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))

## Corner analysis
segs <- map.to.segments(map)
## Some curation is required here
segs4 <- segs[[4]][23:1,]
edge.path <- create.path(list(segs[[1]], segs[[3]], segs4), close=TRUE)

## Plot the outline of the flattened retina
plot(edge.path[,"X"], edge.path[,"Y"], xlab="x", ylab="y")
lines(edge.path[,"X"], edge.path[,"Y"],lwd=2)

## Plot the entire grid
## segments(G[Cu[,1],"x"], G[Cu[,1],"y"],
##         G[Cu[,2],"x"], G[Cu[,2],"y"])

## Now remove points that are outside the retina
## To do this, the intersections of each of the horizontal
## grid lines with the retina are determined
## Points outwith the range of thses connections are discareded
iG.keep <- c()                         # Indicies of rows of G to keep
for (iy in unique(G[,"iy"])) {
  iG <- which(G[,"iy"] == iy)
  y <- unique(G[iG,"y"])
  horiz <- create.path(list(rbind(c(x0 - L*N ,y),
                                  c(x0 + L*M ,y))))
  ## This returns a list of intersection points
  ## There should be an even number of them
  ci <- check.intersection.paths(horiz, edge.path)

  if (!is.null(ci)) {
    for(i in 1:(nrow(ci)/2)) {
      cii <- ci[(i-1)*2+1:2,]
                ## points(cii[,"X"], cii[,"Y"], col="blue")
                iG.keep.new <- which((G[, "x"] > min(cii[,"X"])) &
                                    (G[, "x"] < max(cii[,"X"])) &
                                    (G[,"iy"] == iy))
      ## points(G[iG.keep.new,"x"], G[iG.keep.new,"y"], col="blue")
      iG.keep <- c(iG.keep, iG.keep.new)
    }
  }
}

## Keep only connections within the retina ("any" in place of "all" would
## allow connections to go over the boundary
Cu <- Cu[apply(matrix(Cu %in% iG.keep, nrow(Cu), ncol(Cu)),1,all),]
## Now reorganise G and C so that only the remaining nodes are left in G
## and that the elements of C point to the pruned version of G
G <- G[iG.keep,]
Cu <- cbind(match(Cu[,1], iG.keep), match(Cu[,2], iG.keep))

## Plot the connections that remain after the pruning
segments(G[Cu[,1],"x"], G[Cu[,1],"y"],
         G[Cu[,2],"x"], G[Cu[,2],"y"],col="blue")

## C is the symmetric connectivity matrix
C <- rbind(Cu, Cu[,2:1])
## Matrix to map line segments onto the points they link
B <- matrix(0, nrow(G), nrow(C))
for (i in 1:nrow(C)) {
  B[C[i,1],i] <- 1
}

## Estimate the area. It's roughly equal to the number of remaining points
## times the area of the rhomboid.
area <- length(iG.keep) * L^2 * sqrt(3)/2

## From this we can infer what the radius should be from the formula
## for the area of a sphere which is cut off at a lattitude of phi0
## area = 2 * PI * R^2 * (sin(phi0)+1)
phi0 <- 50*pi/180
R <- sqrt(area/(2*pi*(sin(phi0)+1)))

## Now assign each point to a location in the phi, lambda
## Shift coordinates to rough centre of grid
x <- G[,"x"] - x0 
y <- G[,"y"] - y0 
phi <- -pi/2 + sqrt(x^2 + y^2)/(R)
lambda <- atan2(y, x)

## Initial plot in 3D space
plot.retina(phi, lambda, R, Cu)

## Optimisation and plotting 
opt <- list()
opt$p <- c(phi, lambda)
opt$conv <- 1
while (opt$conv) {
  opt <- optim(opt$p, E, gr=dE,
               method="BFGS", Cu=Cu, L=L, B=B, C=C, R=R, verbose=1)
##               control=list(maxit=200))
  phi    <- opt$p[1:(length(opt$p)/2)]
  lambda <- opt$p[((length(opt$p)/2)+1):length(opt$p)]
  plot.retina(phi, lambda, R, C)
}
## CG min: 288037
##
## BFGS maxit=100 is 261000
## with maxit=10 it is 277226
## with maxit=200 it is 273882
## with maxit=100 and phi0=50 it is 325785






