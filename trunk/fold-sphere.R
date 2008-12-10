source("../../data/Anatomy/cluster-analysis.R")

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


M <- 40                               # Number of grid points in x-dir
N <- 40                               # Number of grid points in y-dir
x0 <- -2000                             # minium x
y0 <- -3600                             # minium y
Delta <- 200                            # increment
G <- expand.grid(ix=1:M, iy=1:N)  # Grid; ix and iy are x & y indicies
G[,"x"] <- x0 + Delta*(G[,"ix"]-G[,"iy"]/2)
G[,"y"] <- y0 + Delta*G[,"iy"]*sqrt(3)/2

## A is the connectivity matrix: One row for each connection,
## which contains indicies of pairs of connected points
## Connect (x,y) to (x,y+1)
## Connect (x,y) to (x+1,y+1)
## Connect (x,y) to (x,y+1)
m <- M*N
## Start points & endpoints
sp <- (1:m)[(1:m %% M) != 0]
ep <- sp + 1
A <- cbind(sp, ep)
sp <- 1:(m-M)
ep <- sp + M
A <- rbind(A, cbind(sp, ep))
sp <- (1:(m-M))[(1:(m-M) %% M) != 0]
ep <- sp + M + 1
A <- rbind(A, cbind(sp, ep))

## Read in data
map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))

## Corner analysis
segs <- map.to.segments(map)
## Some curation is required here
segs4 <- segs[[4]][23:1,]
edge.path <- create.path(list(segs[[1]], segs[[3]], segs4), close=TRUE)

## Plotting
plot(edge.path[,"X"], edge.path[,"Y"])
lines(edge.path[,"X"], edge.path[,"Y"],lwd=2)
# plot.map(map, seginfo=FALSE)

segments(G[A[,1],"x"], G[A[,1],"y"],
         G[A[,2],"x"], G[A[,2],"y"])

G.keep <- c()
## Now remove points that are outside the retina
for (iy in unique(G[,"iy"])) {
  iG <- which(G[,"iy"] == iy)
  y <- unique(G[iG,"y"])
  horiz <- create.path(list(rbind(c(x0 - Delta*N ,y),
                                  c(x0 + Delta*M ,y))))
  lines(horiz[,"X"], horiz[,"Y"], col="red")
  ci <- check.intersection.paths(horiz, edge.path)

  if (!is.null(ci)) {
    if (nrow(ci) > 1) {
    for(i in 1:(nrow(ci)/2)) {
      cii <- ci[(i-1)*2+1:2,]
                points(cii[,"X"], cii[,"Y"], col="blue")
                G.keep.new <- which((G[, "x"] > min(cii[,"X"])) &
                                    (G[, "x"] < max(cii[,"X"])) &
                                    (G[,"iy"] == iy))
                points(G[G.keep.new,"x"], G[G.keep.new,"y"], col="blue")
      if((length(G.keep.new) < 10)) { print(cii) }
    }
    }
  }
}


