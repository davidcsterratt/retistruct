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
create.path <- function(segs) {
  p <- matrix(0, 0, 4)                  # the path matrix
  colnames(p) <- c("X", "Y", "s", "l")
  s.cum <- 0
  for (i in 1:length(segs)) {
    seg <- segs[[i]]
    v <- diff(seg)
    l <- sqrt(apply(v^2, 1, sum))
    s <- s.cum + c(0, cumsum(l))
    s.cum <- s[length(s)]
    if (i > 1) {
      p <- rbind(p, rep(NA, 4))
    }
    p <- rbind(p,
               cbind(seg,
                     s,
                     c(l, NA)))
  }
  return(p)
}

## Find a point a fractional distance f along a path p
find.points.in.path <- function(f, p) {
  s <- p[,"s"]
  F <- s/s[length(s)]                  # fractional distance along path
  # remove NAs. Have to replce with distances in sequence
  F[is.na(F)] <- F[which(is.na(F))-1]          
  ## Find intervals in which points occur
  is <- findInterval(f, F, rightmost.closed=TRUE)
  ## Find fractional distance *within* interval
  f <- (f - F[is])/(F[is+1] - F[is])
  ## Interpolate to find the points
  P <- (1-f)*p[is,c("X","Y")] + f*p[is+1,c("X","Y")]
  return(P)
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
  if(length(out) > 0) {
    return(out)
  } else {
    return(NULL)
  }
}

## Test code:        
## p1 = create.path(list(rbind(c(0, 0), c(1, 1))))
## p2 = create.path(list(rbind(c(0, 1), c(1, 0))))
## check.intersection.paths(p1, p2)
## Output should be
##        X   Y        s1        s2  f1  f2
## [1,] 0.5 0.5 0.7071068 0.7071068 0.5 0.5

###
### Start of code
### 

## Read in data
map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))

## Corner analysis
segs <- map.to.segments(map)
## Some curation is required here
segs4 <- segs[[4]][23:1,]
edge <- join.segments(list(segs[[1]], segs[[3]], segs4))
##map <- segs[[1]]

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

## Hand pick corners
## c(16, 31)
## c(53, 65)
## c(83, 102)
## Define rim as list of line segments
rim <- list(edge[16:31,],
            edge[53:65,],
            edge[83:94,])

## Find distance along edges
rim.path <- create.path(rim)

## Distribute points equally along edges
P <- find.points.in.path(seq(0, by=1/M, len=M), rim.path)

## Experimental code: add a hand-curated tear pair
tear <- list()
tear[[1]]  = create.path(list(edge[42:31,]))
tear[[2]]  = create.path(list(edge[42:53,]))


## Fitting of mesh to data
m <- (N-1)*M+1
A <- 1/L[1:m,1:m]
B <- 1/L[1:m,(m+1):(m+M)]
C <- matrix(0, m, 0)
A[is.na(A)] <- 0
B[is.na(B)] <- 0
V1 <- matrix(0, 0, 2)
colnames(V1) <- c("X", "Y")
V2 <- V1
R1 <- V1
R2 <- V1
S <- matrix(0, 0, 2)
f <- c()                                # sliding points
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
  dev.print(pdf, file="demo1.pdf")
  
  lines(rim.path[,"X"], rim.path[,"Y"], col="red", lwd=2)
  dev.print(pdf, file="demo2.pdf")
  
  points(P[,1], P[,2], pch=16)
  dev.print(pdf, file="demo3.pdf")
  ## text(P[,1], P[,2], labels=1:M)
  
  R <- rbind(Q, P)
  plot.mesh2(R, L)
  dev.print(pdf, file="demo4.pdf")

  ## Next step: try to find points at which the mesh intersects the tears
  ## Work through the connectivity matrix L
  for (i in 1:m) {
    js <- which(!is.na(L[i,]))
    js <- js[js>i]
    for(j in js) {
      path <- create.path(list(rbind(R[i,], R[j,])))
      lines(path[,"X"], path[,"Y"], col="green")
      ci <- check.intersection.paths(path, tear[[1]])
      if (!is.null(ci)) {
        ##        print(ci)
        points(ci[,"X"], ci[,"Y"])
        ## Remove connection from i to j
        if (j <= m) {
          A[i, j] <- 0
          A[j, i] <- 0
          ## Create a new sliding point
          f.new <- ci[1,"f2"]
          f <- c(f, f.new)
          n.f <- length(f)
          C.new <- matrix(0, m, 2)
          C.new[i, 1] <- 1/L[i,j]
          C.new[j, 2] <- 1/L[i,j]          
          C <- cbind(C, C.new)

          ## Find the V and R of this tear
          V1 <- rbind(V1, ci[1,c("V2X", "V2Y")])
          R1 <- rbind(R1, ci[1,c("R2X", "R2Y")])

          ## Find the V and R of the corresponding tear
          ci <- check.intersection.paths(path, tear[[1]])
          V2 <- rbind(V2, ci[1,c("V2X", "V2Y")])
          R2 <- rbind(R2, ci[1,c("R2X", "R2Y")])

          ## Plot the intersection points
          S <- rbind(S,
                     find.points.in.path(f.new, tear[[1]]),
                     find.points.in.path(f.new, tear[[2]]))
          points(inf[,1], inf[,2], col="blue")
          
        } else {
          B[i, j-m] <- 0
        }
      }
    }
  }

  plot.mesh3(rbind(Q, P, S), cbind(A, B, C))
  ## Next step: what happens if the longest links are removed?
  ## Need to find the length of all the links
  ## l <- link.lengths(R, L[1:((N-1)*M)+1,])
  ## m.i <- which.max.matrix(l)
  ## L[m.i[1], m.i[2]] <- NA
  ## L[m.i[2], m.i[1]] <- NA
  ## C[which.max(l)] <- 0
}

## Now try improving on this initial guess by using the energy function
## (l - L)^2
E <- function(p, P, A, B, C, L, tear) {
  ## Construct Q and f from p
  Q <- matrix(p[1:(2*m)], m, 2)
  f <- p[(2*m+1):length(p)]

  ## Get some useful constants
  m <- ncol(A)
  M <- ncol(B) 
  n.s <- ncol(C)
  
  ## Find the intersection points from the sliding points
  S <- matrix(0, 0, 2)
  for (f.new in f) {
    S <- rbind(S,
               find.points.in.path(f.new, tear[[1]]),
               find.points.in.path(f.new, tear[[2]]))
  }

  l    <-   link.lengths(Q, L[1:m,1:m])
  lhat <-   link.lengths2(rbind(Q, P), B)
  ltilde <- link.lengths2(rbind(Q, S), C)
  E <-
    sum(sum((A*L[1:m,1:m]     - l)^2   , na.rm=TRUE)) +
    sum(sum((B*L[1:m,(m+1):(m+M)] - lhat)^2, na.rm=TRUE))
  print(E)
  return(E)
}

dE <- function(p, P, A, B, C, L, tear) {
  ## Construct Q and f from p
  Q <- matrix(p[1:(2*m)], m, 2)
  f <- p[(2*m+1):length(p)]

  print(f)
  
  ## Get some useful constants
  m <- ncol(A)
  M <- ncol(B) 
  n.s <- ncol(C)

  ## Find the intersection points from the sliding points
  S <- matrix(0, 0, 2)
  for (f.new in f) {
    S <- rbind(S,
               find.points.in.path(f.new, tear[[1]]),
               find.points.in.path(f.new, tear[[2]]))
  }

  l    <-   link.lengths(Q, L[1:m,1:m])
  lhat <-   link.lengths2(rbind(Q, P), B)
  ltilde <- link.lengths2(rbind(Q, S), C)

  Ahat <- A * (1 - l/L[1:m,1:m])
  Ahat[is.na(Ahat)] <- 0

  Bhat <- B * (1 - lhat/L[1:m,(m+1):(m+M)])
  Bhat[is.na(Bhat)] <- 0
  
  ## Diagonal matrix of row sums of [ A B ]
  D <- diag(apply(cbind(Ahat, Bhat), 1, sum))
  
  dEdQ <- D %*% Q - cbind(Ahat, Bhat) %*% rbind(Q, P)
  dEdQ <- c(as.vector(dEdQ),  rep(0, length(f)))
  
  return(dEdQ)
}

dev.print(pdf, file="demo5.pdf")
opt <- optim(c(as.vector(Q), f), E, gr=dE, method="BFGS", L=L, P=P,
             A=A, B=B, C=C, tear=tear, control=list(maxit=100))
##plot.map(map)
##plot.mesh2(rbind(Q, P), L)
