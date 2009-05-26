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
  ## Close the path by putting the initial point at the end of the list of
  ## segs
  if (close) {
    segs[[length(segs)]] <- rbind(segs[[length(segs)]],
                                  segs[[1]][1,])
  }
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
## Test code:        
## p1 = create.path(list(rbind(c(0, 0), c(1, 1))))
## p2 = create.path(list(rbind(c(0, 1), c(1, 0))))
## check.intersection.paths(p1, p2)
## Output should be
##        X   Y        s1        s2  f1  f2
## [1,] 0.5 0.5 0.7071068 0.7071068 0.5 0.5

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

