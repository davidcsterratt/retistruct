source("common.R")

## Create points that lie on path p2 at the fraction distances of
## the points along p1
create.corresponding.points <- function(p1, p2) {
  s1 <- p1[,"s"]
  F1 <- s1/s1[length(s1)]                  # fractional distance along path 1
  F1 <- F1[(F1!=0) & (F1!=1)]
  s2 <- p2[,"s"]
  F2 <- s2/s2[length(s2)]                  # fractional distance along path 2

  P <- find.points.in.path(F1, p2)
  P <- rbind(p2[,c("X", "Y")], P)
  s <- sort(c(F2, F1), index.return=TRUE)
  P <- P[s$ix,]
  p <- create.path(list(P))
  return(p)
}

## Stich together a pair of paths to create a new path and a set of correpondances
stitch <- function(tear) {
  par(pch=16)
  p1 <- create.corresponding.points(tear[[1]], tear[[2]])
  p2 <- create.corresponding.points(tear[[2]], tear[[1]])

  p <- create.path(list(rbind(p1[nrow(p1):1,c("X", "Y")], p2[-1,c("X", "Y")])))
  corrs <- cbind((nrow(p1)-1):1, nrow(p1) + (1:(nrow(p1)-1)))
  return(list(p=p, corrs=corrs))
}

for (tear in tears) {
  s <- stitch(tear)
  points(s$p[,c("X", "Y")])
  segments(s$p[s$corrs[,1],"X"], s$p[s$corrs[,1], "Y"],
           s$p[s$corrs[,2],"X"], s$p[s$corrs[,2], "Y"], col="red")
}


