## Return next index in path
path.next <- function(i, g, h) {
  return(ifelse(h[i]==i, g[i], h[i]))
}

## Return sequence of indices in path between i and j, governed by
## pointer vector p
path <- function(i, j, g, h) {
  ## if (length(g) != length(h)) {
  ##   stop("g and h must be the same length")
  ## }
  if (i == j) {
    return(i)
  } else {
    return(c(i, path(path.next(i, g, h), j, g, h)))
  }
}

## Return sequence of indices in path between i and j, governed by
## pointer vector p
path.length <- function(i, j, g, h, P) {
  if (any(is.na(c(i, j)))) {
    stop("i or j contains NA")
  }
  if (i == j) {
    return(0)
  } else {
    if (h[i] == i) {
      return(sqrt(sum((P[i,] - P[g[i],])^2)) + path.length(g[i], j, g, h, P))
    } else {
      return(path.length(h[i], j, g, h, P))
    }
  }
}

new.point <- function(P, i0, i1, f) {
  return(rbind(P,
  (1 - f)*P[i0,] +
  f*     P[i1,]))
}
## ##' @importFrom graphics arrows
## plot.path <- function(P, gf, gb, hf, hb, h, add=FALSE, ...) {
##   if (!add) {
##     plot(P)
##   }  
##   s <- which(!is.na(gf))
##   d <- gf[s]
##   arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="red", lwd=2, ...)

##   s <- which(!is.na(gb))
##   d <- gb[s]
##   arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="blue", ...)

##   s <- which(hf != 1:length(hf))
##   d <- hf[s]
##   arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="red", lty=2, lwd=2, ...)

##   s <- which(hb != 1:length(hb))
##   d <- hb[s]
##   arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="blue", lty=2, lwd=1, ...)

  
##   s <- which(h!=1:nrow(P))
##   d <- h[s]
##   arrows(P[s,1], P[s,2], P[d,1], P[d,2], ...)
  
##   points(P[,1], P[,2], cex=2.2, pch=21, bg="white") 
##   text(P[,1], P[,2], 1:nrow(P), cex=0.7)
## }

##' @importFrom graphics arrows
##' @export
plot.path <- function(x, xaxt="n",
                      xlim=NULL,
                      ylim=NULL,
                      arrows=FALSE,
                      ...) {

  plot(P)
  P <- x$P
  gf <- x$gf
  gb <- x$gb
  hf <- x$hf
  hb <- x$hb
  h  <- x$h
  
  s <- which(!is.na(gf))
  d <- gf[s]
  arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="red", lwd=2, ...)

  s <- which(!is.na(gb))
  d <- gb[s]
  arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="blue", ...)

  s <- which(hf != 1:length(hf))
  d <- hf[s]
  arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="red", lty=2, lwd=2, ...)

  s <- which(hb != 1:length(hb))
  d <- hb[s]
  arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="blue", lty=2, lwd=1, ...)
  
  s <- which(h!=1:nrow(P))
  d <- h[s]
  arrows(P[s,1], P[s,2], P[d,1], P[d,2], ...)
  
  points(P[,1], P[,2], cex=2.2, pch=21, bg="white") 
  text(P[,1], P[,2], 1:nrow(P), cex=0.7)
}

stitch.paths <- function(P,
                         VF0, VF1, VB0, VB1,
                         gf, gb,
                         hf, hb,
                         h,
                         epsilon,
                         report = function(x) {}) {
  ## Compute the total path length along each side of the tear
  Sf <- path.length(VF0, VF1, gf, hf, P)
  Sb <- path.length(VB0, VB1, gb, hb, P)

  ## Initialise forward and backward indicies to first points in each
  ## path
  i0 <- VF0
  j0 <- VB0
  i <- path.next(VF0, gf, hf)
  j <- path.next(VB0, gb, hb)
  sf0 <- 0
  sb0 <- 0
  
  ## Start stitching
  while(!((i == VF1) && (j == VB1))) {
    ## Get distance to current points in either path
    sf <- path.length(VF0, i, gf, hf, P)
    sb <- path.length(VB0, j, gb, hb, P)
    report(paste("i =", i, "i0 =", i0, "j =", j, "j0 =", j0,
                 "sf =", sf, "sf0 =", sf0, "sb =", sb, "sb0 =", sb0))
    if (sf/Sf <= sb/Sb) {
      ## If forward point is behind backward point, project forward to
      ## backward path

      if (any(h[-i] == i)) {
        ## Point is already pointed to
        report(paste("Point", i, "at", sf, "along forward path already pointed to"))
      } else {
        
        if (hf[i0] == i) {
          ## Point is pointed to by its predecessor, so make it point to
          ## the same point
          report(paste("Point", i, "in forward path pointed to by", i0))
          h[i] <- h[i0]
        } else {
          if (abs(sf/Sf*Sb - sb) < epsilon) {
            ## If projection of forward point is within tolerance of
            ## backward point, don't create a new point
            report(paste("Point", j, "at", sb, "along backward path within",
                         epsilon, "of projection from", i, "in forward path")) 
            h[i] <- j
          } else {
            ## Insert a point in the backward path. Point j is ahead,
            f <- (sf/Sf*Sb - sb0)/(sb - sb0)
            P <- new.point(P, j0, j, f)
            k <- nrow(P)
            report(paste("Insert point", k, "at", sf/Sf*Sb,
                         "along backward path; projection from", i,
                         "in forward path"))
            h[k] <- i
            gb[k] <- j
            gf[j] <- k
            gf[k] <- j0
            gb[j0] <- k
            hf[k] <- k
            hb[k] <- k
            j0 <- k

            ## gf[i0] <- k
            ## gf[k] <- i
            ## gb[i] <- k
            ## gb[k] <- i0
          }
        }
      }
      ## Move on to next point
      i0 <- i
      sf0 <- sf
      i <- path.next(i0, gf, hf)
    } else {
      ## If forward point is behind backward point, project forward to
      ## backward path

      if (any(h[-j] == j)) {
        ## Point is already pointed to
        report(paste("Point", j, "at", sb, "along backward path already pointed to"))
      } else {
        if (hb[j0] == j) {
          ## Point is pointed to by its predecessor, so make it point to
          ## the same point
          h[j] <- h[j0]
        } else {
          if (abs(sb/Sb*Sf - sf) < epsilon) {
            ## If projection of forward point is within tolerance of
            ## backward point, don't create a new point
            report(paste("Point", i, "at", sf, "along forward path within",
                         epsilon, "of projection from", j, "in backward path")) 
            h[j] <- i
          } else {
            ## Insert a point in the backward path. Point j is ahead,
            report(paste("Insert point at", sb/Sb*Sf,
                         "along forward path; projection from", j,
                         "in backward path"))
            f <- (sb/Sb*Sf - sf0)/(sf - sf0)
            P <- new.point(P, i0, i, f)
            k <- nrow(P)
            h[k] <- j
            gf[k] <- i
            gb[i] <- k
            gb[k] <- i0
            gf[i0] <- k
            hf[k] <- k
            hb[k] <- k
            i0 <- k
            ## gf[j] <- k
            ## gf[k] <- j0
          }
        }
      }
      ## Move on
      j0 <- j
      sb0 <- sb
      j <- path.next(j0, gb, hb)
    }
  }
  out <- list(P=P, gf=gf, gb=gb, hf=hf, hb=hb, h=h)
  class(out) <- "path"
  return(out)
}
