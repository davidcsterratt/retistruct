##' The member function \code{stitchSubpaths()} stitches together two
##' subpaths of the outline. One subpath is stitched in the forward
##' direction from the point indexed by \code{VF0} to the point
##' indexed by \code{VF1}. The other is stitched in the backward
##' direction from \code{VB0} to \code{VB1}. Each point in the subpath
##' is linked to points in the opposing pathway at an equal or
##' near-equal fraction along. If a point exists in the opposing
##' pathway within a distance \code{epsilon} of the projection, this
##' point is connected. If no point exists within this tolerance, a
##' new point is created.
##'
##' @title Add point fullcuts to the outline
##' @return To the \code{\link{Outline}} object this adds
##' \item{\code{hf}}{point cut mapping in forward direction for
##' points on boundary}
##' \item{\code{hb}}{point cut mapping in backward direction for
##' points on boundary}
##' @export
PathOutline <- R6Class("PathOutline",
  inherit = Outline,
  public = list(
    ##' @field hf Forward fullcuts
    hf = NULL,
    ##' @field hb Backward fullcuts
    hb = NULL,
    ##' @description Add points to the outline register of points
    ##' @param P 2 column matrix of points to add
    ##' @param fid fragment id of the points
    ##' @return The ID of each added point in the register. If points already
    ##'   exist a point will not be created in the register,
    ##'   but an ID will be returned
    addPoints = function(P, fid) {
      pids <- super$addPoints(P, fid)
      ## For *new* points set forward and backward pointers
      newpids <- pids
      if (length(self$hf) > 0) {
        newpids <- setdiff(pids, 1:length(self$hf))
      }
      if (length(newpids) > 0) {
        self$hf[newpids] <- newpids
        self$hb[newpids] <- newpids
      }
      return(pids)
    },
    ##' @description Get next point in path for
    ##' @param pids Point IDs of points to get next position
    nextPoint = function(pids) {
      return(sapply(pids, function(i) {path.next(i, self$gf, self$hf)}))
    },
    ##' @description Insert point at a fractional distance between points
    ##' @param i0 Point ID of first point
    ##' @param i1 Point ID of second point
    ##' @param f Fraction of distance  between points \code{i0} and
    ##'   \code{i1} at which to insert point
    insertPoint = function(i0, i1, f) {
      if (!((self$gf[i0] == i1) || (self$gf[i1] == i0))) {
        stop("Points", i0, "and", i1, "are not connected by an edge")
      }
      if ((f >= 1) | (f <= 0)) {
        stop("f argument should be between 0 and 1. f = ", f)
      }
      fid  <- self$getFragmentIDsFromPointIDs(i0)
      fid1  <- self$getFragmentIDsFromPointIDs(i1)
      if (fid != fid1) {
        stop("Fragment IDs of points differ")
      }
      p <-
        (1 - f)*self$getPoints(i0) +
        f*      self$getPoints(i1)
      ## Find the index of any row of P that matches p
      n <- anyDuplicated(rbind(self$getPointsXY(c(i0, i1)), p[c("X", "Y")]),
                         fromLast=TRUE)
      if (n == 0) {
        ## If the point p doesn't exist
        n <- self$addPoints(p, fid) # n is Index of new point
        ## Update forward and backward pointers
        if (!is.na(self$gf[i0] == i1) && self$gf[i0] == i1) {
          self$gb[n]  <- i0
          self$gf[n]  <- i1
          self$gb[i1] <- n
          self$gf[i0] <- n
        } else {
          self$gb[n]  <- i1
          self$gf[n]  <- i0
          self$gb[i0] <- n
          self$gf[i1] <- n
        }
        ## Update fullcuts
        self$hf[n] <- n
        self$hb[n] <- n
      }
      return(n)
    },
    ##' @description Stitch subpaths
    ##' @param VF0 First vertex of \dQuote{forward} subpath
    ##' @param VF1 Second vertex of \dQuote{forward} subpath
    ##' @param VB0 First vertex of \dQuote{backward} subpath
    ##' @param VB1 Second vertex of \dQuote{backward} subpath
    ##' @param epsilon Minimum distance between points
    stitchSubpaths = function(VF0, VF1, VB0, VB1,
                              epsilon) {
      ## Compute the total path length along each side of the subpath
      Sf <- path.length(VF0, VF1, self$gf, self$hf, self$getPointsScaled())
      Sb <- path.length(VB0, VB1, self$gb, self$hb, self$getPointsScaled())

      report(paste0("\nstitchSubpaths(", VF0, ", ", VF1, ", ", VB0, ", ", VB1, ", ", epsilon, ")\n"))

      report(paste("Sf =", Sf, ", Sb =", Sb, "\n"))

      ## Initialise forward and backward indicies to first points in each
      ## path
      i0 <- VF0
      j0 <- VB0
      i <- path.next(VF0, self$gf, self$hf)
      j <- path.next(VB0, self$gb, self$hb)

      sf0 <- 0
      sb0 <- 0

      ## Start stitching
      while(!((i == VF1) && (j == VB1))) {
        ## Get distance to current points in either path
        sf <- path.length(VF0, i, self$gf, self$hf, self$getPointsScaled())
        sb <- path.length(VB0, j, self$gb, self$hb, self$getPointsScaled())
        report(paste("i =", i, "i0 =", i0, "j =", j, "j0 =", j0,
                     "sf =", sf, "sf0 =", sf0, "sb =", sb, "sb0 =", sb0))
        ## Defensive programming - a test might be better, but at
        ## least this will highlight one class of error
        if (sf - Sf > 1E-10*Sf) {
          stop("Distance along forward path, ", sf, ", is greater than length of forward path, ", Sf, " ; Sf - sf=", Sf - sf)
        }
        if (sb - Sb > 1E-10*Sb) {
          stop("Distance along forward path, ", sb, ", is greater than length of forward path, ", Sb, " ; Sb - sb=", Sb - sb)
        }
        if (sf/Sf <= sb/Sb) {
          ## If forward point is behind backward point, project forward to
          ## backward path

          if (any(self$h[-i] == i)) {
            ## Point is already pointed to
            report(paste("Point", i, "at", sf, "along forward path already pointed to"))
          } else {
            if (self$hf[i0] == i) {
              ## Point is pointed to by its predecessor, so make it point to
              ## the same point
              report(paste("Point", i, "in forward path pointed to by", i0))
              self$h[i] <- self$h[i0]
            } else {
              if (abs(sf/Sf*Sb - sb) < epsilon) {
                ## If projection of forward point is within tolerance of
                ## backward point, don't create a new point
                report(paste("Point", j, "at", sb, "along backward path within",
                             epsilon, "of projection from", i, "in forward path"))
                report(paste("Not creating new point, but setting point", i, "to link to", j))
                self$h[i] <- j
              } else {
                ## Insert a point in the backward path. Point j is ahead.
                f <- (sf/Sf*Sb - sb0)/(sb - sb0)
                k <- self$insertPoint(j0, j, f)
                report(paste("Insert point", k, "at", sf/Sf*Sb,
                                  "along backward path; projection from", i,
                                  "in forward path"))
                self$h[k] <- i
                j0 <- k
                sb0 <- sf/Sf*Sb
              }
            }
          }
          ## Move on to next point
          i0 <- i
          sf0 <- sf
          i <- path.next(i0, self$gf, self$hf)
        } else {
          ## If forward point is behind backward point, project forward to
          ## backward path

          if (any(self$h[-j] == j)) {
            ## Point is already pointed to
            report(paste("Point", j, "at", sb, "along backward path already pointed to"))
          } else {
            if (self$hb[j0] == j) {
              ## Point is pointed to by its predecessor, so make it point to
              ## the same point
              self$h[j] <- self$h[j0]
            } else {
              if (abs(sb/Sb*Sf - sf) < epsilon) {
                ## If projection of forward point is within tolerance of
                ## backward point, don't create a new point
                report(paste("Point", i, "at", sf, "along forward path within",
                                  epsilon, "of projection from", j, "in backward path"))
                self$h[j] <- i
              } else {
                ## Insert a point in the forward path. Point i is ahead.
                f <- (sb/Sb*Sf - sf0)/(sf - sf0)
                k <- self$insertPoint(i0, i, f)
                report(paste("Insert point", k, "at", sb/Sb*Sf,
                                  "along forward path; projection from", j,
                                  "in backward path"))
                self$h[k] <- j
                i0 <- k
                sf0 <- sb/Sb*Sf
              }
            }
          }
          ## Move on
          j0 <- j
          sb0 <- sb
          j <- path.next(j0, self$gb, self$hb)
        }
      }
    }
  )
)

##' @importFrom graphics arrows
##' @export
flatplot.PathOutline <- function(x, axt="n",
                                 xlim=NULL,
                                 ylim=NULL,
                                 paths=FALSE,
                                 ...) {

  NextMethod()

  if (paths) {
    P <- x$P
    gf <- x$gf
    gb <- x$gb
    hf <- x$hf
    hb <- x$hb
    h  <- x$h

    s <- which(!is.na(gf))
    d <- gf[s]
    arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="red", lwd=2)

    s <- which(!is.na(gb))
    d <- gb[s]
    arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="blue")

    s <- which(hf != 1:length(hf))
    d <- hf[s]
    arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="red", lty=2, lwd=2)

    s <- which(hb != 1:length(hb))
    d <- hb[s]
    arrows(P[s,1], P[s,2], P[d,1], P[d,2], col="blue", lty=2, lwd=1)

    s <- which(h!=1:nrow(P))
    d <- h[s]
    arrows(P[s,1], P[s,2], P[d,1], P[d,2])

    points(P[,1], P[,2], cex=2.2, pch=21, bg="white")
    text(P[,1], P[,2], 1:nrow(P), cex=0.7)
  }
}
