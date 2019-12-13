##' Class to triangulate \link{Fragment}s
##'
##' @description A TriangulatedFragment contains a function to create a
##'   triangulated mesh over an fragment, and fields to hold the mesh
##'   information.
##' @import ttutils
##' @importFrom geometry Unique 
##' @author David Sterratt
##' @export
TriangulatedFragment <- R6Class("TriangulatedFragment",
  inherit = Fragment,
  public = list(
    ##' @field T 3 column matrix in which each row contains IDs of
    ##'   points of each triangle
    T = NULL,
    ##' @field A Area of each triangle in the  mesh - has same number of
    ##'   elements as there are rows of \code{T}
    A = NULL,
    ##' @field Cu 2 column matrix in which each row contains IDs of
    ##'   points of edge in mesh
    Cu = NULL,
    ##' @field L Length of each edge in the mesh - has same number of
    ##'   elements as there are rows of \code{Cu}
    L = NULL,
    ##' @field A.signed Signed area of each triangle generated using
    ##'  \code{\link{tri.area.signed}}. Positive sign indicates points are
    ##'  anticlockwise direction; negative indicates clockwise.
    A.signed = NULL,
    ##' @description Constructor
    ##' @param fragment \link{Fragment} to triangulate
    ##' @param n Minimum number of points in the triangulation
    ##' @param suppress.external.steiner If \code{TRUE} prevent the
    ##'   addition of points in the outline. This happens to maintain
    ##'   triangle quality.
    ##' @param report Function to report progress
    initialize=function(fragment, n=200,
                        suppress.external.steiner=FALSE,
                        report=message) {
      P <- fragment$P[,1:2]
      g <- fragment$gf
      self$h <- fragment$h

      if (is.null(self$h)) {
        self$h <- 1:nrow(P)
      }
      ## By default, segments are outline of points in order
      S <- cbind(1:nrow(P), c(2:nrow(P), 1))
      if (!is.null(g)) {
         S <- pointers2segments(g)
      }
      ## Make initial triangulation
      out <- RTriangle::triangulate(RTriangle::pslg(P=P, S=S),
                                    Y=TRUE, j=TRUE, Q=TRUE)

      ## It can be that there are crossovers in the segments. The
      ## triangulate() routine will reveal this as segments that are not
      ## on a boundary. We get rid of these segments by re-triangulating,
      ## only using boundary segments
      out <- RTriangle::triangulate(RTriangle::pslg(P=out$P, S=out$S[out$SB==1,]),
                                    Y=TRUE, j=TRUE, Q=TRUE)
      
      ## Sometimes a point exists which only belongs to one segment. The
      ## point to which it is connected is itself connected by three
      ## segments. We want to get rid of these points, and the easiest way
      ## is to triangulate without the naughty points.
      i.bad <- which(table(out$S)==1)
      if (length(i.bad) > 0) {
        warning(paste("Bad points:", paste(i.bad, collapse=" ")))
        out <- RTriangle::triangulate(RTriangle::pslg(P=P[-i.bad,], S=S),
                                      Y=TRUE, j=TRUE, Q=TRUE)
      }

      ## Now determine the area
      self$A.tot <- sum(tri.area(cbind(out$P, 0), out$T))

      ## Produce refined triangulation
      P <- out$P
      S <- out$S
      if (!is.na(n)) {
        out <- RTriangle::triangulate(RTriangle::pslg(P=P, S=S), a=self$A.tot/n, q=20,
                                      Y=suppress.external.steiner, j=TRUE, Q=TRUE)
      }
      if (any(P != P[1:nrow(P),])) {
        stop("Points changed in triangulation")
      }
      P <- out$P
      T <- out$T

      ## Create pointers from segments

      ## To ensure the correct orientaion, we use the fact that the
      ## triangles are all anticlockwise in orinentation, and that the
      ## orientation of the first row of the segment matrix determines the
      ## orientation of all the other rows.

      ## We therefore find the triangle which contains the first segment
      S <- out$S
      T1 <- which(apply(T, 1, function(x) {all(S[1,] %in% x)}))

      ## Then find out which of the vertices in the triangle is not the
      ## one we need
      i <- which((T[T1,] %in% S[1,]) == FALSE)
      if (i == 3) S[1,] <- T[T1,c(1,2)]
      if (i == 2) S[1,] <- T[T1,c(3,1)]
      if (i == 1) S[1,] <- T[T1,c(2,3)]

      ## Now create the pointers from the segments
      gf <- segments2pointers(S)
      Rset <- na.omit(gf)
      
      ## Derive edge matrix from triangulation
      Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
      Cu <- Unique(Cu, TRUE)

      ## If we are in the business of refining triangles (i.e. specifying
      ## n), remove lines which join non-ajancent parts of the outline
      if (!is.na(n)) {
        for (i in 1:nrow(Cu)) {
          C1 <- Cu[i,1]
          C2 <- Cu[i,2]
          if (all(Cu[i,] %in% Rset)) {
            if (!((C1 == gf[C2]) ||
                  (C2 == gf[C1]))) {
              ## Find triangles containing the line
              ## segments(P[C1,1], P[C1,2], P[C2,1], P[C2,2], col="yellow")
              Tind <- which(apply(T, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
              report(paste("Non-adjacent points in rim connected by line:", C1, C2))
              report(paste("In triangle:", paste(Tind, collapse=", ")))
              ## Find points T1 & T2 in the two triangles which are
              ## not common with the edge
              T1 <- setdiff(T[Tind[1],], Cu[i,])
              T2 <- setdiff(T[Tind[2],], Cu[i,])
              report(paste("Other points in triangles:", T1, T2))
              ## Create a new point at the centroid of the four vertices
              ## C1, C2, T1, T2
              p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
              P <- rbind(P, p)
              n <- nrow(P)
              ## Remove the two old triangles, and create the four new ones
              T[Tind[1],] <- c(n, C1, T1)
              T[Tind[2],] <- c(n, C1, T2)
              T <- rbind(T,
                         c(n, C2, T1),
                         c(n, C2, T2))
            }
          }
        }

        ## Add the new points to the correspondances vector
        self$h <- c(self$h, (length(self$h)+1):nrow(P))

        ## Create the edge matrix from the triangulation
        Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
        Cu <- Unique(Cu, TRUE)
      }

      ## Swap orientation of triangles which have clockwise orientation
      self$A.signed <- tri.area.signed(cbind(P, 0), T)
      T[self$A.signed<0,c(2,3)] <- T[self$A.signed<0,c(3,2)]
      self$A <- abs(self$A.signed)
      self$A.tot <- sum(self$A)

      ## Find lengths of connections
      self$L <- vecnorm(P[Cu[,1],] - P[Cu[,2],])

      ## Check there are no zero-length lines
      if (any(self$L==0)) {
        warning("zero-length lines")
      }

      ## Pad gf to number of points
      if (length(gf) != nrow(P)) {
        gf[(length(gf)+1):nrow(P)] <- NA
      }

      ## Backward pointer
      gb <- gf
      gb[na.omit(gf)] <- which(!is.na(gf))
      
      self$P <- P
      self$T <- T
      self$gf <- gf
      self$gb <- gb
      self$Cu <- Cu
    }
  )
  )
  

## Convert a matrix containing on each line the indices of the points
## forming a segment, and convert this to two sets of ordered pointers
segments2pointers <- function(S) {
  g <- c()
  j <- 1                                # Row of S
  k <- 1                                # Column of S
  while(nrow(S) > 0) {
    i <- S[j,3-k]                       # i is index of the next point
    g[S[j,k]] <- i                      # Set the pointer to i
    S <- S[-j,,drop=FALSE]              # We have used this row of S
    if (nrow(S) == 0) {
      return(g)
    }
    j <- which(S[,1] == i)            # Is i in the first column of S?
    if (length(j) > 1) {
      stop("The segment list is not valid as it contains an element more than twice.")
    }
    if (length(j)) {              # If so, set the current column to 1
      k <- 1
    } else {
      j <- which(S[,2] == i) # Otherwise, look for i in the second column
      k <- 2
      if (!length(j)) {
        stop(paste("No matching index for point", i, "in S."))
        return(NULL)
      }
    }
  }
  return(g)
}

## Convert a set of ordered pointers to a matrix containing on each
## line the indices of the points forming a segment
pointers2segments <- function(g) {
  S1 <- which(!is.na(g))
  S2 <- g[S1]
  return(cbind(S1, S2))
}
