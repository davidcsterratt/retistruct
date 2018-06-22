##' Construct an outline object. This sanitises the input points
##' \code{P}, as described below.
##'
##' @title Fragment constructor
##' @return An \code{Fragment} object containing the following:
##' \item{\code{P}}{A N-by-2 matrix of points of the \code{Outline} arranged in anticlockwise order}
##' \item{\code{gf}}{For each row of \code{P}, the index of \code{P} that is next in the outline travelling anticlockwise (forwards)}
##' \item{\code{gb}}{For each row of \code{P}, the index of \code{P} that is next in the outline travelling clockwise (backwards)}
##' \item{\code{h}}{For each row of \code{P}, the correspondence of that point (which will be to itself initially)}
##' @author David Sterratt
Fragment <- R6Class("Fragment",
  public = list(
    P = NULL,
    gf = NULL,
    gb = NULL,
    h = NULL,
    A.tot = NULL,
    initialize = function(P = NULL, gf = NULL, gb=NULL, h=NULL) {
      self$P = P
      self$gf = gf
      self$gb = gb
      self$h = h
    },
    initializeFromPoints = function(P) {
      if (is.null(P)) {
        stop("P is NULL; it should be a matrix")
      }
      if (!is.matrix(P)) {
        stop("P is not a matrix; it should be")
      }
      if (ncol(P) != 2) {
        stop("P should have two columns")
      }
      P <- remove.identical.consecutive.rows(P)
      self$P <- simplifyFragment(P)
      t <- TriangulatedFragment$new(self, n=NA)
      if (ncol(t$P) != 2) {
        stop("P should have two columns")
      }
      self$P <- t$P
      self$gf <- t$gf
      self$gb <- t$gb
      self$h <- t$h
      self$A.tot <- t$A.tot
    }
  )
)

##' Simplify a fragment object by removing vertices bordering short
##' edges while not encroaching on any of the outline. At present,
##' this is done by finding concave vertices. It is safe to remove
##' these, at the expense of increasing the area a bit.
##'
##' @title Simplify an outline object by removing short edges
##' @param P points to simplify
##' @param min.frac.length the minimum length as a fraction of the
##' total length of the outline. 
##' @param plot whether to display plotting or not during simplification
##' @return Simplified \code{outline} object
##' @author David Sterratt
simplifyFragment <- function(P, min.frac.length=0.001, plot=FALSE) {
  N <- nrow(P)                        # Number of vertices
  Q <- rbind(P, P[1,])                # Convenience variable
  v <- diff(Q)                         # Vector of each edge
  l <- vecnorm(v)                     # Length of each edge
  ## Compute outer products at each vertex
  e <- extprod3d(cbind(v[c(N, 1:(N-1)),], 0), cbind(v, 0))[,3]
  
  ## Find short edges
  S <- l/sum(l) < min.frac.length

  ## Find indices of points that can be removed.
  ## They have to be concave or colinear (e<=0). And they need to border short edges
  i.rem <- which((e <= 0) & (S | (S[c(N, 1:(N-1))])))

  if (plot) {
    ## Plot short edges...
    plot(P, col="white")
    if (any(S)) {
      segments(P[S,1], P[S,2], P[S,1]+v[S,1], P[S,2] + v[S,2], col="red")
    }
    ## and longer ones.
    if (any(!S)) {
      segments(P[!S,1], P[!S,2], P[!S,1]+v[!S,1], P[!S,2] + v[!S,2], col="black")
    }
    ## Plot colinear, convex and concave points
    points(P[e>0,1], P[e>0, 2], col="green")
    points(P[e==0,1], P[e==0, 2], col="orange")
    points(P[e<0,1], P[e<0, 2], col="blue")
    ## Plot points to remove
    points(P[i.rem,1], P[i.rem, 2], pch="X", col="brown")
  }

  ## If there are any points to remove, remove the first one
  if (length(i.rem) > 0) {
    message(paste("simplifyFragment: Removing vertex", i.rem[1]))
    return(simplifyFragment(P[-i.rem[1],],
                            min.frac.length=min.frac.length, plot=plot))
  } else {
    return(P)
  }
}
