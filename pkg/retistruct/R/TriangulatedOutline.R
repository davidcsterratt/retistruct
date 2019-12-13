##' Class containing functions and data relating to Triangulation
##'
##' @description A TriangulatedOutline contains a function to create a
##'   triangulated mesh over an outline, and fields to hold the mesh
##'   information. Note that areas and lengths are all scaled using
##'   the value of the \code{scale} field.
##'
##' @author David Sterratt
##' @export
##' @importFrom R6 R6Class
##' @examples
##' P <- rbind(c(1,1),   c(2,1),  c(2,-1),
##'            c(1,-1),  c(1,-2), c(-1,-2),
##'            c(-1,-1), c(-2,-1),c(-2,1),
##'            c(-1,1),  c(-1,2), c(1,2))
##' o <- TriangulatedOutline$new(P)
##' o$addTear(c(3, 4, 5))
##' o$addTear(c(6, 7, 8))
##' o$addTear(c(9, 10, 11))
##' o$addTear(c(12, 1, 2))
##' flatplot(o)
TriangulatedOutline <- R6Class("TriangulatedOutline",
  inherit = AnnotatedOutline,
  public = list(
    ##' @field T 3 column matrix in which each row contains IDs of
    ##'   points of each triangle
    T =  matrix(NA, 0, 3),
    ##' @field A Area of each triangle in the mesh - has same number of
    ##'   elements as there are rows of \code{T}
    A = NULL,
    ##' @field A.tot Total area of the mesh
    A.tot = NULL,
    ##' @field Cu 2 column matrix in which each row contains IDs of
    ##    points of edge in mesh
    Cu = matrix(NA, 0, 2),
    ##' @field L Length of each edge in the mesh - has same number of
    ##'   elements as there are rows of \code{Cu}
    L = NULL,
    ##' @description Triangulate (mesh) outline
    ##' @param n Desired number of points in mesh
    ##' @param suppress.external.steiner Boolean variable describing
    ##'   whether to insert external Steiner points - see
    ##'   \link{TriangulatedFragment}
    triangulate = function(n=200, suppress.external.steiner=FALSE) {
      self$T <- matrix(NA, 0, 3)
      self$Cu <- matrix(NA, 0, 2)
      t <- TriangulatedFragment$new(self,
                                    n=n,
                                    suppress.external.steiner=suppress.external.steiner,
                                    report=report)
      pids <- self$addPoints(t$P)
      if (length(t$gf) != length(pids)) {
        stop("Number of indices is not equal to number of pids supplied")
      }
      self$mapTriangulatedFragment(t, pids)

      ## Find areas and lengths of connections
      P <- self$getPointsScaled()
      self$A <- tri.area(cbind(P, 0), self$T)
      self$L <- vecnorm(P[self$Cu[,1],] - P[self$Cu[,2],])
      self$A.tot <- sum(self$A)
    },
    ##' @description Map the point IDs of a \link{TriangulatedFragment} on the
    ##'   point IDs of this Outline
    ##' @param fragment \link{TriangulatedFragment} to map
    ##' @param pids Point IDs in TriangulatedOutline of points in \link{TriangulatedFragment}
    mapTriangulatedFragment = function(fragment, pids) {
      self$mapFragment(fragment, pids)
      if (!is.null(fragment$T)) {
        self$T <-  rbind(self$T,  matrix(pids[fragment$T], ncol=3))
        self$Cu <- rbind(self$Cu, matrix(pids[fragment$Cu], ncol=2))
      }
    }
  )
)

##' Plot flat \code{\link{TriangulatedOutline}}.
##' @param x \code{\link{TriangulatedOutline}} object
##' @param axt whether to plot axes
##' @param xlim x-limits
##' @param ylim y-limits
##' @param mesh If \code{TRUE}, plot mesh
##' @param ... Other plotting parameters
##' @method flatplot TriangulatedOutline
##' @importFrom geometry trimesh
##' @author David Sterratt
##' @export
flatplot.TriangulatedOutline <- function(x, axt="n",
                                         xlim=NULL, ylim=NULL,
                                         mesh=TRUE,
                                         ...) {
  NextMethod()

  if (mesh) 
    trimesh(x$T, x$P, col="grey", add=TRUE)
}

##' @import rgl
##' @method depthplot3D TriangulatedOutline
##' @export
depthplot3D.TriangulatedOutline <- function(r, ...) {
  rgl.clear()
  if (nrow(r$T) == 0) {
    warning("Outline not yet triangulated - no depthplot will show")
  }
  P <- r$getPointsScaled()
  triangles3d(matrix(P[t(r$T),"X"], nrow=3),
              matrix(P[t(r$T),"Y"], nrow=3),
              matrix(P[t(r$T),"Z"], nrow=3),
              color="red", alpha=1)
}

