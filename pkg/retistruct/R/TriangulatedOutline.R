##' Constructor for TriangulatedOutline object.
##'
##' @title TriangulatedOutline object
##' @return TriangulatedOutline object, with extra fields for tears
##'   latitude of rim \code{phi0} and index of fixed point \code{i0}.
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
    T =  matrix(NA, 0, 3),
    Cu = matrix(NA, 0, 2),
    L = NULL,
    A = NULL,
    A.tot = NULL,
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

