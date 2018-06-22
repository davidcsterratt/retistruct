##' Plot "flat" (unreconstructed) representation of outline
##' @param x \code{\link{Outline}}, \code{\link{AnnotatedOutline}},  \code{\link{StitchedOutline}} &c object
##' @param axt whether to plot axes
##' @param xlim x limits
##' @param ylim y limits
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
flatplot <- function(x, axt="n", xlim=NULL, ylim=NULL, ...) {
  UseMethod("flatplot")
}

##' @export
flatplot.default <- function(x, axt="n", xlim=NULL, ylim=NULL, ...) {
}

##' Plot projection of a reconstructed outline
##' @param r Object such as a \code{\link{ReconstructedOutline}}
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
projection <- function(r, ...) {
  UseMethod("projection")
}

##' @export
projection.default <- function(r, ...) {
  graphics::plot.new()
}

##' Spherical plot of reconstructed outline
##' @param r Object inheriting \code{\link{ReconstructedOutline}}
##' @param ... Parameters depending on class of \code{r}
##' @author David Sterratt
##' @export
sphericalplot <- function(r, ...) {
  UseMethod("sphericalplot")
}

##' @export
sphericalplot.default <- function(r, ...) {
  rgl.clear()
  rgl.bg(color="white")
}

##' Draw the "flat" outline in 3D with depth information 
##' @param r \code{\link{TriangulatedOutline}} object
##' @param ... Parameters depending on class of \code{r}
##' @author David Sterratt
##' @export
depthplot3D <- function(r, ...) {
  UseMethod("depthplot3D")
}

##' @export
depthplot3D.default <- function(r, ...) {
}


##' Plot the fractional change in length of mesh edges
##' 
##' Plot the fractional change in length of mesh edges. The length of
##' each edge in the mesh in the reconstructed object is plotted
##' against each edge in the spherical object. The points are
##' colour-coded according to the amount of log strain each edge is
##' under.
##' @param r \code{\link{ReconstructedOutline}} object
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
lvsLplot <- function(r, ...) {
  UseMethod("lvsLplot")
}

##' @export
lvsLplot.default <- function(r, ...) {
}

