##' This stitches together the incisions and tears by inserting new
##' points in the tears and creating correspondences between new
##' points.
##'
##' @title Stitch together incisions and tears in an AnnotatedOutline
##' @return
##' \item{\code{Rset}}{the set of points on the rim}
##' \item{\code{i0}}{the index of the landmark}
##' \item{\code{P}}{a new set of mesh points}
##' \item{\code{V0}}{indices of the apex of each tear}
##' \item{\code{VF}}{indices of the forward vertex of each tear}
##' \item{\code{VB}}{indices of the backward vertex of each tear}
##' \item{\code{TFset}}{list containing indices of points in each forward tear}
##' \item{\code{TBset}}{list containing indices of points in each backward tear}
##' \item{\code{gf}}{new forward pointer list}
##' \item{\code{gb}}{new backward pointer list}
##' \item{\code{h}}{correspondence mapping}
##' \item{\code{hf}}{correspondence mapping in forward direction for
##' points on boundary}
##' \item{\code{hb}}{correspondence mapping in backward direction for
##' points on boundary}
##' @author David Sterratt
##' @export
StitchedOutline <- R6Class("StitchedOutline",
  inherit = TriangulatedOutline,
  public = list(
    Rset = NULL,
    TFset = NULL,
    epsilon = NA,
    initialize = function(...) {
      super$initialize(...)
      if (!is.null(self$P) & (nrow(self$getPointsScaled()) > 0)) {
        rs <- self$getRimSet()
        self$hf <- rep(NA, nrow(self$P))
        self$hf[rs] <- rs
        self$hb <- rep(NA, nrow(self$P))
        self$hb[rs] <- rs
        ## Theoretically the maximum tolerance should be half of the
        ## minimum distance between points. We'll make it this, or 0.01%
        ## of the total outline length, whichever is smaller.
        self$epsilon <- min(min(self$getOutlineLengths())/4,
                            sum(self$getRimLengths())*0.01/100)
      }
    },
    stitchTears = function() {
      r <- self$computeTearRelationships(self$tears)

      if (length(r$TFset) == 0) {
        return(NULL)
      }
      ## If not set, set the landmark marker index. Otherwise
      ## check it
      self$Rset <- r$Rset
      if (!(self$i0 %in% self$Rset)) {
        stop(paste("Fixed Point", self$i0, "is not in rim points:",
                   paste(self$Rset, collapse=", ")))
      }
      
      self$hf <- r$hf
      self$hb <- r$hb
      self$h <- r$h
      for (i in 1:nrow(self$tears)) {
        self$stitchSubpaths(self$tears[i,"V0"], self$tears[i,"VF"],
                            self$tears[i,"V0"], self$tears[i,"VB"],
                            epsilon=self$epsilon)
      }
      
      ## Link up points on rim
      self$h[self$Rset] <- self$hf[self$Rset]
      
      ## Make sure that there are no chains of correspondences
      while (!all(self$h==self$h[self$h])) {
        self$h <- self$h[self$h]
      }
      ## self$h <- h
      self$TFset <- r$TFset
    }
  )
)

##' Plot flat \code{\link{StitchedOutline}}. If the optional argument
##' \code{stitch} is \code{TRUE} the user markup is displayed.
##'
##' @title Flat plot of AnnotatedOutline
##' @param x \code{\link{AnnotatedOutline}} object
##' @param axt whether to plot axes
##' @param xlim x-limits
##' @param ylim y-limits
##' @param stitch If \code{TRUE}, plot stitch
##' @param lwd Line width
##' @param ... Other parameters
##' @method flatplot StitchedOutline
##' @author David Sterratt
##' @export
flatplot.StitchedOutline <- function(x, axt="n",
                                     xlim=NULL, ylim=NULL,
                                     stitch=TRUE, lwd=1,
                                     ...) {
  NextMethod()

  if (stitch) {
    for (TF in x$TFset) {
      lines(x$P[TF,], col=getOption("TF.col"), lwd=lwd)
    }
    for (TB in x$TBset) {
      lines(x$P[TB,], col=getOption("TB.col"), lwd=lwd)
    }
    for (j in 1:length(x$h)) {
      if (x$h[j] != j) {
        lines(x$P[c(j, x$h[j]),], col=getOption("stitch.col"), lwd=lwd)
      }
    }
    ## if (!is.null(x$hf)) {
    ##   for (j in 1:length(x$hf)) {
    ##     if (x$hf[j] != j) {
    ##       lines(x$P[c(j, x$hf[j]),], col=getOption("stitch.col"), lwd=lwd)
    ##     }
    ##   }
    ##   for (j in 1:length(x$hb)) {
    ##     if (x$hb[j] != j) {
    ##       lines(x$P[c(j, x$hb[j]),], col=getOption("V.stitch.col"), lwd=lwd)
    ##     }
    ##   }
    ## }
  }
}
