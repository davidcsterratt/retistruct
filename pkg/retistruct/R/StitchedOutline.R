##' Class containing functions and data relating to Stitching outlines
##'
##' @description A StitchedOutline contains a function to stitch the
##'   tears and fullcuts, setting the correspondences \code{hf}, \code{hb} and
##'   \code{h}
##'
##' @author David Sterratt
##' @export
StitchedOutline <- R6Class("StitchedOutline",
  inherit = TriangulatedOutline,
  public = list(
    ##' @field Rset the set of points on the rim
    Rset = NULL,
    ##' @field TFset list containing indices of points in each forward tear
    TFset = NULL,
    ##' @field CFset list containing indices of points in each forward cut
    CFset = NULL,
    ##' @field epsilon the minimum distance between points, set
    ##'   automatically
    epsilon = NA,
    ##' @field tearsStitched Boolean indicating if tears have been stitched
    tearsStitched = FALSE,
    ##' @field fullCutsStitched Boolean indicating if full cuts have been stitched
    fullCutsStitched = FALSE,
    ##' @description Constructor
    ##' @param ... Parameters to superclass constructors
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
    ##' @description Stitch together the incisions and tears by inserting new
    ##'   points in the tears and creating correspondences between new
    ##'   points.
    stitchTears = function() {
      self$tearsStitched = TRUE
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
    },
    ##' @description Stitch together the fullcuts by inserting new
    ##'   points in the tears and creating correspondences between new
    ##'   points.
    stitchFullCuts = function() {
      self$fullCutsStitched = TRUE
      r <- self$computeFullCutRelationships(self$fullcuts)

      if (length(r$CFset) == 0) {
        return(NULL)
      }

      ## If not set, set the landmark marker index. Otherwise
      ## check it
      self$Rset <- r$Rset
      if (!(self$i0 %in% self$Rset)) {
        stop(paste("Fixed Point", self$i0, "is not in rim points:",
                   paste(self$Rset, collapse=", ")))
      }

      for (i in 1:nrow(self$fullcuts)) {
        self$stitchSubpaths(self$fullcuts[i,"VF0"], self$fullcuts[i,"VF1"],
                            self$fullcuts[i,"VB0"], self$fullcuts[i,"VB1"],
                            epsilon=self$epsilon)
      }

      ## Link up points on rim
      self$hf[self$Rset] <- r$hf[self$Rset]
      self$hb[self$Rset] <- r$hb[self$Rset]
      ## self$h <- r$h

      self$h[self$Rset] <- self$hf[self$Rset]

      ## Make sure that there are no chains of correspondences
      while (!all(self$h==self$h[self$h])) {
        self$h <- self$h[self$h]
      }
      self$CFset <- r$CFset
    },
    ##' @description Test if the outline has been stitched
    ##' @return Boolean, indicating if the outline has been stitched or not
    isStitched = function() {
      return(self$tearsStitched & self$fullCutsStitched)
    },
    ##' @description Get point IDs of points on boundaries
    ##' @return List of Point IDs of points on the boundaries.
    ##' If the outline has been stitched,
    ##' the point IDs in each
    ##' element of the list will be ordered in the direction of the
    ##' forward pointer, and the boundary that is longest will be
    ##' named as \code{Rim}. If the outline has not been stitched,
    ##' the list will have one element named \code{Rim}.
    getBoundarySets = function() {
      Bsets <- super$getBoundarySets()
      if (!self$isStitched()) {
        return(Bsets)
      }
      UBset <- Bsets[["Rim"]]
      ## Now separate out the Bsets
      Bsets  <- list()
      Blengths <- c()
      P <- self$getPointsScaled()
      while(length(UBset) > 0) {
        i  <- UBset[1]
        j  <- path.next(i, self$gf, self$hf)
        B12 <- path(i, j, self$gf, self$hf)
        B21 <- path(j, i, self$gf, self$hf)
        BL12 <- path.length(i, j,
                            self$gf, self$hf, P)
        BL21 <- path.length(j, i,
                            self$gf, self$hf, P)
        Bset <- c(B12[-1], B21[-1])
        Bsets <- c(Bsets, list(Bset))
        Blengths <- c(Blengths, BL12 + BL21)
        UBset <- setdiff(UBset, Bset)
      }
      Bsets <- name.list(Bsets)
      names(Bsets)[which.max(Blengths)] <- "Rim"
      return(Bsets)
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
