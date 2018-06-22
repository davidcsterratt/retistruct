##' This stitches together the incisions and tears by inserting new
##' points in the tears and creating correspondences between new
##' points.
##'
##' @title Stitch together incisions and tears in an AnnotatedOutline
##' @return
##' \item{\code{Rset}}{the set of points on the rim}
##' \item{\code{i0}}{the index of the landmark}
##' \item{\code{P}}{a new set of meshpoints}
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
    initialize = function(...) {
      super$initialize(...)
      rs <- self$getRimSet()
      self$hf <- rep(NA, nrow(self$P))
      self$hf[rs] <- rs
      self$hb <- rep(NA, nrow(self$P))
      self$hb[rs] <- rs
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

      V0 <- self$tears[,"V0"]
      VF <- self$tears[,"VF"]
      VB <- self$tears[,"VB"]
      TFset <- r$TFset
      TBset <- r$TBset
      gf <- self$gf
      gb <- self$gb
      hf <- r$hf
      hb <- r$hb
      h <- r$h
      
      ## Insert points on the backward tears corresponding to points on
      ## the forward tears
      sF <-     private$stitchInsertPoints(V0, VF, V0, VB, TFset, TBset,
                                           gf, gb, hf, hb, h,
                                           "Forwards")

      ## Insert points on the forward tears corresponding to points on
      ## the backward tears
      sB <- with(sF,
                 private$stitchInsertPoints(V0, VB, V0, VF, TBset, TFset,
                                            gb, gf, hb, hf, h,
                                            "Backwards"))
      ## Extract data from object
      self$gf <- sB$gb
      self$gb <- sB$gf
      self$hf <- sB$hb
      self$hb <- sB$hf
      h <- sB$h

      ## Link up points on rim
      h[self$Rset] <- hf[self$Rset]
      
      ## Make sure that there are no chains of correspondences
      while (!all(h==h[h])) {
        h <- h[h]
      }
      self$h <- h
      self$TFset <- TFset
    }
  ),
  private = list(
    ## Inner function responsible for inserting the points
    stitchInsertPoints = function(VF0, VF1, VB0, VB1,
                                  TFset, TBset,
                                  gf, gb,
                                  hf, hb, h,
                                  dir,
                                  stitchType="Tear") {
      M <- length(VF0)                       # Number of tears
      ## Iterate through tears to insert new points
      for (j in 1:M) {
        ## Compute the total path length along each side of the tear
        Sf <- path.length(VF0[j], VF1[j], gf, hf, self$getPointsScaled())
        Sb <- path.length(VB0[j], VB1[j], gb, hb, self$getPointsScaled())
        message(paste(stitchType, j, ": Sf =", Sf, "; Sb =", Sb))

        ## For each point in the forward path, create one in the backwards
        ## path at the same fractional location
        message(paste("  ", dir, " path", sep=""))
        for (i in setdiff(TFset[[j]], c(VF0[j], VF1[j]))) {
          sf <- path.length(VF0[j], i, gf, hf, self$getPointsScaled())
          ## If the point isn't at the apex, insert a point
          if (sf > 0) {
            message(paste("    i =", i,
                          "; sf/Sf =", sf/Sf,
                          "; sf =", sf))
            for (k in TBset[[j]]) {
              sb <- path.length(VB0[j], k, gb, hb, self$getPointsScaled())
              message(paste("      k =", format(k, width=4),
                            "; sb/Sb =", sb/Sb,
                            "; sb =", sb))
              if (sb/Sb > sf/Sf) {
                break;
              }
              k0 <- k
              sb0 <- sb
            }

            ## If this point does already not point to another, create
            ## a new point and link to it
            if ((hf[i] == i)) {
              f <- (sf/Sf*Sb - sb0)/(sb - sb0)
              message(paste("      Creating new point: f =", f))
              PXY <- self$getPoints()
              p <-
                (1 - f)*PXY[k0,] +
                f*      PXY[k,]
              ## Find the index of any row of P that matches p
              n <- anyDuplicated(rbind(PXY, p),
                                 fromLast=TRUE)
              if (n == 0) {
                ## If the point p doesn't exist
                n <- self$addPoints(p) # n is Index of new point
                ## Update forward and backward pointers
                gb[n]     <- k
                gf[n]     <- gf[k]
                gb[gf[k]] <- n
                gf[k]     <- n

                ## Update correspondences
                hf[n] <- n
                hb[n] <- n
                h[i] <- n
                h[n] <- n
                message(paste("      Point", i, "points to point", n))
              } else {
                message(paste("      Point", n, "already exists"))
                h[i] <- n
                h[n] <- n
                message(paste("      Point", i, "points to point", n))
              }
            } else {
              ## If not creating a point, set the point to point to the forward pointer 
              h[i] <- hf[i]
            }
          }
        } 
      }
      return(list(hf=hf, hb=hb, gf=gf, gb=gb, h=h))
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
