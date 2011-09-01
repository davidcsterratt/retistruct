plot.flat.stitchedOutline <- function(s, axt="n", ylim=NULL, ...) {
  NextMethod()
  args <- list(...)
  plot.stitch <- is.null(args$stitch) || args$stitch

  if (plot.stitch) {
    with(s, {
      points(P[VF,], col=TFcol, pch="+")
      points(P[VB,], col=TBcol, pch="+")
      points(P[V0,], col=Vcol, pch="+")
      for (TF in TFset) {
        suppressWarnings(lines(P[TF,], col=TFcol, ...))
      }
      for (TB in TBset) {
        suppressWarnings(lines(P[TB,], col=TBcol, ...))
      }
      for (j in 1:length(h)) {
        if (h[j] != j) {
          suppressWarnings(lines(P[c(j, h[j]),], col="blue", ...))
        }
      }
      
      for (j in 1:length(hf)) {
        if (hf[j] != j) {
          suppressWarnings(lines(P[c(j, hf[j]),], col="green", ...))
        }
      }
      for (j in 1:length(hb)) {
        if (hb[j] != j) {
          suppressWarnings(lines(P[c(j, hb[j]),], col="green", ...))
        }
      }
    })
  }
}
