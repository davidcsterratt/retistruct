plot.flat.stitchedOutline <- function(s, axt="n", ylim=NULL, ...) {
  NextMethod()
  args <- list(...)
  plot.stitch <- is.null(args$stitch) || args$stitch

  if (plot.stitch) {
    with(s, {
      points(P[VF,], col="red", pch="+")
      points(P[VB,], col="orange", pch="+")
      points(P[V0,], col="cyan", pch="+")
      for (TF in TFset) {
        suppressWarnings(lines(P[TF,], col="red", ...))
      }
      for (TB in TBset) {
        suppressWarnings(lines(P[TB,], col="orange", ...))
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
