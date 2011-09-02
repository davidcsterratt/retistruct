plot.flat.stitchedOutline <- function(s, axt="n", ylim=NULL, ...) {
  NextMethod()
  args <- list(...)
  plot.stitch <- is.null(args$stitch) || args$stitch

  if (plot.stitch) {
    with(s, {
      for (TF in TFset) {
        suppressWarnings(lines(P[TF,], col=getOption("TF.col"), ...))
      }
      for (TB in TBset) {
        suppressWarnings(lines(P[TB,], col=getOption("TB.col"), ...))
      }
      for (j in 1:length(h)) {
        if (h[j] != j) {
          suppressWarnings(lines(P[c(j, h[j]),], col=getOption("stitch.col"), ...))
        }
      }
      
      for (j in 1:length(hf)) {
        if (hf[j] != j) {
          suppressWarnings(lines(P[c(j, hf[j]),], col=getOption("stitch.col"), ...))
        }
      }
      for (j in 1:length(hb)) {
        if (hb[j] != j) {
          suppressWarnings(lines(P[c(j, hb[j]),], col=getOption("V.stitch.col"), ...))
        }
      }
    })
  }
}
