Outline <- function(P) {
  t <- triangulate.outline(P, n=NA)
  o <- list()
  class(o) <- "outline"
  o$P <- t$P
  o$gf <- t$gf
  o$gb <- t$gb
  return(o)
}

plot.flat.outline <- function(o, axt="n", ylim=NULL, ...) {
  args <- list(...)
  with(o, {
    s <- which(!is.na(gb))                # source index
    d <- na.omit(gb)                      # destination index
    par(mar=c(1.4, 1.4, 1, 1), mgp=c(2, 0.2, 0), tcl=-0.2)
    
    suppressWarnings(plot(P[s,1], P[s,2],
                          pch=".", xaxt=axt, yaxt=axt, xlab="", ylab="",
                          bty="n", ylim=ylim, ...))
    suppressWarnings(segments(P[s,1], P[s,2], P[d,1], P[d,2], ...))
  })
}
