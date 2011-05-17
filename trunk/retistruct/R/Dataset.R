Dataset <- function(o, dataset, Ds, Ss, cols, raw) {
  d <- o
  class(d) <- c("dataset", class(o))
  d$dataset <- dataset
  d$Ds <- Ds
  d$Ss <- Ss
  d$cols <- cols
  d$raw <- raw
  return(d)
}

plot.flat.dataset <- function(d, axt="n", ylim=NULL, ...) {
  args <- list(...)
  plot.datapoints <- is.null(args$datapoints) || args$datapoints
  plot.landmarks <- is.null(args$landmarks) || args$landmarks

  plot.flat.outline(d, axt=axt, ylim=ylim, ...)
  if (plot.datapoints) {
    with(d, {
      for(col in names(Ds)) {
        suppressWarnings(points(Ds[[col]][,1], Ds[[col]][,2], col=cols[[col]], pch=20,cex=0.5, ...))
      }
    })
  }
  if (plot.landmarks) {
    with(d, {
      if (length(Ss) > 0) {
        for(i in 1:length(Ss)) {
          suppressWarnings(lines(Ss[[i]][,1], Ss[[i]][,2], ...))
        }
      }
    })
  }
}
