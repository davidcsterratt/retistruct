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
  args$datapoints <- NULL
  plot.flat.outline(d, axt=axt, ylim=ylim, ...=as.pairlist(args))
  if (plot.datapoints) {
    with(d, {
      for(col in names(Ds)) {
        points(Ds[[col]][,1], Ds[[col]][,2], col=cols[[col]], pch=20,cex=0.5, ...=as.pairlist(args))
      }
    })
  }
}
