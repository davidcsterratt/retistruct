plot.flat <- function(x, axt="n", ylim=NULL, ...) {
  UseMethod("plot.flat")
}

plot.flat.default <- function(x, axt="n", ylim=NULL, ...) {
}

plot.flat.list <- function(x, axt="n", ylim=NULL, ...) {
  plot.flat.annotatedDataset(x, axt, ylim, ...)
}


