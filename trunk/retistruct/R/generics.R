plot.flat <- function(x, axt="n", ylim=NULL, ...) {
  UseMethod("plot.flat")
}

plot.flat.default <- function(x, axt="n", ylim=NULL, ...) {
}

plot.polar <- function(r, show.grid=TRUE,
                       grid.col="gray", grid.bg="transparent", 
                       grid.int.minor=15, grid.int.major=45, ...) {
  UseMethod("plot.polar")
}

plot.polar.default <- function(r, show.grid=TRUE,
                               grid.col="gray", grid.bg="transparent", 
                               grid.int.minor=15, grid.int.major=45, ...) {
  plot.new()
}


