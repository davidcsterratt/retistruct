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

nameLandmark <- function(d, i, name) {
  new.names <- rep("", length(d$Ss))
  new.names[i] <- name
  names(d$Ss) <- new.names
  return(d)
}

getLandmarkID <- function(d, name) {
  id <- which(names(d$Ss) == name)
  if (length(id) == 1) {
    return(id)
  } else {
    return(NA)
  }
}

plot.flat.dataset <- function(d, axt="n", ylim=NULL, ...) {
  args <- list(...)
  plot.datapoints <- is.null(args$datapoints) || args$datapoints
  plot.landmarks <- is.null(args$landmarks) || args$landmarks

  NextMethod()
  ## plot.flat.outline(d, axt=axt, ylim=ylim, ...)
  if (plot.datapoints) {
    with(d, {
      for (col in names(Ds)) {
        suppressWarnings(points(Ds[[col]][,1], Ds[[col]][,2],
                                col=cols[[col]], pch=20, ...))
      }
    })
  }
  if (plot.landmarks) {
    with(d, {
      if (length(Ss) > 0) {
        for (i in 1:length(Ss)) {
          name <- names(Ss)[i]
          col <- ifelse(is.null(name) || (name==""), "default", name)
          suppressWarnings(lines(Ss[[i]][,1], Ss[[i]][,2], col=cols[[col]], ...))
        }
      }
    })
  }
}


