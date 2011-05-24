plot.polar.reconstructedDataset <- function(r, show.grid=TRUE,
                                            grid.col="gray",
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45, ...) {
  NextMethod()

  args <- list(...)
  plot.datapoints <- is.null(args$datapoints) || args$datapoints
  plot.landmarks <- is.null(args$landmarks) || args$landmarks

  phi0d <- r$phi0*180/pi
  grid.pos <- c(seq(-90, phi0d, by=grid.int.minor), phi0d)
  maxlength <- diff(range(grid.pos))

  ## Radial Labels
  angles <- c(0, 90, 180, 270)*pi/180
  xpos <- cos(angles) * maxlength * 1.05
  ypos <- sin(angles) * maxlength * 1.05
  if (r$side=="Right") {
    text(xpos, ypos, c("N", "D", "T", "V"))
  } else {
    text(xpos, ypos, c("T", "D", "N", "V"))
  }

  ## Datapoints
  if (plot.datapoints) {
    with(r, {
      for (i in 1:length(Dss)) {
        phis    <- Dss[[i]][,"phi"]
        lambdas <- Dss[[i]][,"lambda"]
        xpos <- cos(lambdas) * ((phis * 180/pi) + 90)
        ypos <- sin(lambdas) * ((phis * 180/pi) + 90)
        if (r$DVflip)
          ypos <- -ypos
        suppressWarnings(points(xpos, ypos, col=cols[[names(Dss)[i]]],
                                pch=20, ...))
      }
    })
  }

  ## Landmarks
  if (plot.landmarks) {
    with(r, {
      if (length(Sss) > 0) {
        for (i in 1:length(Sss)) {
          name <- names(Ss)[i]
          col <- ifelse(is.null(name) || (name==""), "default", name)
          phi    <- Sss[[i]][,"phi"]
          lambda <- Sss[[i]][,"lambda"]
          x <- cos(lambda) * ((phi * 180/pi) + 90)
          y <- sin(lambda) * ((phi * 180/pi) + 90)
          if (r$DVflip)
            y <- -y
          suppressWarnings(lines(x, y, col=cols[[col]], ...))
        }
      }
    })
  }

  ## Outline
  with(r, {
    for (TF in r$TFset) {
      ## Convert indicies to the spherical frame of reference
      j <- r$ht[TF]
      ## Plot
      x <- with(r, cos(lambda[j]) * ((phi[j] * 180/pi) + 90))
      y <- with(r, sin(lambda[j]) * ((phi[j] * 180/pi) + 90))
      if (r$DVflip)
          y <- -y
      suppressWarnings(lines(x, y, ...))
    }
  })
}
