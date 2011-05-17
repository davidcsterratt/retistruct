plot.polar.reconstructedOutline <- function(r, show.grid=TRUE,
                                            grid.col="gray",
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45, ...) {

  phi0d <- r$phi0*180/pi
  par(mar=c(1, 1, 1, 1))
  grid.pos <- c(seq(-90, phi0d, by=grid.int.minor), phi0d)
  maxlength <- diff(range(grid.pos))
  plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
       type = "n", axes = FALSE, xlab = "", ylab = "", asp=1)

  ## Plot the grid
  ## Tangnential lines
  angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
  if (show.grid) {
    for (i in seq(length(grid.pos), 1, by = -1)) {
      xpos <- cos(angles) * (grid.pos[i] - grid.pos[1])
      ypos <- sin(angles) * (grid.pos[i] - grid.pos[1])
      if (((grid.pos[i] %% grid.int.major) == 0) || (i == length(grid.pos))) {
        col <- "black"
      } else {
        col <- grid.col
      }
      polygon(xpos, ypos, border = col, col = grid.bg)
    }
  }

  ## Radial lines
  angles <- seq(0, 180-grid.int.minor, by = grid.int.minor)
  col <- rep(grid.col, length(angles))
  col[(angles %% grid.int.major) == 0] <- "black"
  angles <- angles * pi/180
  xpos <- cos(angles) * maxlength
  ypos <- sin(angles) * maxlength
  segments(xpos, ypos, -xpos, -ypos, col=col)

  ## Radial Labels
  labels <- c(seq(-90, phi0d, by=grid.int.major), phi0d)
  label.pos <- labels - min(grid.pos)
  text(label.pos, -maxlength/15, labels)
}
